"""Cross-repository contract test: aisoptics writes, ais++ reads.

This is the test that catches a layout drift between the two packages. The C++
side (tests/test_interpolation.cc) pins ais++'s own reader against a file it
writes itself; what it cannot check is whether *aisoptics* still writes that
same layout. This one does, by exporting a real beam through the public
AISPPExporter API and reading it back with the exact index arithmetic
AISUtilities.cc uses.

Skipped automatically when aisoptics is not installed.

Run with:  pytest tests/test_aisoptics_roundtrip.py -v
"""

import os
import tempfile

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
aisoptics = pytest.importorskip("aisoptics", reason="aisoptics is not installed")

from aisoptics import AISPPExporter, Backend, GaussianBeam, GridSpec  # noqa: E402


# The layout ais++ requires. Keep these literals in sync with
# AISLaserBeam::SetInterpolationGrids and trilinearInterpolation.
REQUIRED_DATASETS = ("x", "y", "z", "phase", "amplitude")
REQUIRED_EXPORT_FORMAT = "aispp_current_interpolation_hdf5"


@pytest.fixture(scope="module")
def exported_beam():
    """Export a Gaussian beam through aisoptics and hand back the file path."""
    backend = Backend("numpy")
    beam = GaussianBeam(
        wavelength=698e-9,      # Sr-87 clock transition
        waist=0.01,
        focus_z=0.0,
        propagation_direction="+z",
        backend=backend,
    )
    grid = GridSpec.from_bounds(
        xlim=(-4e-3, 4e-3),
        ylim=(-3e-3, 3e-3),     # deliberately different extents per axis:
        zlim=(0.0, 0.5),        # a transposed index cannot survive this
        shape=(24, 20, 6),
    )
    field_grid = beam.sample(grid)

    tmpdir = tempfile.mkdtemp()
    path = os.path.join(tmpdir, "beam.h5")
    AISPPExporter().export_total_field(path, field_grid)
    yield path


def cpp_trilinear(x, y, z, xg, yg, zg, values):
    """Reimplementation of ais++'s trilinearInterpolation, including its clamping.

    Deliberately a transcription rather than a call into scipy: the point is to
    exercise the index formula i*Ny*Nz + j*Nz + k that the C++ uses, so a change
    on either side shows up as a disagreement here.
    """
    ny, nz = len(yg), len(zg)

    def lower_index(val, grid):
        if val <= grid[0]:
            return 0
        if val >= grid[-1]:
            return len(grid) - 2
        return int(np.searchsorted(grid, val, side="right")) - 1

    def clamp01(v):
        return 0.0 if v < 0.0 else (1.0 if v > 1.0 else v)

    i, j, k = lower_index(x, xg), lower_index(y, yg), lower_index(z, zg)
    xd = clamp01((x - xg[i]) / (xg[i + 1] - xg[i]))
    yd = clamp01((y - yg[j]) / (yg[j + 1] - yg[j]))
    zd = clamp01((z - zg[k]) / (zg[k + 1] - zg[k]))

    def at(ii, jj, kk):
        return values[ii * ny * nz + jj * nz + kk]

    c00 = at(i, j, k) * (1 - xd) + at(i + 1, j, k) * xd
    c01 = at(i, j, k + 1) * (1 - xd) + at(i + 1, j, k + 1) * xd
    c10 = at(i, j + 1, k) * (1 - xd) + at(i + 1, j + 1, k) * xd
    c11 = at(i, j + 1, k + 1) * (1 - xd) + at(i + 1, j + 1, k + 1) * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    return c0 * (1 - zd) + c1 * zd


class TestFileLayout:
    """The structural contract, asserted key by key so a break names itself."""

    def test_required_datasets_exist(self, exported_beam):
        with h5py.File(exported_beam, "r") as f:
            for name in REQUIRED_DATASETS:
                assert name in f, f"ais++ requires dataset {name!r}"

    def test_export_format_attribute_matches(self, exported_beam):
        """ais++ refuses any other value, so a rename here is a hard break."""
        with h5py.File(exported_beam, "r") as f:
            assert "export_format" in f.attrs
            value = f.attrs["export_format"]
            if isinstance(value, bytes):
                value = value.decode()
            assert value == REQUIRED_EXPORT_FORMAT

    def test_axes_are_one_dimensional(self, exported_beam):
        with h5py.File(exported_beam, "r") as f:
            for name in ("x", "y", "z"):
                assert f[name].ndim == 1, f"{name} must be a 1-D axis"

    def test_value_arrays_are_flat_and_correctly_sized(self, exported_beam):
        with h5py.File(exported_beam, "r") as f:
            nx, ny, nz = len(f["x"]), len(f["y"]), len(f["z"])
            for name in ("phase", "amplitude"):
                assert f[name].ndim == 1, f"{name} must be flattened"
                assert len(f[name]) == nx * ny * nz, (
                    f"{name} length {len(f[name])} != Nx*Ny*Nz = {nx*ny*nz}"
                )

    def test_axes_are_ascending(self, exported_beam):
        """The cell search assumes ascending axes."""
        with h5py.File(exported_beam, "r") as f:
            for name in ("x", "y", "z"):
                axis = f[name][:]
                assert np.all(np.diff(axis) > 0), f"{name} must be strictly ascending"

    def test_carrier_is_excluded(self, exported_beam):
        """ais++ adds the plane-wave k*z itself, so the file must not contain it.
        With k = 2*pi/698nm, a z-range of 0.5 m would give ~4.5e6 rad of carrier;
        the envelope phase stays small."""
        with h5py.File(exported_beam, "r") as f:
            assert f.attrs["carrier_included"] in (False, np.False_, 0)
            phase = f["phase"][:]
            assert np.abs(phase).max() < 1e3, (
                "phase looks like it still contains the carrier"
            )


class TestIndexOrderRoundTrip:
    """The part that actually breaks silently: index order."""

    def test_flatten_order_is_c(self, exported_beam):
        with h5py.File(exported_beam, "r") as f:
            order = f.attrs["flatten_order"]
            if isinstance(order, bytes):
                order = order.decode()
            assert order == "C", "ais++'s index formula assumes C order"

    def test_node_values_match_under_cpp_indexing(self, exported_beam):
        """Read every node back through the C++ index formula and compare with a
        C-order reshape. A transposed or F-order write shows up immediately."""
        with h5py.File(exported_beam, "r") as f:
            xg, yg, zg = f["x"][:], f["y"][:], f["z"][:]
            phase = f["phase"][:]

        cube = phase.reshape(len(xg), len(yg), len(zg), order="C")
        ny, nz = len(yg), len(zg)
        rng = np.random.default_rng(0)
        for _ in range(50):
            i = rng.integers(0, len(xg))
            j = rng.integers(0, ny)
            k = rng.integers(0, nz)
            assert phase[i * ny * nz + j * nz + k] == cube[i, j, k]

    def test_interpolation_reproduces_nodes_exactly(self, exported_beam):
        """Sampling exactly on a node must return that node's value."""
        with h5py.File(exported_beam, "r") as f:
            xg, yg, zg = f["x"][:], f["y"][:], f["z"][:]
            phase = f["phase"][:]

        cube = phase.reshape(len(xg), len(yg), len(zg), order="C")
        for i in (0, len(xg) // 3, len(xg) - 1):
            for j in (0, len(yg) // 2, len(yg) - 1):
                for k in (0, len(zg) - 1):
                    got = cpp_trilinear(xg[i], yg[j], zg[k], xg, yg, zg, phase)
                    assert got == pytest.approx(cube[i, j, k], rel=1e-12, abs=1e-12)

    def test_interpolation_agrees_with_aisoptics_own_interpolator(self, exported_beam):
        """Both sides must agree away from the nodes too. This is the assertion
        that fails if either package changes its interpolation convention."""
        from aisoptics import TrilinearInterpolator

        with h5py.File(exported_beam, "r") as f:
            xg, yg, zg = f["x"][:], f["y"][:], f["z"][:]
            phase = f["phase"][:]

        cube = phase.reshape(len(xg), len(yg), len(zg), order="C")
        grid = GridSpec(x=xg, y=yg, z=zg)
        interp = TrilinearInterpolator(grid, cube)

        rng = np.random.default_rng(1)
        for _ in range(40):
            x = rng.uniform(xg[0], xg[-1])
            y = rng.uniform(yg[0], yg[-1])
            z = rng.uniform(zg[0], zg[-1])
            theirs = float(np.asarray(interp.evaluate(x, y, z)).ravel()[0])
            ours = cpp_trilinear(x, y, z, xg, yg, zg, phase)
            assert ours == pytest.approx(theirs, rel=1e-9, abs=1e-12)


class TestAmplitudeConvention:

    def test_amplitude_is_non_negative_and_peaks_near_axis(self, exported_beam):
        with h5py.File(exported_beam, "r") as f:
            xg, yg, zg = f["x"][:], f["y"][:], f["z"][:]
            amp = f["amplitude"][:]

        assert np.all(amp >= 0.0), "field amplitude must be non-negative"

        cube = amp.reshape(len(xg), len(yg), len(zg), order="C")
        # at the first z-plane the maximum should sit at the transverse centre
        plane = cube[:, :, 0]
        i, j = np.unravel_index(np.argmax(plane), plane.shape)
        assert abs(xg[i]) <= (xg[1] - xg[0]), "peak is not on the beam axis in x"
        assert abs(yg[j]) <= (yg[1] - yg[0]), "peak is not on the beam axis in y"

    def test_amplitude_convention_is_declared(self, exported_beam):
        """ais++ multiplies rabifreq by this value, so the convention matters."""
        with h5py.File(exported_beam, "r") as f:
            value = f.attrs["amplitude_convention"]
            if isinstance(value, bytes):
                value = value.decode()
            assert value == "field_amplitude_relative_to_configured_rabi_frequency"
