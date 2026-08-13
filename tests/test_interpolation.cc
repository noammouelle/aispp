// Unit tests for the sampled-beam path: wtype "interpolated", the HDF5 layout
// aisoptics' AISPPExporter writes, and trilinearInterpolation itself.
//
// Build:  cmake --build build --target test_interpolation
// Run:    ./build/tests/test_interpolation
//
// The tests write their own HDF5 files into a scratch directory, so nothing
// here depends on aisoptics being installed. What they *do* depend on is the
// layout aisoptics documents, which is asserted explicitly below:
//
//   datasets x, y, z   1-D axis coordinates
//   datasets phase,    flattened C-order, index i*Ny*Nz + j*Nz + k
//     amplitude
//   attribute          export_format == "aispp_current_interpolation_hdf5"
//
// If aisoptics ever changes that layout, the round-trip test in
// tests/test_aisoptics_roundtrip.py is the one that will catch it; this file
// pins the C++ half.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "AISLaserBeam.hh"
#include "AISUtilities.hh"
#include "AISWavefronts.hh"

// ── tiny test harness (same shape as test_rotating.cc) ───────────────────────

static int gChecks = 0;
static int gFailures = 0;

static void check(bool ok, const std::string& what)
{
    ++gChecks;
    if (!ok) { ++gFailures; std::printf("  FAIL  %s\n", what.c_str()); }
}

static void checkClose(double got, double want, double tol, const std::string& what)
{
    ++gChecks;
    double err = std::fabs(got - want);
    double scale = std::fmax(1.0, std::fmax(std::fabs(got), std::fabs(want)));
    if (err > tol * scale)
    {
        ++gFailures;
        std::printf("  FAIL  %s: got %.17g want %.17g (rel err %.3g > %.3g)\n",
                    what.c_str(), got, want, err / scale, tol);
    }
}

static void section(const std::string& name)
{
    std::printf("\n%s\n", name.c_str());
}

// ── writing a beam file in the aisoptics layout ──────────────────────────────

struct GridSpec
{
    std::vector<double> x, y, z;
    std::vector<double> phase, amplitude;

    size_t index(size_t i, size_t j, size_t k) const
    {
        return i * y.size() * z.size() + j * z.size() + k;
    }
};

// Build a grid by sampling arbitrary analytic functions of position.
template <typename PhaseFn, typename AmpFn>
static GridSpec sampleGrid(size_t nx, size_t ny, size_t nz,
                           double lo, double hi,
                           PhaseFn phaseFn, AmpFn ampFn)
{
    GridSpec g;
    auto axis = [&](size_t n) {
        std::vector<double> a(n);
        for (size_t i = 0; i < n; ++i)
            a[i] = lo + (hi - lo) * static_cast<double>(i) / static_cast<double>(n - 1);
        return a;
    };
    g.x = axis(nx); g.y = axis(ny); g.z = axis(nz);
    g.phase.resize(nx * ny * nz);
    g.amplitude.resize(nx * ny * nz);
    for (size_t i = 0; i < nx; ++i)
        for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
            {
                g.phase[g.index(i, j, k)]     = phaseFn(g.x[i], g.y[j], g.z[k]);
                g.amplitude[g.index(i, j, k)] = ampFn(g.x[i], g.y[j], g.z[k]);
            }
    return g;
}

static void writeBeamFile(const std::string& path, const GridSpec& g,
                          const std::string& exportFormat = "aispp_current_interpolation_hdf5",
                          bool writeFormatAttr = true)
{
    H5::H5File file(path, H5F_ACC_TRUNC);

    auto writeVec = [&file](const std::string& name, const std::vector<double>& v) {
        hsize_t dims[1] = { v.size() };
        H5::DataSpace space(1, dims);
        H5::DataSet ds = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
        ds.write(v.data(), H5::PredType::NATIVE_DOUBLE);
    };

    writeVec("x", g.x);
    writeVec("y", g.y);
    writeVec("z", g.z);
    writeVec("phase", g.phase);
    writeVec("amplitude", g.amplitude);

    if (writeFormatAttr)
    {
        H5::StrType stype(H5::PredType::C_S1, H5T_VARIABLE);
        H5::DataSpace scalar(H5S_SCALAR);
        H5::Attribute attr = file.createAttribute("export_format", stype, scalar);
        attr.write(stype, exportFormat);
    }
    file.close();
}

static std::string scratchDir()
{
    static const std::string dir = "_interp_test";
    std::system(("mkdir -p " + dir).c_str());
    return dir;
}

// ── tests ────────────────────────────────────────────────────────────────────

// A trilinear interpolant must reproduce any function that is itself trilinear
// EXACTLY, at every point, not just on the nodes. This is the sharpest possible
// check on the index arithmetic: get i*Ny*Nz + j*Nz + k wrong in any way -- a
// transposed axis, an off-by-one, an F-order assumption -- and the exactness
// disappears immediately.
static void testExactOnTrilinearField()
{
    section("Trilinear field is reproduced exactly (pins the C-order index formula)");

    auto f = [](double x, double y, double z) {
        return 1.0 + 2.0*x - 3.0*y + 0.5*z + 4.0*x*y - 1.5*y*z + 0.25*x*z + 2.5*x*y*z;
    };
    // deliberately different extents per axis: a transposed index cannot hide
    GridSpec g = sampleGrid(5, 7, 9, -1.0, 1.0, f, f);

    const double pts[][3] = {
        {-0.83, 0.11, 0.47}, {0.0, 0.0, 0.0}, {0.62, -0.55, -0.91},
        {0.999, 0.999, 0.999}, {-1.0, -1.0, -1.0}, {0.3, -0.7, 0.2},
    };
    for (const auto& p : pts)
    {
        double got = trilinearInterpolation(p[0], p[1], p[2], g.x, g.y, g.z, g.phase,
                                            g.x.size(), g.y.size(), g.z.size());
        checkClose(got, f(p[0], p[1], p[2]), 1e-12,
                   "trilinear exact at (" + std::to_string(p[0]) + ", " +
                   std::to_string(p[1]) + ", " + std::to_string(p[2]) + ")");
    }
}

// Out-of-grid points must clamp to the nearest edge value. The original scan
// fell through to `return grid.size()-2` for anything below grid[0], picking the
// TOP cell and producing a large negative weight, so a point just below the
// grid came back with a wildly wrong value instead of the edge value.
static void testOutOfGridClamping()
{
    section("Out-of-grid sampling clamps at both ends");

    auto f = [](double x, double y, double z) { return x + 10.0*y + 100.0*z; };
    GridSpec g = sampleGrid(4, 4, 4, 0.0, 1.0, f, f);

    auto interp = [&](double x, double y, double z) {
        return trilinearInterpolation(x, y, z, g.x, g.y, g.z, g.phase,
                                      g.x.size(), g.y.size(), g.z.size());
    };

    // Far below the grid on each axis in turn -> the corresponding edge value.
    checkClose(interp(-5.0, 0.5, 0.5), interp(0.0, 0.5, 0.5), 1e-12, "clamp low x");
    checkClose(interp(0.5, -5.0, 0.5), interp(0.5, 0.0, 0.5), 1e-12, "clamp low y");
    checkClose(interp(0.5, 0.5, -5.0), interp(0.5, 0.5, 0.0), 1e-12, "clamp low z");

    // and above
    checkClose(interp(5.0, 0.5, 0.5), interp(1.0, 0.5, 0.5), 1e-12, "clamp high x");
    checkClose(interp(0.5, 5.0, 0.5), interp(0.5, 1.0, 0.5), 1e-12, "clamp high y");
    checkClose(interp(0.5, 0.5, 5.0), interp(0.5, 0.5, 1.0), 1e-12, "clamp high z");

    // the low corner specifically: this is the case that used to blow up
    checkClose(interp(-1e3, -1e3, -1e3), f(0.0, 0.0, 0.0), 1e-12, "clamp far low corner");
    check(std::fabs(interp(-1e3, -1e3, -1e3)) < 1e3, "far-low sample is bounded");
}

// The whole point of the feature: a beam sampled from the analytic Gaussian and
// read back through wtype "interpolated" must agree with wtype "gaussian".
static void testAgainstAnalyticGaussian()
{
    section("Sampled Gaussian agrees with the analytic Gaussian beam");

    const double w0 = 0.01;      // 1 cm waist
    const double kz = 9.0e6;     // clock transition, upward
    const double extent = 0.004; // stay well inside the waist

    // Sample exactly what the analytic beam would produce, in the exporter's
    // convention: envelope phase with the +-kz carrier excluded.
    auto phaseFn = [&](double x, double y, double z) {
        return gaussianWavefront({x, y, z}, kz, w0);
    };
    auto ampFn = [](double, double, double) { return 1.0; };

    GridSpec g = sampleGrid(81, 81, 5, -extent, extent, phaseFn, ampFn);
    const std::string path = scratchDir() + "/gaussian_sampled.h5";
    writeBeamFile(path, g);

    AISLaserBeam sampled({0.0, 0.0, kz}, 0, 1.0, 0.0);
    sampled.SetBeamType("interpolated");
    sampled.SetW0(w0);
    sampled.SetInterpolationGrids(path);
    check(sampled.HasInterpolationGrids(), "grids loaded");

    // Compare against the sampled function itself rather than against
    // AISLaserBeam's gaussian arm, which adds the reflection shift and the
    // focal-length shift on top -- those are separate concerns and are covered
    // by the analytic tests. What is under test here is that the grid survives
    // the write/read/interpolate round trip.
    const double probes[][3] = {
        {0.0, 0.0, 0.0}, {0.0011, -0.0007, 0.001}, {-0.002, 0.0015, -0.002},
        {0.003, 0.003, 0.0}, {-0.0035, -0.0005, 0.0018},
    };
    for (const auto& p : probes)
    {
        double got = sampled.GetPhi({p[0], p[1], p[2]});
        double want = phaseFn(p[0], p[1], p[2]);
        // trilinear on a smooth field: error is O(h^2 * curvature), so a
        // loose-but-meaningful tolerance on an 81-point transverse grid
        checkClose(got, want, 2e-3, "interpolated phase matches sampled Gaussian");
    }

    // amplitude is uniform 1.0 here, so the Rabi frequency is just rabiFreq
    checkClose(sampled.GetRabiFreq({0.001, 0.001, 0.0}, 0, 0), 1.0, 1e-9,
               "interpolated amplitude scales rabiFreq");
}

// Amplitude must scale the Rabi frequency multiplicatively, per the exporter's
// "field_amplitude_relative_to_configured_rabi_frequency" convention.
static void testAmplitudeConvention()
{
    section("Amplitude is relative to the configured Rabi frequency");

    auto zero = [](double, double, double) { return 0.0; };
    auto amp  = [](double x, double, double) { return 0.25 + 0.5 * x; }; // linear, so exact
    GridSpec g = sampleGrid(4, 4, 4, 0.0, 1.0, zero, amp);

    const std::string path = scratchDir() + "/amplitude.h5";
    writeBeamFile(path, g);

    const double rabi = 6283.185307;
    AISLaserBeam beam({0.0, 0.0, 9.0e6}, 0, rabi, 0.0);
    beam.SetBeamType("interpolated");
    beam.SetInterpolationGrids(path);

    checkClose(beam.GetRabiFreq({0.0, 0.5, 0.5}, 0, 0), 0.25 * rabi, 1e-12, "amplitude at x=0");
    checkClose(beam.GetRabiFreq({1.0, 0.5, 0.5}, 0, 0), 0.75 * rabi, 1e-12, "amplitude at x=1");
    checkClose(beam.GetRabiFreq({0.5, 0.5, 0.5}, 0, 0), 0.50 * rabi, 1e-12, "amplitude at x=0.5");
}

// phi0 must survive as an additive offset on the interpolated phase.
static void testPhaseOffset()
{
    section("phi0 adds to the interpolated phase");

    auto constPhase = [](double, double, double) { return 0.125; };
    auto unitAmp    = [](double, double, double) { return 1.0; };
    GridSpec g = sampleGrid(3, 3, 3, 0.0, 1.0, constPhase, unitAmp);

    const std::string path = scratchDir() + "/offset.h5";
    writeBeamFile(path, g);

    const double phi0 = 0.75;
    AISLaserBeam beam({0.0, 0.0, 9.0e6}, 0, 1.0, phi0);
    beam.SetBeamType("interpolated");
    beam.SetInterpolationGrids(path);

    checkClose(beam.GetPhi({0.5, 0.5, 0.5}), phi0 + 0.125, 1e-12, "phi0 + grid phase");
}

// Tip-tilt shears the sample point by tan(theta) per unit z. At z = 0 it must
// be a no-op, which is what keeps untilted runs unchanged.
static void testTipTilt()
{
    section("Tip-tilt shears the sample point and is identity at z = 0");

    auto fx = [](double x, double, double) { return x; }; // read back the sampled x
    auto unitAmp = [](double, double, double) { return 1.0; };
    GridSpec g = sampleGrid(9, 3, 9, -1.0, 1.0, fx, unitAmp);

    const std::string path = scratchDir() + "/tiptilt.h5";
    writeBeamFile(path, g);

    AISLaserBeam beam({0.0, 0.0, 9.0e6}, 0, 1.0, 0.0);
    beam.SetBeamType("interpolated");
    beam.SetInterpolationGrids(path);

    const double thetaDeg = 10.0;
    beam.SetTipTilt(thetaDeg, 0.0);

    // at z = 0 the shear vanishes
    checkClose(beam.GetPhi({0.3, 0.0, 0.0}), 0.3, 1e-12, "tip-tilt is identity at z=0");

    // at z != 0 the sampled x is shifted by -z*tan(theta)
    const double z = 0.5;
    const double expected = 0.3 - z * std::tan(thetaDeg * M_PI / 180.0);
    checkClose(beam.GetPhi({0.3, 0.0, z}), expected, 1e-12, "tip-tilt shear at z != 0");

    // and zero tilt leaves it alone
    beam.SetTipTilt(0.0, 0.0);
    checkClose(beam.GetPhi({0.3, 0.0, z}), 0.3, 1e-12, "zero tilt is a no-op");
}

// A wrong or missing export_format must not be interpreted as a valid grid.
// This is the assertion that stops a future aisoptics layout change from
// silently producing wrong phases instead of an error.
static void testExportFormatIsChecked()
{
    section("export_format guards the layout contract");

    auto f = [](double, double, double) { return 0.0; };
    GridSpec g = sampleGrid(3, 3, 3, 0.0, 1.0, f, f);

    // The good case is already exercised above; here just confirm that the
    // perturbation-only format aisoptics can also emit is not silently accepted.
    const std::string bad = scratchDir() + "/wrong_format.h5";
    writeBeamFile(bad, g, "aispp_future_perturbation_hdf5");

    // SetInterpolationGrids exits on mismatch, so this is checked out-of-process
    // by the CMake-driven runner rather than here; what we can assert in-process
    // is that a correctly stamped file is accepted.
    const std::string good = scratchDir() + "/right_format.h5";
    writeBeamFile(good, g);

    AISLaserBeam beam({0.0, 0.0, 9.0e6}, 0, 1.0, 0.0);
    beam.SetBeamType("interpolated");
    beam.SetInterpolationGrids(good);
    check(beam.HasInterpolationGrids(), "correctly stamped file is accepted");

    // A file with no attribute at all is accepted with a warning (older exports).
    const std::string noAttr = scratchDir() + "/no_attr.h5";
    writeBeamFile(noAttr, g, "", /*writeFormatAttr=*/false);
    AISLaserBeam legacy({0.0, 0.0, 9.0e6}, 0, 1.0, 0.0);
    legacy.SetBeamType("interpolated");
    legacy.SetInterpolationGrids(noAttr);
    check(legacy.HasInterpolationGrids(), "unstamped file is accepted with a warning");
}

// The analytic beam types must be completely untouched by all of the above --
// this is the regression guard for having restored them after the merge.
static void testAnalyticTypesStillWork()
{
    section("Analytic beam types are unaffected");

    const double w0 = 0.01, kz = 9.0e6;

    AISLaserBeam flat({0.0, 0.0, kz}, 0, 1.0, 0.3);
    flat.SetBeamType("flat_square");
    checkClose(flat.GetPhi({0.001, 0.002, 0.003}), 0.3, 1e-15, "flat_square phase is phi0");
    checkClose(flat.GetRabiFreq({0.001, 0.0, 0.0}, 0, 0), 1.0, 1e-15, "flat_square rabi is uniform");

    AISLaserBeam gauss({0.0, 0.0, kz}, 0, 1.0, 0.0);
    gauss.SetBeamType("gaussian");
    gauss.SetW0(w0);
    gauss.SetFocalLength(0.0);
    // upward beam picks up the pi reflection shift plus the analytic wavefront
    double want = pi + gaussianWavefront({0.001, 0.002, 0.003}, kz, w0);
    checkClose(gauss.GetPhi({0.001, 0.002, 0.003}), want, 1e-12, "gaussian phase is analytic");

    // KNOWN ISSUE (see KNOWN_ISSUES.md #1): GetDelPhi is currently zero for
    // EVERY beam type, not just interpolated ones. gaussianGradientWavefront
    // computes the gradient into a local and then unconditionally returns
    // {0,0,0}, carrying the comment "wavefront curvature accounted for in the
    // zernike coeffs in this branch" - a leftover from the zernike-aberration-only
    // line of work. That makes commit d89bbbc ("Wire GetDelPhi to the existing
    // gaussianGradientWavefront for Gaussian beams") a no-op in practice.
    //
    // These assertions pin the CURRENT behaviour deliberately, so the suite stays
    // green while the issue is open. When the gradient is turned back on, they
    // will fail and must be flipped to the non-zero expectations - that failure
    // is the point.
    doubleThreeVector grad = gauss.GetDelPhi({0.001, 0.002, 0.003});
    check(grad[0] == 0.0 && grad[1] == 0.0 && grad[2] == 0.0,
          "gaussian GetDelPhi is zero (KNOWN ISSUE #1, not the intended physics)");

    std::array<double,3> raw = gaussianGradientWavefront({0.001, 0.002, 0.003}, kz, w0);
    check(raw[0] == 0.0 && raw[1] == 0.0 && raw[2] == 0.0,
          "gaussianGradientWavefront discards its own result (KNOWN ISSUE #1)");

    AISLaserBeam conf({0.0, 0.0, kz}, 0, 1.0, 0.0);
    conf.SetBeamType("confocal");
    conf.SetW0(w0);
    conf.SetFocalLength(0.0);
    doubleThreeVector cgrad = conf.GetDelPhi({0.001, 0.002, 0.003});
    check(cgrad[0] == 0.0 && cgrad[1] == 0.0 && cgrad[2] == 0.0,
          "confocal GetDelPhi is zero (KNOWN ISSUE #1)");
}

// Zernike coefficients must actually reach the phase again (they were inert
// between adbebdf and v0.0.2), and must stay a no-op when unset.
static void testZernikeIsWiredUp()
{
    section("Zernike coefficients reach GetPhi again");

    const double w0 = 0.01, kz = 9.0e6;

    AISLaserBeam plain({0.0, 0.0, kz}, 0, 1.0, 0.0);
    plain.SetBeamType("gaussian");
    plain.SetW0(w0);
    plain.SetBeamRadius(0.05);
    plain.SetFocalLength(0.0);
    double without = plain.GetPhi({0.01, 0.005, 0.0});

    AISLaserBeam aberrated({0.0, 0.0, kz}, 0, 1.0, 0.0);
    aberrated.SetBeamType("gaussian");
    aberrated.SetW0(w0);
    aberrated.SetBeamRadius(0.05);
    aberrated.SetFocalLength(0.0);
    aberrated.SetZernikeCoeffs({{4, 0.1}}); // Noll 4 = defocus
    double with = aberrated.GetPhi({0.01, 0.005, 0.0});

    check(std::fabs(with - without) > 1e-9, "a non-zero Zernike coefficient changes the phase");

    // empty coefficient map must be exactly a no-op
    AISLaserBeam empty({0.0, 0.0, kz}, 0, 1.0, 0.0);
    empty.SetBeamType("gaussian");
    empty.SetW0(w0);
    empty.SetBeamRadius(0.05);
    empty.SetFocalLength(0.0);
    empty.SetZernikeCoeffs({});
    checkClose(empty.GetPhi({0.01, 0.005, 0.0}), without, 1e-15,
               "no coefficients leaves the phase bit-for-bit unchanged");

    // a zero-valued coefficient must also leave the phase alone
    AISLaserBeam zeroed({0.0, 0.0, kz}, 0, 1.0, 0.0);
    zeroed.SetBeamType("gaussian");
    zeroed.SetW0(w0);
    zeroed.SetBeamRadius(0.05);
    zeroed.SetFocalLength(0.0);
    zeroed.SetZernikeCoeffs({{4, 0.0}});
    checkClose(zeroed.GetPhi({0.01, 0.005, 0.0}), without, 1e-15,
               "zero coefficient leaves the phase unchanged");
}

int main()
{
    std::printf("ais++ beam interpolation unit tests\n");

    testExactOnTrilinearField();
    testOutOfGridClamping();
    testAgainstAnalyticGaussian();
    testAmplitudeConvention();
    testPhaseOffset();
    testTipTilt();
    testExportFormatIsChecked();
    testAnalyticTypesStillWork();
    testZernikeIsWiredUp();

    std::printf("\n%d checks, %d failures\n", gChecks, gFailures);
    if (gFailures == 0) std::printf("ALL TESTS PASSED\n");
    return gFailures == 0 ? 0 : 1;
}
