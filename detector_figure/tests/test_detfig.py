import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from detfig import ConfigError, load, standalone, write  # noqa: E402
from detfig.config import parse_override  # noqa: E402
from detfig.labels import spread  # noqa: E402
from detfig.panels.trajectory import clouds_from_instrument, landing_time  # noqa: E402

CONFIGS = sorted((ROOT / "configs").glob("*.yaml"))
FIGURES = [p for p in CONFIGS if p.name != "defaults.yaml"]


@pytest.mark.parametrize("path", FIGURES, ids=lambda p: p.stem)
def test_every_config_generates_tex(path):
    tex = standalone(load(path))
    assert tex.count("\\begin{tikzpicture}") == 1
    assert "\\end{document}" in tex


def test_extends_merges_dicts_and_replaces_lists():
    cfg = load(ROOT / "configs" / "magis100.yaml")
    assert [n["id"] for n in cfg.instrument.nodes] == ["bottom", "middle", "top"]
    assert cfg.instrument.shield.show is True          # inherited from defaults
    assert cfg.node.source.width == 1.2                 # inherited from generic.yaml
    assert cfg.node.arm == 0.8                          # overridden in magis100.yaml


def test_typo_is_reported_with_suggestion(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("instrument:\n  shiled: {show: false}\n")
    with pytest.raises(ConfigError, match="did you mean 'shield'"):
        load(bad)


def test_override_from_command_line():
    cfg = load(ROOT / "configs" / "generic.yaml", overrides=parse_override("overview.height=12"))
    assert cfg.overview.height == 12


def test_unknown_anchor_lists_known_ones(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("extends: %s\noverview:\n  labels:\n    - {text: x, anchor: nosuch}\n"
                        % (ROOT / "configs" / "generic.yaml"))
    with pytest.raises(ValueError, match="telescope"):
        standalone(load(cfg_file))


def test_unquoted_comma_in_flow_label_is_caught(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("extends: %s\noverview:\n  labels:\n    - {text: a\\,b, anchor: beam}\n"
                        % (ROOT / "configs" / "generic.yaml"))
    with pytest.raises(ConfigError, match="quote it"):
        load(cfg_file)


@pytest.mark.parametrize("yaml_text, message", [
    ("overview:\n  labels:\n    - {text: x, anchor: beam, sid: left}\n", "did you mean 'side'"),
    ("instrument:\n  laser: {side: rigth}\n", "did you mean 'right'"),
    ("instrument:\n  nodes:\n    - {id: a, z: 0, source: Left}\n", "not one of"),
    ("instrument:\n  baseline: 10\n", "outside 0..baseline"),
    ("overview:\n  scale_bar: {step: -5}\n", "step"),
    ("layout:\n  panels: [overview, overview]\n", "twice"),
    ("instrument:\n  nodes:\n    - {id: a, z: 0}\n    - {id: a, z: 100}\n", "unique id"),
    ("style:\n  colors:\n    beam: #B2182B\n", "unquoted '#'"),
])
def test_validation_errors(tmp_path, yaml_text, message):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text(f"extends: {ROOT / 'configs' / 'generic.yaml'}\n" + yaml_text)
    with pytest.raises(ConfigError, match=message):
        load(cfg_file)


def test_trajectory_clouds_match_panel_a():
    cfg = load(ROOT / "configs" / "aion10.yaml")
    clouds = clouds_from_instrument(cfg)
    assert len(clouds) == 2
    name, z0, v0, z_end = clouds[0]
    apex = z0 + v0 ** 2 / (2 * cfg.trajectory.g)
    assert apex == pytest.approx(cfg.instrument.atoms.trajectories[0]["apex"])


def test_landing_time_of_a_drop():
    # dropped from 100 m, detected at 50 m: t = sqrt(2 * 50 / g)
    assert landing_time(100, 0, 50, 9.81) == pytest.approx((100 / 9.81) ** 0.5)


def test_baseline_break_and_source_stages_render(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text(f"extends: {ROOT / 'configs' / 'generic.yaml'}\n"
                        "overview: {baseline_break: {show: true, z: 30, gap: 1.0}}\n"
                        "node: {stages: [oven, zeeman, mot2d]}\n")
    tex = standalone(load(cfg_file))
    assert tex.count("fill=white, draw=none") >= 1      # the break mask


def test_unknown_source_stage_is_rejected(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text(f"extends: {ROOT / 'configs' / 'generic.yaml'}\nnode: {{stages: [oven, zeman]}}\n")
    with pytest.raises(ConfigError, match="did you mean 'zeeman'"):
        load(cfg_file)


def _mzplots_available():
    import importlib.util
    return importlib.util.find_spec("lmt_sequence_diagram") is not None


@pytest.mark.skipif(not _mzplots_available(), reason="mz-plots not installed")
def test_mzplots_engine_embeds_lmt_diagram():
    cfg = load(ROOT / "configs" / "generic.yaml")
    assert cfg.trajectory.engine == "mzplots"
    tex = standalone(cfg)
    assert "\\usepackage{braket}" in tex                 # mz-plots' packages added
    assert "xscale=" in tex and "port" in tex              # its picture, scaled in


def test_mzplots_engine_falls_back_when_missing(tmp_path, monkeypatch):
    import detfig.panels.trajectory as trajectory
    monkeypatch.setattr(trajectory, "load_mzplots", lambda cfg: None)
    cfg = load(ROOT / "configs" / "generic.yaml")
    tex = standalone(cfg)
    assert "braket" not in tex and "\\pi/2" in tex       # builtin panel instead
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text(f"extends: {ROOT / 'configs' / 'generic.yaml'}\n"
                        "trajectory: {mzplots: {fallback: false}}\n")
    with pytest.raises(ValueError, match="mz-plots|lmt_sequence_diagram"):
        standalone(load(cfg_file))


def test_spread_keeps_order_and_gap():
    ys = spread([0.0, 0.05, 0.1, 3.0], [1, 2, 1, 1], gap=0.4)
    assert ys[0] < ys[1] < ys[2] < ys[3]
    assert ys[1] - ys[0] >= 0.4 * 1.5 - 1e-9
    assert ys[3] == pytest.approx(3.0)


def test_section_rejects_unknown_kind(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("section:\n  layers:\n    - {kind: hexagon, radius: 1}\n")
    with pytest.raises(ValueError, match="unknown kind"):
        standalone(load(cfg_file))


@pytest.mark.skipif(not shutil.which("pdflatex"), reason="pdflatex not installed")
@pytest.mark.parametrize("path", FIGURES, ids=lambda p: p.stem)
def test_every_config_compiles(path, tmp_path):
    cfg = load(path, overrides={"output": {"png_dpi": 0}})
    files = write(cfg, tmp_path)
    assert files["pdf"].stat().st_size > 10_000
