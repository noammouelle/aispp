import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from detfig import ConfigError, load, standalone, write  # noqa: E402
from detfig.config import parse_override  # noqa: E402
from detfig.labels import spread  # noqa: E402

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
    assert cfg.node.arm == 1.0                          # inherited from generic.yaml


def test_typo_is_reported_with_suggestion(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("instrument:\n  shiled: {show: false}\n")
    with pytest.raises(ConfigError, match="did you mean 'shield'"):
        load(bad)


def test_override_from_command_line():
    cfg = load(ROOT / "configs" / "generic.yaml", overrides=parse_override("instrument.baseline=10"))
    assert cfg.instrument.baseline == 10


def test_unknown_anchor_lists_known_ones(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("extends: %s\noverview:\n  labels:\n    - {text: x, anchor: nosuch}\n"
                        % (ROOT / "configs" / "generic.yaml"))
    with pytest.raises(KeyError, match="telescope"):
        standalone(load(cfg_file))


def test_unquoted_comma_in_flow_label_is_caught(tmp_path):
    cfg_file = tmp_path / "c.yaml"
    cfg_file.write_text("extends: %s\noverview:\n  labels:\n    - {text: a\\,b, anchor: beam}\n"
                        % (ROOT / "configs" / "generic.yaml"))
    with pytest.raises(ValueError, match="quote it"):
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
