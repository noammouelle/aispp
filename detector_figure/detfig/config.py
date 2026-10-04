"""YAML configuration: loading, inheritance and typo checking.

A config file may start with ``extends: other.yaml`` (path relative to the
file).  Files are deep-merged on top of ``configs/defaults.yaml``:
dictionaries merge key by key, everything else (numbers, strings, *lists*)
is replaced wholesale.  So to change one colour you write only that colour,
but to change the list of nodes you write the whole list.

Any key that does not exist in ``defaults.yaml`` raises a ``ConfigError``
naming the closest valid key.  This catches typos such as ``shiled`` that
would otherwise be ignored silently.  Keys under a dictionary whose default
is empty (``{}``) are free-form and not checked.
"""

import difflib
from pathlib import Path

import yaml

DEFAULTS = Path(__file__).resolve().parent.parent / "configs" / "defaults.yaml"


class ConfigError(ValueError):
    pass


class Cfg(dict):
    """A dict with attribute access, recursively: ``cfg.style.colors.beam``."""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"config has no key '{key}' (have: {', '.join(self)})") from None

    def __setattr__(self, key, value):
        self[key] = value


def _wrap(obj):
    if isinstance(obj, dict):
        return Cfg({k: _wrap(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [_wrap(v) for v in obj]
    return obj


def deep_merge(base, override):
    """Return ``base`` updated recursively with ``override`` (neither is modified)."""
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def _read(path, seen=()):
    path = Path(path).resolve()
    if path in seen:
        raise ConfigError(f"circular 'extends' involving {path}")
    try:
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
    except yaml.YAMLError as err:
        raise ConfigError(f"{path}: not valid YAML (text with commas, braces or a leading "
                          f"special character must be quoted):\n{err}") from None
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: top level must be a mapping")
    parent = data.pop("extends", None)
    if parent:
        data = deep_merge(_read(path.parent / parent, seen + (path,)), data)
    return data


def check_keys(user, reference, where="config"):
    """Raise ConfigError for keys in ``user`` that ``reference`` does not have."""
    for key, value in user.items():
        if key not in reference:
            hint = difflib.get_close_matches(str(key), [str(k) for k in reference], n=1)
            msg = f"unknown key '{where}.{key}'"
            raise ConfigError(msg + (f" (did you mean '{hint[0]}'?)" if hint else ""))
        ref = reference[key]
        if isinstance(value, dict) and isinstance(ref, dict) and ref:
            check_keys(value, ref, f"{where}.{key}")


def load(*paths, overrides=None):
    """Load defaults, then each YAML file in turn, then a dict of overrides."""
    with open(DEFAULTS) as fh:
        defaults = yaml.safe_load(fh)
    cfg = defaults
    for p in paths:
        user = _read(p)
        check_keys(user, defaults)
        cfg = deep_merge(cfg, user)
    if overrides:
        check_keys(overrides, defaults)
        cfg = deep_merge(cfg, overrides)
    cfg = _wrap(cfg)
    validate(cfg)
    return cfg


# Allowed keys of the items of each list in the config (lists replace
# wholesale, so they cannot be checked against defaults.yaml).
LABEL_KEYS = {"text", "anchor", "side", "dx", "dy", "shift", "line"}
MARKER_KEYS = {"n", "anchor", "dx", "dy", "leader"}
LIST_ITEM_KEYS = {
    "instrument.nodes": {"id", "z", "source", "camera"},
    "instrument.atoms.trajectories": {"id", "start", "apex", "end", "dx"},
    "overview.labels": LABEL_KEYS,
    "overview.markers": MARKER_KEYS,
    "node.labels": LABEL_KEYS,
    "node.markers": MARKER_KEYS,
    "trajectory.clouds": {"z0", "v0", "z_end", "label"},
    "section.layers": {"kind", "radius", "size", "thickness", "sides", "color", "label", "side",
                       "angle", "layers", "spacing", "count", "wire", "pair_spread"},
}

# Keys whose value must be one of a fixed set.
CHOICES = {
    "instrument.laser.side": ("left", "right"),
    "instrument.laser.route": ("direct", "periscope"),
    "instrument.laser.tip_tilt_on": ("axis", "transfer"),
    "instrument.coils.style": ("bars", "turns"),
    "overview.scale_bar.style": ("dimension", "axis"),
    "node.source_side": (None, "left", "right"),
    "node.below": ("tube", "stub"),
    "node.launch_lattice": ("crossed", "standing"),
}


def _get(cfg, dotted):
    for part in dotted.split("."):
        cfg = cfg[part]
    return cfg


def _choice(value, allowed, where):
    if value not in allowed:
        hint = difflib.get_close_matches(str(value), [str(a) for a in allowed if a], n=1)
        raise ConfigError(f"{where}: {value!r} is not one of {', '.join(map(str, allowed))}"
                          + (f" (did you mean '{hint[0]}'?)" if hint else ""))


def validate(cfg):
    """Check what defaults.yaml cannot: list items, enumerations, consistency."""
    for dotted, allowed in LIST_ITEM_KEYS.items():
        items = _get(cfg, dotted) or []
        if not isinstance(items, list):
            raise ConfigError(f"{dotted} must be a list")
        for i, item in enumerate(items):
            if not isinstance(item, dict):
                raise ConfigError(f"{dotted}[{i}] must be a mapping like {{key: value, ...}}")
            text = str(item.get("text", item.get("label", "")))
            if text.endswith("\\") and not text.endswith("\\\\"):
                raise ConfigError(f"{dotted}[{i}]: text {text!r} ends in a backslash - in YAML flow "
                                  "style ({text: ..., ...}) a comma ends the text, so quote it: "
                                  "text: 'a\\,b'")
            for key in item:
                if key not in allowed:
                    hint = difflib.get_close_matches(str(key), sorted(allowed), n=1)
                    raise ConfigError(f"unknown key '{key}' in {dotted}[{i}]"
                                      + (f" (did you mean '{hint[0]}'?)" if hint else "")
                                      + f"; allowed: {', '.join(sorted(allowed))}")
    for dotted, allowed in CHOICES.items():
        _choice(_get(cfg, dotted), allowed, dotted)

    ins = cfg.instrument
    if not ins.nodes:
        raise ConfigError("instrument.nodes must contain at least one node")
    ids = [n.get("id") for n in ins.nodes]
    if None in ids or len(set(ids)) != len(ids):
        raise ConfigError(f"instrument.nodes: every node needs a unique id (got {ids})")
    if not ins.baseline or ins.baseline <= 0:
        raise ConfigError("instrument.baseline must be a positive length in metres")
    for n in ins.nodes:
        _choice(n.get("source", "none"), ("left", "right", "none"), f"instrument.nodes[{n['id']}].source")
        if not -1e-9 <= n["z"] <= ins.baseline * (1 + 1e-9):
            raise ConfigError(f"node '{n['id']}' at z = {n['z']} m lies outside 0..baseline "
                              f"({ins.baseline} m); heights are metres above the bottom")
    for stage in cfg.node.stages:
        _choice(stage, ("oven", "zeeman", "mot2d"), "node.stages")
    for side in cfg.node.cameras or []:
        _choice(side, ("left", "right"), "node.cameras")
    for name, value in cfg.style.colors.items():
        if value is None:
            raise ConfigError(f"style.colors.{name} is empty - in YAML an unquoted '#' starts a "
                              "comment, so write colours as \"B2182B\" or '#B2182B'")
    for lay in cfg.section.layers:
        if str(lay.get("color", "")).startswith("#"):
            raise ConfigError(f"section layer colour {lay['color']!r}: define it under style.colors "
                              "and refer to it by name")
        if lay.get("kind") == "coils" and lay.get("count", 16) < 1:
            raise ConfigError("section coils layer needs count >= 1")
    if cfg.overview.scale_bar.step < 0:
        raise ConfigError("overview.scale_bar.step must be >= 0 (0 = no ticks)")
    panels = cfg.layout.panels
    if len(set(panels)) != len(panels):
        raise ConfigError(f"layout.panels lists a panel twice: {panels}")
    if "section" in panels and not cfg.section.layers:
        raise ConfigError("section.layers is empty; add layers or remove 'section' from layout.panels")


def parse_override(text):
    """Turn ``'instrument.baseline=10'`` into ``{'instrument': {'baseline': 10}}``."""
    key, _, value = text.partition("=")
    if not _:
        raise ConfigError(f"override '{text}' must look like key.sub=value")
    result = yaml.safe_load(value)
    for part in reversed(key.strip().split(".")):
        result = {part: result}
    return result
