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
    with open(path) as fh:
        data = yaml.safe_load(fh) or {}
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
    return _wrap(cfg)


def parse_override(text):
    """Turn ``'instrument.baseline=10'`` into ``{'instrument': {'baseline': 10}}``."""
    key, _, value = text.partition("=")
    if not _:
        raise ConfigError(f"override '{text}' must look like key.sub=value")
    result = yaml.safe_load(value)
    for part in reversed(key.strip().split(".")):
        result = {part: result}
    return result
