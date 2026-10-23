#!/usr/bin/env python3
"""Compatibility loader for chains written by the official no-qPCR runner.

The manuscript runner is executed with ``runpy(..., run_name="__main__")``.
Consequently, two runner-local classes are recorded in ``mcmc.pkl`` as
``__main__.FilteringWithoutQPCR`` and
``__main__.SubjectLogTrajectorySetWithoutQPCR``.  A later Python process must
restore those names before unpickling the trusted, locally generated chain.

This module centralises that narrowly scoped compatibility shim.  It does not
make arbitrary pickle files safe; callers must only load trusted MDSINE2
artifacts produced by this workflow.
"""

from __future__ import annotations

from contextlib import contextmanager
import importlib.util
from pathlib import Path
import sys
from types import ModuleType
from typing import Any, Callable, Iterator


DEFAULT_OFFICIAL_RUNNER = Path(
    "/home/qqt/MDSINE2_Server_Package/MDSINE2_Paper/scripts/"
    "semisynthetic2/inference/mdsine2/helpers/mdsine2_ra.py"
)
SHIM_CLASS_NAMES = (
    "FilteringWithoutQPCR",
    "SubjectLogTrajectorySetWithoutQPCR",
)
_RUNNER_MODULE_NAME = "_fig3_mdsine2_official_noqpcr_runner"


def _load_runner_module(runner: Path) -> ModuleType:
    runner = Path(runner).resolve()
    if not runner.is_file():
        raise FileNotFoundError(f"Official no-qPCR runner not found: {runner}")

    cached = sys.modules.get(_RUNNER_MODULE_NAME)
    if cached is not None and Path(getattr(cached, "__file__", "")).resolve() == runner:
        return cached

    spec = importlib.util.spec_from_file_location(_RUNNER_MODULE_NAME, runner)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import official no-qPCR runner: {runner}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[_RUNNER_MODULE_NAME] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(_RUNNER_MODULE_NAME, None)
        raise
    return module


@contextmanager
def noqpcr_main_shim(
    runner: Path = DEFAULT_OFFICIAL_RUNNER,
) -> Iterator[dict[str, Any]]:
    """Temporarily restore the historical ``__main__`` pickle class names."""

    module = _load_runner_module(Path(runner))
    main_module = sys.modules.get("__main__")
    if main_module is None:  # pragma: no cover - Python always supplies this
        raise RuntimeError("sys.modules['__main__'] is unavailable")

    sentinel = object()
    previous = {name: getattr(main_module, name, sentinel) for name in SHIM_CLASS_NAMES}
    installed: dict[str, str] = {}
    try:
        for name in SHIM_CLASS_NAMES:
            if not hasattr(module, name):
                raise AttributeError(
                    f"Official no-qPCR runner lacks required class {name!r}"
                )
            cls = getattr(module, name)
            setattr(main_module, name, cls)
            installed[name] = f"{cls.__module__}.{cls.__qualname__}"
        yield {
            "runner": str(Path(runner).resolve()),
            "runner_module": module.__name__,
            "installed_into": "__main__",
            "classes": installed,
        }
    finally:
        for name, value in previous.items():
            if value is sentinel:
                if hasattr(main_module, name):
                    delattr(main_module, name)
            else:
                setattr(main_module, name, value)


def load_noqpcr_mcmc(
    path: Path,
    runner: Path = DEFAULT_OFFICIAL_RUNNER,
    loader: Callable[[str], Any] | None = None,
) -> Any:
    """Load a trusted no-qPCR ``mcmc.pkl`` after installing its class shim.

    ``loader`` is injectable for focused tests.  Production callers omit it,
    in which case ``mdsine2.BaseMCMC.load`` is used.
    """

    source = Path(path).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    with noqpcr_main_shim(Path(runner)):
        if loader is None:
            try:
                import mdsine2 as md2
            except ImportError as exc:  # pragma: no cover - environment-dependent
                raise RuntimeError(
                    "MDSINE2 is required to load mcmc.pkl; use the mdsine2 conda "
                    "environment and a writable NUMBA_CACHE_DIR."
                ) from exc
            loader = md2.BaseMCMC.load
        return loader(str(source))
