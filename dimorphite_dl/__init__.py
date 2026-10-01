"""Adds hydrogen atoms to molecular representations as specified by pH"""

import contextlib
import os
import sys

from loguru import logger

from .protonate.run import protonate_smiles

__all__ = ["protonate_smiles"]

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"

logger.disable("dimorphite_dl")

LOG_FORMAT = (
    "<green>{time:HH:mm:ss}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
)

# Loguru ids of the sinks enable_logging added, so a later call can replace
# them without touching sinks that belong to the host application.
_handler_ids: list[int] = []


def enable_logging(
    level_set: int,
    stdout_set: bool = True,
    file_path: str | None = None,
    log_format: str = LOG_FORMAT,
    colorize: bool = True,
) -> None:
    r"""Enable logging.

    Calling again replaces the sinks added by the previous call. Sinks added
    by anything else, including loguru's default stderr sink, are left in
    place; an application that wants only these sinks should call
    `logger.remove()` first.

    Args:
        level: Requested log level: `10` is debug, `20` is info.
        stdout_set: Write logs to the console. They go to stderr, despite the
            name, so that stdout carries only protonated SMILES.
        file_path: Also write logs to files here.
        colorize: Color the console output. File output is never colored.
    """
    # Replace only this package's sinks. logger.configure(handlers=...) would
    # also remove every sink the host application added.
    for handler_id in _handler_ids:
        # The host may already have removed it, e.g. with logger.remove().
        with contextlib.suppress(ValueError):
            logger.remove(handler_id)
    _handler_ids.clear()
    if stdout_set:
        _handler_ids.append(
            logger.add(
                sys.stderr, level=level_set, format=log_format, colorize=colorize
            )
        )
    if isinstance(file_path, str):
        _handler_ids.append(
            logger.add(
                file_path,
                level=level_set,
                format=log_format,
                # Color codes would be written into the file verbatim.
                colorize=False,
            )
        )

    logger.enable("dimorphite_dl")


def _env_flag(name: str, default: bool) -> bool:
    """Parse a boolean environment variable leniently, because values like
    "true" or "yes" made literal_eval raise and broke the import.

    Args:
        name: Environment variable name.
        default: Value when the variable is unset.

    Returns:
        True for 1/true/yes/on (any case), False for anything else.
    """
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _env_log_level(name: str, default: int) -> int:
    """Parse a log level given as a number or a name such as "INFO", because
    int() on a name raised and broke the import.

    Args:
        name: Environment variable name.
        default: Level when the variable is unset or unrecognized.

    Returns:
        The numeric log level.
    """
    value = os.environ.get(name)
    if value is None:
        return default
    value = value.strip()
    if value.isdigit():
        return int(value)
    try:
        return logger.level(value.upper()).no
    except ValueError:
        return default


if _env_flag("DIMORPHITE_DL_LOG", False):
    level = _env_log_level("DIMORPHITE_DL_LOG_LEVEL", 20)
    stdout = _env_flag("DIMORPHITE_DL_STDOUT", True)
    log_file_path = os.environ.get("DIMORPHITE_DL_LOG_FILE_PATH", None)
    enable_logging(level, stdout, log_file_path)
