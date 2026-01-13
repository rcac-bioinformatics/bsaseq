"""Logging configuration for bsaseq.

This module provides a consistent logging setup using rich for
formatted console output with timestamps and progress bars.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path
    from typing import TypeVar

    T = TypeVar("T")

# Global console for stderr output
console = Console(stderr=True)

# Track if logging has been set up
_logging_configured = False


def setup_logging(
    level: int = logging.INFO,
    show_time: bool = True,
    show_path: bool = False,
) -> None:
    """Configure logging with rich handler.

    This function sets up the root logger to use rich formatting
    for console output. Call this once at application startup.

    Args:
        level: Logging level (default: INFO).
        show_time: Whether to show timestamps (default: True).
        show_path: Whether to show file paths in log messages.
    """
    global _logging_configured

    if _logging_configured:
        return

    # Configure rich handler
    handler = RichHandler(
        console=console,
        show_time=show_time,
        show_path=show_path,
        rich_tracebacks=True,
        tracebacks_show_locals=False,
        markup=True,
    )
    handler.setFormatter(logging.Formatter("%(message)s"))

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.addHandler(handler)

    # Configure bsaseq logger
    bsaseq_logger = logging.getLogger("bsaseq")
    bsaseq_logger.setLevel(level)

    _logging_configured = True


def get_logger(name: str) -> logging.Logger:
    """Get a logger with the given name.

    The logger will be a child of the bsaseq logger, ensuring
    consistent formatting and configuration.

    Args:
        name: Logger name (typically __name__ from the calling module).

    Returns:
        Configured logger instance.

    Example:
        >>> logger = get_logger(__name__)
        >>> logger.info("Processing started")
    """
    # Ensure logging is configured
    if not _logging_configured:
        setup_logging()

    # If name starts with 'bsaseq.', use it directly
    # Otherwise, prefix with 'bsaseq.'
    if name.startswith("bsaseq."):
        return logging.getLogger(name)
    return logging.getLogger(f"bsaseq.{name}")


def create_progress() -> Progress:
    """Create a rich progress bar for long operations.

    Returns a Progress context manager configured with appropriate
    columns for tracking file processing operations.

    Returns:
        Configured Progress instance.

    Example:
        >>> with create_progress() as progress:
        ...     task = progress.add_task("Processing", total=100)
        ...     for i in range(100):
        ...         progress.advance(task)
    """
    return Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
        transient=True,
    )


def progress_iterator(
    iterator: Iterator[T],
    total: int | None = None,
    description: str = "Processing",
) -> Iterator[T]:
    """Wrap an iterator with a progress bar.

    Args:
        iterator: Iterator to wrap.
        total: Total number of items (if known).
        description: Description to show in progress bar.

    Yields:
        Items from the iterator.

    Example:
        >>> for item in progress_iterator(items, total=len(items)):
        ...     process(item)
    """
    with create_progress() as progress:
        task = progress.add_task(description, total=total)
        for item in iterator:
            yield item
            progress.advance(task)


def print_info(message: str) -> None:
    """Print an info message to stderr.

    Args:
        message: Message to print.
    """
    console.print(f"[blue]INFO:[/blue] {message}")


def print_warning(message: str) -> None:
    """Print a warning message to stderr.

    Args:
        message: Message to print.
    """
    console.print(f"[yellow]WARNING:[/yellow] {message}")


def print_error(message: str) -> None:
    """Print an error message to stderr.

    Args:
        message: Message to print.
    """
    console.print(f"[red]ERROR:[/red] {message}")


def print_success(message: str) -> None:
    """Print a success message to stderr.

    Args:
        message: Message to print.
    """
    console.print(f"[green]SUCCESS:[/green] {message}")


def print_stats(stats: dict[str, int | float | str], title: str = "Statistics") -> None:
    """Print statistics in a formatted table.

    Args:
        stats: Dictionary of statistic names and values.
        title: Title for the statistics block.
    """
    from rich.table import Table

    table = Table(title=title, show_header=False)
    table.add_column("Metric", style="cyan")
    table.add_column("Value", justify="right")

    for key, value in stats.items():
        if isinstance(value, float):
            table.add_row(key, f"{value:,.2f}")
        elif isinstance(value, int):
            table.add_row(key, f"{value:,}")
        else:
            table.add_row(key, str(value))

    console.print(table)


# Convenience aliases for clearer semantic naming
log_info = print_info
log_warning = print_warning
log_error = print_error
log_success = print_success


def log_step(step: int, total: int, description: str) -> None:
    """Log a step in a multi-step process.

    Displays step progress in a consistent format, useful for showing
    progress through major analysis phases.

    Args:
        step: Current step number (1-indexed).
        total: Total number of steps.
        description: Description of the current step.

    Example:
        >>> log_step(1, 5, "Loading VCF file")
        [1/5] Loading VCF file
    """
    console.print(f"[bold cyan][{step}/{total}][/bold cyan] {description}")


def track_progress(
    iterable: Iterator[T],
    total: int | None = None,
    description: str = "Processing",
) -> Iterator[T]:
    """Wrap an iterable with a progress bar.

    This is an alias for progress_iterator with a clearer name,
    matching common usage patterns.

    Args:
        iterable: Iterable to wrap.
        total: Total number of items (if known).
        description: Description to show in progress bar.

    Yields:
        Items from the iterable.

    Example:
        >>> for item in track_progress(items, total=len(items), description="Filtering"):
        ...     process(item)
    """
    yield from progress_iterator(iterable, total=total, description=description)


def print_section(title: str) -> None:
    """Print a section header.

    Args:
        title: Section title to display.
    """
    from rich.rule import Rule

    console.print(Rule(title, style="blue"))


def print_file_created(path: str | Path) -> None:
    """Print a message indicating a file was created.

    Args:
        path: Path to the created file.
    """
    from pathlib import Path as PathlibPath

    path = PathlibPath(path)
    console.print(f"  [dim]Created:[/dim] {path.name}")
