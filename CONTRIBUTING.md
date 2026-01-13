# Contributing to bsaseq

Thank you for your interest in contributing to bsaseq! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Code Style](#code-style)
- [Testing](#testing)
- [Documentation](#documentation)
- [Pull Request Process](#pull-request-process)
- [Reporting Issues](#reporting-issues)
- [Feature Requests](#feature-requests)

## Code of Conduct

This project follows a standard code of conduct. Please be respectful and constructive in all interactions.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Set up the development environment
4. Create a feature branch
5. Make your changes
6. Submit a pull request

## Development Setup

### Prerequisites

- Python 3.9-3.12
- Git
- (Optional) Docker for container testing

### Installation

```bash
# Clone your fork
git clone https://github.com/rcac-bioinformatics/bsaseq.git
cd bsaseq

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/macOS
# or: venv\Scripts\activate  # Windows

# Install in development mode with all dependencies
pip install -e ".[dev]"

# Verify installation
bsaseq --version
pytest tests/ -v
```

### IDE Setup

For VS Code, recommended extensions:
- Python
- Pylance
- Ruff

For PyCharm:
- Enable ruff integration in Settings > Tools > Ruff

## Making Changes

### Branch Naming

Use descriptive branch names:
- `feature/add-new-statistic` - New features
- `fix/vcf-parsing-error` - Bug fixes
- `docs/update-readme` - Documentation
- `refactor/improve-window-calculation` - Code refactoring

### Commit Messages

Write clear, concise commit messages:

```
Add tricube smoothing option for window analysis

- Implement tricube weight function
- Add --smoothing CLI option
- Update tests for new functionality
```

## Code Style

We use [ruff](https://github.com/astral-sh/ruff) for linting and formatting.

### Running the Linter

```bash
# Check for issues
ruff check src/ tests/

# Auto-fix issues
ruff check --fix src/ tests/

# Format code
ruff format src/ tests/
```

### Style Guidelines

- Follow PEP 8
- Maximum line length: 100 characters
- Use type hints for function signatures
- Write docstrings for public functions (Google style)

Example:

```python
def calculate_delta_af(
    af_high: float,
    af_low: float,
) -> float:
    """Calculate delta allele frequency between bulks.

    Args:
        af_high: Allele frequency in high bulk.
        af_low: Allele frequency in low bulk.

    Returns:
        Difference in allele frequency (high - low).

    Raises:
        ValueError: If frequencies are outside [0, 1].
    """
    if not (0 <= af_high <= 1 and 0 <= af_low <= 1):
        raise ValueError("Allele frequencies must be between 0 and 1")
    return af_high - af_low
```

## Testing

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ -v --cov=bsaseq --cov-report=term-missing

# Run specific test file
pytest tests/test_windows.py -v

# Run specific test
pytest tests/test_windows.py::TestWindowCalculation::test_window_calculation_basic -v
```

### Writing Tests

- Place tests in `tests/` directory
- Name test files `test_*.py`
- Name test functions `test_*`
- Use fixtures from `conftest.py`
- Aim for >80% code coverage

Example:

```python
import pytest
from bsaseq.analysis.windows import calculate_windows_to_list


class TestWindowCalculation:
    """Tests for window calculation functions."""

    def test_window_calculation_basic(self, sample_variants):
        """Test basic window calculation."""
        windows = calculate_windows_to_list(
            sample_variants,
            window_size=1000000,
            step_size=250000,
        )

        assert len(windows) > 0
        assert all(w.n_variants > 0 for w in windows)

    def test_empty_variants_returns_empty(self):
        """Test that empty input returns empty output."""
        windows = calculate_windows_to_list(
            [],
            window_size=1000000,
            step_size=250000,
        )

        assert windows == []
```

### Test Fixtures

Common fixtures are defined in `tests/conftest.py`:
- `sample_vcf_path` - Temporary VCF file
- `sample_variants` - List of Variant objects
- `sample_windows` - List of Window objects
- `sample_regions` - List of CandidateRegion objects

## Documentation

### Docstrings

All public functions, classes, and modules should have docstrings following Google style.

### README Updates

Update README.md when:
- Adding new features
- Changing CLI options
- Modifying output formats

### Changelog

Add entries to CHANGELOG.md under `[Unreleased]` for:
- New features (Added)
- Bug fixes (Fixed)
- Breaking changes (Changed)
- Deprecations (Deprecated)
- Removals (Removed)

## Pull Request Process

1. **Update your branch**
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

2. **Run checks locally**
   ```bash
   ruff check src/ tests/
   pytest tests/ -v
   ```

3. **Create pull request**
   - Use a descriptive title
   - Reference related issues
   - Describe changes and motivation
   - Include test results

4. **Address review feedback**
   - Respond to all comments
   - Push additional commits for changes
   - Request re-review when ready

5. **Merge requirements**
   - All CI checks passing
   - At least one approving review
   - No merge conflicts

## Reporting Issues

### Before Reporting

1. Check existing issues
2. Verify with latest version: `pip install --upgrade bsaseq`
3. Reproduce with minimal example

### Issue Template

```markdown
**Description**
Brief description of the issue.

**Steps to Reproduce**
1. Run command: `bsaseq run --vcf test.vcf ...`
2. Observe error

**Expected Behavior**
What should happen.

**Actual Behavior**
What actually happens.

**Environment**
- bsaseq version: `bsaseq --version`
- Python version: `python --version`
- Operating system: Linux/macOS/Windows
- Installation method: pip/conda/docker

**Error Output**
```
Paste full error traceback here
```

**Sample Data**
If possible, provide minimal VCF that reproduces the issue.
```

## Feature Requests

### Before Requesting

1. Check existing issues and PRs
2. Consider if it fits the project scope
3. Think about implementation approach

### Feature Request Template

```markdown
**Feature Description**
What feature would you like?

**Use Case**
Why is this feature needed? What problem does it solve?

**Proposed Implementation**
How might this be implemented? (Optional)

**Alternatives Considered**
What alternatives have you considered?
```

## Questions?

- Open a GitHub Discussion
- Check existing documentation
- Review closed issues

Thank you for contributing to bsaseq!
