.PHONY: help install dev test lint format type-check docs docker docker-snpeff clean publish-test publish

PYTHON := python
PIP := pip
PYTEST := pytest
RUFF := ruff

# Default target
help:
	@echo "bsaseq development commands"
	@echo ""
	@echo "Setup:"
	@echo "  install       Install package"
	@echo "  dev           Install with dev dependencies"
	@echo ""
	@echo "Quality:"
	@echo "  test          Run tests"
	@echo "  test-cov      Run tests with coverage"
	@echo "  lint          Run linter"
	@echo "  format        Format code"
	@echo "  type-check    Run type checker"
	@echo "  check         Run all checks (lint + type-check + test)"
	@echo ""
	@echo "Build:"
	@echo "  build         Build package"
	@echo "  docker        Build Docker image"
	@echo "  docker-snpeff Build Docker image with snpEff"
	@echo ""
	@echo "Publish:"
	@echo "  publish-test  Publish to TestPyPI"
	@echo "  publish       Publish to PyPI"
	@echo ""
	@echo "Maintenance:"
	@echo "  clean         Remove build artifacts"
	@echo "  clean-all     Remove all generated files"

# Installation
install:
	$(PIP) install -e .

dev:
	$(PIP) install -e ".[dev]"

# Testing
test:
	$(PYTEST) tests/ -v

test-cov:
	$(PYTEST) tests/ -v --cov=bsaseq --cov-report=term-missing --cov-report=html

test-fast:
	$(PYTEST) tests/ -v -x --tb=short

# Code quality
lint:
	$(RUFF) check src/ tests/

lint-fix:
	$(RUFF) check --fix src/ tests/

format:
	$(RUFF) format src/ tests/
	$(RUFF) check --fix src/ tests/

type-check:
	mypy src/bsaseq --ignore-missing-imports

check: lint type-check test

# Build
build:
	$(PYTHON) -m build

build-check:
	$(PYTHON) -m build
	twine check dist/*

# Docker
docker:
	docker build -t bsaseq:latest .

docker-snpeff:
	docker build --build-arg INSTALL_SNPEFF=true -t bsaseq:snpeff .

docker-test:
	docker run --rm bsaseq:latest --help
	docker run --rm bsaseq:latest --version

# Publishing
publish-test: clean build
	twine upload --repository testpypi dist/*

publish: clean build
	twine upload dist/*

# Cleanup
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf src/*.egg-info/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true

clean-test:
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf coverage.xml
	rm -rf htmlcov/

clean-lint:
	rm -rf .ruff_cache/
	rm -rf .mypy_cache/

clean-all: clean clean-test clean-lint
	rm -rf .venv/
	rm -rf venv/

# Development utilities
version:
	@$(PYTHON) -c "import bsaseq; print(bsaseq.__version__)"

deps:
	$(PIP) list --outdated

update-deps:
	$(PIP) install --upgrade pip
	$(PIP) install -e ".[dev]" --upgrade

# Documentation (placeholder for future sphinx docs)
docs:
	@echo "Documentation build not yet configured"
	@echo "See README.md for current documentation"

# Shortcuts
t: test
l: lint
f: format
c: check
