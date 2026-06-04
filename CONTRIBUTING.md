# Contributing to Valyte

Thank you for your interest in contributing to Valyte! This document provides guidelines and instructions for contributing.

## How to Contribute

### Reporting Bugs

If you find a bug, please [open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new) with:

- A clear, descriptive title
- Steps to reproduce the problem
- Expected vs. actual behavior
- Your environment — Python version, Valyte version (`valyte --version`), and OS
- Relevant VASP output files (or minimal reproducing examples) if applicable

### Requesting Features

Feature requests are welcome. [Open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new) with:

- A clear description of the feature and why it would be useful
- Example usage or expected output, if possible
- Any relevant references (papers, existing tools, etc.)

### Submitting Pull Requests

1. **Fork** the repository and create a new branch from `main`:

   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Install** in development mode:

   ```bash
   git clone https://github.com/<your-username>/Valyte-Project
   cd Valyte-Project
   pip install -e .
   ```

3. **Make your changes.** Follow the existing code style:
   - Use clear, descriptive variable and function names
   - Add docstrings for new functions
   - Keep commits focused — one logical change per commit

4. **Test** your changes with real VASP output files to ensure correctness.

5. **Push** your branch and [open a pull request](https://github.com/nikyadav002/Valyte-Project/pulls):
   - Reference any related issues (e.g., "Fixes #12")
   - Describe what changed and why
   - Include sample plots or output if the change is visual

## Development Setup

```bash
# Clone your fork
git clone https://github.com/<your-username>/Valyte-Project
cd Valyte-Project

# Install in editable mode
pip install -e .

# Verify the CLI works
valyte --help
```

### Dependencies

Valyte depends on:

| Package | Purpose |
|---|---|
| `numpy` | Numerical operations |
| `matplotlib` | Plotting |
| `pymatgen` | Crystal structure parsing, POTCAR generation |
| `scipy` | Curve fitting (effective mass) |
| `click` | CLI framework |
| `seekpath` | K-path generation |

## Code Style

- **Python**: Follow PEP 8 conventions
- **CLI**: Use [Click](https://click.palletsprojects.com/) decorators consistently
- **Plots**: Match the existing Valyte aesthetic — clean typography, gradient fills, publication-quality defaults

## Documentation

If your change adds a new command, flag, or modifies existing behavior:

1. Update the relevant page in `docs/`
2. Update `docs/cli-reference.md` with the new flags
3. Update `README.md` feature tables if a new command is added

Documentation is built with [MkDocs Material](https://squidfunk.github.io/mkdocs-material/). To preview locally:

```bash
pip install mkdocs-material
mkdocs serve
```

## Questions?

If you're unsure about anything, feel free to [open an issue](https://github.com/nikyadav002/Valyte-Project/issues/new) to discuss before starting work. We're happy to help!
