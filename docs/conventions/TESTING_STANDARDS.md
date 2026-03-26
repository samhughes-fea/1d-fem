# Testing standards

Pytest standards for this repository for writing and organizing automated tests with `pytest`.

---

## 1. Goals

- Make tests **readable** and self-documenting.
- Make failures **easy to diagnose**.
- Keep tests **fast**, **reliable**, and **deterministic**.
- Enable **automation** in CI without special per-developer knowledge.

---

## 2. Test layout and discovery

### 2.1 Directory structure

- Tests live either under a **single top-level `tests/`** directory, or under **per-package `tests/`** directories (see §2.3).
- Mirror the application's package structure where practical.

**Single-package example:**

```text
src/
  myproject/
    core/
      loads.py
    api/
      routes.py

tests/
  core/
    test_loads.py
  api/
    test_routes.py
```

### 2.2 pytest configuration

- **Repo root:** Add a root `pytest.ini` for project-wide settings (e.g. shared markers such as `slow`, `integration`, and common `addopts`). This is the single place for markers used across multiple packages.
- **Discovery:** Either set `testpaths = tests` in the root (single top-level `tests/`) or run pytest with explicit paths per package (see §10).
- **Per-package config (multi-package repos):** Each package may have its own `pytest.ini` with `testpaths = tests`, `python_files`, `python_classes`, `python_functions`, and package-specific `addopts` or markers.

**Minimal root `pytest.ini` (shared markers):**

```ini
[pytest]
markers =
    slow: Slow-running tests (e.g. optimization, > 1s)
    integration: Integration tests requiring full pipeline
```

**Full discovery from root (single top-level tests/):**

```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
```

**Per-package `pytest.ini` (inside e.g. `blade_structure_fatigue/`):**

```ini
[pytest]
testpaths = tests
python_files = test_*.py
python_classes = Test*
python_functions = test_*
addopts = -v --strict-markers --tb=short --disable-warnings
markers =
    unit: Unit tests
    integration: Integration tests
    slow: Slow running tests
```

All tests must be discoverable without special per-developer knowledge (either `pytest` from root with suitable `testpaths`, or `pytest <package>/tests` from root, or `pytest` from within each package directory).

Avoid custom discovery rules unless strictly necessary.

### 2.3 Multi-package (monorepo) layout

When a monorepo contains **multiple installable packages** (not this repo) (e.g. `blade_structure_fatigue`, `blade_structure_precompute`, `blade_structure_utilities`), use **per-package** test directories:

- Each package has its own `tests/` directory: `<package>/tests/`.
- Under `<package>/tests/`, mirror that package's module layout where practical (e.g. `tests/section_buckling/`, `tests/stress/unit/`).
- Shared fixtures for that package live in `<package>/tests/conftest.py`; subdirectories may have their own `conftest.py` for scoped fixtures.
- Nested `pytest.ini` files may define markers or options for that subtree (e.g. `tests/section_selection/pytest.ini`); root and package-level `pytest.ini` remain the main config.

**Example (blade-structure style):**

```text
blade_structure_fatigue/
  tests/
    conftest.py
    internal_force_history/
      test_internal_force_history.py
    stress_history/
      test_stress_history.py

blade_structure_utilities/
  tests/
    section_buckling/
      conftest.py
      test_plate_buckling.py
    stress/
      unit/
        test_analytical_validation.py
```

Ensure the repo root is on `PYTHONPATH` or packages are installed in editable mode (`pip install -e .`) so imports like `from blade_structure_utilities.stress import ...` resolve when running tests from any directory.

---

## 3. Naming conventions

### 3.1 Test files

Test files must be named: `test_*.py`

Examples: `test_loads.py`, `test_routes.py`, `test_integration_api.py`.

**Suffix (Option A):** The part after `test_` should be the **module or subpackage** under test (snake_case), matching the directory that mirrors the package. When one module has multiple test files, use `test_<module>_<topic>.py`.

**Allowed exceptions:**

- `test_<package>_import.py` – import smoke tests.
- `test_<module>_integration.py` – dedicated integration tests for that module.
- `test_suite_<N>_<topic>.py` – numbered suite batches (e.g. direct_stiffness_sampling, stiffness_selection).
- `test_<area>_debug.py` – debug/exploratory tests for that area.

**Avoid:** Generic suffixes like `test_utilities.py` or `test_helpers.py` for files that are the main tests for a module; prefer `test_<module>.py` or `test_<module>_utilities.py`. Shared test infrastructure (helper classes/functions with no test cases) should not be named `test_*.py` so pytest does not collect them; use e.g. `<module>_helpers.py` or `_test_helpers.py`.

### 3.2 Test classes

Use classes only to group related tests, not to share state.

Class names:

- Start with `Test`.
- Use PascalCase.

Examples:

```python
class TestComputeLoads:
    ...

class TestBladeApi:
    ...
```

### 3.3 Test functions

All test functions must:

- Start with `test_`.
- Be written in snake_case.
- Follow the pattern:

  **`test_<unitUnderTest>_<conditionOrInput>_<expectedOutcome>`**

Examples:

- `test_compute_tip_zero_bill_returns_zero`
- `test_validate_blade_length_negative_value_raises_valueerror`
- `test_api_get_loads_unauthenticated_returns_401`

**Guideline:** The name alone should make the test's purpose clear without reading the body.

---

## 4. Test structure (AAA / Given–When–Then)

Each test must follow a three-part structure:

1. **Arrange (Given)** – set up data, dependencies, and environment.
2. **Act (When)** – perform exactly one action under test.
3. **Assert (Then)** – verify behavior at the end of the test.

**Example (AAA):**

```python
def test_compute_tip_zero_bill_returns_zero():
    # Arrange
    bill = 0.0
    rate = 0.15

    # Act
    result = compute_tip(bill, rate)

    # Assert
    assert result == 0.0
```

**Example (Given–When–Then comments):**

```python
def test_validate_blade_length_negative_value_raises_valueerror():
    # Given a negative blade length
    blade_length = -5.0

    # When validating the blade
    with pytest.raises(ValueError):
        validate_blade_length(blade_length)

    # Then a ValueError is raised
```

**Rules:**

- Prefer one logical behavior per test.
- Avoid conditionals and loops inside tests; split into multiple tests or use parametrization.
- Assertions should be at the end of the test.

---

## 5. Fixtures

### 5.1 General rules

- Use pytest fixtures for reusable setup instead of class-level setup or global state.
- Place shared fixtures in `tests/conftest.py` (in multi-package repos, each package has its own `tests/conftest.py`).
- Name fixtures descriptively: `load_model`, `db_session`, `api_client`.

**Example:**

```python
# tests/conftest.py
import pytest
from myproject.core.loads import LoadModel

@pytest.fixture
def load_model():
    return LoadModel()
```

**Usage:**

```python
# tests/core/test_loads.py
def test_compute_loads_zero_wind_returns_zero(load_model):
    # Arrange
    wind_speed = 0.0

    # Act
    load = load_model.compute(wind_speed)

    # Assert
    assert load == 0.0
```

### 5.2 Fixture guidelines

- Fixtures should be explicit dependencies: prefer function parameters over global imports.
- Avoid hidden side effects (e.g., mutating shared global state).
- Use fixture scopes (`function`, `module`, `session`) deliberately; default to `function`.

---

## 6. Parametrization

Use parametrization for systematic input variation instead of loops in tests.

**Example:**

```python
import pytest

@pytest.mark.parametrize(
    "wind_speed, expected_load",
    [
        (0.0, 0.0),
        (5.0, 100.0),
        (10.0, 400.0),
    ],
)
def test_compute_loads_various_speeds(load_model, wind_speed, expected_load):
    # Act
    load = load_model.compute(wind_speed)

    # Assert
    assert load == expected_load
```

**Rules:**

- Keep parameter sets small and meaningful.
- Name parameters clearly; avoid generic names like `x`, `y` unless trivial.

---

## 7. Test types and placement

| Type | Description | Placement |
|------|-------------|-----------|
| **Unit tests** | Test a single function or method in isolation. Use lightweight fixtures and no real external services. | Next to the relevant module under `tests/` subdirectories. |
| **Integration tests** | Test interactions between components (e.g., DB + service, API + service). | Can be grouped into folders like `tests/integration/`. |
| **End-to-end / system tests** | Exercise the full application stack where needed. Kept minimal due to cost. | May live in `tests/e2e/`. |

---

## 8. BDD and higher-level behaviors (optional)

For user-facing or cross-functional behaviors, BDD-style tests are optional.

If adopted:

- Use `.feature` files with Gherkin syntax under `tests/features/`.
- Keep BDD scenarios for high-level behavior, not low-level unit details.
- Step definitions live in Python files under `tests/steps/`.

**Example feature:**

```gherkin
Feature: Blade load validation

  Scenario: Reject loads above ultimate factor
    Given a blade design with ultimate factor 1.5
    And a load model predicting factor 1.6
    When I validate the blade loads
    Then the design is marked as unsafe
```

---

## 9. Style and quality guidelines

- Tests must be **deterministic** (no dependence on time, randomness, or external state without control).
- Avoid real external services; use:
  - Fakes, mocks, or stubs.
  - In-memory databases where suitable.
- Keep tests fast; slow tests should be:
  - Marked with `@pytest.mark.slow`.
  - Excluded from default runs if necessary.

**Example:**

```python
import pytest

@pytest.mark.slow
def test_full_pipeline_with_real_files(...):
    ...
```

---

## 10. Running tests

- **Local development:** Run `pytest` at minimum before pushing.
- **Single top-level `tests/`:** From repo root, `pytest` (with root `testpaths = tests`).
- **Multi-package (per-package `tests/`):**
  - From repo root, run one package: `pytest blade_structure_fatigue/tests` or `pytest blade_structure_utilities/tests`.
  - Run all packages: `pytest blade_structure_fatigue/tests blade_structure_utilities/tests blade_structure_precompute/tests blade_structure_vibration/tests` (adjust to actual package list).
  - From inside a package: `cd blade_structure_fatigue && pytest` (uses that package's `pytest.ini` and `testpaths = tests`).
- **Imports:** Ensure repo root is on `PYTHONPATH` or install packages in editable mode so test imports resolve from any working directory.
- **CI:** Must run `pytest` on every PR (with appropriate paths in multi-package repos). Consider `pytest --maxfail=1 --disable-warnings` for faster feedback.


---

## 11. Code review checklist (for tests)

Reviewers should ensure:

- [ ] File, class, and function names follow this standard.
- [ ] Each test follows AAA / Given–When–Then.
- [ ] Tests are readable without needing to inspect implementation.
- [ ] Fixtures are used appropriately; no unnecessary global state.
- [ ] Parametrization is used instead of loops where appropriate.
- [ ] Slow or integration tests are clearly marked.

---

*End of testing standards*
