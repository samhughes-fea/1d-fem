# Finish Option A: processing_OOP → processing

Run these in order from the repo root (`c:\Users\s1834431\Programs\fem_model`).

## 1. Remove empty/partial processing (if present)

```cmd
rmdir /s /q processing
```

## 2. Rename the package directory

```cmd
git mv processing_OOP processing
```

## 3. Verify

```cmd
python -c "from processing.static.operations.solver import SolveCondensedSystem; print('ok')"
python -m pytest tests/test_formulation_cache_shape_functions.py tests/test_nodal_section_forces_projector.py -v
```

## 4. Remove helper scripts and commit

```cmd
del _rename_processing.py _copy_remaining.py run_option_a.bat REFACTOR_OPTION_A_FINISH.md
git add -A
git commit -m "refactor: rename package processing_OOP to processing"
```
