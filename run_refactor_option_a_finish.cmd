@echo off
REM Finish Option A: processing_OOP -> processing (run from repo root)
cd /d "%~dp0"

if exist processing rmdir /s /q processing
git mv processing_OOP processing
python -c "from processing.static.operations.solver import SolveCondensedSystem; print('ok')"
python -m pytest tests/test_formulation_cache_shape_functions.py tests/test_nodal_section_forces_projector.py -v

del _rename_processing.py _copy_remaining.py run_option_a.bat REFACTOR_OPTION_A_FINISH.md run_refactor_option_a_finish.cmd 2>nul
git add -A
git commit -m "refactor: rename package processing_OOP to processing"
echo Done.
