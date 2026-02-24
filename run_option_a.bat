@echo off
REM Option A: Remove partial processing/, then git mv processing_OOP processing
cd /d "%~dp0"
if exist processing rmdir /s /q processing
git mv processing_OOP processing
echo Done. Run: python -m pytest tests/test_formulation_cache_shape_functions.py tests/test_nodal_section_forces_projector.py -v
echo Then remove _rename_processing.py and _copy_remaining.py, and commit.
