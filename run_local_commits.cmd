@echo off
REM Modular local commits (conventional commits). Run from repo root.
REM Run "git status" first and skip any step that has nothing to commit.
cd /d "%~dp0"

echo === 1. docs(processing): fix file path comments ===
git add processing/static/results/compute_secondary/energy.py processing/static/results/compute_secondary/strain.py processing/static/results/compute_secondary/stress.py processing/static/results/map_reader.py
git diff --cached --quiet || git commit -m "docs(processing): fix file path comments in four modules"

echo === 2. chore: remove refactor helper scripts ===
git add REFACTOR_OPTION_A_FINISH.md _copy_remaining.py _rename_processing.py run_option_a.bat run_refactor_option_a_finish.cmd 2>nul
git diff --cached --quiet || git commit -m "chore: remove refactor helper scripts and finish note"

echo === 3. docs(plans): implementation plans ===
git add docs/plans/ 2>nul
git diff --cached --quiet || git commit -m "docs(plans): add implementation plans for B2, formulation cache, section forces"

echo === 4. refactor(jobs): mesh-variant job layout ===
git add jobs/ 2>nul
git diff --cached --quiet || git commit -m "refactor(jobs): switch to mesh-variant job layout and add distributed mesh variants"

echo === 5. feat(pre_processing): mesh variants and mesh_generator ===
git add pre_processing/mesh_library/ pre_processing/mesh_library/schemes/mesh_generator.py 2>nul
git diff --cached --quiet || git commit -m "feat(pre_processing): add distributed mesh variants and update mesh_generator"

echo === 6. feat(post_processing): verification visualisers ===
git add post_processing/ 2>nul
git diff --cached --quiet || git commit -m "feat(post_processing): add verification visualisers and results hierarchy docs"

echo === 7. refactor: rename processing_OOP to processing (only if still uncommitted) ===
git add processing/ 2>nul
git add -u processing_OOP/ 2>nul
git diff --cached --quiet || git commit -m "refactor: rename package processing_OOP to processing"

echo.
echo === Done. Remaining ===
git status --short
