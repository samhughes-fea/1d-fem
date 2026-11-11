from __future__ import annotations

from pathlib import Path


def build_runner_result(
    *,
    family: str,
    summary_file: str | Path,
    reference_dir: str | Path,
    expected_files: list[str],
    present_files: list[str],
    ready: bool,
    passed: bool,
    status: str,
    tolerances: dict[str, float | str],
    metrics: dict[str, float | int | str | bool],
    artifact_paths: dict[str, str | Path],
) -> dict[str, str | bool | float | int]:
    out: dict[str, str | bool | float | int] = {
        "family": family,
        "ready": bool(ready),
        "passed": bool(passed),
        "status": str(status),
        "summary_file": str(summary_file),
        "reference_dir": str(reference_dir),
        "expected_files": ",".join(expected_files),
        "present_files": ",".join(present_files),
        "tolerances": ",".join(f"{key}={value}" for key, value in tolerances.items()),
    }
    for key, value in metrics.items():
        out[key] = value
    for key, value in artifact_paths.items():
        out[key] = str(value)
    return out


def write_runner_summary(
    summary_file: str | Path,
    *,
    family: str,
    ready: bool,
    passed: bool,
    status: str,
    reference_dir: str | Path,
    expected_files: list[str],
    present_files: list[str],
    tolerances: dict[str, float | str],
    metrics: dict[str, float | int | str | bool],
    artifact_paths: dict[str, str | Path],
) -> None:
    lines = [
        f"family={family}",
        f"status={status}",
        f"ready={ready}",
        f"passed={passed}",
        f"reference_dir={reference_dir}",
        f"expected_files={','.join(expected_files)}",
        f"present_files={','.join(present_files)}",
        f"tolerance_keys={','.join(tolerances.keys())}",
    ]
    for key, value in tolerances.items():
        lines.append(f"tolerance_{key}={value}")
    for key, value in artifact_paths.items():
        lines.append(f"{key}={value}")
    for key, value in metrics.items():
        lines.append(f"{key}={value}")
    Path(summary_file).write_text("\n".join(lines) + "\n", encoding="utf-8")
