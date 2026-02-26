"""
One-off script: scan all job result directories under post_processing/results,
read their .log files, and report any lines containing errors.
Run from repo root: python workflow_orchestrator/review_job_logs.py
"""
import os
import re
import sys

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.abspath(os.path.join(script_dir, ".."))
results_base = os.path.join(repo_root, "post_processing", "results")

# Patterns that indicate something went wrong
ERROR_PATTERNS = [
    re.compile(r"ERROR", re.I),
    re.compile(r"CRITICAL", re.I),
    re.compile(r"Traceback", re.I),
    re.compile(r"Exception:", re.I),
    re.compile(r"Error:", re.I),
    re.compile(r"\bFailed\b", re.I),
    re.compile(r"failed\b", re.I),
    re.compile(r"Invalid\b", re.I),
    re.compile(r"ValueError|RuntimeError|FileNotFoundError|OSError", re.I),
]

def main():
    if not os.path.isdir(results_base):
        print(f"No results directory at {results_base}")
        return 0

    dirs = [d for d in os.listdir(results_base) if os.path.isdir(os.path.join(results_base, d))]
    # Job result dirs have a 'logs' subdir and look like job_XXXX_nN_...
    job_result_dirs = []
    for d in dirs:
        logs_path = os.path.join(results_base, d, "logs")
        if os.path.isdir(logs_path):
            job_result_dirs.append(d)

    if not job_result_dirs:
        print(f"No job result directories (with logs/) found under {results_base}")
        return 0

    print(f"Found {len(job_result_dirs)} job result director(y/ies). Scanning .log files for errors.\n")
    any_errors = False
    for job_dir in sorted(job_result_dirs):
        logs_path = os.path.join(results_base, job_dir, "logs")
        log_files = [f for f in os.listdir(logs_path) if f.endswith(".log")]
        for log_file in sorted(log_files):
            path = os.path.join(logs_path, log_file)
            try:
                with open(path, "r", encoding="utf-8", errors="replace") as f:
                    for i, line in enumerate(f, 1):
                        for pat in ERROR_PATTERNS:
                            if pat.search(line):
                                if not any_errors:
                                    any_errors = True
                                print(f"{job_dir}/logs/{log_file}:{i}: {line.rstrip()}")
                                break
            except Exception as e:
                print(f"{job_dir}/logs/{log_file}: (read error) {e}")
                any_errors = True

    if not any_errors:
        print("No error-like lines found in any of the scanned log files.")
    return 1 if any_errors else 0

if __name__ == "__main__":
    sys.exit(main())
