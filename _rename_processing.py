"""One-off script to rename processing_OOP to processing. Run from repo root."""
import os
import shutil
import sys

root = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(root, "processing_OOP")
dst = os.path.join(root, "processing")

if not os.path.isdir(src):
    print("processing_OOP not found", file=sys.stderr)
    sys.exit(1)

if os.path.exists(dst):
    print("Removing existing processing/ ...", file=sys.stderr)
    shutil.rmtree(dst)

print("Copying processing_OOP -> processing ...", flush=True)
shutil.copytree(src, dst)
print("Removing processing_OOP ...", flush=True)
shutil.rmtree(src)
print("Done: processing_OOP renamed to processing", flush=True)
