"""Copy remaining files from processing_OOP to processing (only missing files)."""
import os
import shutil

root = os.path.dirname(os.path.abspath(__file__))
src = os.path.join(root, "processing_OOP")
dst = os.path.join(root, "processing")
n = 0
for r, dirs, files in os.walk(src):
    for f in files:
        rel = os.path.relpath(os.path.join(r, f), src)
        p1 = os.path.join(src, rel)
        p2 = os.path.join(dst, rel)
        if not os.path.exists(p2) or os.path.getmtime(p1) > os.path.getmtime(p2):
            os.makedirs(os.path.dirname(p2), exist_ok=True)
            shutil.copy2(p1, p2)
            n += 1
            print(rel)
print("Copied", n, "files")
