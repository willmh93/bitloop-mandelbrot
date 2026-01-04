from pathlib import Path
import os
import shutil
import subprocess
import sys

def find_marker_upwards(start: Path | str | None = None, marker: str = ".vcpkg") -> Path | None:
    p = Path(start or Path.cwd()).resolve()
    for parent in (p, *p.parents):
        cand = parent / marker
        if cand.is_dir():
            return cand
    return None

def resolve_vcpkg_exe(marker_dir: Path | None) -> str:
    import os, shutil
    from pathlib import Path

    exe_name = "vcpkg.exe" if os.name == "nt" else "vcpkg"

    # 0) Explicit override
    p = os.environ.get("VCPKG_EXE")
    if p and Path(p).exists():
        return str(Path(p))

    # 1) VCPKG_ROOT
    root = os.environ.get("VCPKG_ROOT")
    if root:
        cand = Path(root) / exe_name
        if cand.exists():
            return str(cand)

    # 2) **INSIDE** the found .vcpkg directory (your case)
    if marker_dir:
        cand = marker_dir / exe_name
        if cand.exists():
            return str(cand)
        # a couple of common subpaths just in case
        for sub in ("downloads", "bin"):
            cand2 = marker_dir / sub / exe_name
            if cand2.exists():
                return str(cand2)

    # 3) A checked-out sibling 'vcpkg' folder up the tree
    start = Path.cwd().resolve()
    for parent in (start, *start.parents):
        cand = parent / "vcpkg" / exe_name
        if cand.exists():
            return str(cand)

    # 4) PATH fallback
    which = shutil.which(exe_name)
    if which:
        return which

    raise FileNotFoundError("Could not find 'vcpkg' (tried VCPKG_EXE, VCPKG_ROOT, .vcpkg/, vcpkg/ up-tree, PATH).")

def update_baseline(cwd: Path | str | None = None) -> int:
    cwd = Path(cwd or Path.cwd()).resolve()
    marker = find_marker_upwards(cwd, ".vcpkg")
    vcpkg = resolve_vcpkg_exe(marker)
    # Run in the *current project dir* so it updates that vcpkg.json's baseline
    print(f"Running: {vcpkg} x-update-baseline (cwd={cwd})")
    proc = subprocess.run([vcpkg, "x-update-baseline"], cwd=cwd)
    return proc.returncode

if __name__ == "__main__":
    code = update_baseline()
    if code != 0:
        sys.exit(code)
