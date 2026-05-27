from pathlib import Path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]
CPP_BUILD_DIR = REPO_ROOT / "cpp" / "build"


def ensure_cpp_binaries():
    subprocess.run(["make", "-C", "cpp"], cwd=REPO_ROOT, check=True)
    logsort = CPP_BUILD_DIR / "logsort"
    logthin = CPP_BUILD_DIR / "logthin"
    if not (logsort.exists() and logthin.exists()):
        raise RuntimeError("Expected C++ binaries were not produced by the build.")
    return logsort, logthin
