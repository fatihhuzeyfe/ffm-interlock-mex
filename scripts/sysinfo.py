"""Memory / CPU information without Task Manager (standard library only, Windows API via ctypes).

    python scripts/sysinfo.py            # one snapshot
    python scripts/sysinfo.py --watch 30 # every 30 s (Ctrl+C to stop)
"""
import os, sys, time, ctypes, subprocess


def memory_gb():
    """(total, available) physical memory in GB, or (None, None) if unknown."""
    if os.name == "nt":
        class MEMSTAT(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
        st = MEMSTAT()
        st.dwLength = ctypes.sizeof(MEMSTAT)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
            return st.ullTotalPhys / 1e9, st.ullAvailPhys / 1e9
        return None, None
    try:                                         # Linux
        info = dict(l.split(":", 1) for l in open("/proc/meminfo"))
        kb = lambda k: float(info[k].split()[0]) / 1e6
        return kb("MemTotal"), kb("MemAvailable")
    except Exception:
        return None, None


ABQ_IMAGES = ("abqcaek.exe", "standard.exe", "pre.exe", "package.exe", "explicit.exe",
              "abqsmastddriver.exe")


def abaqus_process_list(mine_only=False):
    """[(image name, PID, owner)] of running Abaqus executables (exact image names).
    The owner is only visible for the current user's processes ('N/A' for other users)."""
    if os.name != "nt":
        return []
    try:
        out = subprocess.run("tasklist /V /FO CSV /NH", shell=True, capture_output=True, text=True,
                             errors="ignore").stdout
    except Exception:
        return []
    me = os.environ.get("USERNAME", "").lower()
    res = []
    for line in out.splitlines():
        f = [x.strip('"') for x in line.strip().split('","')]
        if f and f[0].lower() in ABQ_IMAGES:
            owner = f[6] if len(f) > 6 else ""
            mine = bool(me) and owner.lower().rsplit(chr(92), 1)[-1] == me   # DOMAIN\\user
            if mine or not mine_only:
                res.append((f[0], f[1] if len(f) > 1 else "?", owner if owner else "?"))
    return res


def abaqus_processes(mine_only=False):
    """Number of running Abaqus executables (CAE kernels and solvers)."""
    if os.name != "nt":
        return None
    return len(abaqus_process_list(mine_only))


def snapshot():
    tot, av = memory_gb()
    n = abaqus_processes()
    s = time.strftime("%H:%M:%S")
    if tot:
        s += f"  RAM: {av:.1f} GB free of {tot:.1f} GB ({100 * (1 - av / tot):.0f} % used)"
    if n is not None:
        s += f"  Abaqus processes: {n}"
        if 0 < n <= 10:
            s += "  [" + ", ".join(f"{a} pid {b} owner {c}" for a, b, c in abaqus_process_list()) + "]"
    return s


if __name__ == "__main__":
    if "--watch" in sys.argv:
        dt = float(sys.argv[sys.argv.index("--watch") + 1])
        while True:
            print(snapshot(), flush=True)
            time.sleep(dt)
    print(snapshot())
