"""Run the compiled-kernel validation suite and assert only the known a2 conflict fails."""
import re
import subprocess
import sys
from pathlib import Path

SUITE = Path(__file__).with_name("kernel_suite.py")


def test_kernel_suite():
    proc = subprocess.run([sys.executable, "-B", str(SUITE)],
                          capture_output=True, text=True)
    summary = re.search(r"SUMMARY: (\d+)/(\d+) checks passed", proc.stdout)
    assert summary, proc.stdout[-3000:] + proc.stderr[-1000:]
    passed, total = int(summary.group(1)), int(summary.group(2))
    assert proc.returncode == 0, proc.stdout[-3000:]
    assert "FAILED: a2 " in proc.stdout
    assert passed == total - 1, (passed, total)


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
