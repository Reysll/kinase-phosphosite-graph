"""
PostToolUse hook: run ruff format on any Python file Claude just wrote or edited.
Receives Claude's tool response JSON on stdin.
"""
import json
import subprocess
import sys
from pathlib import Path

d = json.load(sys.stdin)
file_path = d.get("tool_input", {}).get("file_path", "")

if not file_path or not file_path.endswith(".py"):
    sys.exit(0)

project_root = Path(__file__).resolve().parent.parent.parent
ruff = project_root / ".venv312" / "Scripts" / "ruff.exe"

if not ruff.exists():
    sys.exit(0)

subprocess.run([str(ruff), "format", file_path], check=False)
