import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from parmmosahn.utils import ToolExecutionError, run_command


def test_run_command_success():
    """Test that run_command returns stdout when the subprocess succeeds."""
    cmd = ["echo", "hello"]

    # We patch subprocess.run so it doesn't actually execute on the system
    with patch("subprocess.run") as mock_run:
        # Create a mock CompletedProcess object
        mock_run.return_value = subprocess.CompletedProcess(
            args=cmd, returncode=0, stdout="hello\n"
        )

        output = run_command(cmd)

        assert output == "hello\n"
        mock_run.assert_called_once_with(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
        )


def test_run_command_failure():
    """Test that run_command raises ToolExecutionError when the subprocess fails."""
    cmd = ["false"]

    with patch("subprocess.run") as mock_run:
        # Simulate a tool failure
        mock_run.side_effect = subprocess.CalledProcessError(
            returncode=1, cmd=cmd, stderr="Command not found"
        )

        with pytest.raises(ToolExecutionError) as exc_info:
            run_command(cmd)

        assert "Tool failed with exit code 1" in str(exc_info.value)


def test_run_command_with_file_output(tmp_path: Path):
    """Test run_command when redirecting stdout and stderr to files."""
    cmd = ["some", "tool"]
    stdout_file = tmp_path / "out.log"
    stderr_file = tmp_path / "err.log"

    with patch("subprocess.run") as mock_run:
        # No return value needed when outputs are redirected to files
        mock_run.return_value = subprocess.CompletedProcess(args=cmd, returncode=0)

        output = run_command(cmd, stdout_file=stdout_file, stderr_file=stderr_file)

        assert output is None  # Should return None when writing to files
        mock_run.assert_called_once()

        # Verify it used the file handles
        called_kwargs = mock_run.call_args.kwargs
        assert called_kwargs["stdout"].name == str(stdout_file)
        assert called_kwargs["stderr"].name == str(stderr_file)
