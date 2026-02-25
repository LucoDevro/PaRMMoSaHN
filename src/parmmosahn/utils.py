import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


class ToolExecutionError(Exception):
    """Custom exception raised when an external bio-tool fails."""

    pass


def run_command(
    cmd: list[str], stdout_file: Path | None = None, stderr_file: Path | None = None
) -> str | None:
    """
    Runs an external command robustly.
    If stdout_file and stderr_file are provided, outputs are written to them instead of captured.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    try:
        # Check if BOTH files are provided for file-based output
        if stdout_file is not None and stderr_file is not None:
            with stdout_file.open("a") as out, stderr_file.open("a") as err:
                subprocess.run(cmd, stdout=out, stderr=err, text=True, check=True)
            return None

        # In all other cases (including if only one file is provided), capture to string
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        return result.stdout

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {' '.join(cmd)}")
        if not (stdout_file and stderr_file):
            logger.error(f"Standard Error output:\n{e.stderr}")
        else:
            logger.error(f"Check {stderr_file} for detailed error logs.")
        raise ToolExecutionError(f"Tool failed with exit code {e.returncode}") from e
