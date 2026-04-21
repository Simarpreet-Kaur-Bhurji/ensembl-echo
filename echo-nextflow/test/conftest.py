"""
conftest.py — pytest configuration for ECHO test suite.

Logging
-------
  Console : INFO level — summary lines only (set via --log-cli-level in main.py)
  File    : two session-level files written to echo-nextflow/test/logs/:
              info_<timestamp>.log  — INFO and above (mirrors console)
              debug_<timestamp>.log — full DEBUG trace for all tests

To see console logs when running pytest directly:
  pytest echo-nextflow/test/ -v --log-cli-level=INFO
"""

import logging
import os
from datetime import datetime

import pytest

LOG_FORMAT      = "%(asctime)s [%(levelname)-8s] %(name)s: %(message)s"
LOG_DATE_FORMAT = "%H:%M:%S"
LOGS_DIR        = os.path.join(os.path.dirname(__file__), "logs")


def pytest_addoption(parser):
    parser.addoption(
        "--cluster-id-mode",
        default="same",
        choices=["same", "random"],
        help="Cluster ID mode for selection tests: "
             "'same' = all proteins share one cluster ID (default), "
             "'random' = each protein gets a random cluster ID.",
    )
    parser.addoption(
        "--run-pipeline",
        action="store_true",
        default=False,
        help="(test_e2e) Execute the Nextflow pipeline before validating outputs. "
             "Requires 'nextflow' on PATH and the singularity image in params.test.yaml.",
    )
    parser.addoption(
        "--outdir",
        default=None,
        help="(test_e2e) Path to an existing pipeline output directory to validate. "
             "Overrides the default nextflow_test/ location.",
    )


def pytest_configure(config):
    """Set root logger level so DEBUG records reach file handlers."""
    logging.basicConfig(
        format=LOG_FORMAT,
        datefmt=LOG_DATE_FORMAT,
        level=logging.DEBUG,
    )


@pytest.fixture(scope="session")
def cluster_id_mode(request):
    """Returns 'same' or 'random' depending on --cluster-id-mode CLI option."""
    return request.config.getoption("--cluster-id-mode")


@pytest.fixture(scope="session", autouse=True)
def log_dir():
    """
    Session-scoped: creates echo-nextflow/test/logs/ and attaches two FileHandlers
    for the entire session:
      info_<timestamp>.log  — INFO and above (mirrors console output)
      debug_<timestamp>.log — full DEBUG trace across all tests
    Yields the logs/ directory path so tests can store artifacts in subdirectories.
    """
    os.makedirs(LOGS_DIR, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    formatter = logging.Formatter(fmt=LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    info_handler = logging.FileHandler(os.path.join(LOGS_DIR, f"info_{timestamp}.log"))
    info_handler.setLevel(logging.INFO)
    info_handler.setFormatter(formatter)

    debug_handler = logging.FileHandler(os.path.join(LOGS_DIR, f"debug_{timestamp}.log"))
    debug_handler.setLevel(logging.DEBUG)
    debug_handler.setFormatter(formatter)

    root = logging.getLogger()
    root.addHandler(info_handler)
    root.addHandler(debug_handler)

    logging.getLogger("echo.test").info("session log dir: %s", LOGS_DIR)

    yield LOGS_DIR

    logging.getLogger("echo.test").info("session complete — logs saved to %s", LOGS_DIR)
    for h in (info_handler, debug_handler):
        h.close()
        root.removeHandler(h)
