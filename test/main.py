#!/usr/bin/env python3
"""
main.py — test runner entry point for the ECHO pipeline test suite.

Usage examples:
  # Run all tests
  python test/main.py

  # Run a specific test file
  python test/main.py --test test_selection

  # Run a specific test function
  python test/main.py --test test_selection::test_more_taxa_than_n_selects_exactly_n

  # Run with extra pytest flags (e.g. stop on first failure)
  python test/main.py --test test_smoke -x

  # Run quietly (no verbose)
  python test/main.py --verbose false
"""

import argparse
import logging
import sys
import os

import pytest


TEST_DIR = os.path.dirname(os.path.abspath(__file__))

log = logging.getLogger("echo.test")


def setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=level,
        stream=sys.stdout,
    )

AVAILABLE_TESTS = [
    "test_selection",
    "test_dedup",
    "test_smoke",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="ECHO pipeline test runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"Available test files: {', '.join(AVAILABLE_TESTS)}",
    )
    parser.add_argument(
        "--test", "-t",
        default=None,
        help=(
            "Test to run. Can be a file name (e.g. test_selection), "
            "or file::function (e.g. test_selection::test_more_taxa_than_n_selects_exactly_n). "
            "Omit to run all tests."
        ),
    )
    parser.add_argument(
        "--verbose", "-v",
        default="true",
        choices=["true", "false"],
        help="Verbose output (default: true)",
    )
    parser.add_argument(
        "-x",
        action="store_true",
        help="Stop after first failure",
    )
    parser.add_argument(
        "--cluster-id-mode",
        default="same",
        choices=["same", "random"],
        help="Cluster ID mode for selection tests: "
             "'same' = all proteins share one cluster ID (default), "
             "'random' = each protein gets a random cluster ID.",
    )
    parser.add_argument(
        "--tb",
        default="short",
        choices=["short", "long", "no", "line"],
        help="Traceback style on failure (default: short)",
    )
    return parser.parse_args()


def build_pytest_args(args):
    pytest_args = []

    # test target
    if args.test:
        # support file-only (test_selection) or file::function
        if "::" in args.test:
            file_part, func_part = args.test.split("::", 1)
            pytest_args.append(os.path.join(TEST_DIR, f"{file_part}.py::{func_part}"))
        else:
            pytest_args.append(os.path.join(TEST_DIR, f"{args.test}.py"))
    else:
        pytest_args.append(TEST_DIR)

    pytest_args += ["--cluster-id-mode", args.cluster_id_mode]

    if args.verbose == "true":
        pytest_args.append("-v")
    # console always shows INFO only; DEBUG goes to per-test debug.log files
    pytest_args += ["--log-cli-level", "INFO"]

    # show log format consistently with setup_logging
    pytest_args += [
        "--log-cli-format", "%(asctime)s [%(levelname)-8s] %(name)s: %(message)s",
        "--log-cli-date-format", "%H:%M:%S",
    ]

    if args.x:
        pytest_args.append("-x")

    pytest_args += ["--tb", args.tb]

    return pytest_args


def main():
    args = parse_args()
    setup_logging(verbose=args.verbose == "true")

    log.info("ECHO pipeline test runner starting")
    log.debug("Test directory: %s", TEST_DIR)

    if args.test:
        log.info("Target test: %s", args.test)
    else:
        log.info("Target test: all (%s)", ", ".join(AVAILABLE_TESTS))

    pytest_args = build_pytest_args(args)
    log.debug("pytest args: %s", pytest_args)

    log.info("Running: pytest %s", " ".join(pytest_args))
    exit_code = pytest.main(pytest_args)

    if exit_code == 0:
        log.info("All tests passed")
    else:
        log.error("Tests failed (exit code %d)", exit_code)

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
