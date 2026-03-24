# -*- coding: utf-8 -*-
"""
CLI entry point for AutoPoly agent API.

All output is JSON to stdout, errors to stderr. Exit 0 on success, 1 on failure.

Usage:
    autopoly info
    autopoly validate <config.json>
    autopoly validate --stdin
    autopoly generate <config.json>
    autopoly generate --stdin
    autopoly describe <smiles_string>
"""
import argparse
import json
import sys

from . import agent


def main():
    parser = argparse.ArgumentParser(
        prog="autopoly",
        description="AutoPoly: config-driven polymer system generation for LAMMPS",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # info
    subparsers.add_parser("info", help="Show available options, limits, and examples")

    # validate
    val_parser = subparsers.add_parser("validate", help="Validate a config without generating")
    val_parser.add_argument("config_file", nargs="?", help="Path to config JSON file")
    val_parser.add_argument("--stdin", action="store_true", help="Read config from stdin")

    # generate
    gen_parser = subparsers.add_parser("generate", help="Generate a polymer system")
    gen_parser.add_argument("config_file", nargs="?", help="Path to config JSON file")
    gen_parser.add_argument("--stdin", action="store_true", help="Read config from stdin")

    # describe
    desc_parser = subparsers.add_parser("describe", help="Describe a SMILES string")
    desc_parser.add_argument("smiles", help="SMILES string to describe")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "info":
        _output(agent.info())

    elif args.command == "validate":
        config = _load_config(args)
        result = agent.validate(config)
        _output(result.to_dict(), exit_code=0 if result.success else 1)

    elif args.command == "generate":
        config = _load_config(args)
        result = agent.generate(config)
        _output(result.to_dict(), exit_code=0 if result.success else 1)

    elif args.command == "describe":
        _output(agent.describe_smiles(args.smiles))


def _load_config(args) -> dict:
    """Load config from file or stdin."""
    if args.stdin:
        raw = sys.stdin.read()
    elif args.config_file:
        try:
            with open(args.config_file) as f:
                raw = f.read()
        except FileNotFoundError:
            _error(f"File not found: {args.config_file}")
        except PermissionError:
            _error(f"Permission denied: {args.config_file}")
    else:
        _error("Provide a config file path or use --stdin")
        return {}  # unreachable

    try:
        return json.loads(raw)
    except json.JSONDecodeError as e:
        _error(f"Invalid JSON: {e}")
        return {}  # unreachable


def _output(data: dict, exit_code: int = 0):
    """Print JSON to stdout and exit."""
    print(json.dumps(data, indent=2))
    sys.exit(exit_code)


def _error(msg: str):
    """Print error to stderr and exit 1."""
    print(json.dumps({"error": msg}), file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
