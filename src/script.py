"""Compatibility entry point. Running without --publish is a read-only dry run."""
from tracker.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
