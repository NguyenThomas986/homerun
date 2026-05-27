"""
HomerunState.py
───────────────
Decides which pipeline steps still need to execute for a given sample row.

Logic per step column:
  ""          → PENDING   → needs to run
  "PENDING"   →           → needs to run
  "FAILED"    →           → needs to run if --retry-failed, else skip
  "DONE"      →           → skip
  any other   → treat as an output path → step is DONE

The step ordering dependency is enforced: if an upstream step is not DONE,
all downstream steps are also skipped (they cannot proceed without inputs).
"""

from __future__ import annotations

from typing import Dict, Any, List, Optional

STEP_ORDER = ["trim", "star", "tagdir", "qc", "deseq2"]

# Statuses that mean "this step has not yet succeeded"
NOT_DONE = {"", "PENDING", "FAILED", None}


def _is_done(value: Any) -> bool:
    """A step is done when its column holds any non-empty non-pending value."""
    if value is None:
        return False
    s = str(value).strip()
    return s not in NOT_DONE


def _is_failed(value: Any) -> bool:
    return str(value).strip().upper() == "FAILED"


def resolve_steps(
    row: Dict[str, Any],
    allowed_steps: Optional[List[str]] = None,
    retry_failed: bool = False,
) -> List[str]:
    """
    Return the ordered list of step names that should execute for this sample.

    Parameters
    ----------
    row           : dict of CSV column → value for one sample
    allowed_steps : if given, only steps in this list are eligible
    retry_failed  : if True, FAILED steps are reset to PENDING and re-run
    """
    to_run: List[str] = []
    upstream_blocked = False   # once an upstream step is not-done, block the rest

    for step in STEP_ORDER:
        val = row.get(step, "")

        # If an earlier step hasn't finished, downstream steps can't run
        if upstream_blocked:
            break

        # Limit to explicitly requested steps
        if allowed_steps and step not in allowed_steps:
            # If the step is not done, block downstream
            if not _is_done(val):
                upstream_blocked = True
            continue

        if _is_done(val):
            # Already complete — do nothing, don't block
            continue

        if _is_failed(val) and not retry_failed:
            # Failed but not retrying — block downstream
            upstream_blocked = True
            continue

        # Needs to run (PENDING / empty / FAILED+retry_failed)
        to_run.append(step)

        # A step that needs to run blocks all subsequent steps
        # (they will be scheduled in the same run, not blocked)
        # Only block if we're NOT going to run this step right now.
        # We handle the actual blocking inside process_sample on failure.

    return to_run
