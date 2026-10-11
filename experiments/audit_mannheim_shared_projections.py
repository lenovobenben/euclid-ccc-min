"""25 E 接触核心后，以一圆同时提供两个圆心线截点 J 的静态检查。

一圆与已有 ell 的两交点若是 J_i,J_j，可替代两条 O3 K' 辅助线，
令两套 6 E 后缀共用该圆，净省 1 E。只作数值发现，不证明下界。
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from itertools import combinations
from pathlib import Path
from time import monotonic

from search_mannheim_center_locus_2e import point_value
from search_mannheim_joint_root_kp import PROFILES, build_carrier_state
from search_mannheim_shared_circle_replacement import available, dual_circle_candidates, eligible_centers


CORE_AFTER_SHARED = ("P0_K_core_line_1", "P0_K_core_line_2", "P0_tau",
                     "P2_K_core_line_1", "P2_K_core_line_2", "P2_tau",
                     "P1_tau_from_S", "P3_tau_from_S")


def screen(selected=(), mode="radial"):
    started = monotonic()
    planned = tuple(node for profile in selected for node in (
        (f"{profile}_O3_Kp",) if mode == "radial"
        else (f"{profile}_O3_Kp", f"{profile}_center_locus_from_Kp")))
    state, replays, audit = build_carrier_state(extra_ids=CORE_AFTER_SHARED + planned)
    root = tuple(point_value(replay.objects, "Mannheim_S_center_locus") for replay in replays)
    if len(state.program) != 25 + len(planned) or not available(state, root):
        raise AssertionError("必须合法重放 25 E 接触核心，且已有根心")
    targets = {profile: tuple(point_value(replay.objects, f"{profile}_center_locus_J") for replay in replays)
               for profile in PROFILES}
    pairs = []
    remaining = tuple(profile for profile in PROFILES if profile not in selected)
    for first, second in combinations(remaining, 2):
        candidates = dual_circle_candidates(state, targets[first], targets[second], len(state.program) + 1)
        pairs.append({"profiles": [first, second],
                      "eligible_centers": len(eligible_centers(state, targets[first], targets[second])),
                      "candidates": [candidate.describe() for candidate in candidates]})
    report = {"schema": "shared-projections-static/v1", "base_e": len(state.program),
              "planned_profiles": selected, "planned_mode": mode, "planned_objects": planned,
              "policy": "known-first", "base_points": len(state.points), "base_objects": len(state.drawables),
              "already_available": [profile for profile, values in targets.items() if available(state, values)],
              "pairs": pairs, "elapsed_seconds": round(monotonic() - started, 3),
              "binding_audit": asdict(audit),
              "complete": True, "scope": "static circles through two J points; no dynamic helpers"}
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--planned-helpers", action="store_true")
    args = parser.parse_args()
    report = screen()
    if args.planned_helpers:
        states = [report]
        for count in (1, 2):
            for selected in combinations(PROFILES, count):
                for mode in ("radial", "locus"):
                    states.append(screen(selected, mode))
        report = {"schema": "shared-projections-planned/v1", "states": states}
    Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
