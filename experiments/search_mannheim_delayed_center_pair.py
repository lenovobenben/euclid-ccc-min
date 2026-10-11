"""五步双目标后缀：非圆心首对象，随后两个双圆心载体，再画两目标圆。

用未绘制的圆心连线、垂直平分线和半径轨迹筛选首对象；虚拟交点
只供检验，不加入可用点集。目标是 25+5+3*6=48 E，数值命中需精确验证。
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from time import monotonic

from mannheim_search_calibration import BindingAudit, apply_calibrated_candidate, close_point
from mannheim_virtual_witnesses import bisector, witness_helpers
from replay_mannheim_three_block_dependencies import collapse_point
from replay_mannheim_ordered_branches import branch_data
from search_mannheim_center_locus_2e import object_value
from search_mannheim_double_kp_global_sequential_2e import float_scalar
from search_mannheim_joint_root_kp import PROFILES, build_carrier_state, fast_target_candidates
from audit_mannheim_shared_projections import CORE_AFTER_SHARED
from search_mannheim_root_center_2e import drawable_contains_targets
from search_mannheim_shared_circle_replacement import available, dual_circle_candidates, eligible_centers
from search_parallel_3e import DrawableBundle, intersections


def pair_carriers(state, first_targets, second_targets, move):
    ids = frozenset(point.point_id for point in state.points)
    lines = tuple(candidate for candidate in fast_target_candidates(state, first_targets, ids, move)
                  if candidate.kind == "line" and drawable_contains_targets(candidate.drawable, second_targets))
    return lines + dual_circle_candidates(state, first_targets, second_targets, move)


def pair_data(state, replays, profile):
    sample_targets = tuple(replay.verify_pair(profile, branch_data(replay, profile)["tau"],
                                            allow_repeated_physical_signs=True) for replay in replays)
    keys = tuple(sorted(sample_targets[0]))
    if len(keys) != 2 or any(tuple(sorted(row)) != keys for row in sample_targets):
        raise AssertionError("配对目标键必须一致且恰有两个")
    targets = tuple(tuple(tuple(float_scalar(value) for value in collapse_point(row[key]["center"]))
                          for row in sample_targets) for key in keys)
    contacts = tuple(tuple(tuple(float_scalar(value) for value in collapse_point(row[key]["contact_3"]))
                           for row in sample_targets) for key in keys)
    if any(available(state, target) for target in targets) or not all(available(state, contact) for contact in contacts):
        raise AssertionError("后缀初态应有两个接触点，但没有目标圆心")
    return keys, targets


def run_case(state, replays, profile, seconds):
    started = monotonic()
    deadline = started + seconds
    keys, targets = pair_data(state, replays, profile)
    first_targets, second_targets = targets
    if any(drawable_contains_targets(item, values)
           for item in state.drawables for values in targets):
        raise AssertionError("延迟载体初态不能已有经过任一目标圆心的对象")
    if pair_carriers(state, first_targets, second_targets, 26):
        raise AssertionError("初态已有可画共同载体，需保留无需新定义点的首步形状")
    axis = DrawableBundle("virtual_axis", "line", tuple(object_value(replay.objects, f"{profile}_center_locus_from_Kp")[1]
                           for replay in replays))
    perpendicular = DrawableBundle("virtual_bisector", "line", tuple(bisector(p, q)
                                    for p, q in zip(first_targets, second_targets, strict=True)))
    centers = eligible_centers(state, first_targets, second_targets)
    goals = [axis, perpendicular]
    for center in centers:
        values = tuple((p[0], p[1], (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2)
                       for p, q in zip(center.values, first_targets, strict=True))
        goals.append(DrawableBundle(f"virtual_radius_{center.point_id}", "circle", values))
    helpers, witness_report = witness_helpers(state, goals, 26, deadline)
    report = {"profile": profile, "target_keys": keys, "policy": "known-first",
              "base_e": len(state.program), "base_points": len(state.points), "base_objects": len(state.drawables),
              "static_eligible_centers": len(centers), "witness_screen": witness_report,
              "initial_center_objects": [sum(drawable_contains_targets(item, target) for item in state.drawables)
                                         for target in targets],
              "helper_examples": [candidate.describe() for candidate in helpers[:10]],
              "candidate_fingerprint": hashlib.sha256("\n".join(c.describe() for c in helpers).encode()).hexdigest(),
              "checked_helpers": 0, "skipped_early_center_carriers": 0,
              "second_carriers": 0, "third_carriers": 0, "carrier_examples": [], "hits": [], "budget_seconds": seconds}
    audit = BindingAudit()
    for first in helpers:
        if monotonic() >= deadline:
            break
        report["checked_helpers"] += 1
        if any(drawable_contains_targets(first.drawable, values) for values in targets):
            report["skipped_early_center_carriers"] += 1
            continue
        one = apply_calibrated_candidate(state, first, 26, audit=audit)
        seconds_candidates = pair_carriers(one, first_targets, second_targets, 27)
        report["second_carriers"] += len(seconds_candidates)
        for second in seconds_candidates:
            two = apply_calibrated_candidate(one, second, 27, audit=audit)
            thirds = pair_carriers(two, first_targets, second_targets, 28)
            report["third_carriers"] += len(thirds)
            if len(report["carrier_examples"]) < 10:
                report["carrier_examples"].append({"first": first.describe(), "second": second.describe(),
                                                   "third_count": len(thirds)})
            for third in thirds:
                rows = tuple(intersections(second.kind, p, third.kind, q)
                             for p, q in zip(second.drawable.values, third.drawable.values, strict=True))
                if all(all(any(close_point(point, target) for point in row) for target in pair)
                       for row, pair in zip(rows, zip(*targets, strict=True), strict=True)):
                    report["hits"].append([first.describe(), second.describe(), third.describe()])
        if report["hits"]:
            break
    report["complete"] = witness_report["complete"] and report["checked_helpers"] == len(helpers)
    report["binding_audit"] = asdict(audit)
    report["elapsed_seconds"] = round(monotonic() - started, 3)
    print(json.dumps(report, ensure_ascii=False), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds-per-case", type=float, default=30)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.seconds_per_case <= 0:
        parser.error("预算必须为正")
    state, replays, audit = build_carrier_state(extra_ids=CORE_AFTER_SHARED)
    report = {"schema": "delayed-center-pair/v1", "base_binding_audit": asdict(audit), "cases": []}
    for profile in PROFILES:
        report["cases"].append(run_case(state, replays, profile, args.seconds_per_case))
        Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        if report["cases"][-1]["hits"]:
            break


if __name__ == "__main__":
    main()
