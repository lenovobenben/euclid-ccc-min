"""保留三条共享弦，用一个辅助对象和一个圆恢复缺少的两个 K'。

若再用至多 2 E 取得 S，则 13+3+2+2+4+24=48 E。
数值候选发现；固定现有后缀接口，不覆盖两个独立载体等其它替换。
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from time import monotonic

from mannheim_search_calibration import apply_calibrated_candidate, close_point
from mannheim_virtual_witnesses import bisector, witness_helpers
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_joint_root_kp import (
    PROFILES, SHARED_IDS, build_carrier_state, fast_target_candidates,
)
from search_mannheim_root_center_2e import drawable_contains_targets
from search_mannheim_root_center_3e import equidistant_from_targets
from search_parallel_3e import (
    Candidate, DrawableBundle, circle_through, drawable_signature, generate_candidates, intersections,
)


def available(state, targets):
    return any(all(close_point(p, q) for p, q in zip(point.values, targets, strict=True))
               for point in state.points)


def eligible_centers(state, first_targets, second_targets):
    """圆通过两目标的必要条件：圆心到它们等距。无需绘制垂直平分线。"""
    result = []
    for point in state.points:
        for center, first, second in zip(point.values, first_targets, second_targets, strict=True):
            left = (first[0] - center[0]) ** 2 + (first[1] - center[1]) ** 2
            right = (second[0] - center[0]) ** 2 + (second[1] - center[1]) ** 2
            if abs(left - right) > 1e-8 * max(1.0, left, right):
                break
        else:
            result.append(point)
    return tuple(result)


def dual_circle_candidates(state, first_targets, second_targets, move):
    existing = {drawable_signature(item) for item in state.drawables}
    seen, result = set(), []
    for center in eligible_centers(state, first_targets, second_targets):
        for through in state.points:
            if center.point_id == through.point_id:
                continue
            if not equidistant_from_targets(center, through, first_targets):
                continue
            if not equidistant_from_targets(center, through, second_targets):
                continue
            values = tuple(circle_through(p, q) for p, q in zip(center.values, through.values, strict=True))
            if any(value is None for value in values):
                continue
            drawable = DrawableBundle(f"move_{move}", "circle", values)
            signature = drawable_signature(drawable)
            if signature in existing or signature in seen:
                continue
            seen.add(signature)
            result.append(Candidate("circle", center.point_id, through.point_id, drawable))
    return tuple(result)


def two_step_root_paths(state, targets, move, audit):
    if available(state, targets):
        return [()]
    old_carriers = tuple(item for item in state.drawables if drawable_contains_targets(item, targets))
    ids = frozenset(point.point_id for point in state.points)
    firsts = fast_target_candidates(state, targets, ids, move)
    paths = []
    for first in firsts:
        for old in old_carriers:
            rows = tuple(intersections(first.kind, p, old.kind, q)
                         for p, q in zip(first.drawable.values, old.values, strict=True))
            if all(any(close_point(p, target) for p in row) for row, target in zip(rows, targets, strict=True)):
                return [(first.describe(),)]
        one = apply_calibrated_candidate(state, first, move, audit=audit)
        seconds = fast_target_candidates(one, targets, frozenset(p.point_id for p in one.points), move + 1)
        for second in seconds:
            rows = tuple(intersections(first.kind, p, second.kind, q)
                         for p, q in zip(first.drawable.values, second.drawable.values, strict=True))
            if all(any(close_point(p, target) for p in row) for row, target in zip(rows, targets, strict=True)):
                paths.append((first.describe(), second.describe()))
    return paths


def run_case(omitted, seconds, helper_kind="both", witness_pruning=False):
    started = monotonic()
    kept = tuple(node for node in SHARED_IDS if node != omitted)
    state, replays, audit = build_carrier_state(kept)
    graph = replays[0].objects.graph
    root = tuple(point_value(replay.objects, "Mannheim_S_center_locus") for replay in replays)
    kps = {profile: tuple(point_value(replay.objects, graph.dependencies[f"{profile}_O3_Kp"][1])
                         for replay in replays) for profile in PROFILES}
    missing = [profile for profile, values in kps.items() if not available(state, values)]
    if len(missing) != 2 or len(state.program) != 16:
        raise AssertionError("三个共享弦应保留两个 K'，另两个尚未出现")
    first_targets, second_targets = (kps[profile] for profile in missing)
    static = fast_target_candidates(state, first_targets, frozenset(p.point_id for p in state.points), 17)
    static = tuple(candidate for candidate in static if drawable_contains_targets(candidate.drawable, second_targets))
    omitted_drawable = DrawableBundle(omitted, "line", tuple(object_value(replay.objects, omitted)[1]
                                     for replay in replays))
    omitted_signature = drawable_signature(omitted_drawable)
    all_firsts = generate_candidates(state, 17)
    first_pool = all_firsts
    witness_report = None
    if witness_pruning:
        goals = [DrawableBundle("virtual_bisector", "line", tuple(bisector(p, q)
                                for p, q in zip(first_targets, second_targets, strict=True)))]
        for center in eligible_centers(state, first_targets, second_targets):
            values = tuple((p[0], p[1], (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2)
                           for p, q in zip(center.values, first_targets, strict=True))
            goals.append(DrawableBundle("virtual_radius", "circle", values))
        first_pool, witness_report = witness_helpers(state, goals, 17, started + seconds)
    firsts = tuple(sorted((candidate for candidate in first_pool
                          if helper_kind == "both" or candidate.kind == helper_kind),
                         key=lambda candidate: candidate.describe()))
    report = {"omitted": omitted, "kept": kept, "policy": "known-first", "missing_profiles": missing,
              "base_e": 16, "base_points": len(state.points), "base_objects": len(state.drawables),
              "static_dual_carriers": [candidate.describe() for candidate in static],
              "static_alternatives": sum(drawable_signature(c.drawable) != omitted_signature for c in static),
              "static_eligible_circle_centers": len(eligible_centers(state, first_targets, second_targets)),
              "first_root_carriers": sum(drawable_contains_targets(c.drawable, root) for c in all_firsts),
              "helper_kind": helper_kind, "all_first_candidates": len(all_firsts),
              "witness_pruning": witness_pruning, "witness_screen": witness_report,
              "first_candidates": len(firsts), "checked_firsts": 0, "dual_circles": 0,
              "center_helpers": 0, "restored_missing_points": 0, "binding_failures": 0,
              "skipped_standard_line": 0, "hits": [], "complete": False, "budget_seconds": seconds,
              "candidate_fingerprint": hashlib.sha256("\n".join(c.describe() for c in firsts).encode()).hexdigest(),
              "base_binding_audit": asdict(audit)}
    print(json.dumps({"initial": report}, ensure_ascii=False), flush=True)
    next_progress = monotonic() + 15
    for first in firsts:
        if monotonic() - started >= seconds:
            break
        report["checked_firsts"] += 1
        if drawable_signature(first.drawable) == omitted_signature:
            # 本轮只研究替代载体，排除还原标准弦的辅助操作。
            report["skipped_standard_line"] += 1
            continue
        one = apply_calibrated_candidate(state, first, 17, audit=audit)
        if eligible_centers(one, first_targets, second_targets):
            report["center_helpers"] += 1
        circles = dual_circle_candidates(one, first_targets, second_targets, 18)
        report["dual_circles"] += len(circles)
        for circle in circles:
            two = apply_calibrated_candidate(one, circle, 18, audit=audit)
            if not all(available(two, kps[profile]) for profile in missing):
                report["binding_failures"] += 1
                continue
            report["restored_missing_points"] += 1
            for root_path in two_step_root_paths(two, root, 19, audit):
                report["hits"].append([first.describe(), circle.describe(), *root_path])
        if report["hits"]:
            break
        if monotonic() >= next_progress:
            print(json.dumps({"omitted": omitted, "checked": report["checked_firsts"],
                              "center_helpers": report["center_helpers"], "dual_circles": report["dual_circles"]}), flush=True)
            next_progress = monotonic() + 15
    report["complete"] = report["checked_firsts"] == len(firsts) and (witness_report is None or witness_report["complete"])
    report["next_index"] = report["checked_firsts"]
    report["elapsed_seconds"] = round(monotonic() - started, 3)
    report["binding_audit"] = asdict(audit)
    print(json.dumps(report, ensure_ascii=False), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds-per-case", type=float, default=15)
    parser.add_argument("--helper-kind", choices=("both", "line", "circle"), default="both")
    parser.add_argument("--witness-pruning", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.seconds_per_case <= 0:
        parser.error("时间预算必须为正")
    report = {"schema": "shared-circle-replacement/v1", "cases": []}
    for omitted in SHARED_IDS:
        report["cases"].append(run_case(omitted, args.seconds_per_case, args.helper_kind, args.witness_pruning))
        # 各子例正常完成后即保存；预算包含该例的准备，检查发生在首步之间。
        Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        if report["cases"][-1]["hits"]:
            break


if __name__ == "__main__":
    main()
