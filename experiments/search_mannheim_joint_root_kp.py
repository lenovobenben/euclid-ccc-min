"""有时间预算的联合核心筛查：四条 K' 共享弦后用至多 3 E 取得根心。

13 E 公共准备 + 4 E 共享弦 + 3 E 根心 + 4 E 接触弦 + 24 E 后缀
若成功并保留所需接口，给出 48 E 候选。仅覆盖共享弦先行子类。
浮点束、已知根排除和全部根配对均用于发现候选；不是下界证明。
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
from time import monotonic

import numpy as np

from mannheim_search_calibration import BindingAudit, apply_calibrated_candidate, close_point
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_root_center_2e import FIXTURES, drawable_contains_targets
from search_mannheim_root_center_3e import collinear_with_targets, equidistant_from_targets
from search_parallel_3e import (
    Candidate, DrawableBundle, PointBundle, State, circle_through, drawable_signature,
    generate_candidates, intersections, line_through, point_signature,
)


SHARED_IDS = ("P0_Kp_core_line_1", "P0_Kp_core_line_2",
              "P2_Kp_core_line_1", "P1_Kp_locus_line_2")
PROFILES = ("P0", "P2", "P1", "P3")


def bind_named_intersections(state, replays, object_ids):
    """权威重放中由现有对象求交得到的命名点，保留已证明的几何身份。

    只加入父对象已经存在、且各样本都实际位于父对象交集的点。
    不把根心/目标数值注入搜索状态。
    """

    seen = {point_signature(point) for point in state.points}
    additions = []
    graph = replays[0].objects.graph
    for node_id in sorted(set(replays[0].objects.point_registry.values())):
        dependencies = graph.dependencies[node_id]
        if len(dependencies) != 2 or not set(dependencies) <= object_ids:
            continue
        values = tuple(point_value(replay.objects, node_id) for replay in replays)
        for replay, value in zip(replays, values, strict=True):
            parents = [object_value(replay.objects, parent) for parent in dependencies]
            roots = intersections(parents[0][0], parents[0][1], parents[1][0], parents[1][1])
            if not any(close_point(root, value) for root in roots):
                raise AssertionError(f"{node_id} 不是已画父对象的有限交点")
        point = PointBundle(node_id, values)
        signature = point_signature(point)
        if signature not in seen:
            seen.add(signature)
            additions.append(point)
    return State(state.points + tuple(additions), state.drawables, state.program)


def build_joint_state(policy="known-first"):
    replays = []
    for centers, radii in FIXTURES:
        replay = KpCenterLocusReplay("joint_root_kp", centers=centers, radii=radii, emit=False)
        report = replay.run()
        if report["all_targets"] != 49 or set(report["branches"].values()) != {"regular"}:
            raise AssertionError("校准夹具必须属于 49 E 正规分支")
        replays.append(replay)
    graph = replays[0].objects.graph
    paid_ids = tuple(graph.paid_order[:13]) + SHARED_IDS
    initial = tuple(PointBundle(f"O{i+1}", tuple(tuple(map(float, replay.centers[i]))
                    for replay in replays)) for i in range(3))
    initial_objects = tuple(DrawableBundle(node, "circle", tuple(object_value(replay.objects, node)[1]
                            for replay in replays)) for node in ("Gamma1", "Gamma2", "Gamma3"))
    state = State(initial, initial_objects)
    object_ids = {"Gamma1", "Gamma2", "Gamma3"}
    audit = BindingAudit()
    for move, node in enumerate(paid_ids, 1):
        references = graph.dependencies[node]
        # 使用重放的点身份；检查每个定位点确已存在于各样本交点闭包。
        for reference in references:
            values = tuple(point_value(replay.objects, reference) for replay in replays)
            if not any(all(close_point(value, point.values[i]) for i, value in enumerate(values))
                       for point in state.points):
                raise AssertionError(f"{node} 使用了尚未构造的 {reference}")
        rows = tuple(object_value(replay.objects, node) for replay in replays)
        drawable = DrawableBundle(node, rows[0][0], tuple(value for _, value in rows))
        candidate = Candidate(drawable.kind, *references, drawable)
        state = apply_calibrated_candidate(state, candidate, move, policy=policy, audit=audit)
        # 保留权威对象 ID，以核对免费命名点的父对象。
        renamed = DrawableBundle(node, drawable.kind, drawable.values)
        state = State(state.points, state.drawables[:-1] + (renamed,), state.program)
        object_ids.add(node)
        state = bind_named_intersections(state, replays, object_ids)
    targets = tuple(point_value(replay.objects, "Mannheim_S_center_locus") for replay in replays)
    kps = {}
    for profile in PROFILES:
        kp_id = graph.dependencies[f"{profile}_O3_Kp"][1]
        values = tuple(point_value(replay.objects, kp_id) for replay in replays)
        if not any(all(close_point(p, q) for p, q in zip(point.values, values, strict=True))
                   for point in state.points):
            raise AssertionError(f"{profile} 的 K' 尚未出现")
        kps[profile] = kp_id
    if any(all(close_point(p, q) for p, q in zip(point.values, targets, strict=True))
           for point in state.points):
        raise AssertionError("17 E 状态意外已有根心")
    if any(drawable_contains_targets(drawable, targets) for drawable in state.drawables):
        raise AssertionError("17 E 状态意外已有经过根心的对象")
    return state, targets, kps, audit


def fast_target_candidates(state, targets, required_ids, move):
    """NumPy 先筛首样本的必要关系，再逐项核对全部样本。

    首样本预筛容差较宽；最终核对沿用既有目标对象谓词。
    """

    points = state.points
    coordinates = np.asarray([point.values[0] for point in points])
    relative = coordinates - np.asarray(targets[0])
    required = np.asarray([point.point_id in required_ids for point in points])
    allowed = required[:, None] | required[None, :]
    lengths = np.linalg.norm(relative, axis=1)
    determinant = relative[:, None, 0] * relative[None, :, 1] - relative[:, None, 1] * relative[None, :, 0]
    line_mask = np.triu(allowed & (np.abs(determinant) <= 1e-7 * np.maximum(1.0, lengths[:, None] * lengths[None, :])), 1)
    delta = coordinates[:, None, :] - coordinates[None, :, :]
    through_squared = np.sum(delta * delta, axis=2)
    target_squared = np.sum(relative * relative, axis=1)[:, None]
    circle_mask = allowed & (np.abs(through_squared - target_squared) <= 1e-7 * np.maximum(1.0, np.maximum(through_squared, target_squared)))
    np.fill_diagonal(circle_mask, False)
    existing = {drawable_signature(drawable) for drawable in state.drawables}
    seen, result = set(), []
    for kind, mask in (("line", line_mask), ("circle", circle_mask)):
        for first_index, second_index in zip(*np.nonzero(mask), strict=True):
            first, second = points[first_index], points[second_index]
            predicate = collinear_with_targets if kind == "line" else equidistant_from_targets
            if not predicate(first, second, targets):
                continue
            operation = line_through if kind == "line" else circle_through
            values = tuple(operation(p, q) for p, q in zip(first.values, second.values, strict=True))
            if any(value is None for value in values):
                continue
            drawable = DrawableBundle(f"move_{move}", kind, values)
            signature = drawable_signature(drawable)
            if signature in seen or signature in existing:
                continue
            seen.add(signature)
            result.append(Candidate(kind, first.point_id, second.point_id, drawable))
    return tuple(result)


def third_target_candidates(base, two, targets, move):
    # 第三步也可能仅使用第一步产生的点；无需依赖第二步的新交点。
    required = frozenset(point.point_id for point in two.points[len(base.points):])
    return fast_target_candidates(two, targets, required, move)


def run_search(args):
    start = monotonic()
    state, targets, kps, audit = build_joint_state(args.policy)
    firsts = tuple(sorted(generate_candidates(state, 18), key=lambda candidate: candidate.describe()))
    fingerprint = hashlib.sha256("\n".join(candidate.describe() for candidate in firsts).encode()).hexdigest()
    positive_firsts = [first for first in firsts if drawable_contains_targets(first.drawable, targets)]
    report = {"schema": "joint-root-kp-screen/v1", "policy": args.policy, "samples": len(FIXTURES),
              "base_e": 17, "base_points": len(state.points), "base_objects": len(state.drawables),
              "kp_ids": kps, "candidate_fingerprint": fingerprint, "first_candidates": len(firsts),
              "first_target_objects": len(positive_firsts), "start_index": args.start_index,
              "checked_firsts": 0, "second_target_objects": 0, "third_target_objects": 0,
              "second_examples": [],
              "hits": [], "budget_seconds": args.seconds, "complete": False,
              "base_binding_audit": asdict(audit), "maximum_first_state_points": len(state.points)}
    if positive_firsts:
        raise RuntimeError("发现更短形状；先精确审计首步根心对象再扩展")
    print(json.dumps({"initial": report, "binding_audit": asdict(audit)}, ensure_ascii=False), flush=True)
    next_progress = monotonic() + 15
    for index in range(args.start_index, len(firsts)):
        if monotonic() - start >= args.seconds or (args.limit is not None and report["checked_firsts"] >= args.limit):
            break
        first = firsts[index]
        one = apply_calibrated_candidate(state, first, 18, policy=args.policy, audit=audit)
        report["maximum_first_state_points"] = max(report["maximum_first_state_points"], len(one.points))
        required = frozenset(point.point_id for point in one.points[len(state.points):])
        if required:
            seconds = fast_target_candidates(one, targets, required, 19)
            report["second_target_objects"] += len(seconds)
            for second in seconds:
                two = apply_calibrated_candidate(one, second, 19, policy=args.policy, audit=audit)
                thirds = third_target_candidates(state, two, targets, 20)
                report["third_target_objects"] += len(thirds)
                if len(report["second_examples"]) < 10:
                    report["second_examples"].append({"first_index": index, "first": first.describe(),
                                                       "second": second.describe(), "thirds": len(thirds)})
                for third in thirds:
                    roots = tuple(intersections(second.kind, p, third.kind, q)
                                  for p, q in zip(second.drawable.values, third.drawable.values, strict=True))
                    if all(any(close_point(point, target) for point in row)
                           for row, target in zip(roots, targets, strict=True)):
                        report["hits"].append([first.describe(), second.describe(), third.describe()])
        report["checked_firsts"] += 1
        if report["hits"]:
            break
        if monotonic() >= next_progress:
            print(json.dumps({"checked": report["checked_firsts"], "next_index": index + 1,
                              "seconds_through_S": report["second_target_objects"],
                              "elapsed_seconds": round(monotonic() - start, 3)}, ensure_ascii=False), flush=True)
            next_progress = monotonic() + 15
    report["next_index"] = args.start_index + report["checked_firsts"]
    report["covered_range"] = [args.start_index, report["next_index"]]
    # complete 沿用 v1 的“已到枚举尾部”含义；coverage_complete 才表示本次覆盖全部首步。
    report["complete"] = report["next_index"] == len(firsts)
    report["coverage_complete"] = args.start_index == 0 and report["complete"]
    report["elapsed_seconds"] = round(monotonic() - start, 3)
    report["binding_audit"] = asdict(audit)
    if args.output:
        Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", choices=("known-first", "all"), default="known-first")
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--start-index", type=int, default=0)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--output")
    args = parser.parse_args()
    if args.seconds <= 0 or args.start_index < 0 or (args.limit is not None and args.limit < 1):
        parser.error("时间预算、起始索引或数量限制不合法")
    run_search(args)


if __name__ == "__main__":
    main()
