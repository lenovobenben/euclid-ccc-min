"""Gergonne 平行四边形接口的有界入口检查；不提供 E 下界。

重放 13、17、25 E 真实状态，检查中心 C、根心 G 和四极点的已有载体
及一步可画载体。仅对 C 展开两步；没有旧载体时首步必须经过 C。
解析目标只用于检验，绝不作为定义点或对象加入状态。
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from time import monotonic

from audit_gergonne_polar_interface import similarity_center
from audit_mannheim_shared_projections import CORE_AFTER_SHARED
from check_mannheim_degenerate_fixture import line_through as exact_line
from mannheim_search_calibration import BindingAudit, apply_calibrated_candidate, close_point
from search_mannheim_center_locus_2e import point_value
from search_mannheim_joint_root_kp import (
    PROFILES, SHARED_IDS, build_carrier_state, fast_target_candidates,
)
from search_mannheim_root_center_2e import FIXTURES, drawable_contains_targets
from search_mannheim_shared_circle_replacement import available
from search_parallel_3e import circle_through, drawable_signature, intersections, line_through
from verify_mannheim_contact_chord_identity import PROFILES as SIGN_PROFILES


def pole_of_line(line, center, radius):
    a, b, c = line
    denominator = a * center[0] + b * center[1] + c
    if denominator == 0:
        raise ValueError("入口夹具不能使用无穷远极点")
    return (center[0] - radius**2 * a / denominator,
            center[1] - radius**2 * b / denominator)


def diagnostic_targets(replays):
    targets = {"C": [], **{profile: [] for profile in PROFILES}}
    for replay in replays:
        centers, radii = replay.centers, replay.radii
        origin, radius = centers[2], radii[2]
        targets["C"].append(pole_of_line(exact_line(*centers[:2]), origin, radius))
        for profile in PROFILES:
            _, sigma2, sigma3 = SIGN_PROFILES[profile]["sigma"]
            first = similarity_center(centers, radii, 0, sigma3)
            second = similarity_center(centers, radii, 1, sigma3 * sigma2)
            targets[profile].append(pole_of_line(exact_line(first, second), origin, radius))
    # 目标坐标是 oracle；不通过 PointBundle 注入构造状态。
    targets = {key: tuple(tuple(map(float, point)) for point in values)
               for key, values in targets.items()}
    targets["G"] = tuple(point_value(replay.objects, "Mannheim_S_center_locus")
                         for replay in replays)
    return targets


def check_legal_candidate(state, candidate):
    points = {point.point_id: point for point in state.points}
    first, second = points[candidate.first_id], points[candidate.second_id]
    operation = line_through if candidate.kind == "line" else circle_through
    values = tuple(operation(p, q) for p, q in zip(first.values, second.values, strict=True))
    if any(value is None for value in values):
        raise AssertionError("候选使用重合定义点")
    rebuilt = type(candidate.drawable)(candidate.drawable.drawable_id, candidate.kind, values)
    if drawable_signature(rebuilt) != drawable_signature(candidate.drawable):
        raise AssertionError("候选不由真实已有定义点产生")


def finite_target_intersection(first, second, targets):
    rows = tuple(intersections(first.kind, p, second.kind, q)
                 for p, q in zip(first.values, second.values, strict=True))
    return all(any(close_point(root, target) for root in roots)
               for roots, target in zip(rows, targets, strict=True))


def screen_state(base_e):
    started = monotonic()
    state, replays, audit = build_carrier_state(
        () if base_e == 13 else SHARED_IDS,
        extra_ids=CORE_AFTER_SHARED if base_e == 25 else (),
    )
    if len(state.program) != base_e:
        raise AssertionError("入口实际计费与状态标签不符")
    targets = diagnostic_targets(replays)
    rows, firsts_by_target = {}, {}
    for name, values in targets.items():
        is_available = available(state, values)
        old = [obj.drawable_id for obj in state.drawables if drawable_contains_targets(obj, values)]
        firsts = () if is_available else fast_target_candidates(
            state, values, frozenset(point.point_id for point in state.points), base_e + 1,
        )
        for first in firsts:
            check_legal_candidate(state, first)
        firsts_by_target[name] = firsts
        rows[name] = {"available": is_available, "old_carriers": old,
                      "first_screen_skipped": is_available,
                      "first_carriers": [first.describe() for first in firsts]}

    # 只展开 C。G、各极点若无旧载体且无首步载体，则无需展开两步。
    # 25 E 中极点已有接触弦载体，此时不能据零首步排除两步路径。
    c_values = targets["C"]
    if rows["C"]["available"] or rows["C"]["old_carriers"]:
        raise AssertionError("本轮两步剪枝要求 C 不可用且没有旧载体")
    expansions = []
    for first in firsts_by_target["C"]:
        extra_audit = BindingAudit()
        one = apply_calibrated_candidate(state, first, base_e + 1, audit=extra_audit)
        seconds = fast_target_candidates(
            one, c_values, frozenset(point.point_id for point in one.points), base_e + 2,
        )
        hits = []
        for second in seconds:
            check_legal_candidate(one, second)
            if finite_target_intersection(first.drawable, second.drawable, c_values):
                two = apply_calibrated_candidate(one, second, base_e + 2, audit=extra_audit)
                hits.append({"second": second.describe(), "available_in_binding": available(two, c_values)})
        expansions.append({"first": first.describe(), "points_after_first": len(one.points),
                           "available_after_first": available(one, c_values),
                           "second_carriers": [second.describe() for second in seconds],
                           "finite_intersection_hits": hits, "binding_audit": asdict(extra_audit)})
    return {"base_e": base_e, "points": len(state.points), "objects": len(state.drawables),
            "targets": rows, "c_two_step_expansions": expansions,
            "base_binding_audit": asdict(audit),
            "elapsed_seconds": round(monotonic() - started, 3)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    cases = [screen_state(base_e) for base_e in (13, 17, 25)]
    report = {"schema": "gergonne-parallelogram-entry/v1", "policy": "known-first",
              "samples": len(FIXTURES), "score_certificate": False, "lower_bound": False,
              "fixtures": [{"centers": [[str(value) for value in point] for point in centers],
                            "radii": [str(radius) for radius in radii]}
                           for centers, radii in FIXTURES],
              "model_scope_complete": True, "construction_coverage_complete": False,
              "scope": "static target carriers; C two-step completions only",
              "limitations": ["multi-fixture floating screening", "known-first root pairing only",
                              "unequal real-root counts skipped", "no dynamic G or pole search",
                              "no joint core/locus or reflection search"], "cases": cases}
    Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"states": [{"base_e": case["base_e"], "points": case["points"],
                                 "c_first_carriers": len(case["targets"]["C"]["first_carriers"]),
                                 "c_second_carriers": sum(len(row["second_carriers"])
                                                          for row in case["c_two_step_expansions"]),
                                 "c_hits": sum(len(row["finite_intersection_hits"])
                                               for row in case["c_two_step_expansions"]),
                                 "elapsed_seconds": case["elapsed_seconds"]}
                                for case in cases]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
