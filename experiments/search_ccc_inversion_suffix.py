"""公共反演圆后的有界 5 E 双目标后缀筛查。

早期单圆心载体、两半径后共同载体、以及首输出圆先行均显式检查。
非圆心首对象的延迟共同载体复用已有见证剪枝。不是一般参数下界。
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from time import monotonic

from audit_ccc_shared_inversion import omega_candidates, verify_inversion
from audit_mannheim_shared_projections import CORE_AFTER_SHARED
from mannheim_search_calibration import BindingAudit, apply_calibrated_candidate, close_point
from search_mannheim_center_locus_2e import point_value
from search_mannheim_delayed_center_pair import pair_carriers, pair_data, run_case
from search_mannheim_joint_root_kp import PROFILES, build_carrier_state, fast_target_candidates
from search_mannheim_root_center_2e import drawable_contains_targets
from search_mannheim_shared_circle_replacement import available
from search_parallel_3e import Candidate, DrawableBundle, circle_through


def matching_point(state, values):
    return next((point for point in state.points
                 if all(close_point(p, q) for p, q in zip(point.values, values, strict=True))), None)


def draw_output(state, center_values, contact_values, audit):
    center, contact = matching_point(state, center_values), matching_point(state, contact_values)
    if center is None or contact is None:
        return None
    values = tuple(circle_through(p, q) for p, q in zip(center.values, contact.values, strict=True))
    if any(value is None for value in values):
        raise AssertionError("目标输出圆不能是零半径圆")
    candidate = Candidate("circle", center.point_id, contact.point_id,
                          DrawableBundle("output", "circle", values))
    return apply_calibrated_candidate(state, candidate, len(state.program) + 1, audit=audit)


def finish_pair(state, targets, contacts, audit):
    for center, contact in zip(targets, contacts, strict=True):
        state = draw_output(state, center, contact, audit)
        if state is None:
            return None
    return state


def early_case(state, replays, profile, deadline):
    keys, targets = pair_data(state, replays, profile)
    contacts = tuple(tuple(point_value(replay.objects, replay.objects.graph.dependencies[
        replay.targets[f"{profile}:{key}"]["output_id"]][1]) for replay in replays) for key in keys)
    return early_target_search(state, keys, targets, contacts, deadline)


def early_target_search(state, keys, targets, contacts, deadline):
    """与 CCC 参数提取分开，便于用已知几何阳性校准时间顺序。"""
    started = monotonic()
    if any(drawable_contains_targets(obj, target) for obj in state.drawables for target in targets):
        raise AssertionError("早期载体子类要求初态无圆心旧载体")
    if pair_carriers(state, *targets, len(state.program) + 1):
        raise AssertionError("已有一步共同载体，应先审计更短接口")
    audit = BindingAudit()
    report = {"target_keys": keys, "first_carriers": [], "checked_firsts": 0,
              "second_carriers": 0, "second_shared_carriers": 0, "third_shared_carriers": 0,
              "first_center_available_after_two": 0, "output_first_branches": 0,
              "output_first_third_carriers": 0, "three_auxiliary_hits": 0,
              "output_first_hits": 0, "examples": [], "hits": [], "complete": True}
    ids = frozenset(point.point_id for point in state.points)
    for index in (0, 1):
        firsts = fast_target_candidates(state, targets[index], ids, len(state.program) + 1)
        report["first_carriers"].append({"key": keys[index], "candidates": [c.describe() for c in firsts]})
        for first in firsts:
            if monotonic() >= deadline:
                report["complete"] = False
                break
            report["checked_firsts"] += 1
            one = apply_calibrated_candidate(state, first, len(state.program) + 1, audit=audit)
            other = 1 - index
            seconds = fast_target_candidates(one, targets[other], frozenset(p.point_id for p in one.points), len(one.program) + 1)
            report["second_carriers"] += len(seconds)
            for second in seconds:
                if monotonic() >= deadline:
                    report["complete"] = False
                    break
                shared = drawable_contains_targets(second.drawable, targets[index])
                report["second_shared_carriers"] += shared
                two = apply_calibrated_candidate(one, second, len(one.program) + 1, audit=audit)
                first_ready = available(two, targets[index])
                report["first_center_available_after_two"] += first_ready
                # 未定第一圆心时，第三对象必须同时经过两圆心。
                thirds = (fast_target_candidates(two, targets[other], frozenset(p.point_id for p in two.points), len(two.program) + 1)
                          if first_ready else pair_carriers(two, *targets, len(two.program) + 1))
                report["third_shared_carriers"] += sum(drawable_contains_targets(c.drawable, targets[index]) for c in thirds)
                if len(report["examples"]) < 8:
                    report["examples"].append({"first": first.describe(), "second": second.describe(), "third_count": len(thirds)})
                for third in thirds:
                    if monotonic() >= deadline:
                        report["complete"] = False
                        break
                    three = apply_calibrated_candidate(two, third, len(two.program) + 1, audit=audit)
                    done = finish_pair(three, targets, contacts, audit)
                    if done is not None:
                        report["three_auxiliary_hits"] += 1
                        report["hits"].append(list(done.program[len(state.program):]))
                # 输出圆先行：只有第一圆心在第二辅助对象后真实可用才画它。
                if first_ready:
                    output = draw_output(two, targets[index], contacts[index], audit)
                    report["output_first_branches"] += 1
                    thirds = fast_target_candidates(output, targets[other], frozenset(p.point_id for p in output.points), len(output.program) + 1)
                    report["output_first_third_carriers"] += len(thirds)
                    for third in thirds:
                        if monotonic() >= deadline:
                            report["complete"] = False
                            break
                        three = apply_calibrated_candidate(output, third, len(output.program) + 1, audit=audit)
                        done = draw_output(three, targets[other], contacts[other], audit)
                        if done is not None:
                            report["output_first_hits"] += 1
                            report["hits"].append(list(done.program[len(state.program):]))
            if report["hits"] or not report["complete"]:
                break
        if report["hits"] or not report["complete"]:
            break
    if any(len(program) != 5 for program in report["hits"]):
        raise AssertionError("命中须展开为恰好五个实际绘制动作")
    if report["hits"]:
        report["complete"] = False  # 阳性后停止，不把尚未展开的首步计为完成。
    report["binding_audit"] = asdict(audit)
    report["elapsed_seconds"] = round(monotonic() - started, 3)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds-per-case", type=float, default=15)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.seconds_per_case <= 0:
        parser.error("每例预算必须为正")
    _, replays, _ = build_carrier_state(extra_ids=CORE_AFTER_SHARED)
    fixtures, values = verify_inversion(replays)
    omega = DrawableBundle("virtual_omega", "circle", values)
    suffix_ids = {profile: tuple(node for node in replays[0].objects.graph.paid_order[25:]
                                if node.startswith(profile + "_") or node.startswith("target_" + profile + ":"))
                  for profile in PROFILES}
    report = {"schema": "inversion-suffix-screen/v1", "policy": "known-first",
              "entry_domain": "selected seed pair has two finite real intersections",
              "score_certificate": False, "construction_coverage_complete": False,
              "scope": "three one-pair seed states; early-center firsts and delayed shared-carrier helpers; five-draw suffixes only",
              "limitations": ["known-first floating closure skips unequal real-root counts",
                              "virtual witness roots use all pairings; virtual points are not inputs",
                              "no cross-profile shared auxiliary programs",
                              "no changed contact core or changed seed-pair program",
                              "no general-parameter E certificate"],
              "fixtures": fixtures, "cases": [], "seconds_per_case": args.seconds_per_case}
    for seed in ("P2", "P1", "P3"):
        base, rows, base_audit = build_carrier_state(extra_ids=CORE_AFTER_SHARED + suffix_ids[seed])
        _, candidates = omega_candidates(base, omega)
        if not candidates:
            raise AssertionError("所选种子在多夹具闭包中必须有公共圆周点")
        entry_audit = BindingAudit()
        state = apply_calibrated_candidate(base, candidates[0], len(base.program) + 1, audit=entry_audit)
        if len(state.program) != 32:
            raise AssertionError("种子对完成并画公共圆后应为 32 E")
        for profile in PROFILES:
            if profile == seed:
                continue
            started = monotonic()
            deadline = started + args.seconds_per_case
            early = early_case(state, rows, profile, deadline)
            remaining = deadline - monotonic()
            delayed = None if remaining <= 0 or early["hits"] else run_case(state, rows, profile, remaining)
            case = {"seed": seed, "profile": profile, "base_e": 32, "base_points": len(state.points),
                    "base_binding_audit": asdict(base_audit), "omega_entry": candidates[0].describe(),
                    "entry_binding_audit": asdict(entry_audit), "early": early, "delayed": delayed,
                    "complete": early["complete"] and delayed is not None and delayed["complete"],
                    "elapsed_seconds": round(monotonic() - started, 3)}
            report["cases"].append(case)
            Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({"seed": seed, "profile": profile, "early_firsts": early["checked_firsts"],
                              "early_second_shared": early["second_shared_carriers"],
                              "early_third_shared": early["third_shared_carriers"],
                              "delayed_helpers": None if delayed is None else delayed["checked_helpers"],
                              "hits": len(early["hits"]) + (0 if delayed is None else len(delayed["hits"])),
                              "complete": case["complete"], "elapsed_seconds": case["elapsed_seconds"]}, ensure_ascii=False), flush=True)
            if early["hits"] or (delayed is not None and delayed["hits"]):
                return


if __name__ == "__main__":
    main()
