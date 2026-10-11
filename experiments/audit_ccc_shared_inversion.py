"""CCC 判别式乘积与公共反演接口；代数共享不等于 E 节省。

符号恒等式用整数多项式展开，配对圆像用实二次域精确核对。
圆的可画性另在 known-first 多夹具浮点闭包中筛查；不注入解析点。
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import asdict
from itertools import combinations
from pathlib import Path
from time import monotonic

from audit_mannheim_shared_projections import CORE_AFTER_SHARED
from mannheim_search_calibration import BindingAudit, apply_calibrated_candidate, close_point
from mannheim_virtual_witnesses import coincides
from search_mannheim_joint_root_kp import PROFILES, build_carrier_state, fast_target_candidates
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_root_center_2e import drawable_contains_targets
from search_parallel_3e import Candidate, DrawableBundle, circle_through
from verify_mannheim_contact_chord_identity import (
    PROFILES as SIGN_PROFILES, c, d, h, p, q,
)
from verify_mannheim_tau_concurrence import MANNHEIM_POINT_S


def dot(first, second):
    return sum(a * b for a, b in zip(first, second, strict=True))


def verify_discriminants():
    gx, gy, w = MANNHEIM_POINT_S
    v = d * h
    k = gx**2 + gy**2 - w**2
    factors = {
        "12_same": d**2 - (p - q)**2,
        "12_opposite": d**2 - (p + q)**2,
        "13_same": c**2 + h**2 - (p - 1)**2,
        "13_opposite": c**2 + h**2 - (p + 1)**2,
        "23_same": (c + d)**2 + h**2 - (q - 1)**2,
        "23_opposite": (c + d)**2 + h**2 - (q + 1)**2,
    }
    occurrences, rows = Counter(), []
    for profile in PROFILES:
        _, sigma2, sigma3 = SIGN_PROFILES[profile]["sigma"]
        ux = h * (p - sigma2 * q)
        uy = -c * (p - sigma2 * q) - d * (p - sigma3)
        a = ux**2 + uy**2 - v**2
        b = gx * ux + gy * uy - sigma3 * w * v
        keys = ("12_" + ("same" if sigma2 == 1 else "opposite"),
                "13_" + ("same" if sigma3 == 1 else "opposite"),
                "23_" + ("same" if sigma2 == sigma3 else "opposite"))
        product = factors[keys[0]] * factors[keys[1]] * factors[keys[2]]
        # B^2-A*K = F12*F13*F23/(4*d^2*h^2)。清分母逐项比较。
        if not (b**2 - a * k - v**2 * product).is_zero:
            raise AssertionError(f"{profile} 判别式因式分解失败")
        occurrences.update(keys)
        rows.append({"profile": profile, "factors": keys, "identity": True})
    # 已逐类精确证明因式分解；每个因子恰出现两次即可证明乘积恒等式。
    if occurrences != Counter({key: 2 for key in factors}):
        raise AssertionError("四类根式乘积不是六因子乘积的平方")
    return {"profiles": rows, "four_radicand_product_is_square": True,
            "factor_occurrences": dict(occurrences),
            "ring": "Z[c,d,h,p,q]", "score_certificate": False}


def exact_radical_center(centers, radii):
    rows = []
    for index in (1, 2):
        delta = tuple(a - b for a, b in zip(centers[index], centers[0], strict=True))
        constant = dot(centers[index], centers[index]) - dot(centers[0], centers[0])
        constant += radii[0]**2 - radii[index]**2
        rows.append((2 * delta[0], 2 * delta[1], constant))
    a, b, e = rows[0]
    c0, d0, f = rows[1]
    determinant = a * d0 - b * c0
    return ((e * d0 - b * f) / determinant, (a * f - e * c0) / determinant)


def verify_inversion(replays):
    fixtures, omega_values = [], []
    for replay in replays:
        g = exact_radical_center(replay.centers, replay.radii)
        delta = tuple(a - b for a, b in zip(g, replay.centers[2], strict=True))
        power = dot(delta, delta) - replay.radii[2]**2
        if power <= 0:
            raise AssertionError("互相离圆夹具的公共幂必须为正")
        for center, radius in zip(replay.centers, replay.radii, strict=True):
            relative = tuple(a - b for a, b in zip(center, g, strict=True))
            if dot(relative, relative) - radius**2 != power:
                raise AssertionError("输入圆不是公共反演的不变圆")
        checks = []
        for profile in PROFILES:
            signs = SIGN_PROFILES[profile]["sigma"]
            product = 1
            for first, second in combinations(range(3), 2):
                edge = tuple(a - b for a, b in zip(replay.centers[first], replay.centers[second], strict=True))
                product *= dot(edge, edge) - (signs[first] * replay.radii[first] - signs[second] * replay.radii[second])**2
            edge1 = tuple(a - b for a, b in zip(replay.centers[1], replay.centers[0], strict=True))
            edge2 = tuple(a - b for a, b in zip(replay.centers[2], replay.centers[0], strict=True))
            area = edge1[0] * edge2[1] - edge1[1] * edge2[0]
            discriminant_quarter = product / (4 * area**2)
            pair = [target["circle"] for target in replay.targets.values() if target["profile"] == profile]
            if len(pair) != 2:
                raise AssertionError("正规方向类须恰有两个目标圆")
            for source, target in ((pair[0], pair[1]), (pair[1], pair[0])):
                relative = tuple(a - b for a, b in zip(source[0], g, strict=True))
                denominator = dot(relative, relative) - source[1]
                inverse_center = tuple(origin + power * value / denominator
                                       for origin, value in zip(g, relative, strict=True))
                inverse_radius_squared = power**2 * source[1] / (denominator * denominator)
                if inverse_center != target[0] or inverse_radius_squared != target[1]:
                    raise AssertionError(f"{profile} 的公共反演圆像不匹配另一解")
            checks.append({"profile": profile, "both_directions_exact": True,
                           "discriminant_quarter": str(discriminant_quarter),
                           "pair_intersection_margin": str(power - discriminant_quarter),
                           "pair_has_two_real_intersections": discriminant_quarter < power})
        fixtures.append({"centers": [[str(value) for value in point] for point in replay.centers],
                         "radii": [str(value) for value in replay.radii],
                         "radical_center": [str(value) for value in g],
                         "inversion_power": str(power), "pairs": checks})
        omega_values.append((*map(float, g), float(power)))
    return fixtures, tuple(omega_values)


def screen_omega(selected, suffix_ids, omega):
    started = monotonic()
    planned = tuple(node for profile in selected for node in suffix_ids[profile])
    state, replays, audit = build_carrier_state(extra_ids=CORE_AFTER_SHARED + planned)
    if len(state.program) != 25 + 6 * len(selected):
        raise AssertionError("实际前缀费用不符合完成目标对的状态")
    center = next(point for point in state.points
                  if all(close_point(value, row[:2]) for value, row in zip(point.values, omega.values, strict=True)))
    candidates = []
    for point in state.points:
        values = tuple(circle_through(g, through) for g, through in zip(center.values, point.values, strict=True))
        if any(value is None for value in values):
            continue
        drawable = DrawableBundle("omega_candidate", "circle", values)
        if coincides(drawable, omega):
            candidates.append(Candidate("circle", center.point_id, point.point_id, drawable))
    completion = None
    if candidates:
        extra_audit = BindingAudit()
        # 必须使用真实已有点定义的圆，不能把 virtual_omega 加入状态。
        one = apply_calibrated_candidate(state, candidates[0], len(state.program) + 1, audit=extra_audit)
        new_points = one.points[len(state.points):]
        loci, centers, shared_carriers = {}, {}, {}
        for profile in PROFILES:
            if profile in selected:
                continue
            lines = tuple(object_value(replay.objects, f"{profile}_center_locus_from_Kp")[1]
                          for replay in replays)
            loci[profile] = [point.point_id for point in new_points
                             if not all(close_point(p, g) for p, g in zip(point.values, center.values, strict=True))
                             and all(abs(a * x + b * y + c0) <= 1e-8 * max(1.0, abs(x), abs(y), abs(c0))
                                     for (a, b, c0), (x, y) in zip(lines, point.values, strict=True))]
            target_keys = [key for key, target in replays[0].targets.items() if target["profile"] == profile]
            center_values = []
            for key in target_keys:
                values = tuple(point_value(replay.objects, replay.objects.graph.dependencies[
                    replay.targets[key]["output_id"]][0]) for replay in replays)
                centers[key] = [point.point_id for point in new_points
                                if all(close_point(p, target) for p, target in zip(point.values, values, strict=True))]
                center_values.append(values)
            firsts = fast_target_candidates(
                one, center_values[0], frozenset(point.point_id for point in one.points), len(one.program) + 1,
            )
            shared_carriers[profile] = [candidate.describe() for candidate in firsts
                                       if drawable_contains_targets(candidate.drawable, center_values[1])]
        completion = {"draw": candidates[0].describe(), "new_points": len(new_points),
                      "new_points_on_remaining_loci": loci, "new_remaining_centers": centers,
                      "static_shared_center_carriers": shared_carriers,
                      "binding_audit": asdict(extra_audit)}
    return {"completed_profiles": selected, "base_e": len(state.program),
            "points": len(state.points), "objects": len(state.drawables),
            "omega_existing": [obj.drawable_id for obj in state.drawables if coincides(obj, omega)],
            "omega_drawable_through_points": [candidate.describe() for candidate in candidates],
            "omega_completion": completion, "binding_audit": asdict(audit),
            "elapsed_seconds": round(monotonic() - started, 3)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    started = monotonic()
    symbolic = verify_discriminants()
    _, replays, _ = build_carrier_state(extra_ids=CORE_AFTER_SHARED)
    fixtures, omega_values = verify_inversion(replays)
    omega = DrawableBundle("virtual_omega", "circle", omega_values)
    suffix_ids = {profile: tuple(node for node in replays[0].objects.graph.paid_order[25:]
                                if node.startswith(profile + "_") or node.startswith("target_" + profile + ":"))
                  for profile in PROFILES}
    if any(len(nodes) != 6 for nodes in suffix_ids.values()):
        raise AssertionError("正规后缀必须每方向类六个对象")
    states = [screen_omega(selected, suffix_ids, omega)
              for count in range(4) for selected in combinations(PROFILES, count)]
    report = {"schema": "ccc-shared-inversion/v1", "symbolic": symbolic,
              "exact_inversion_fixtures": fixtures, "omega_entry_states": states,
              "policy": "known-first", "virtual_omega_is_not_an_input": True,
              "score_certificate": False, "construction_coverage_complete": False,
              "scope": "15 completed-pair states; static omega entry; one omega draw; free loci and centers; static shared-center carriers",
              "limitations": ["known-first floating closure skips unequal root counts",
                              "no dynamic inverse-circle construction",
                              "no general-parameter E coverage certificate"],
              "elapsed_seconds": round(monotonic() - started, 3)}
    Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"symbolic": symbolic, "exact_pairs": len(fixtures) * len(PROFILES),
                      "states": [{"completed": row["completed_profiles"],
                                  "omega_through_points": len(row["omega_drawable_through_points"]),
                                  "remaining_locus_hits": ([] if row["omega_completion"] is None else
                                      [key for key, values in row["omega_completion"]["new_points_on_remaining_loci"].items() if values]),
                                  "shared_center_hits": (0 if row["omega_completion"] is None else
                                      sum(len(values) for values in row["omega_completion"]["static_shared_center_carriers"].values()))}
                                 for row in states], "elapsed_seconds": report["elapsed_seconds"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
