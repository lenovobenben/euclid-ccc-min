"""未绘制轨迹的有限见证剪枝；虚拟点只作目标，不成为可用点。"""

from __future__ import annotations

from dataclasses import asdict
from math import hypot
from time import monotonic

from mannheim_search_calibration import BindingAudit, root_bundles
from search_mannheim_joint_root_kp import fast_target_candidates
from search_parallel_3e import drawable_signature, generate_candidates, intersections, point_signature, PointBundle


def bisector(first, second):
    dx, dy = second[0] - first[0], second[1] - first[1]
    norm = hypot(dx, dy)
    if norm <= 1e-9:
        raise ValueError("目标点必须不同")
    a, b = dx / norm, dy / norm
    return a, b, -a * (first[0] + second[0]) / 2 - b * (first[1] + second[1]) / 2


def coincides(first, second):
    if first.kind != second.kind:
        return False
    for p, q in zip(first.values, second.values, strict=True):
        residual = max(abs(a - b) for a, b in zip(p, q))
        if first.kind == "line":
            residual = min(residual, max(abs(a + b) for a, b in zip(p, q)))
        scale = max(1.0, *(abs(value) for value in p), *(abs(value) for value in q))
        if residual > 1e-8 * scale:
            return False
    return True


def witness_helpers(state, goals, move, deadline):
    """必要见证：新定义点在虚拟轨迹与某旧对象的有限交集。

    虚拟见证剩余根用全部配对；首对象仍只由 state 中的实际点定义。
    轨迹与旧对象重合时不使用有限剪枝，而保留全部首对象。
    """
    audit = BindingAudit()
    helpers, witness_signatures = {}, set()
    point_ids = frozenset(point.point_id for point in state.points)
    pairs_checked, coincident_pairs = 0, 0
    complete = True
    for goal in goals:
        for old in state.drawables:
            if monotonic() >= deadline:
                complete = False
                break
            pairs_checked += 1
            if coincides(goal, old):
                coincident_pairs += 1
                for candidate in generate_candidates(state, move):
                    helpers.setdefault(drawable_signature(candidate.drawable), candidate)
                continue
            rows = tuple(intersections(goal.kind, p, old.kind, q)
                         for p, q in zip(goal.values, old.values, strict=True))
            for _, values in root_bundles(rows, state.points, policy="all", audit=audit):
                signature = point_signature(PointBundle("virtual", values))
                if signature in witness_signatures:
                    continue
                witness_signatures.add(signature)
                for candidate in fast_target_candidates(state, values, point_ids, move):
                    helpers.setdefault(drawable_signature(candidate.drawable), candidate)
        if not complete:
            break
    report = {"goal_count": len(goals), "goal_old_pairs": len(goals) * len(state.drawables),
              "checked_goal_old_pairs": pairs_checked, "witnesses": len(witness_signatures),
              "helper_candidates": len(helpers), "complete": complete,
              "witness_binding_audit": asdict(audit), "coincident_goal_old_pairs": coincident_pairs}
    return tuple(sorted(helpers.values(), key=lambda candidate: candidate.describe())), report
