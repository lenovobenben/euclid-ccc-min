"""隔离校准：排除已知交点后绑定新点，并可诊断全部跨样本根配对。

不改写旧搜索器；all 配对只供发现候选，不提供合法选根证书。
任何命中仍须以现行确定交点 E 模型精确重放。
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import dist

from search_parallel_3e import (
    Candidate, DrawableBundle, PointBundle, State, intersections, point_signature,
)


@dataclass
class BindingAudit:
    object_pairs: int = 0
    removed_known_roots: int = 0
    known_root_realignments: int = 0
    root_count_mismatches: int = 0
    mixed_root_bundles: int = 0


def close_point(first, second, tolerance=1e-8):
    return dist(first, second) <= tolerance * max(1.0, dist(first, (0, 0)), dist(second, (0, 0)))


def root_bundles(rows, known_points=(), *, policy="known-first", audit=None):
    """已知交点用已有点身份排除；其余根按排序或全部配对。

    不同样本实交点数量不同仍不合并；这样的候选需要分支覆盖验证。
    """

    if policy not in ("known-first", "all"):
        raise ValueError(policy)
    if audit is None:
        audit = BindingAudit()
    audit.object_pairs += 1
    remaining = [list(enumerate(row)) for row in rows]
    for known in known_points:
        matches = [
            next((index for index, point in row if close_point(point, value)), None)
            for row, value in zip(remaining, known.values, strict=True)
        ]
        if all(index is not None for index in matches):
            audit.removed_known_roots += 1
            audit.known_root_realignments += len(set(matches)) > 1
            for row, index in zip(remaining, matches, strict=True):
                row[:] = [(i, point) for i, point in row if i != index]
    counts = {len(row) for row in remaining}
    if len(counts) != 1:
        audit.root_count_mismatches += 1
        return ()
    if not remaining or not remaining[0]:
        return ()
    if policy == "all":
        combinations = product(*remaining)
    else:
        combinations = zip(*remaining, strict=True)
    result = []
    for combination in combinations:
        indices = tuple(item[0] for item in combination)
        values = tuple(item[1] for item in combination)
        audit.mixed_root_bundles += len(set(indices)) > 1
        result.append((indices, values))
    return tuple(result)


def apply_calibrated_candidate(state, candidate, move, *, policy="known-first", audit=None):
    new_drawable = DrawableBundle(f"object_{move}", candidate.kind, candidate.drawable.values)
    seen = {point_signature(point) for point in state.points}
    new_points = []
    for existing in state.drawables:
        rows = tuple(
            intersections(new_drawable.kind, new, existing.kind, old)
            for new, old in zip(new_drawable.values, existing.values, strict=True)
        )
        for indices, values in root_bundles(rows, state.points, policy=policy, audit=audit):
            suffix = ",".join(map(str, indices))
            point = PointBundle(f"object_{move}&{existing.drawable_id}[{suffix}]", values)
            signature = point_signature(point)
            if signature not in seen:
                seen.add(signature)
                new_points.append(point)
    return State(state.points + tuple(new_points), state.drawables + (new_drawable,),
                 state.program + (candidate.describe(),))
