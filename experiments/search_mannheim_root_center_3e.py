"""在 Mannheim 13 E 公共前缀后筛查三步根心构造。

``search_mannheim_root_center_2e.py`` 已确认：初始对象和所有可直接画出的
首步对象都不经过三个样本的根心。故三步程序若存在，第二步必须使用
第一步新交点画出第一条经过根心的对象，第三步再画另一条经过根心的
对象。这里按这个必要形状枚举，并只生成至少使用一个新增点的目标对象。

搜索使用三个严格正规 ``D8`` 夹具的浮点束。命中需要精确重放，零命中
不是三步不可能性的形式证明。
"""

from __future__ import annotations

import argparse
from itertools import combinations
from math import hypot
from time import monotonic

from search_mannheim_root_center_2e import (
    build_state,
    drawable_contains_targets,
)
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    State,
    apply_candidate,
    circle_through,
    drawable_signature,
    generate_candidates,
    line_through,
)


def collinear_with_targets(first, second, targets) -> bool:
    for p, q, target in zip(first.values, second.values, targets, strict=True):
        first_x = p[0] - target[0]
        first_y = p[1] - target[1]
        second_x = q[0] - target[0]
        second_y = q[1] - target[1]
        determinant = first_x * second_y - first_y * second_x
        scale = max(
            1.0,
            hypot(first_x, first_y) * hypot(second_x, second_y),
        )
        if abs(determinant) > 1e-8 * scale:
            return False
    return True


def equidistant_from_targets(center, through, targets) -> bool:
    for c, p, target in zip(center.values, through.values, targets, strict=True):
        through_squared = (p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2
        target_squared = (target[0] - c[0]) ** 2 + (target[1] - c[1]) ** 2
        scale = max(1.0, through_squared, target_squared)
        if abs(through_squared - target_squared) > 1e-8 * scale:
            return False
    return True


def target_candidates(
    state: State,
    targets,
    required_ids: frozenset[str],
    move: int,
):
    existing = {drawable_signature(item) for item in state.drawables}
    seen = set()
    candidates = []
    for first, second in combinations(state.points, 2):
        if (
            first.point_id not in required_ids
            and second.point_id not in required_ids
        ):
            continue
        if not collinear_with_targets(first, second, targets):
            continue
        values = tuple(
            line_through(first_value, second_value)
            for first_value, second_value in zip(
                first.values, second.values, strict=True
            )
        )
        if any(value is None for value in values):
            continue
        drawable = DrawableBundle(f"move_{move}", "line", values)
        signature = drawable_signature(drawable)
        if signature in existing or signature in seen:
            continue
        seen.add(signature)
        candidates.append(
            Candidate("line", first.point_id, second.point_id, drawable)
        )

    for center in state.points:
        for through in state.points:
            if center.point_id == through.point_id:
                continue
            if (
                center.point_id not in required_ids
                and through.point_id not in required_ids
            ):
                continue
            if not equidistant_from_targets(center, through, targets):
                continue
            values = tuple(
                circle_through(center_value, through_value)
                for center_value, through_value in zip(
                    center.values, through.values, strict=True
                )
            )
            if any(value is None for value in values):
                continue
            drawable = DrawableBundle(f"move_{move}", "circle", values)
            signature = drawable_signature(drawable)
            if signature in existing or signature in seen:
                continue
            seen.add(signature)
            candidates.append(
                Candidate("circle", center.point_id, through.point_id, drawable)
            )
    if any(
        not drawable_contains_targets(candidate.drawable, targets)
        for candidate in candidates
    ):
        raise AssertionError("目标对象的快速必要条件产生了伪命中")
    return tuple(candidates)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--first-kind",
        choices=("line", "circle", "both"),
        default="both",
    )
    args = parser.parse_args()

    initial, targets = build_state()
    initial_point_count = len(initial.points)
    first_candidates = generate_candidates(initial, 14)
    if args.first_kind != "both":
        first_candidates = tuple(
            candidate
            for candidate in first_candidates
            if candidate.kind == args.first_kind
        )

    start = monotonic()
    second_target_candidates = 0
    third_target_candidates = 0
    firsts_with_target_second = 0
    hits = []
    for first_index, first in enumerate(first_candidates, start=1):
        state_one = apply_candidate(initial, first, 14)
        first_new_ids = frozenset(
            point.point_id for point in state_one.points[initial_point_count:]
        )
        if not first_new_ids:
            continue
        seconds = target_candidates(
            state_one,
            targets,
            first_new_ids,
            15,
        )
        if seconds:
            firsts_with_target_second += 1
        second_target_candidates += len(seconds)
        for second in seconds:
            state_two = apply_candidate(state_one, second, 15)
            all_new_ids = frozenset(
                point.point_id for point in state_two.points[initial_point_count:]
            )
            thirds = target_candidates(
                state_two,
                targets,
                all_new_ids,
                16,
            )
            third_target_candidates += len(thirds)
            for third in thirds:
                hits.append(
                    (first.describe(), second.describe(), third.describe())
                )
        if first_index % 1000 == 0:
            print(
                "progress",
                {
                    "checked": first_index,
                    "total": len(first_candidates),
                    "firsts_with_target_second": firsts_with_target_second,
                    "second_target_candidates": second_target_candidates,
                    "hits": len(hits),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    print(
        "root_center_3e_search",
        {
            "samples": len(targets),
            "first_kind": args.first_kind,
            "initial_points": initial_point_count,
            "initial_drawables": len(initial.drawables),
            "first_candidates": len(first_candidates),
            "firsts_with_target_second": firsts_with_target_second,
            "second_target_candidates": second_target_candidates,
            "third_target_candidates": third_target_candidates,
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for hit in hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
