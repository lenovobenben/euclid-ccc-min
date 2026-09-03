"""搜索双 ``K'`` 平行方向类的共享定心圆能否在 2 E 内构造。

双平行截点恒等式给出：``P2/P3`` 的圆心线与 ``ell=O1O2`` 的交点
``J2,J3`` 关于一点 ``H`` 对称。若以 ``H`` 为圆心、过 ``J2,J3`` 的圆
能在 2 E 内画出，再连接 ``S J2``、``S J3``，两个方向类的四个目标只
需 12 E 后缀，比当前两套 7 E 后缀少 2 E。

脚本重放所有接触弦核心及免费交点。第一步枚举任意已有两点可画的线或
圆，检查它是否同时使圆心 ``H`` 和共享定心圆上的一个非圆心点可用；
若是，第二步即可画共享定心圆。``--merge-fixtures`` 切换到先行 ``P0``
简单合并的另一组 55 E 饱和夹具。

这是三个不同半径夹具上的同步浮点筛查。命中必须精确重放；零命中只
排除当前对象闭包与上述两步依赖形状。
"""

from __future__ import annotations

import argparse
from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    FIXTURES,
    MERGE_FIXTURES,
    build_state,
    candidate_has_witness,
    candidates,
    on_drawable,
    same_drawable,
    same_point,
)
from search_parallel_3e import intersections


TOLERANCE = 1e-7


def line_intersection(first, second):
    roots = intersections("line", first, "line", second)
    if len(roots) != 1:
        raise AssertionError("目标圆心线没有与 ell 产生唯一有限交点")
    return roots[0]


def target_circle_data(state, target_lines):
    ell = state.drawables[3]
    if ell.kind != "line":
        raise AssertionError("第一个计费对象不是 ell")
    centers = []
    circles = []
    j_rows = []
    for sample_index, ell_value in enumerate(ell.values):
        j2 = line_intersection(
            target_lines["P2"][sample_index], ell_value
        )
        j3 = line_intersection(
            target_lines["P3"][sample_index], ell_value
        )
        center = ((j2[0] + j3[0]) / 2, (j2[1] + j3[1]) / 2)
        radius_squared = (
            (j2[0] - j3[0]) ** 2 + (j2[1] - j3[1]) ** 2
        ) / 4
        if radius_squared <= TOLERANCE:
            raise AssertionError("双平行截点没有严格分离")
        centers.append(center)
        circles.append((center[0], center[1], radius_squared))
        j_rows.append((j2, j3))
    return tuple(centers), tuple(circles), tuple(j_rows)


def point_on_circle(point, circle) -> bool:
    center_x, center_y, radius_squared = circle
    residual = abs(
        (point[0] - center_x) ** 2
        + (point[1] - center_y) ** 2
        - radius_squared
    )
    return residual <= TOLERANCE * max(1.0, radius_squared)


def build_circle_witnesses(state, target_circles):
    witnesses = []
    for drawable in state.drawables:
        sample_rows = []
        for circle, old_value in zip(
            target_circles,
            drawable.values,
            strict=True,
        ):
            roots = intersections(
                "circle",
                circle,
                drawable.kind,
                old_value,
            )
            if not roots:
                break
            sample_rows.append(roots)
        else:
            witnesses.append(
                (
                    drawable.drawable_id,
                    drawable.kind,
                    drawable.values,
                    tuple(sample_rows),
                )
            )
    return tuple(witnesses)


def candidate_yields_center(candidate, centers, support_drawables) -> bool:
    if not all(
        on_drawable(center, candidate.kind, value)
        for center, value in zip(
            centers,
            candidate.drawable.values,
            strict=True,
        )
    ):
        return False
    return any(
        not (
            candidate.kind == old.kind
            and all(
                same_drawable(candidate.kind, value, old_value)
                for value, old_value in zip(
                    candidate.drawable.values,
                    old.values,
                    strict=True,
                )
            )
        )
        for old in support_drawables
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merge-fixtures", action="store_true")
    parser.add_argument("--include-finite-targets", action="store_true")
    args = parser.parse_args()
    fixtures = MERGE_FIXTURES if args.merge_fixtures else FIXTURES

    start = monotonic()
    state, target_lines, _, paid_sequence = build_state(
        fixtures=fixtures,
        include_finite_targets=args.include_finite_targets,
        merge_branch=args.merge_fixtures,
    )
    centers, target_circles, j_rows = target_circle_data(
        state,
        target_lines,
    )
    existing_center_points = tuple(
        point.point_id
        for point in state.points
        if all(
            same_point(value, center)
            for value, center in zip(
                point.values,
                centers,
                strict=True,
            )
        )
    )
    existing_circle_points = tuple(
        point.point_id
        for point in state.points
        if all(
            point_on_circle(value, circle)
            for value, circle in zip(
                point.values,
                target_circles,
                strict=True,
            )
        )
    )
    center_support = tuple(
        drawable
        for drawable in state.drawables
        if all(
            on_drawable(center, drawable.kind, value)
            for center, value in zip(
                centers,
                drawable.values,
                strict=True,
            )
        )
    )
    circle_witnesses = build_circle_witnesses(state, target_circles)

    hits = []
    checked = 0
    for candidate in candidates(state):
        checked += 1
        yields_center = candidate_yields_center(
            candidate,
            centers,
            center_support,
        )
        yields_circle_point = candidate_has_witness(
            candidate,
            circle_witnesses,
        )
        if (
            (existing_center_points or yields_center)
            and (existing_circle_points or yields_circle_point)
        ):
            hits.append(
                {
                    "candidate": candidate.describe(),
                    "yields_center": yields_center,
                    "circle_point_witness": yields_circle_point,
                }
            )
        if checked % 100_000 == 0:
            print(
                "progress",
                {
                    "checked": checked,
                    "hits": len(hits),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    print(
        "mannheim_double_kp_center_circle_2e_search",
        {
            "samples": len(fixtures),
            "merge_fixtures": args.merge_fixtures,
            "include_finite_targets": args.include_finite_targets,
            "paid_objects": len(paid_sequence),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "existing_center_points": existing_center_points,
            "existing_circle_points": existing_circle_points,
            "center_support_objects": len(center_support),
            "circle_witness_objects": len(circle_witnesses),
            "candidates": checked,
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    print("target_intercepts", j_rows)
    for hit in hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
