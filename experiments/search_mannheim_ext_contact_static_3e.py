"""筛查正规 Mannheim 单根的静态 3 E 直接接触点程序。

初态是合法 5 E 平行前缀与 ``P0`` 四条批量线组成的 9 E 状态。目标
接触点 ``M3`` 已在输入圆 ``Gamma3`` 上，因此只须再画一个经过它的
直线或圆。这里限定前两步为两个都能从初态直接画出的不同对象；它们
产生交点 ``Q``，第三步使用一个初态点 ``P`` 和 ``Q`` 画目标轨迹。

逐个 ``P`` 检查三种完整的定义方向：``Line(P,Q)``、以 ``P`` 为圆心
经过 ``Q`` 的圆，以及以 ``Q`` 为圆心经过 ``P`` 的圆。前两种要求
``Q`` 分别位于 ``P M3`` 或圆 ``Circle(P,M3)`` 上；第三种要求 ``Q``
位于线段 ``P M3`` 的垂直平分线上。程序在三个严格正规夹具上同步
绑定交点。

这是多实例浮点筛查，不是 3 E 下界。它不覆盖第二步依赖首步新点的
动态形状，不覆盖第三步使用两个派生点的形状，也不覆盖需要四个新对象
的直接接触点程序。任何命中仍须用精确算术独立重放。
"""

from __future__ import annotations

from collections import Counter, defaultdict

from replay_mannheim_ordered_branches import (
    OrderedBranchReplay,
    branch_data,
)
from replay_mannheim_three_block_dependencies import collapse_point
from search_mannheim_ext_tau_4e import (
    FIXTURES,
    build_state,
    point_signature,
)
from search_parallel_3e import (
    DrawableBundle,
    PointBundle,
    circle_through,
    drawable_signature,
    generate_candidates,
    intersections,
    line_through,
)


EPSILON = 1e-9
TARGET_TOLERANCE = 2e-7
INITIAL_COST = 9
MODES = ("line", "center_known", "center_new")


def scalar_value(value) -> float:
    if hasattr(value, "approximate"):
        return value.approximate()
    return float(value)


def target_points() -> tuple[tuple[float, float], ...]:
    rows = []
    for index, fixture in enumerate(FIXTURES):
        replay = OrderedBranchReplay(
            f"ext_contact_static_{index}",
            *fixture,
            emit=False,
        )
        data = branch_data(replay, "P0")
        if data["kind"] != "regular":
            raise AssertionError("直接接触点搜索夹具的 P0 必须严格正规")
        pair = replay.verify_pair(
            "P0",
            data["tau"],
            allow_repeated_physical_signs=True,
        )
        key = next(key for key in pair if key.startswith("+++"))
        point = collapse_point(pair[key]["contact_3"])
        rows.append(tuple(scalar_value(value) for value in point))
    return tuple(rows)


def perpendicular_bisector(first, second):
    midpoint_x = (first[0] + second[0]) / 2
    midpoint_y = (first[1] + second[1]) / 2
    delta_x = second[0] - first[0]
    delta_y = second[1] - first[1]
    norm = (delta_x * delta_x + delta_y * delta_y) ** 0.5
    if norm <= EPSILON:
        return None
    a = delta_x / norm
    b = delta_y / norm
    c = -(a * midpoint_x + b * midpoint_y)
    if a < -EPSILON or (abs(a) <= EPSILON and b < 0):
        a, b, c = -a, -b, -c
    return (a, b, c)


def support_locus(mode, support, target):
    if mode == "line":
        return "line", line_through(support, target)
    if mode == "center_known":
        return "circle", circle_through(support, target)
    return "line", perpendicular_bisector(support, target)


def final_drawable(
    mode: str,
    support: PointBundle,
    values: tuple[tuple[float, float], ...],
) -> DrawableBundle | None:
    if mode == "line":
        rows = tuple(
            line_through(support_value, value)
            for support_value, value in zip(
                support.values,
                values,
                strict=True,
            )
        )
        kind = "line"
    elif mode == "center_known":
        rows = tuple(
            circle_through(support_value, value)
            for support_value, value in zip(
                support.values,
                values,
                strict=True,
            )
        )
        kind = "circle"
    else:
        rows = tuple(
            circle_through(value, support_value)
            for support_value, value in zip(
                support.values,
                values,
                strict=True,
            )
        )
        kind = "circle"
    if any(row is None for row in rows):
        return None
    return DrawableBundle("target_locus", kind, rows)  # type: ignore[arg-type]


def contains_targets(drawable, targets) -> bool:
    for value, target in zip(drawable.values, targets, strict=True):
        if drawable.kind == "line":
            residual = abs(value[0] * target[0] + value[1] * target[1] + value[2])
            scale = max(1.0, abs(value[2]), abs(target[0]), abs(target[1]))
        else:
            distance_squared = (
                (target[0] - value[0]) ** 2
                + (target[1] - value[1]) ** 2
            )
            residual = abs(distance_squared - value[2])
            scale = max(1.0, distance_squared, abs(value[2]))
        if residual > TARGET_TOLERANCE * scale:
            return False
    return True


def main() -> None:
    state, _ = build_state()
    targets = target_points()
    candidates = generate_candidates(state, INITIAL_COST + 1)
    existing = {
        drawable_signature(item): item.drawable_id
        for item in state.drawables
    }

    first_target_hits = tuple(
        candidate.describe()
        for candidate in candidates
        if contains_targets(candidate.drawable, targets)
    )
    if first_target_hits:
        raise AssertionError(
            "发现一步直接接触轨迹，需要精确重放："
            f"{first_target_hits[:3]}"
        )

    locus_count = 0
    shared_groups = 0
    coincident_support_groups = 0
    existing_final_groups = 0
    existing_final_ids = Counter()
    hits = []

    for mode in MODES:
        for support in state.points:
            locus_rows = tuple(
                support_locus(mode, support_value, target)
                for support_value, target in zip(
                    support.values,
                    targets,
                    strict=True,
                )
            )
            if any(value is None for _, value in locus_rows):
                continue
            kinds = {kind for kind, _ in locus_rows}
            if len(kinds) != 1:
                continue
            locus_count += 1
            locus_kind = locus_rows[0][0]
            locus_values = tuple(value for _, value in locus_rows)
            groups = defaultdict(list)

            for candidate in candidates:
                rows = tuple(
                    intersections(
                        locus_kind,
                        locus_value,
                        candidate.kind,
                        candidate_value,
                    )
                    for locus_value, candidate_value in zip(
                        locus_values,
                        candidate.drawable.values,
                        strict=True,
                    )
                )
                root_counts = {len(row) for row in rows}
                if len(root_counts) != 1 or not rows[0]:
                    continue
                for root_index in range(len(rows[0])):
                    values = tuple(row[root_index] for row in rows)
                    groups[point_signature(values)].append(
                        (candidate, values)
                    )

            for group in groups.values():
                descriptions = {
                    candidate.describe() for candidate, _ in group
                }
                if len(descriptions) < 2:
                    continue
                shared_groups += 1
                _, values = group[0]
                final = final_drawable(mode, support, values)
                if final is None:
                    coincident_support_groups += 1
                    continue
                if not contains_targets(final, targets):
                    raise AssertionError("支撑轨迹没有产生目标接触轨迹")
                final_signature = drawable_signature(final)
                if final_signature in existing:
                    existing_final_groups += 1
                    existing_final_ids[existing[final_signature]] += 1
                    continue
                hits.append(
                    (
                        mode,
                        support.point_id,
                        tuple(sorted(descriptions))[:3],
                    )
                )

    if hits:
        raise AssertionError(
            "发现静态 3 E 直接接触轨迹，需要精确重放："
            f"{hits[:3]}"
        )

    print(
        "mannheim_ext_contact_static_3e_search",
        {
            "samples": len(FIXTURES),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "one_step_candidates": len(candidates),
            "first_target_hits": len(first_target_hits),
            "support_loci": locus_count,
            "shared_locus_intersection_groups": shared_groups,
            "coincident_support_groups": coincident_support_groups,
            "existing_final_groups": existing_final_groups,
            "existing_final_ids": dict(sorted(existing_final_ids.items())),
            "three_e_contact_locus_hits": 0,
        },
    )


if __name__ == "__main__":
    main()
