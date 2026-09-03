"""筛查正规 Mannheim 单根接触弦能否在 4 E 内显式画出。

初态是 5 E 合法平行前缀与 ``P0`` 的四条批量线，共 9 E。目标是
Mannheim 接触弦 ``tau``。本轮限定程序形状为：第二步产生第一个
``tau`` 上点，第三步产生第二个不同点，第四步连接二点。初态没有这种
点，任何首步候选也不产生这种点。

搜索在三个严格正规夹具上同步执行。为避免枚举重复的旧点对候选，第二
步被拆成两类：初态已可画的候选由预计算的 ``tau`` 交点签名匹配；依赖
首步新点的候选单独增量生成。所有命中仍须另作精确重放。

这是多实例浮点筛查，不是 4 E 下界。它不覆盖前两步都没有弦上点、
第三步由一圆同时产生两个弦上点的延迟形状，也不排除第四步画出与
第三圆具有同一公共弦的圆，或绕过完整接触弦而直接取得 ``+++`` 接触点
或目标圆心的程序。
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from time import monotonic

from replay_mannheim_ordered_branches import OrderedBranchReplay
from search_mannheim_center_locus_2e import object_value
from search_mannheim_root_center_2e import (
    FIXTURES,
    add_known_drawable,
    make_sample,
)
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
    circle_through,
    drawable_signature,
    generate_candidates,
    intersections,
    line_through,
)


ROUND_POINT = 6
ROUND_DRAWABLE = 7
INITIAL_COST = 9
PROGRESS_INTERVAL = 250
P0_BATCH = ("alphaA", "aB", "a1A", "alpha1B")


@dataclass(frozen=True, slots=True)
class TauPointPath:
    first: Candidate
    second: Candidate
    signatures: frozenset[tuple[float, ...]]


def point_signature(values) -> tuple[float, ...]:
    return tuple(
        round(coordinate, ROUND_POINT)
        for point in values
        for coordinate in point
    )


def tolerant_drawable_signature(drawable: DrawableBundle):
    return (
        drawable.kind,
        *(
            round(coordinate, ROUND_DRAWABLE)
            for value in drawable.values
            for coordinate in value
        ),
    )


def tau_intersection_signatures(drawable, tau_lines):
    rows = tuple(
        intersections("line", tau, drawable.kind, value)
        for tau, value in zip(tau_lines, drawable.values, strict=True)
    )
    root_counts = {len(row) for row in rows}
    if len(root_counts) != 1 or not rows[0]:
        return frozenset()
    return frozenset(
        point_signature(tuple(row[root_index] for row in rows))
        for root_index in range(len(rows[0]))
    )


def line_is_tau(drawable, tau_lines) -> bool:
    if drawable.kind != "line":
        return False
    return all(
        max(
            abs(left - right)
            for left, right in zip(value, tau, strict=True)
        )
        <= 2e-7
        for value, tau in zip(drawable.values, tau_lines, strict=True)
    )


def build_state():
    samples = tuple(make_sample(*fixture) for fixture in FIXTURES)
    points = tuple(
        PointBundle(
            point_id,
            tuple(sample.points[point_id] for sample in samples),
        )
        for point_id in ("O1", "O2", "O3")
    )
    drawables = tuple(
        DrawableBundle(
            drawable_id,
            samples[0].drawables[drawable_id][0],
            tuple(sample.drawables[drawable_id][1] for sample in samples),
        )
        for drawable_id in ("Gamma1", "Gamma2", "Gamma3")
    )
    state = State(points, drawables)
    paid_ids = (
        "ell",
        "parallel_c0",
        "parallel_cX",
        "parallel_diameter",
        "ell3",
        *(f"batch_{key}" for key in P0_BATCH),
    )
    for move, drawable_id in enumerate(paid_ids, start=1):
        state = add_known_drawable(
            state,
            drawable_id,
            samples[0].drawables[drawable_id][0],
            tuple(sample.drawables[drawable_id][1] for sample in samples),
            move,
        )

    tau_lines = []
    for index, fixture in enumerate(FIXTURES):
        replay = OrderedBranchReplay(
            f"ext_tau_{index}",
            *fixture,
            emit=False,
        )
        report = replay.run()
        if report["branches"]["P0"] != "regular":
            raise AssertionError("接触弦搜索夹具的 P0 必须严格正规")
        tau_lines.append(object_value(replay.objects, "P0_tau")[1])
    return state, tuple(tau_lines)


def make_candidate(kind, first, second, move):
    values = []
    for first_value, second_value in zip(
        first.values,
        second.values,
        strict=True,
    ):
        value = (
            line_through(first_value, second_value)
            if kind == "line"
            else circle_through(first_value, second_value)
        )
        if value is None:
            return None
        values.append(value)
    drawable = DrawableBundle(f"move_{move}", kind, tuple(values))
    return Candidate(kind, first.point_id, second.point_id, drawable)


def generate_incremental_candidates(
    old_state: State,
    state: State,
    move: int,
) -> tuple[Candidate, ...]:
    """只生成至少使用一个新增点的候选。"""

    old_count = len(old_state.points)
    existing = {drawable_signature(item) for item in state.drawables}
    seen = set()
    candidates = []

    for first_index, first in enumerate(state.points):
        for second_index in range(first_index + 1, len(state.points)):
            if first_index < old_count and second_index < old_count:
                continue
            second = state.points[second_index]
            candidate = make_candidate("line", first, second, move)
            if candidate is None:
                continue
            signature = drawable_signature(candidate.drawable)
            if signature not in existing and signature not in seen:
                seen.add(signature)
                candidates.append(candidate)

    for center_index, center in enumerate(state.points):
        for through_index, through in enumerate(state.points):
            if center_index == through_index:
                continue
            if center_index < old_count and through_index < old_count:
                continue
            candidate = make_candidate("circle", center, through, move)
            if candidate is None:
                continue
            signature = drawable_signature(candidate.drawable)
            if signature not in existing and signature not in seen:
                seen.add(signature)
                candidates.append(candidate)
    return tuple(candidates)


def state_tau_intersection_signatures(state, tau_lines):
    return frozenset().union(
        *(
            tau_intersection_signatures(drawable, tau_lines)
            for drawable in state.drawables
        )
    )


def candidate_creates_tau_point(candidate, existing_signatures, tau_lines):
    return (
        tau_intersection_signatures(candidate.drawable, tau_lines)
        & existing_signatures
    )


def main() -> None:
    start = monotonic()
    state, tau_lines = build_state()
    initial_tau_points = tuple(
        point.point_id
        for point in state.points
        if all(
            abs(tau[0] * value[0] + tau[1] * value[1] + tau[2]) <= 2e-7
            for value, tau in zip(point.values, tau_lines, strict=True)
        )
    )
    if initial_tau_points:
        raise AssertionError("9 E 初态意外已有接触弦上的点")

    first_candidates = generate_candidates(state, INITIAL_COST + 1)
    initial_tau_signatures = state_tau_intersection_signatures(
        state,
        tau_lines,
    )
    if any(
        candidate_creates_tau_point(
            candidate,
            initial_tau_signatures,
            tau_lines,
        )
        for candidate in first_candidates
    ):
        raise AssertionError("首步候选意外产生接触弦上的点")

    # 初态候选之间共享同一 tau 交点时，两步即可产生该点。先建立静态
    # 索引，随后每个首步只增量生成真正依赖新点的第二步候选。
    static_by_signature = defaultdict(list)
    for candidate in first_candidates:
        for signature in tau_intersection_signatures(
            candidate.drawable,
            tau_lines,
        ):
            static_by_signature[signature].append(candidate)

    paths: dict[tuple, TauPointPath] = {}
    incremental_path_keys = set()
    incremental_second_candidates = 0
    for first_index, first in enumerate(first_candidates, start=1):
        first_signatures = tau_intersection_signatures(
            first.drawable,
            tau_lines,
        )
        for signature in first_signatures:
            for second in static_by_signature[signature]:
                if second.describe() == first.describe():
                    continue
                key = tuple(sorted((first.describe(), second.describe())))
                paths[key] = TauPointPath(
                    first,
                    second,
                    frozenset({signature}),
                )

        state_one = apply_candidate(state, first, INITIAL_COST + 1)
        existing_drawables = {
            tolerant_drawable_signature(drawable)
            for drawable in state_one.drawables
        }
        state_one_tau_signatures = state_tau_intersection_signatures(
            state_one,
            tau_lines,
        )
        for second in generate_incremental_candidates(
            state,
            state_one,
            INITIAL_COST + 2,
        ):
            incremental_second_candidates += 1
            if (
                tolerant_drawable_signature(second.drawable)
                in existing_drawables
            ):
                continue
            signatures = candidate_creates_tau_point(
                second,
                state_one_tau_signatures,
                tau_lines,
            )
            if not signatures:
                continue
            key = (first.describe(), second.describe())
            paths[key] = TauPointPath(first, second, signatures)
            incremental_path_keys.add(key)

        if first_index % PROGRESS_INTERVAL == 0:
            print(
                "progress",
                {
                    "first": first_index,
                    "first_total": len(first_candidates),
                    "incremental_second_candidates": (
                        incremental_second_candidates
                    ),
                    "two_step_tau_point_paths": len(paths),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    four_step_hits = []
    third_candidates = 0
    for path in paths.values():
        state_one = apply_candidate(state, path.first, INITIAL_COST + 1)
        second = next(
            (
                candidate
                for candidate in generate_candidates(
                    state_one,
                    INITIAL_COST + 2,
                )
                if candidate.describe() == path.second.describe()
            ),
            None,
        )
        if second is None:
            continue
        state_two = apply_candidate(state_one, second, INITIAL_COST + 2)
        existing_drawables = {
            tolerant_drawable_signature(drawable)
            for drawable in state_two.drawables
        }
        state_two_tau_signatures = state_tau_intersection_signatures(
            state_two,
            tau_lines,
        )
        for third in generate_candidates(state_two, INITIAL_COST + 3):
            third_candidates += 1
            if (
                tolerant_drawable_signature(third.drawable)
                in existing_drawables
            ):
                continue
            if line_is_tau(third.drawable, tau_lines):
                four_step_hits.append(
                    (path.first.describe(), second.describe(), third.describe())
                )
                continue
            new_signatures = candidate_creates_tau_point(
                third,
                state_two_tau_signatures,
                tau_lines,
            ) - path.signatures
            if new_signatures:
                four_step_hits.append(
                    (
                        path.first.describe(),
                        second.describe(),
                        third.describe(),
                        "Line(first_tau_point, second_tau_point)",
                    )
                )

    if four_step_hits:
        raise AssertionError(
            "发现候选 4 E 接触弦程序，需要逐项精确重放："
            f"{four_step_hits[:3]}"
        )

    repeated_static_groups = tuple(
        candidates
        for candidates in static_by_signature.values()
        if len(candidates) > 1
    )
    print(
        "mannheim_ext_tau_4e_search",
        {
            "samples": len(FIXTURES),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "first_candidates": len(first_candidates),
            "static_tau_intersection_groups": len(repeated_static_groups),
            "incremental_second_candidates": incremental_second_candidates,
            "two_step_tau_point_paths": len(paths),
            "incremental_tau_point_paths": len(incremental_path_keys),
            "third_candidates": third_candidates,
            "four_step_tau_hits": 0,
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for index, candidates in enumerate(repeated_static_groups, start=1):
        print(
            "static_tau_group",
            index,
            tuple(candidate.describe() for candidate in candidates),
        )


if __name__ == "__main__":
    main()
