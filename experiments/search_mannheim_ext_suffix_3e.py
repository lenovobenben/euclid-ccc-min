"""筛查正规 Mannheim 单根的 4 E 定心后缀能否降为 3 E。

初始状态精确重放 ``P0`` 的 5 E 合法前缀、四条批量线和五线接触弦
核心，共 14 E；此时 ``+++`` 在第三圆上的接触点已经出现。若最终三步
能够画出目标圆，则目标圆心必须在前两步后出现。初始状态没有经过目标
圆心的已画对象，所以第一、第二个新对象都必须经过该圆心。

程序在三个严格正规 ``D8`` 夹具上同步枚举。首步候选中只有已有方案的
第三圆接触半径经过目标圆心；画出它以后，全部第二步候选均不再经过
目标圆心。因此没有发现 3 E 后缀。

这是多实例浮点筛查，不是 3 E 下界。零命中只排除当前夹具套件、点
绑定和基础 ``Line/ Circle`` 候选宇宙；任何未来命中都必须精确重放。
"""

from __future__ import annotations

from fractions import Fraction
from time import monotonic

from replay_mannheim_ordered_branches import OrderedBranchReplay
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_root_center_2e import (
    add_known_drawable,
    drawable_contains_targets,
    make_sample,
)
from search_parallel_3e import (
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
    generate_candidates,
)


F = Fraction
FIXTURES = (
    (
        ((F(0), F(0)), (F(13), F(0)), (F(4), F(15))),
        (F(4), F(2), F(1)),
    ),
    (
        ((F(0), F(0)), (F(17), F(0)), (F(5), F(12))),
        (F(5), F(3), F(1)),
    ),
    (
        ((F(0), F(0)), (F(12), F(0)), (F(1), F(11))),
        (F(7), F(3), F(2)),
    ),
)

PRE_SUFFIX_COST = 14
EXPECTED_FIRST = "Line(O3, object_14&Gamma3[1])"


def same_point_bundle(values, targets) -> bool:
    return all(
        abs(point[0] - target[0]) + abs(point[1] - target[1]) <= 1e-7
        for point, target in zip(values, targets, strict=True)
    )


def build_state():
    samples = tuple(make_sample(*fixture) for fixture in FIXTURES)
    replays = []
    for index, (centers, radii) in enumerate(FIXTURES):
        replay = OrderedBranchReplay(
            f"ext_suffix_{index}",
            centers,
            radii,
            emit=False,
        )
        report = replay.run()
        if report["branches"]["P0"] != "regular":
            raise AssertionError("单根后缀搜索夹具的 P0 必须严格正规")
        if report["first_ext"] != 18:
            raise AssertionError("单根后缀搜索夹具没有校准到 18 E")
        replays.append(replay)

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

    paid_sequences = tuple(
        tuple(replay.objects.graph.paid_order[:PRE_SUFFIX_COST])
        for replay in replays
    )
    if len(set(paid_sequences)) != 1:
        raise AssertionError("三个夹具的 14 E 前缀对象顺序不一致")
    paid_sequence = paid_sequences[0]
    if paid_sequence[-1] != "P0_tau":
        raise AssertionError("14 E 前缀没有结束于 P0 接触弦")

    for move, node_id in enumerate(paid_sequence, start=1):
        rows = tuple(object_value(replay.objects, node_id) for replay in replays)
        kinds = {kind for kind, _ in rows}
        if len(kinds) != 1:
            raise AssertionError(f"{node_id} 的对象种类不一致")
        state = add_known_drawable(
            state,
            node_id,
            rows[0][0],
            tuple(value for _, value in rows),
            move,
        )

    target_centers = []
    for replay in replays:
        target_key = next(
            key for key in replay.targets if key.startswith("P0:+++@")
        )
        target_centers.append(
            point_value(replay.objects, f"P0_{target_key}_center")
        )
    return state, tuple(target_centers)


def main() -> None:
    start = monotonic()
    state, target_centers = build_state()
    existing_target_points = tuple(
        point.point_id
        for point in state.points
        if same_point_bundle(point.values, target_centers)
    )
    existing_objects_through_center = tuple(
        drawable.drawable_id
        for drawable in state.drawables
        if drawable_contains_targets(drawable, target_centers)
    )
    if existing_target_points or existing_objects_through_center:
        raise AssertionError("14 E 前缀意外已经定出目标圆心")

    first_candidates = generate_candidates(state, PRE_SUFFIX_COST + 1)
    first_hits = tuple(
        candidate
        for candidate in first_candidates
        if drawable_contains_targets(candidate.drawable, target_centers)
    )
    if tuple(candidate.describe() for candidate in first_hits) != (
        EXPECTED_FIRST,
    ):
        raise AssertionError("经过目标圆心的首步候选发生变化")

    state_one = apply_candidate(
        state,
        first_hits[0],
        PRE_SUFFIX_COST + 1,
    )
    second_candidates = generate_candidates(state_one, PRE_SUFFIX_COST + 2)
    second_hits = tuple(
        candidate.describe()
        for candidate in second_candidates
        if drawable_contains_targets(candidate.drawable, target_centers)
    )
    if second_hits:
        raise AssertionError("发现候选 3 E 定心后缀，需要精确重放")

    print(
        "mannheim_ext_suffix_3e_search",
        {
            "samples": len(FIXTURES),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "first_candidates": len(first_candidates),
            "first_center_hits": len(first_hits),
            "first_hit": first_hits[0].describe(),
            "points_after_radius": len(state_one.points),
            "second_candidates": len(second_candidates),
            "second_center_hits": len(second_hits),
            "three_e_suffix_hits": 0,
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )


if __name__ == "__main__":
    main()
