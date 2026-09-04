"""搜索完成一个双 ``K'`` 平行块后，另一块能否再用 2 E 完成。

初始状态与 ``search_mannheim_double_kp_joint_5e`` 相同，随后明确加入
``P2`` 或 ``P3`` 的标准两圆一线修复及其全部免费交点。枚举一个新辅助
对象；若它与已有对象产生一个严格不同于剩余有限对角点、且位于目标
平行线上的交点，最后再画目标线即可，总联合修复为 5 E。

这是三样本浮点筛查。命中必须精确重放；无命中只排除这种顺序与依赖
形状，不构成 2 E 下界。
"""

from __future__ import annotations

from itertools import combinations
from math import hypot
from time import monotonic

from search_mannheim_double_kp_joint_5e import (
    FIXTURES,
    build_state,
    candidate_target_witness,
    target_point_pair,
)
from replay_mannheim_three_block_dependencies import role_batch_keys
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
    circle_through,
    generate_candidates,
    intersections,
    line_through,
)


def point(state, point_id: str) -> PointBundle:
    return next(item for item in state.points if item.point_id == point_id)


def role_point(state, profile: str, role: str) -> PointBundle:
    return point(state, role_batch_keys(profile)[role])


def ensure_free_point(
    state: State,
    point_id: str,
    expected_values,
) -> tuple[State, PointBundle]:
    for item in state.points:
        if all(
            hypot(actual[0] - expected[0], actual[1] - expected[1]) <= 1e-6
            for actual, expected in zip(
                item.values,
                expected_values,
                strict=True,
            )
        ):
            return state, item

    for first, second in combinations(state.drawables, 2):
        for first_value, second_value, expected in zip(
            first.values,
            second.values,
            expected_values,
            strict=True,
        ):
            roots = intersections(
                first.kind,
                first_value,
                second.kind,
                second_value,
            )
            if not any(
                hypot(root[0] - expected[0], root[1] - expected[1]) <= 1e-6
                for root in roots
            ):
                break
        else:
            point_bundle = PointBundle(point_id, tuple(expected_values))
            return (
                State(
                    state.points + (point_bundle,),
                    state.drawables,
                    state.program,
                ),
                point_bundle,
            )
    raise AssertionError("标准修复的预期免费交点没有出现")


def add_candidate(state, kind: str, first: PointBundle, second: PointBundle, move: int):
    constructor = line_through if kind == "line" else circle_through
    values = tuple(
        constructor(first_value, second_value)
        for first_value, second_value in zip(
            first.values,
            second.values,
            strict=True,
        )
    )
    if any(value is None for value in values):
        raise AssertionError("标准修复对象没有定义")
    candidate = Candidate(
        kind,
        first.point_id,
        second.point_id,
        DrawableBundle(f"move_{move}", kind, values),
    )
    return apply_candidate(state, candidate, move)


def build_after_repair(profile: str, initial=None):
    state, targets = build_state() if initial is None else initial
    first = role_point(state, profile, "x")
    second = role_point(state, profile, "y")
    finite = point(state, f"{profile}_K")

    state = add_candidate(state, "circle", first, finite, 1)
    reflected_values = tuple(
        (
            2 * x[0] - k[0],
            2 * x[1] - k[1],
        )
        for x, k in zip(first.values, finite.values, strict=True)
    )
    state, reflected = ensure_free_point(
        state,
        f"{profile}_reflected",
        reflected_values,
    )

    state = add_candidate(state, "circle", second, reflected, 2)
    q_values = tuple(
        (
            k[0] + x[0] - y[0],
            k[1] + x[1] - y[1],
        )
        for k, x, y in zip(
            finite.values,
            first.values,
            second.values,
            strict=True,
        )
    )
    state, q = ensure_free_point(state, f"{profile}_Q", q_values)
    state = add_candidate(state, "line", finite, q, 3)
    return state, targets


def main() -> None:
    for first_profile, remaining_profile in (("P2", "P3"), ("P3", "P2")):
        state, targets = build_after_repair(first_profile)
        finite = point(state, f"{remaining_profile}_K")
        free_pair = target_point_pair(state, targets[remaining_profile])
        if free_pair is not None:
            raise AssertionError(
                f"{remaining_profile} 目标线在加入辅助对象前已经可画"
            )

        start = monotonic()
        candidates = generate_candidates(state, 4)
        hits = []
        for candidate in candidates:
            witness = candidate_target_witness(
                candidate,
                state,
                targets[remaining_profile],
                finite,
            )
            if witness is not None:
                hits.append((candidate.describe(), witness))

        print(
            "summary",
            {
                "samples": len(FIXTURES),
                "first_profile": first_profile,
                "remaining_profile": remaining_profile,
                "initial_points": len(state.points),
                "initial_drawables": len(state.drawables),
                "candidates": len(candidates),
                "hits": len(hits),
                "elapsed_seconds": round(monotonic() - start, 3),
            },
        )
        for hit in hits:
            print("candidate", hit)


if __name__ == "__main__":
    main()
