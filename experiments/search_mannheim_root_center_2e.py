"""筛查 13 E Mannheim 公共前缀后能否在两步内取得根心 ``S``。

搜索在三个严格正规 ``D8`` 夹具上同步执行。初始状态实际重放合法 5 E
平行前缀和八条批量线，并保留每个新对象与既有对象产生的全部同型有限
实交点。若初始对象都不经过 ``S``，任何两步程序的第一个新对象必须
经过 ``S``；否则第二个对象无法与第一步或旧对象在 ``S`` 相交。脚本先
用这一必要条件过滤首步，再枚举可能的第二步。

这是多实例浮点筛查，不是两步不可能性的精确证明。任何命中必须另作
精确重放；零命中只排除当前夹具套件和点绑定策略覆盖的依赖形状。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import hypot, sqrt
from time import monotonic

from scan_mannheim_degeneracies import analyze_fixture, is_d8
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
    circle_through,
    drawable_signature,
    generate_candidates,
    line_through,
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
        ((F(0), F(0)), (F(19), F(0)), (F(10), F(11))),
        (F(6), F(4), F(2)),
    ),
)

BATCH_KEYS = (
    "alphaA",
    "aB",
    "a1A",
    "alpha1B",
    "alphaB",
    "aA",
    "a1B",
    "alpha1A",
)


@dataclass(frozen=True, slots=True)
class PrefixSample:
    points: dict[str, tuple[float, float]]
    drawables: dict[str, tuple[str, tuple[float, ...]]]
    root_center: tuple[float, float]


def normalized_line(first, second):
    value = line_through(first, second)
    if value is None:
        raise AssertionError("公共前缀出现未定义直线")
    return value


def normalized_circle(center, through):
    value = circle_through(center, through)
    if value is None:
        raise AssertionError("公共前缀出现零半径圆")
    return value


def radical_center(centers, radii) -> tuple[float, float]:
    o1, o2, o3 = centers
    r1, r2, r3 = radii
    d = o2[0] - o1[0]
    x = (
        o2[0] ** 2
        + o2[1] ** 2
        - r2**2
        - o1[0] ** 2
        - o1[1] ** 2
        + r1**2
    ) / (2 * d)
    y = (
        o3[0] ** 2
        + o3[1] ** 2
        - r3**2
        - o1[0] ** 2
        - o1[1] ** 2
        + r1**2
        - 2 * (o3[0] - o1[0]) * x
    ) / (2 * (o3[1] - o1[1]))
    return float(x), float(y)


def make_sample(centers, radii) -> PrefixSample:
    if not is_d8(centers, radii) or analyze_fixture(centers, radii):
        raise AssertionError("根心搜索夹具必须是严格正规 D8")
    o1, o2, o3 = (
        tuple(float(coordinate) for coordinate in point) for point in centers
    )
    r1, r2, r3 = (float(radius) for radius in radii)
    distance_13 = hypot(o3[0] - o1[0], o3[1] - o1[1])
    x = (o1[0] + distance_13, o1[1])
    q = (o3[0], 2 * o1[1] - o3[1])
    r = (2 * x[0] - q[0], 2 * x[1] - q[1])
    points = {
        "O1": o1,
        "O2": o2,
        "O3": o3,
        "alpha": (o1[0] - r1, o1[1]),
        "a": (o1[0] + r1, o1[1]),
        "a1": (o2[0] - r2, o2[1]),
        "alpha1": (o2[0] + r2, o2[1]),
        "A": (o3[0] - r3, o3[1]),
        "B": (o3[0] + r3, o3[1]),
        "parallel_X": x,
        "parallel_Q": q,
        "parallel_R": r,
    }
    drawables = {
        "Gamma1": ("circle", (o1[0], o1[1], r1**2)),
        "Gamma2": ("circle", (o2[0], o2[1], r2**2)),
        "Gamma3": ("circle", (o3[0], o3[1], r3**2)),
        "ell": ("line", normalized_line(o1, o2)),
        "parallel_c0": ("circle", normalized_circle(o1, o3)),
        "parallel_cX": ("circle", normalized_circle(x, o3)),
        "parallel_diameter": ("line", normalized_line(x, q)),
        "ell3": ("line", normalized_line(o3, r)),
    }
    for key in BATCH_KEYS:
        if key.startswith("alpha1"):
            source = "alpha1"
        elif key.startswith("alpha"):
            source = "alpha"
        elif key.startswith("a1"):
            source = "a1"
        else:
            source = "a"
        drawables[f"batch_{key}"] = (
            "line",
            normalized_line(points[source], points[key[-1]]),
        )
    return PrefixSample(
        points=points,
        drawables=drawables,
        root_center=radical_center(centers, radii),
    )


def add_known_drawable(
    state: State,
    drawable_id: str,
    kind: str,
    values,
    move: int,
) -> State:
    candidate = Candidate(
        kind,
        drawable_id,
        drawable_id,
        DrawableBundle(drawable_id, kind, tuple(values)),
    )
    return apply_candidate(state, candidate, move)


def build_state() -> tuple[State, tuple[tuple[float, float], ...]]:
    samples = tuple(make_sample(*fixture) for fixture in FIXTURES)
    initial_points = tuple(
        PointBundle(
            point_id,
            tuple(sample.points[point_id] for sample in samples),
        )
        for point_id in ("O1", "O2", "O3")
    )
    initial_drawables = tuple(
        DrawableBundle(
            drawable_id,
            samples[0].drawables[drawable_id][0],
            tuple(sample.drawables[drawable_id][1] for sample in samples),
        )
        for drawable_id in ("Gamma1", "Gamma2", "Gamma3")
    )
    state = State(initial_points, initial_drawables)
    paid_ids = (
        "ell",
        "parallel_c0",
        "parallel_cX",
        "parallel_diameter",
        "ell3",
        *(f"batch_{key}" for key in BATCH_KEYS),
    )
    for move, drawable_id in enumerate(paid_ids, start=1):
        kind = samples[0].drawables[drawable_id][0]
        values = tuple(sample.drawables[drawable_id][1] for sample in samples)
        state = add_known_drawable(state, drawable_id, kind, values, move)
    if len(state.drawables) != 3 + 13:
        raise AssertionError("公共前缀没有恰好加入 13 个计费对象")
    return state, tuple(sample.root_center for sample in samples)


def point_is_target(point: PointBundle, targets) -> bool:
    return all(
        hypot(value[0] - target[0], value[1] - target[1]) <= 1e-7
        for value, target in zip(point.values, targets, strict=True)
    )


def drawable_contains_targets(drawable: DrawableBundle, targets) -> bool:
    for kind_value, target in zip(drawable.values, targets, strict=True):
        if drawable.kind == "line":
            a, b, c = kind_value
            if abs(a * target[0] + b * target[1] + c) > 1e-7:
                return False
        else:
            center_x, center_y, radius_squared = kind_value
            actual = (target[0] - center_x) ** 2 + (target[1] - center_y) ** 2
            if abs(actual - radius_squared) > 1e-6 * max(1.0, radius_squared):
                return False
    return True


def main() -> None:
    state, targets = build_state()
    initial_target_points = tuple(
        point.point_id for point in state.points if point_is_target(point, targets)
    )
    old_objects_through_target = tuple(
        drawable.drawable_id
        for drawable in state.drawables
        if drawable_contains_targets(drawable, targets)
    )
    if initial_target_points or old_objects_through_target:
        raise AssertionError("公共前缀意外已经给出根心或经过根心的对象")

    start = monotonic()
    first_candidates = generate_candidates(state, 14)
    first_through_target = tuple(
        candidate
        for candidate in first_candidates
        if drawable_contains_targets(candidate.drawable, targets)
    )
    hits = []
    second_candidates = 0
    for first in first_through_target:
        state_one = apply_candidate(state, first, 14)
        for second in generate_candidates(state_one, 15):
            second_candidates += 1
            if drawable_contains_targets(second.drawable, targets):
                hits.append((first.describe(), second.describe()))

    report = {
        "samples": len(FIXTURES),
        "initial_points": len(state.points),
        "initial_drawables": len(state.drawables),
        "initial_target_points": len(initial_target_points),
        "old_objects_through_target": len(old_objects_through_target),
        "first_candidates": len(first_candidates),
        "first_through_target": len(first_through_target),
        "second_candidates": second_candidates,
        "hits": len(hits),
        "elapsed_seconds": round(monotonic() - start, 3),
    }
    print("root_center_2e_search", report)
    for hit in hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
