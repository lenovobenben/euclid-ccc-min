"""筛查 Mannheim 平行退化块能否用 3 E 完成目标弦。

``centered`` 的三个样本满足 ``P0`` 的 ``K'`` 位于无穷远点，并额外
满足有限对角点 ``K = O3``；``noncentered`` 的三个样本只满足前一条件。
初始状态保留第三圆、用于构造有限 ``K`` 的两条弦，以及一条作为目标
方向的弦；搜索目标是画出经过 ``K`` 且平行于方向弦的直线，或画出一个
与第三圆的公共弦恰为该目标线的圆。

居中套件的最短候选不依赖初始三条弦，裁剪后确实只有 3 个计费祖先。
非居中套件若命中，还必须把候选实际使用的初始弦计回 E 分数。

这是多实例浮点筛查，不是 3 E 下界证明。
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import hypot
from time import monotonic

from search_parallel_3e import (
    Circle,
    DrawableBundle,
    Line,
    Point,
    PointBundle,
    State,
    apply_candidate,
    generate_candidates,
    intersect_circle_circle,
    intersect_line_circle,
    intersect_line_line,
    line_through,
)


@dataclass(frozen=True, slots=True)
class Fixture:
    radius_1: float
    radius_2: float
    radius_3: float
    distance_12: float
    x3: float
    y3: float


CENTERED_FIXTURES = (
    Fixture(3.0, 2.0, 1.0, 9.0, 6.0, 4.0),
    Fixture(4.0, 2.0, 1.0, 8.0, 6.0, 3.0),
    Fixture(5.0, 4.0, 1.0, 14.0, 8.0, 6.0),
)

NONCENTERED_FIXTURES = (
    Fixture(5.0, 4.0, 1.0, 28.0, 13.0, 13.5),
    Fixture(6.0, 2.0, 1.0, 21.0, 7.0, 9.0),
    Fixture(6.0, 3.0, 1.0, 27.0, 9.0, 12.0),
)


def other_circle_point(
    line: Line,
    circle: Circle,
    known: Point,
) -> Point:
    intersections = intersect_line_circle(line, circle)
    candidates = tuple(
        point
        for point in intersections
        if hypot(point[0] - known[0], point[1] - known[1]) > 1e-7
    )
    if len(candidates) != 1:
        raise AssertionError("不能唯一绑定批量线的第二交点")
    return candidates[0]


def normalized_parallel_through(point: Point, reference: Line) -> Line:
    a, b, _ = reference
    return (a, b, -a * point[0] - b * point[1])


def initial_state(
    fixtures: tuple[Fixture, ...],
    require_centered: bool,
) -> tuple[State, tuple[Line, ...]]:
    point_rows = {name: [] for name in ("O", "x", "y", "z", "w", "K")}
    drawable_rows = {
        name: [] for name in ("Gamma", "K_left", "K_right", "direction")
    }
    targets = []
    for fixture in fixtures:
        o3 = (fixture.x3, fixture.y3)
        gamma = (fixture.x3, fixture.y3, fixture.radius_3**2)
        named = {
            "alpha": (-fixture.radius_1, 0.0),
            "a": (fixture.radius_1, 0.0),
            "a1": (fixture.distance_12 - fixture.radius_2, 0.0),
            "alpha1": (fixture.distance_12 + fixture.radius_2, 0.0),
            "A": (fixture.x3 - fixture.radius_3, fixture.y3),
            "B": (fixture.x3 + fixture.radius_3, fixture.y3),
        }
        batch_points = {}
        for key in ("alphaA", "aB", "a1A", "alpha1B"):
            source_name = key[:-1]
            endpoint_name = key[-1]
            batch_line = line_through(
                named[source_name],
                named[endpoint_name],
            )
            if batch_line is None:
                raise AssertionError("批量线没有定义")
            batch_points[key] = other_circle_point(
                batch_line,
                gamma,
                named[endpoint_name],
            )
        x = batch_points["aB"]
        y = batch_points["alphaA"]
        z = batch_points["a1A"]
        w = batch_points["alpha1B"]
        k_left = line_through(x, w)
        k_right = line_through(y, z)
        direction = line_through(x, y)
        parallel_partner = line_through(z, w)
        if None in (k_left, k_right, direction, parallel_partner):
            raise AssertionError("P0 弦没有定义")
        k_intersections = intersect_line_line(k_left, k_right)  # type: ignore[arg-type]
        if len(k_intersections) != 1:
            raise AssertionError("有限对角点 K 没有定义")
        k = k_intersections[0]
        centered = hypot(k[0] - o3[0], k[1] - o3[1]) <= 1e-7
        if centered != require_centered:
            raise AssertionError("夹具的 K=O3 分支分类错误")
        if require_centered:
            if hypot(
                x[0] + w[0] - 2 * o3[0],
                x[1] + w[1] - 2 * o3[1],
            ) > 1e-7:
                raise AssertionError("夹具的 x,w 不是对径点")
            if hypot(
                y[0] + z[0] - 2 * o3[0],
                y[1] + z[1] - 2 * o3[1],
            ) > 1e-7:
                raise AssertionError("夹具的 y,z 不是对径点")
        if intersect_line_line(direction, parallel_partner):  # type: ignore[arg-type]
            raise AssertionError("夹具的 K' 定义弦不平行")

        for name, point in zip(
            ("O", "x", "y", "z", "w", "K"),
            (o3, x, y, z, w, k),
            strict=True,
        ):
            point_rows[name].append(point)
        drawable_rows["Gamma"].append(gamma)
        drawable_rows["K_left"].append(k_left)
        drawable_rows["K_right"].append(k_right)
        drawable_rows["direction"].append(direction)
        targets.append(
            normalized_parallel_through(k, direction)  # type: ignore[arg-type]
        )

    points = tuple(
        PointBundle(name, tuple(values)) for name, values in point_rows.items()
    )
    drawables = (
        DrawableBundle("Gamma", "circle", tuple(drawable_rows["Gamma"])),
        DrawableBundle("K_left", "line", tuple(drawable_rows["K_left"])),
        DrawableBundle("K_right", "line", tuple(drawable_rows["K_right"])),
        DrawableBundle("direction", "line", tuple(drawable_rows["direction"])),
    )
    return State(points, drawables), tuple(targets)


def is_target_line(candidate: DrawableBundle, targets: tuple[Line, ...]) -> bool:
    if candidate.kind != "line":
        return False
    return all(
        max(abs(left - right) for left, right in zip(value, target, strict=True))
        <= 1e-7
        for value, target in zip(candidate.values, targets, strict=True)
    )


def is_target_chord_circle(
    candidate: DrawableBundle,
    state: State,
    targets: tuple[Line, ...],
) -> bool:
    if candidate.kind != "circle":
        return False
    gamma = next(item for item in state.drawables if item.drawable_id == "Gamma")
    for value, gamma_value, target in zip(
        candidate.values,
        gamma.values,
        targets,
        strict=True,
    ):
        intersections = intersect_circle_circle(
            value,  # type: ignore[arg-type]
            gamma_value,  # type: ignore[arg-type]
        )
        if len(intersections) != 2:
            return False
        for point in intersections:
            if abs(target[0] * point[0] + target[1] * point[1] + target[2]) > 1e-7:
                return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--suite",
        choices=("centered", "noncentered"),
        default="centered",
    )
    args = parser.parse_args()
    fixture_sets = {
        "centered": (CENTERED_FIXTURES, True),
        "noncentered": (NONCENTERED_FIXTURES, False),
    }
    fixtures, require_centered = fixture_sets[args.suite]
    start = monotonic()
    initial, targets = initial_state(fixtures, require_centered)
    first_candidates = generate_candidates(initial, 1)
    states_after_two = 0
    two_step_hits: list[tuple[str, ...]] = []
    third_candidates = 0
    hits: list[tuple[str, ...]] = []
    for first_index, first in enumerate(first_candidates, start=1):
        state_one = apply_candidate(initial, first, 1)
        for second in generate_candidates(state_one, 2):
            state_two = apply_candidate(state_one, second, 2)
            states_after_two += 1
            if is_target_line(second.drawable, targets) or is_target_chord_circle(
                second.drawable,
                state_one,
                targets,
            ):
                two_step_hits.append(state_two.program)
            for third in generate_candidates(state_two, 3):
                third_candidates += 1
                if is_target_line(third.drawable, targets) or is_target_chord_circle(
                    third.drawable,
                    state_two,
                    targets,
                ):
                    hits.append(state_two.program + (third.describe(),))
        if first_index % 10 == 0:
            print(
                "progress",
                {
                    "first": first_index,
                    "first_total": len(first_candidates),
                    "states_after_two": states_after_two,
                    "third_candidates": third_candidates,
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )
    print(
        "summary",
        {
            "suite": args.suite,
            "samples": len(fixtures),
            "first_candidates": len(first_candidates),
            "states_after_two": states_after_two,
            "two_step_hits": len(two_step_hits),
            "third_candidates": third_candidates,
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for hit in two_step_hits:
        print("two_step_candidate")
        for move, description in enumerate(hit, start=1):
            print(f"  {move}. {description}")
    for hit in hits:
        print("candidate")
        for move, description in enumerate(hit, start=1):
            print(f"  {move}. {description}")


if __name__ == "__main__":
    main()
