"""筛查 Mannheim 平行直径是否能在 3 E 内出现。

这是一项启发式、多实例数值筛查，不是下界证明。搜索状态遵循项目的基础
规则：直线由两个已有点确定，圆由已有圆心和已有圆上一点确定；每画一个
对象后，和全部既有对象的有限实交点自动加入点集。

目标有两种等价形式：

1. 实际画出过 O3 且平行于 O1O2 的直线；
2. 画出一个与 Gamma3 的公共弦恰为该直径的圆，从而免费得到 A、B。

程序在三个一般位置样本上同步执行同一组点依赖。只有三个样本都命中的
程序才报告为候选，以过滤单个校准实例中的偶然共线和共圆。
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import hypot, sqrt
from time import monotonic


EPSILON = 1e-9
ROUND_DIGITS = 9


@dataclass(frozen=True, slots=True)
class Sample:
    distance_12: float
    x3: float
    y3: float
    radius_1: float
    radius_2: float
    radius_3: float


SAMPLES = (
    Sample(13.0, 4.3, 15.0, 4.0, 2.0, 1.0),
    Sample(11.2, 2.7, 12.1, 3.2, 1.7, 0.8),
    Sample(15.4, 6.1, 10.8, 2.6, 4.1, 1.2),
)


Point = tuple[float, float]
Line = tuple[float, float, float]
Circle = tuple[float, float, float]
Drawable = Line | Circle


@dataclass(frozen=True, slots=True)
class PointBundle:
    point_id: str
    values: tuple[Point, ...]


@dataclass(frozen=True, slots=True)
class DrawableBundle:
    drawable_id: str
    kind: str
    values: tuple[Drawable, ...]


@dataclass(frozen=True, slots=True)
class Candidate:
    kind: str
    first_id: str
    second_id: str
    drawable: DrawableBundle

    def describe(self) -> str:
        operation = "Line" if self.kind == "line" else "Circle"
        return f"{operation}({self.first_id}, {self.second_id})"


@dataclass(frozen=True, slots=True)
class State:
    points: tuple[PointBundle, ...]
    drawables: tuple[DrawableBundle, ...]
    program: tuple[str, ...] = ()


def line_through(first: Point, second: Point) -> Line | None:
    a = first[1] - second[1]
    b = second[0] - first[0]
    c = first[0] * second[1] - second[0] * first[1]
    norm = hypot(a, b)
    if norm <= EPSILON:
        return None
    a, b, c = a / norm, b / norm, c / norm
    if a < -EPSILON or (abs(a) <= EPSILON and b < 0):
        a, b, c = -a, -b, -c
    return (a, b, c)


def circle_through(center: Point, through: Point) -> Circle | None:
    radius_squared = (
        (center[0] - through[0]) ** 2
        + (center[1] - through[1]) ** 2
    )
    if radius_squared <= EPSILON:
        return None
    return (center[0], center[1], radius_squared)


def intersect_line_line(first: Line, second: Line) -> tuple[Point, ...]:
    a, b, c = first
    d, e, f = second
    denominator = a * e - b * d
    if abs(denominator) <= EPSILON:
        return ()
    return (
        (
            (b * f - c * e) / denominator,
            (c * d - a * f) / denominator,
        ),
    )


def intersect_line_circle(line: Line, circle: Circle) -> tuple[Point, ...]:
    a, b, c = line
    center_x, center_y, radius_squared = circle
    signed_distance = a * center_x + b * center_y + c
    foot_x = center_x - a * signed_distance
    foot_y = center_y - b * signed_distance
    half_chord_squared = radius_squared - signed_distance**2
    if half_chord_squared < -EPSILON:
        return ()
    if abs(half_chord_squared) <= EPSILON:
        return ((foot_x, foot_y),)
    half_chord = sqrt(half_chord_squared)
    return tuple(
        sorted(
            (
                (foot_x + b * half_chord, foot_y - a * half_chord),
                (foot_x - b * half_chord, foot_y + a * half_chord),
            )
        )
    )


def intersect_circle_circle(
    first: Circle,
    second: Circle,
) -> tuple[Point, ...]:
    x0, y0, radius_0_squared = first
    x1, y1, radius_1_squared = second
    delta_x = x1 - x0
    delta_y = y1 - y0
    center_distance_squared = delta_x**2 + delta_y**2
    if center_distance_squared <= EPSILON:
        return ()
    center_distance = sqrt(center_distance_squared)
    radius_0 = sqrt(radius_0_squared)
    radius_1 = sqrt(radius_1_squared)
    if center_distance > radius_0 + radius_1 + EPSILON:
        return ()
    if center_distance < abs(radius_0 - radius_1) - EPSILON:
        return ()
    along = (
        radius_0_squared
        - radius_1_squared
        + center_distance_squared
    ) / (2 * center_distance)
    half_chord_squared = radius_0_squared - along**2
    if half_chord_squared < -EPSILON:
        return ()
    unit_x = delta_x / center_distance
    unit_y = delta_y / center_distance
    base_x = x0 + along * unit_x
    base_y = y0 + along * unit_y
    if abs(half_chord_squared) <= EPSILON:
        return ((base_x, base_y),)
    half_chord = sqrt(half_chord_squared)
    return tuple(
        sorted(
            (
                (
                    base_x - half_chord * unit_y,
                    base_y + half_chord * unit_x,
                ),
                (
                    base_x + half_chord * unit_y,
                    base_y - half_chord * unit_x,
                ),
            )
        )
    )


def intersections(
    first_kind: str,
    first: Drawable,
    second_kind: str,
    second: Drawable,
) -> tuple[Point, ...]:
    if first_kind == "line" and second_kind == "line":
        return intersect_line_line(first, second)  # type: ignore[arg-type]
    if first_kind == "line" and second_kind == "circle":
        return intersect_line_circle(first, second)  # type: ignore[arg-type]
    if first_kind == "circle" and second_kind == "line":
        return intersect_line_circle(second, first)  # type: ignore[arg-type]
    return intersect_circle_circle(first, second)  # type: ignore[arg-type]


def rounded(values: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(round(value, ROUND_DIGITS) for value in values)


def point_signature(point: PointBundle) -> tuple[float, ...]:
    return rounded(
        tuple(coordinate for value in point.values for coordinate in value)
    )


def drawable_signature(drawable: DrawableBundle) -> tuple[str, float, ...]:
    return (
        drawable.kind,
        *rounded(
            tuple(
                coordinate
                for value in drawable.values
                for coordinate in value
            )
        ),
    )


def initial_state() -> State:
    point_rows: dict[str, list[Point]] = {
        "O1": [],
        "O2": [],
        "O3": [],
        "g1_minus": [],
        "g1_plus": [],
        "g2_minus": [],
        "g2_plus": [],
    }
    line_values: list[Line] = []
    gamma_1_values: list[Circle] = []
    gamma_2_values: list[Circle] = []
    gamma_3_values: list[Circle] = []
    for sample in SAMPLES:
        point_rows["O1"].append((0.0, 0.0))
        point_rows["O2"].append((sample.distance_12, 0.0))
        point_rows["O3"].append((sample.x3, sample.y3))
        point_rows["g1_minus"].append((-sample.radius_1, 0.0))
        point_rows["g1_plus"].append((sample.radius_1, 0.0))
        point_rows["g2_minus"].append(
            (sample.distance_12 - sample.radius_2, 0.0)
        )
        point_rows["g2_plus"].append(
            (sample.distance_12 + sample.radius_2, 0.0)
        )
        line_values.append((0.0, 1.0, 0.0))
        gamma_1_values.append((0.0, 0.0, sample.radius_1**2))
        gamma_2_values.append(
            (
                sample.distance_12,
                0.0,
                sample.radius_2**2,
            )
        )
        gamma_3_values.append(
            (sample.x3, sample.y3, sample.radius_3**2)
        )
    points = tuple(
        PointBundle(point_id, tuple(values))
        for point_id, values in point_rows.items()
    )
    drawables = (
        DrawableBundle("ell", "line", tuple(line_values)),
        DrawableBundle("Gamma1", "circle", tuple(gamma_1_values)),
        DrawableBundle("Gamma2", "circle", tuple(gamma_2_values)),
        DrawableBundle("Gamma3", "circle", tuple(gamma_3_values)),
    )
    return State(points, drawables)


def generate_candidates(state: State, move: int) -> tuple[Candidate, ...]:
    existing = {drawable_signature(item) for item in state.drawables}
    seen: set[tuple[str, float, ...]] = set()
    candidates: list[Candidate] = []
    for first_index, first in enumerate(state.points):
        for second in state.points[first_index + 1 :]:
            values: list[Line] = []
            for first_value, second_value in zip(
                first.values,
                second.values,
                strict=True,
            ):
                value = line_through(first_value, second_value)
                if value is None:
                    break
                values.append(value)
            else:
                drawable = DrawableBundle(
                    f"move_{move}",
                    "line",
                    tuple(values),
                )
                signature = drawable_signature(drawable)
                if signature not in existing and signature not in seen:
                    seen.add(signature)
                    candidates.append(
                        Candidate(
                            "line",
                            first.point_id,
                            second.point_id,
                            drawable,
                        )
                    )
    for center in state.points:
        for through in state.points:
            if center.point_id == through.point_id:
                continue
            values: list[Circle] = []
            for center_value, through_value in zip(
                center.values,
                through.values,
                strict=True,
            ):
                value = circle_through(center_value, through_value)
                if value is None:
                    break
                values.append(value)
            else:
                drawable = DrawableBundle(
                    f"move_{move}",
                    "circle",
                    tuple(values),
                )
                signature = drawable_signature(drawable)
                if signature not in existing and signature not in seen:
                    seen.add(signature)
                    candidates.append(
                        Candidate(
                            "circle",
                            center.point_id,
                            through.point_id,
                            drawable,
                        )
                    )
    return tuple(candidates)


def apply_candidate(state: State, candidate: Candidate, move: int) -> State:
    new_drawable = DrawableBundle(
        f"object_{move}",
        candidate.kind,
        candidate.drawable.values,
    )
    known_point_signatures = {point_signature(point) for point in state.points}
    new_points: list[PointBundle] = []
    for existing in state.drawables:
        rows = tuple(
            intersections(
                new_drawable.kind,
                new_value,
                existing.kind,
                existing_value,
            )
            for new_value, existing_value in zip(
                new_drawable.values,
                existing.values,
                strict=True,
            )
        )
        root_counts = {len(row) for row in rows}
        if len(root_counts) != 1 or not rows[0]:
            continue
        for root_index in range(len(rows[0])):
            point = PointBundle(
                f"{new_drawable.drawable_id}&{existing.drawable_id}"
                f"[{root_index}]",
                tuple(row[root_index] for row in rows),
            )
            signature = point_signature(point)
            if signature in known_point_signatures:
                continue
            known_point_signatures.add(signature)
            new_points.append(point)
    return State(
        points=state.points + tuple(new_points),
        drawables=state.drawables + (new_drawable,),
        program=state.program + (candidate.describe(),),
    )


def is_target_line(drawable: DrawableBundle) -> bool:
    if drawable.kind != "line":
        return False
    for value, sample in zip(drawable.values, SAMPLES, strict=True):
        target = (0.0, 1.0, -sample.y3)
        if max(abs(left - right) for left, right in zip(value, target)) > 1e-7:
            return False
    return True


def is_target_chord_circle(drawable: DrawableBundle) -> bool:
    if drawable.kind != "circle":
        return False
    for value, sample in zip(drawable.values, SAMPLES, strict=True):
        center_x, center_y, radius_squared = value
        if abs(center_x - sample.x3) > 1e-7:
            return False
        center_delta_squared = (center_y - sample.y3) ** 2
        if (
            abs(
                center_delta_squared
                - radius_squared
                + sample.radius_3**2
            )
            > 1e-6
        ):
            return False
        if abs(center_y - sample.y3) <= 1e-7:
            return False
    return True


def main() -> None:
    start = monotonic()
    initial = initial_state()
    first_candidates = generate_candidates(initial, 1)
    states_after_one = 0
    states_after_two = 0
    third_candidates = 0
    hits: list[tuple[str, ...]] = []

    for first_index, first in enumerate(first_candidates, start=1):
        state_one = apply_candidate(initial, first, 1)
        states_after_one += 1
        for second in generate_candidates(state_one, 2):
            state_two = apply_candidate(state_one, second, 2)
            states_after_two += 1
            for third in generate_candidates(state_two, 3):
                third_candidates += 1
                if is_target_line(third.drawable):
                    hits.append(state_two.program + (third.describe(),))
                elif is_target_chord_circle(third.drawable):
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
            "samples": len(SAMPLES),
            "first_candidates": len(first_candidates),
            "states_after_one": states_after_one,
            "states_after_two": states_after_two,
            "third_candidates": third_candidates,
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for hit in hits:
        print("candidate")
        for move, description in enumerate(hit, start=1):
            print(f"  {move}. {description}")


if __name__ == "__main__":
    main()
