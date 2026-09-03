"""有界精确扫描 Mannheim 五线块的构造特有退化。

扫描使用水平公共直径 ``O1O2``，枚举整数圆心与整数半径。输入域 ``D8``
以及所有入射、合并、平行和相切判定均为精确运算；主方向失败时再检查
另外两条圆心方向。它用于发现参数分支和反例夹具，不是对连续参数域的
完备证明。
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import gcd, isqrt, lcm

from check_mannheim_degenerate_fixture import (
    Quadratic,
    add,
    determinant,
    dot,
    line_intersection,
    line_through,
    multiply,
    same_line,
    second_circle_point,
    solve_dot_system,
    subtract,
)


F = Fraction
Point = tuple[Fraction, Fraction]
Line = tuple[Fraction, Fraction, Fraction]

ROLE_KEYS = {
    "P0": ("alphaA", "aB", "a1A", "alpha1B"),
    "P1": ("alphaB", "aA", "a1B", "alpha1A"),
    "P2": ("alphaB", "aA", "alpha1B", "a1A"),
    "P3": ("alphaA", "aB", "alpha1A", "a1B"),
}

ROLE_NAMES = ("alpha2", "a2", "a2_prime", "alpha2_prime")

CHORD_ROLES = {
    "K_left": ("a2", "alpha2_prime"),
    "K_right": ("alpha2", "a2_prime"),
    "Kp_left": ("a2", "alpha2"),
    "Kp_right": ("a2_prime", "alpha2_prime"),
}


def canonical_line(line: Line) -> tuple[int, int, int]:
    denominator = lcm(*(item.denominator for item in line))
    values = tuple(int(item * denominator) for item in line)
    divisor = gcd(gcd(abs(values[0]), abs(values[1])), abs(values[2]))
    normalized = tuple(item // divisor for item in values)
    first_nonzero = next(item for item in normalized if item != 0)
    if first_nonzero < 0:
        normalized = tuple(-item for item in normalized)
    return normalized  # type: ignore[return-value]


def is_d8(
    centers: tuple[Point, Point, Point],
    radii: tuple[Fraction, ...],
) -> bool:
    o1, o2, o3 = centers
    r1, r2, r3 = radii
    if not r1 > r2 > r3 > 0:
        return False
    e2 = subtract(o2, o1)
    e3 = subtract(o3, o1)
    if determinant(e2, e3) == 0:
        return False
    distance_squares = (
        dot(e2, e2),
        dot(subtract(o3, o2), subtract(o3, o2)),
        dot(e3, e3),
    )
    if len(set(distance_squares)) != 3:
        return False
    for distance_squared, radius_sum in zip(
        distance_squares,
        (r1 + r2, r2 + r3, r3 + r1),
        strict=True,
    ):
        if distance_squared <= radius_sum**2:
            return False

    b = solve_dot_system(
        e2,
        e3,
        (dot(e2, e2) - r2**2 + r1**2) / 2,
        (dot(e3, e3) - r3**2 + r1**2) / 2,
    )
    constant = dot(b, b) - r1**2
    for sigma2 in (1, -1):
        for sigma3 in (1, -1):
            sigma = (1, sigma2, sigma3)
            u_vector = solve_dot_system(
                e2,
                e3,
                r1 - sigma2 * r2,
                r1 - sigma3 * r3,
            )
            quadratic = dot(u_vector, u_vector) - 1
            linear = 2 * dot(b, u_vector) - 2 * r1
            discriminant = linear**2 - 4 * quadratic * constant
            if quadratic == 0 or discriminant <= 0:
                return False
            for root_sign in (-1, 1):
                root = Quadratic(
                    -linear / (2 * quadratic),
                    F(root_sign, 2) / quadratic,
                    discriminant,
                )
                for direction_sign, radius in zip(sigma, radii, strict=True):
                    if root.sign() * (root + direction_sign * radius).sign() <= 0:
                        return False
    return True


def make_fixture(
    r1: int,
    r2: int,
    r3: int,
    distance: int,
    u: int,
    v: int,
) -> tuple[tuple[Point, Point, Point], tuple[Fraction, ...]]:
    return (
        ((F(0), F(0)), (F(distance), F(0)), (F(u), F(v))),
        (F(r1), F(r2), F(r3)),
    )


def rational_square_root(value: Fraction) -> Fraction | None:
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if (
        numerator * numerator == value.numerator
        and denominator * denominator == value.denominator
    ):
        return F(numerator, denominator)
    return None


def lift(point: Point, discriminant: Fraction):
    return (
        Quadratic(point[0], 0, discriminant),
        Quadratic(point[1], 0, discriminant),
    )


def orientation_failures(
    centers: tuple[Point, Point, Point],
    radii: tuple[Fraction, ...],
    order: tuple[int, int, int],
) -> set[str]:
    first_index, second_index, third_index = order
    delta = subtract(centers[second_index], centers[first_index])
    distance_squared = dot(delta, delta)
    distance = rational_square_root(distance_squared)
    if distance is None:
        direction = (
            Quadratic(0, delta[0] / distance_squared, distance_squared),
            Quadratic(0, delta[1] / distance_squared, distance_squared),
        )
        local_centers = tuple(lift(point, distance_squared) for point in centers)
    else:
        direction = multiply(1 / distance, delta)
        local_centers = centers

    first_center = local_centers[first_index]
    second_center = local_centers[second_index]
    third_center = local_centers[third_index]
    first_radius = radii[first_index]
    second_radius = radii[second_index]
    third_radius = radii[third_index]
    named_points = {
        "alpha": add(first_center, multiply(-first_radius, direction)),
        "a": add(first_center, multiply(first_radius, direction)),
        "a1": add(second_center, multiply(-second_radius, direction)),
        "alpha1": add(second_center, multiply(second_radius, direction)),
        "A": add(third_center, multiply(-third_radius, direction)),
        "B": add(third_center, multiply(third_radius, direction)),
    }
    batch_points = {
        source_name + endpoint_name: second_circle_point(
            named_points[source_name],
            named_points[endpoint_name],
            (third_center, third_radius),
            allow_tangent_alias=True,
        )
        for source_name in ("alpha", "a", "a1", "alpha1")
        for endpoint_name in ("A", "B")
    }

    failures: set[str] = set()
    for class_id, role_keys in ROLE_KEYS.items():
        roles = dict(
            zip(
                ROLE_NAMES,
                (batch_points[key] for key in role_keys),
                strict=True,
            )
        )
        chords = {}
        for chord_name, (first_name, second_name) in CHORD_ROLES.items():
            if roles[first_name] == roles[second_name]:
                failures.add(f"{class_id}:undefined:{chord_name}")
                continue
            chords[chord_name] = line_through(roles[first_name], roles[second_name])

        intersections = {}
        for point_name, left_name, right_name in (
            ("K", "K_left", "K_right"),
            ("Kp", "Kp_left", "Kp_right"),
        ):
            if left_name not in chords or right_name not in chords:
                continue
            left = chords[left_name]
            right = chords[right_name]
            if same_line(left, right):
                failures.add(f"{class_id}:coincident:{point_name}")
            elif determinant((left[0], left[1]), (right[0], right[1])) == 0:
                failures.add(f"{class_id}:parallel:{point_name}")
            else:
                intersections[point_name] = line_intersection(left, right)
        if set(intersections) != {"K", "Kp"}:
            continue
        if intersections["K"] == intersections["Kp"]:
            failures.add(f"{class_id}:collapse:K=Kp")
            continue
        tau = line_through(intersections["K"], intersections["Kp"])
        a, b, c = tau
        center_value = a * third_center[0] + b * third_center[1] + c
        margin = (
            third_radius**2 * (a * a + b * b)
            - center_value * center_value
        )
        margin_sign = (
            margin.sign()
            if isinstance(margin, Quadratic)
            else (margin > 0) - (margin < 0)
        )
        if margin_sign == 0:
            failures.add(f"{class_id}:tau_tangent")
        elif margin_sign < 0:
            failures.add(f"{class_id}:tau_disjoint")
    return failures


def analyze_fixture(
    centers: tuple[Point, Point, Point],
    radii: tuple[Fraction, ...],
) -> set[str]:
    o1, o2, o3 = centers
    r1, r2, r3 = radii
    direction = (F(1), F(0))
    named_points = {
        "alpha": add(o1, multiply(-r1, direction)),
        "a": add(o1, multiply(r1, direction)),
        "a1": add(o2, multiply(-r2, direction)),
        "alpha1": add(o2, multiply(r2, direction)),
        "A": add(o3, multiply(-r3, direction)),
        "B": add(o3, multiply(r3, direction)),
    }
    batch_points: dict[str, Point] = {}
    batch_lines: dict[str, Line] = {}
    events: set[str] = set()
    for source_name in ("alpha", "a", "a1", "alpha1"):
        for endpoint_name in ("A", "B"):
            key = source_name + endpoint_name
            source = named_points[source_name]
            endpoint = named_points[endpoint_name]
            batch_lines[key] = line_through(source, endpoint)
            second = second_circle_point(
                source,
                endpoint,
                (o3, r3),
                allow_tangent_alias=True,
            )
            batch_points[key] = second
            if second == endpoint:
                events.add(f"batch_alias:{key}")

    batch_line_keys = {
        canonical_line(line): key for key, line in batch_lines.items()
    }
    for class_id, role_keys in ROLE_KEYS.items():
        roles = dict(
            zip(
                ROLE_NAMES,
                (batch_points[key] for key in role_keys),
                strict=True,
            )
        )
        for first_index, first_name in enumerate(ROLE_NAMES):
            for second_name in ROLE_NAMES[first_index + 1 :]:
                if roles[first_name] == roles[second_name]:
                    events.add(f"{class_id}:merge:{first_name}={second_name}")

        chords: dict[str, Line] = {}
        for chord_name, (first_name, second_name) in CHORD_ROLES.items():
            if roles[first_name] == roles[second_name]:
                events.add(f"{class_id}:undefined:{chord_name}")
                continue
            chord = line_through(roles[first_name], roles[second_name])
            chords[chord_name] = chord
            canonical = canonical_line(chord)
            if canonical in batch_line_keys:
                events.add(
                    f"{class_id}:reuse:{chord_name}={batch_line_keys[canonical]}"
                )

        intersections: dict[str, Point] = {}
        for point_name, left_name, right_name in (
            ("K", "K_left", "K_right"),
            ("Kp", "Kp_left", "Kp_right"),
        ):
            if left_name not in chords or right_name not in chords:
                continue
            left = chords[left_name]
            right = chords[right_name]
            if same_line(left, right):
                events.add(f"{class_id}:coincident:{point_name}")
            elif determinant((left[0], left[1]), (right[0], right[1])) == 0:
                events.add(f"{class_id}:parallel:{point_name}")
            else:
                intersections[point_name] = line_intersection(left, right)

        if set(intersections) != {"K", "Kp"}:
            continue
        if intersections["K"] == intersections["Kp"]:
            events.add(f"{class_id}:collapse:K=Kp")
            continue
        tau = line_through(intersections["K"], intersections["Kp"])
        tau_key = canonical_line(tau)
        if tau_key in batch_line_keys:
            events.add(f"{class_id}:reuse:tau={batch_line_keys[tau_key]}")
        if tau_key in {canonical_line(line) for line in chords.values()}:
            events.add(f"{class_id}:reuse:tau=chord")
        a, b, c = tau
        margin = (
            r3**2 * (a**2 + b**2)
            - (a * o3[0] + b * o3[1] + c) ** 2
        )
        if margin == 0:
            events.add(f"{class_id}:tau_tangent")
        elif margin < 0:
            events.add(f"{class_id}:tau_disjoint")
    return events


def scan(max_radius: int, max_coordinate: int, show_reuse: bool) -> None:
    tested = 0
    in_domain = 0
    event_counts: Counter[str] = Counter()
    first_fixtures: dict[
        str,
        tuple[tuple[Point, Point, Point], tuple[Fraction, ...]],
    ] = {}
    orientation_failure_counts: Counter[str] = Counter()
    first_orientation_failures: dict[
        str,
        tuple[tuple[Point, Point, Point], tuple[Fraction, ...]],
    ] = {}
    coverage_counts: Counter[str] = Counter()
    first_coverage_fixtures: dict[
        str,
        tuple[tuple[Point, Point, Point], tuple[Fraction, ...]],
    ] = {}
    primary_orientation = (0, 1, 2)
    alternate_orientations = {
        "13->2": (0, 2, 1),
        "23->1": (1, 2, 0),
    }
    for r1 in range(3, max_radius + 1):
        for r2 in range(2, r1):
            for r3 in range(1, r2):
                for distance in range(r1 + r2 + 1, max_coordinate + 1):
                    for u in range(-3, distance + 4):
                        for v in range(1, max_coordinate + 1):
                            tested += 1
                            fixture = make_fixture(r1, r2, r3, distance, u, v)
                            if not is_d8(*fixture):
                                continue
                            in_domain += 1
                            for event in analyze_fixture(*fixture):
                                event_counts[event] += 1
                                first_fixtures.setdefault(event, fixture)
                            primary_failures = orientation_failures(
                                *fixture,
                                primary_orientation,
                            )
                            regular_orientations = [] if primary_failures else ["12->3"]
                            for failure in primary_failures:
                                event = f"12->3:{failure}"
                                orientation_failure_counts[event] += 1
                                first_orientation_failures.setdefault(event, fixture)
                            if primary_failures:
                                for orientation_name, order in (
                                    alternate_orientations.items()
                                ):
                                    failures = orientation_failures(*fixture, order)
                                    if not failures:
                                        regular_orientations.append(orientation_name)
                                    for failure in failures:
                                        event = f"{orientation_name}:{failure}"
                                        orientation_failure_counts[event] += 1
                                        first_orientation_failures.setdefault(event, fixture)
                            coverage = ",".join(regular_orientations) or "none"
                            coverage_counts[coverage] += 1
                            first_coverage_fixtures.setdefault(coverage, fixture)

    print(
        "scan",
        {
            "tested": tested,
            "in_D8": in_domain,
            "events": sum(event_counts.values()),
        },
    )
    for event in sorted(event_counts):
        if not show_reuse and ":reuse:" in event:
            continue
        centers, radii = first_fixtures[event]
        print(
            event,
            {
                "count": event_counts[event],
                "centers": tuple(
                    tuple(str(value) for value in point)
                    for point in centers
                ),
                "radii": tuple(str(value) for value in radii),
            },
        )
    print("orientation_coverage")
    for coverage in sorted(coverage_counts):
        centers, radii = first_coverage_fixtures[coverage]
        print(
            coverage,
            {
                "count": coverage_counts[coverage],
                "centers": tuple(
                    tuple(str(value) for value in point)
                    for point in centers
                ),
                "radii": tuple(str(value) for value in radii),
            },
        )
    print("orientation_failures")
    for event in sorted(orientation_failure_counts):
        centers, radii = first_orientation_failures[event]
        print(
            event,
            {
                "count": orientation_failure_counts[event],
                "centers": tuple(
                    tuple(str(value) for value in point)
                    for point in centers
                ),
                "radii": tuple(str(value) for value in radii),
            },
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-radius", type=int, default=6)
    parser.add_argument("--max-coordinate", type=int, default=18)
    parser.add_argument("--show-reuse", action="store_true")
    arguments = parser.parse_args()
    scan(arguments.max_radius, arguments.max_coordinate, arguments.show_reuse)


if __name__ == "__main__":
    main()
