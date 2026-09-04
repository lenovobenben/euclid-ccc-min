"""精确检查 D8 内 Mannheim 批量线相切和定心线重合的实例。

批量线相切时，反对应点与已知端点重合，但四个五线块仍可继续。某个
目标的两条定心半径重合时，改用另一个输入圆的接触半径即可恢复圆心。
两种分支都不增加 E。所有正确性判定都使用 Fraction 和实二次扩域。
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import isqrt

from replay_mannheim_fixed import Quadratic


F = Fraction
Point = tuple[Fraction, Fraction]
Line = tuple[Fraction, Fraction, Fraction]
Circle = tuple[Point, Fraction]

CENTERS: tuple[Point, Point, Point] = (
    (F(0), F(0)),
    (F(25), F(0)),
    (F(9), F(12)),
)
RADII = (F(10), F(5), F(1))

ROOT_INTERVALS = {
    (1, 1, 1): (
        (F(-2501, 100), F(-2499, 100)),
        (F(699, 100), F(701, 100)),
    ),
    (1, 1, -1): (
        (F(-3451, 100), F(-3450, 100)),
        (F(817, 100), F(818, 100)),
    ),
    (1, -1, 1): (
        (F(-1515, 100), F(-1514, 100)),
        (F(1260, 100), F(1261, 100)),
    ),
    (1, -1, -1): (
        (F(-1626, 100), F(-1625, 100)),
        (F(1530, 100), F(1531, 100)),
    ),
}


def add(first: Point, second: Point) -> Point:
    return (first[0] + second[0], first[1] + second[1])


def subtract(first: Point, second: Point) -> Point:
    return (first[0] - second[0], first[1] - second[1])


def multiply(scalar: Fraction, point: Point) -> Point:
    return (scalar * point[0], scalar * point[1])


def dot(first: Point, second: Point) -> Fraction:
    return first[0] * second[0] + first[1] * second[1]


def determinant(first: Point, second: Point) -> Fraction:
    return first[0] * second[1] - first[1] * second[0]


def line_through(first: Point, second: Point) -> Line:
    line = (
        first[1] - second[1],
        second[0] - first[0],
        first[0] * second[1] - second[0] * first[1],
    )
    if line[0] == 0 and line[1] == 0:
        raise ValueError("不能经过同一点画线")
    return line


def line_intersection(first: Line, second: Line) -> Point:
    a, b, c = first
    d, e, f = second
    denominator = a * e - b * d
    if denominator == 0:
        raise ValueError("两条直线没有唯一有限交点")
    return (
        (b * f - c * e) / denominator,
        (c * d - a * f) / denominator,
    )


def solve_dot_system(
    first: Point,
    second: Point,
    first_value: Fraction,
    second_value: Fraction,
) -> Point:
    denominator = determinant(first, second)
    if denominator == 0:
        raise ValueError("点积方程组奇异")
    return (
        (
            first_value * second[1]
            - first[1] * second_value
        )
        / denominator,
        (
            first[0] * second_value
            - first_value * second[0]
        )
        / denominator,
    )


def assert_d8() -> dict[tuple[int, int, int], tuple[Fraction, Fraction, Fraction]]:
    o1, o2, o3 = CENTERS
    r1, r2, r3 = RADII
    e2 = subtract(o2, o1)
    e3 = subtract(o3, o1)
    distance_squares = (
        dot(e2, e2),
        dot(subtract(o3, o2), subtract(o3, o2)),
        dot(e3, e3),
    )

    if not r1 > r2 > r3 > 0:
        raise AssertionError("半径排序不成立")
    if determinant(e2, e3) == 0:
        raise AssertionError("三个圆心共线")
    if len(set(distance_squares)) != 3:
        raise AssertionError("三个圆心距离不全异")
    for distance_squared, radius_sum in zip(
        distance_squares,
        (r1 + r2, r2 + r3, r3 + r1),
        strict=True,
    ):
        if distance_squared <= radius_sum**2:
            raise AssertionError("输入圆不是两两严格外离")

    b = solve_dot_system(
        e2,
        e3,
        (dot(e2, e2) - r2**2 + r1**2) / 2,
        (dot(e3, e3) - r3**2 + r1**2) / 2,
    )
    constant = dot(b, b) - r1**2
    coefficients = {}
    for sigma2, sigma3 in product((1, -1), repeat=2):
        sigma = (1, sigma2, sigma3)
        u = solve_dot_system(
            e2,
            e3,
            r1 - sigma2 * r2,
            r1 - sigma3 * r3,
        )
        quadratic = dot(u, u) - 1
        linear = 2 * dot(b, u) - 2 * r1
        discriminant = linear**2 - 4 * quadratic * constant
        if quadratic == 0 or discriminant <= 0:
            raise AssertionError(f"{sigma} 没有两个不同有限实根")

        intervals = ROOT_INTERVALS[sigma]
        for left, right in intervals:
            left_value = quadratic * left**2 + linear * left + constant
            right_value = quadratic * right**2 + linear * right + constant
            if left_value * right_value >= 0:
                raise AssertionError(f"{sigma} 的隔离区间没有异号端点")
            for direction_sign, radius in zip(sigma, RADII, strict=True):
                first_zero = min(F(0), -direction_sign * radius)
                second_zero = max(F(0), -direction_sign * radius)
                if not (right < first_zero or left > second_zero):
                    raise AssertionError(
                        f"{sigma} 的根区间不能推出 D8 符号条件"
                    )
        coefficients[sigma] = (quadratic, linear, discriminant)
    return coefficients


def second_circle_point(
    source: Point,
    known: Point,
    circle: Circle,
    *,
    allow_tangent_alias: bool = False,
) -> Point:
    center, radius = circle
    direction = subtract(source, known)
    denominator = dot(direction, direction)
    if denominator == 0:
        raise ValueError("批量线的两个定义点重合")
    second_parameter = (
        -2 * dot(subtract(known, center), direction) / denominator
    )
    if second_parameter == 0:
        if allow_tangent_alias:
            return known
        raise ValueError("批量线与第三圆相切，没有不同的第二交点")
    second = add(known, multiply(second_parameter, direction))
    if dot(subtract(second, center), subtract(second, center)) != radius**2:
        raise AssertionError("第二交点不在第三圆上")
    return second


def assert_tangent_batch(
    points: dict[str, Point],
    third_circle: Circle,
    key: str,
) -> None:
    source_key = key[:-1]
    end_key = key[-1]
    source = points[source_key]
    endpoint = points[end_key]
    tangent = line_through(source, endpoint)
    line_a, line_b, line_c = tangent
    third_center, third_radius = third_circle
    distance_numerator = (
        line_a * third_center[0]
        + line_b * third_center[1]
        + line_c
    ) ** 2
    distance_denominator = line_a**2 + line_b**2
    if distance_numerator != third_radius**2 * distance_denominator:
        raise AssertionError(f"{key} 不是第三圆的切线")
    try:
        second_circle_point(source, endpoint, third_circle)
    except ValueError:
        return
    raise AssertionError(f"退化批量线 {key} 不应产生不同的第二交点")


def same_line(first: Line, second: Line) -> bool:
    a, b, c = first
    d, e, f = second
    return a * e == b * d and a * f == c * d and b * f == c * e


def build_primary_taus() -> dict[str, Line]:
    direction = (F(1), F(0))
    points = {
        "alpha": add(CENTERS[0], multiply(-RADII[0], direction)),
        "a": add(CENTERS[0], multiply(RADII[0], direction)),
        "a1": add(CENTERS[1], multiply(-RADII[1], direction)),
        "alpha1": add(CENTERS[1], multiply(RADII[1], direction)),
        "A": add(CENTERS[2], multiply(-RADII[2], direction)),
        "B": add(CENTERS[2], multiply(RADII[2], direction)),
    }
    third_circle = (CENTERS[2], RADII[2])
    if (points["a"], points["B"]) != (
        (F(10), F(0)),
        (F(10), F(12)),
    ):
        raise AssertionError("O1O2 方向的退化点与夹具声明不符")
    assert_tangent_batch(points, third_circle, "aB")
    batch_points = {
        source + end: second_circle_point(
            points[source],
            points[end],
            third_circle,
            allow_tangent_alias=True,
        )
        for source in ("alpha", "a", "a1", "alpha1")
        for end in ("A", "B")
    }
    if batch_points["aB"] != points["B"]:
        raise AssertionError("相切批量线没有把 a2 绑定为既有端点 B")

    roles = {
        "P0": ("alphaA", "aB", "a1A", "alpha1B"),
        "P1": ("alphaB", "aA", "a1B", "alpha1A"),
        "P2": ("alphaB", "aA", "alpha1B", "a1A"),
        "P3": ("alphaA", "aB", "alpha1A", "a1B"),
    }
    expected_taus = {
        "P0": (F(2), F(1), F(-31)),
        "P1": (F(19), F(12), F(-302)),
        "P2": (F(736), F(483), F(-11753)),
        "P3": (F(93), F(44), F(-1434)),
    }
    taus = {}
    for class_id, keys in roles.items():
        alpha2, a2, a2_prime, alpha2_prime = (
            batch_points[key] for key in keys
        )
        k = line_intersection(
            line_through(a2, alpha2_prime),
            line_through(alpha2, a2_prime),
        )
        k_prime = line_intersection(
            line_through(a2, alpha2),
            line_through(a2_prime, alpha2_prime),
        )
        tau = line_through(k, k_prime)
        if not same_line(tau, expected_taus[class_id]):
            raise AssertionError(f"{class_id} 的 tau 与夹具精确值不符")
        a, b, c = tau
        secant_margin = (
            RADII[2] ** 2 * (a**2 + b**2)
            - (a * CENTERS[2][0] + b * CENTERS[2][1] + c) ** 2
        )
        if secant_margin <= 0:
            raise AssertionError(f"{class_id} 的 tau 不割第三圆于两点")
        taus[class_id] = tau
    return taus


def rational_square_root(value: Fraction) -> Fraction | None:
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if (
        numerator * numerator == value.numerator
        and denominator * denominator == value.denominator
    ):
        return F(numerator, denominator)
    return None


def lift(point: Point, discriminant: Fraction | None):
    if discriminant is None:
        return point
    return (
        Quadratic(point[0], 0, discriminant),
        Quadratic(point[1], 0, discriminant),
    )


def similarity_center(first_index: int, second_index: int, kind: str) -> Point:
    first_center = CENTERS[first_index]
    second_center = CENTERS[second_index]
    first_radius = RADII[first_index]
    second_radius = RADII[second_index]
    if kind == "ext":
        denominator = first_radius - second_radius
        return (
            (
                first_radius * second_center[0]
                - second_radius * first_center[0]
            )
            / denominator,
            (
                first_radius * second_center[1]
                - second_radius * first_center[1]
            )
            / denominator,
        )
    denominator = first_radius + second_radius
    return (
        (
            first_radius * second_center[0]
            + second_radius * first_center[0]
        )
        / denominator,
        (
            first_radius * second_center[1]
            + second_radius * first_center[1]
        )
        / denominator,
    )


def assert_targets(
    coefficients: dict[
        tuple[int, int, int],
        tuple[Fraction, Fraction, Fraction],
    ],
    taus: dict[str, Line],
) -> dict[str, str]:
    e2 = subtract(CENTERS[1], CENTERS[0])
    e3 = subtract(CENTERS[2], CENTERS[0])
    b_vector = solve_dot_system(
        e2,
        e3,
        (dot(e2, e2) - RADII[1] ** 2 + RADII[0] ** 2) / 2,
        (dot(e3, e3) - RADII[2] ** 2 + RADII[0] ** 2) / 2,
    )
    class_by_sigma = {
        (1, 1, 1): "P0",
        (1, 1, -1): "P1",
        (1, -1, -1): "P2",
        (1, -1, 1): "P3",
    }
    recovery_pairs = {}
    physical_signs = set()
    for sigma, (quadratic, linear, discriminant) in coefficients.items():
        u_vector = solve_dot_system(
            e2,
            e3,
            RADII[0] - sigma[1] * RADII[1],
            RADII[0] - sigma[2] * RADII[2],
        )
        tau = taus[class_by_sigma[sigma]]
        for root_sign in (-1, 1):
            discriminant_root = rational_square_root(discriminant)
            if discriminant_root is None:
                t = Quadratic(
                    -linear / (2 * quadratic),
                    F(root_sign, 2) / quadratic,
                    discriminant,
                )
                t_sign = t.sign()
                extension_discriminant = discriminant
            else:
                t = (
                    -linear + root_sign * discriminant_root
                ) / (2 * quadratic)
                t_sign = (t > 0) - (t < 0)
                extension_discriminant = None
            if t_sign == 0:
                raise AssertionError("D8 根不应为零")
            radius = t if t_sign > 0 else -t
            physical = tuple(t_sign * item for item in sigma)
            sign_name = "".join("+" if item > 0 else "-" for item in physical)
            physical_signs.add(sign_name)
            center = add(
                lift(b_vector, extension_discriminant),
                multiply(t, lift(u_vector, extension_discriminant)),
            )
            contacts = []
            for index in range(3):
                input_center = lift(CENTERS[index], extension_discriminant)
                delta = subtract(center, input_center)
                center_distance = radius + physical[index] * RADII[index]
                contact = add(
                    input_center,
                    multiply(
                        physical[index] * RADII[index] / center_distance,
                        delta,
                    ),
                )
                if (
                    dot(
                        subtract(contact, input_center),
                        subtract(contact, input_center),
                    )
                    != RADII[index] ** 2
                ):
                    raise AssertionError(f"{sign_name} 的接触点不在输入圆上")
                expected_distance_squared = center_distance * center_distance
                if dot(delta, delta) != expected_distance_squared:
                    raise AssertionError(f"{sign_name} 的切触等式不成立")
                contacts.append(contact)

            contact_3 = contacts[2]
            if (
                tau[0] * contact_3[0]
                + tau[1] * contact_3[1]
                + tau[2]
                != 0
            ):
                raise AssertionError(f"{sign_name} 的第三圆接触点不在 tau 上")
            if (
                dot(subtract(center, contact_3), subtract(center, contact_3))
                != radius * radius
            ):
                raise AssertionError(f"{sign_name} 的输出圆半径错误")

            chosen_index = None
            for index in (1, 0):
                radial_determinant = determinant(
                    subtract(
                        contacts[index],
                        lift(CENTERS[index], extension_discriminant),
                    ),
                    subtract(
                        contact_3,
                        lift(CENTERS[2], extension_discriminant),
                    ),
                )
                if radial_determinant == 0:
                    continue
                kind = "ext" if physical[index] == physical[2] else "int"
                h = lift(
                    similarity_center(index, 2, kind),
                    extension_discriminant,
                )
                if determinant(
                    subtract(contacts[index], h),
                    subtract(contact_3, h),
                ) != 0:
                    raise AssertionError(
                        f"{sign_name} 的反对应接触点不经过相似中心"
                    )
                recovered = line_intersection(
                    line_through(
                        lift(CENTERS[index], extension_discriminant),
                        contacts[index],
                    ),
                    line_through(
                        lift(CENTERS[2], extension_discriminant),
                        contact_3,
                    ),
                )
                if recovered != center:
                    raise AssertionError(f"{sign_name} 没有恢复出目标圆心")
                chosen_index = index
                break
            if chosen_index is None:
                raise AssertionError(f"{sign_name} 的三条接触半径全部重合")
            recovery_pairs[sign_name] = f"Gamma{chosen_index + 1}/Gamma3"

    if physical_signs != {
        "+++",
        "++-",
        "+-+",
        "+--",
        "-++",
        "-+-",
        "--+",
        "---",
    }:
        raise AssertionError("没有精确恢复全部八种物理切向符号")
    expected_recovery_pairs = {
        "+++": "Gamma1/Gamma3",
        "++-": "Gamma2/Gamma3",
        "+-+": "Gamma2/Gamma3",
        "+--": "Gamma2/Gamma3",
        "-++": "Gamma2/Gamma3",
        "-+-": "Gamma2/Gamma3",
        "--+": "Gamma2/Gamma3",
        "---": "Gamma2/Gamma3",
    }
    if recovery_pairs != expected_recovery_pairs:
        raise AssertionError("定心圆对分支与夹具的精确结果不符")
    return recovery_pairs


def main() -> None:
    coefficients = assert_d8()
    taus = build_primary_taus()
    recovery_pairs = assert_targets(coefficients, taus)
    print(
        "fixture",
        {
            "in_D8": True,
            "batch_aB": "tangent alias B",
            "five_line_blocks": "four valid secants",
            "score": "18/65 E upper bound unchanged",
        },
    )
    print("recovery_pairs", dict(sorted(recovery_pairs.items())))
    for sigma in sorted(coefficients, reverse=True):
        quadratic, linear, discriminant = coefficients[sigma]
        print(
            "D8",
            sigma,
            {
                "A": str(quadratic),
                "B": str(linear),
                "Delta": str(discriminant),
                "root_intervals": tuple(
                    tuple(str(value) for value in interval)
                    for interval in ROOT_INTERVALS[sigma]
                ),
            },
        )


if __name__ == "__main__":
    main()
