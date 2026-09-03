"""精确检查一个位于 D8 内、但固定 O1O2 Mannheim 分支退化的实例。

该夹具说明 D8 只排除了目标圆退化，不能自动保证某一固定 Mannheim
方向程序的全部中间对象非退化。所有判定都使用 Fraction。
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


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
        raise ValueError("批量线与第三圆相切，没有不同的第二交点")
    second = add(known, multiply(second_parameter, direction))
    if dot(subtract(second, center), subtract(second, center)) != radius**2:
        raise AssertionError("第二交点不在第三圆上")
    return second


def direction_points(
    order: tuple[int, int, int],
    direction: Point,
) -> tuple[dict[str, Point], Circle]:
    first_index, second_index, third_index = order
    first_center = CENTERS[first_index]
    second_center = CENTERS[second_index]
    third_center = CENTERS[third_index]
    first_radius = RADII[first_index]
    second_radius = RADII[second_index]
    third_radius = RADII[third_index]
    if dot(direction, direction) != 1:
        raise AssertionError("公共方向不是单位向量")
    return (
        {
            "alpha": add(first_center, multiply(-first_radius, direction)),
            "a": add(first_center, multiply(first_radius, direction)),
            "a1": add(second_center, multiply(-second_radius, direction)),
            "alpha1": add(second_center, multiply(second_radius, direction)),
            "A": add(third_center, multiply(-third_radius, direction)),
            "B": add(third_center, multiply(third_radius, direction)),
        },
        (third_center, third_radius),
    )


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


def assert_all_center_directions_degenerate() -> None:
    first_points, first_third = direction_points(
        (0, 1, 2),
        (F(1), F(0)),
    )
    if (first_points["a"], first_points["B"]) != (
        (F(10), F(0)),
        (F(10), F(12)),
    ):
        raise AssertionError("O1O2 方向的退化点与夹具声明不符")
    assert_tangent_batch(first_points, first_third, "aB")

    second_points, second_third = direction_points(
        (0, 2, 1),
        (F(3, 5), F(4, 5)),
    )
    if (second_points["a"], second_points["A"]) != (
        (F(6), F(8)),
        (F(22), F(-4)),
    ):
        raise AssertionError("O1O3 方向的退化点与夹具声明不符")
    assert_tangent_batch(second_points, second_third, "aA")

    third_points, third_circle = direction_points(
        (1, 2, 0),
        (F(-4, 5), F(3, 5)),
    )
    batch_points = {
        source + end: second_circle_point(
            third_points[source],
            third_points[end],
            third_circle,
        )
        for source in ("alpha", "a", "a1", "alpha1")
        for end in ("A", "B")
    }
    repeated = (F(154, 17), F(72, 17))
    if batch_points["aB"] != repeated or batch_points["a1A"] != repeated:
        raise AssertionError("O2O3 方向的两个第二交点没有按预期重合")

    alpha2 = batch_points["alphaA"]
    a2 = batch_points["aB"]
    a2_prime = batch_points["a1A"]
    alpha2_prime = batch_points["alpha1B"]
    k = line_intersection(
        line_through(a2, alpha2_prime),
        line_through(alpha2, a2_prime),
    )
    k_prime = line_intersection(
        line_through(a2, alpha2),
        line_through(a2_prime, alpha2_prime),
    )
    if k != repeated or k_prime != repeated:
        raise AssertionError("O2O3 方向的 P0 没有按预期退化为 K=K'")
    try:
        line_through(k, k_prime)
    except ValueError:
        return
    raise AssertionError("O2O3 方向不应能画出 P0 的 tau")


def main() -> None:
    coefficients = assert_d8()
    assert_all_center_directions_degenerate()
    print(
        "fixture",
        {
            "in_D8": True,
            "O1O2_branch": "tangent at aB",
            "O1O3_branch": "tangent at aA",
            "O2O3_branch": "P0 has K=K'",
        },
    )
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
