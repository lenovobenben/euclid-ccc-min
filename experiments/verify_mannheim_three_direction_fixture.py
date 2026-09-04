"""精确验证一个三条圆心方向都发生 Mannheim 角色合并的 D8 夹具。

三个圆的半径为 ``(3, 2, 1)``，圆心三角形边长 ``(a, b, c)`` 由三个
整数多项式在有理盒中的唯一公共零点定义。脚本使用有理区间 Krawczyk
算子证明该零点存在且唯一，再用有理区间验证输入域 ``D8`` 的全部严格
条件。没有把打印的小数近似用于任何正确性判定。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


F = Fraction


@dataclass(frozen=True, slots=True)
class Interval:
    low: Fraction
    high: Fraction

    def __post_init__(self) -> None:
        if self.low > self.high:
            raise ValueError("区间端点逆序")

    @classmethod
    def point(cls, value: int | Fraction) -> Interval:
        rational = F(value)
        return cls(rational, rational)

    @staticmethod
    def coerce(value: int | Fraction | Interval) -> Interval:
        return value if isinstance(value, Interval) else Interval.point(value)

    def __add__(self, other: int | Fraction | Interval) -> Interval:
        value = self.coerce(other)
        return Interval(self.low + value.low, self.high + value.high)

    __radd__ = __add__

    def __neg__(self) -> Interval:
        return Interval(-self.high, -self.low)

    def __sub__(self, other: int | Fraction | Interval) -> Interval:
        return self + (-self.coerce(other))

    def __rsub__(self, other: int | Fraction | Interval) -> Interval:
        return self.coerce(other) - self

    def __mul__(self, other: int | Fraction | Interval) -> Interval:
        value = self.coerce(other)
        products = (
            self.low * value.low,
            self.low * value.high,
            self.high * value.low,
            self.high * value.high,
        )
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def reciprocal(self) -> Interval:
        if self.low <= 0 <= self.high:
            raise ZeroDivisionError("区间包含零")
        return Interval(1 / self.high, 1 / self.low)

    def __truediv__(self, other: int | Fraction | Interval) -> Interval:
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: int | Fraction | Interval) -> Interval:
        return self.coerce(other) / self

    def __pow__(self, exponent: int) -> Interval:
        if exponent < 0:
            return (self.reciprocal()) ** (-exponent)
        result = Interval.point(1)
        for _ in range(exponent):
            result *= self
        return result

    def midpoint(self) -> Fraction:
        return (self.low + self.high) / 2

    def strictly_positive(self) -> bool:
        return self.low > 0

    def strictly_negative(self) -> bool:
        return self.high < 0

    def contains_zero(self) -> bool:
        return self.low <= 0 <= self.high

    def inside(self, other: Interval) -> bool:
        return other.low < self.low and self.high < other.high

    def approximate(self) -> tuple[float, float]:
        return (float(self.low), float(self.high))


def square_root_bounds(value: Interval) -> Interval:
    if value.low <= 0:
        raise ValueError("只对严格正区间开平方")
    scale = 10**30

    def lower_bound(rational: Fraction) -> Fraction:
        scaled_square = rational.numerator * scale**2 // rational.denominator
        return F(isqrt(scaled_square), scale)

    low = lower_bound(value.low)
    high_floor = lower_bound(value.high)
    if high_floor * high_floor == value.high:
        high = high_floor
    else:
        high = high_floor + F(1, scale)
    return Interval(low, high)


def fixture_polynomials(a, b, c):
    return (
        a**3 - 3 * a**2 - a * b**2 - a * c**2 + 4 * a + b**2 - c**2,
        (
            a**2 * b
            + 2 * a**2
            - b**3
            - 8 * b**2
            + b * c**2
            - 30 * b
            - 2 * c**2
        ),
        (
            a**2 * c
            + a**2
            + b**2 * c
            - b**2
            - c**3
            - 9 * c**2
            - 40 * c
        ),
    )


def fixture_jacobian(a, b, c):
    return (
        (
            3 * a**2 - 6 * a - b**2 - c**2 + 4,
            2 * b * (1 - a),
            -2 * c * (a + 1),
        ),
        (
            2 * a * (b + 2),
            a**2 - 3 * b**2 - 16 * b + c**2 - 30,
            2 * c * (b - 2),
        ),
        (
            2 * a * (c + 1),
            2 * b * (c - 1),
            a**2 + b**2 - 3 * c**2 - 18 * c - 40,
        ),
    )


def invert_matrix(matrix: tuple[tuple[Fraction, ...], ...]):
    size = len(matrix)
    augmented = [
        [F(value) for value in row]
        + [F(row_index == column_index) for column_index in range(size)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            row
            for row in range(column, size)
            if augmented[row][column] != 0
        )
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        pivot_value = augmented[column][column]
        augmented[column] = [
            value / pivot_value for value in augmented[column]
        ]
        for row in range(size):
            if row == column:
                continue
            multiplier = augmented[row][column]
            augmented[row] = [
                augmented[row][index]
                - multiplier * augmented[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in augmented)


def sum_intervals(values) -> Interval:
    return sum(values, start=Interval.point(0))


def verify_unique_root(box: tuple[Interval, Interval, Interval]):
    midpoint = tuple(interval.midpoint() for interval in box)
    midpoint_value = fixture_polynomials(*midpoint)
    inverse = invert_matrix(fixture_jacobian(*midpoint))
    interval_jacobian = fixture_jacobian(*box)

    center = tuple(
        midpoint[row]
        - sum(
            inverse[row][column] * midpoint_value[column]
            for column in range(3)
        )
        for row in range(3)
    )
    error_matrix = tuple(
        tuple(
            Interval.point(row == column)
            - sum_intervals(
                inverse[row][inner] * interval_jacobian[inner][column]
                for inner in range(3)
            )
            for column in range(3)
        )
        for row in range(3)
    )
    displacement = tuple(
        Interval(interval.low - point, interval.high - point)
        for interval, point in zip(box, midpoint, strict=True)
    )
    image = tuple(
        Interval.point(center[row])
        + sum_intervals(
            error_matrix[row][column] * displacement[column]
            for column in range(3)
        )
        for row in range(3)
    )
    if not all(
        image_interval.inside(box_interval)
        for image_interval, box_interval in zip(image, box, strict=True)
    ):
        raise AssertionError("Krawczyk 像没有严格落在隔离盒内部")
    return image


def dot_norm(coefficients, x, y, cross):
    first, second = coefficients
    return (
        first**2 * x
        + 2 * first * second * cross
        + second**2 * y
    )


def solve_dot_system(x, y, cross, first_value, second_value):
    denominator = x * y - cross**2
    return (
        (first_value * y - second_value * cross) / denominator,
        (second_value * x - first_value * cross) / denominator,
    )


def dot_of_vectors(first, second, x, y, cross):
    return (
        first[0] * (second[0] * x + second[1] * cross)
        + first[1] * (second[0] * cross + second[1] * y)
    )


def verify_d8(box: tuple[Interval, Interval, Interval]):
    a, b, c = box
    radius_1 = F(3)
    radius_2 = F(2)
    radius_3 = F(1)
    if not (a - 5).strictly_positive():
        raise AssertionError("Gamma1 与 Gamma2 不是严格外离")
    if not (b - 4).strictly_positive():
        raise AssertionError("Gamma1 与 Gamma3 不是严格外离")
    if not (c - 3).strictly_positive():
        raise AssertionError("Gamma2 与 Gamma3 不是严格外离")
    for first, second in ((a, b), (a, c), (b, c)):
        if (first - second).contains_zero():
            raise AssertionError("圆心距离的隔离区间发生重叠")

    x = a**2
    y = b**2
    z = c**2
    cross = (x + y - z) / 2
    gram_determinant = x * y - cross**2
    if not gram_determinant.strictly_positive():
        raise AssertionError("圆心三角形退化")

    b_vector = solve_dot_system(
        x,
        y,
        cross,
        (x - radius_2**2 + radius_1**2) / 2,
        (y - radius_3**2 + radius_1**2) / 2,
    )
    constant = dot_norm(b_vector, x, y, cross) - radius_1**2
    reports = {}
    for sigma_2 in (1, -1):
        for sigma_3 in (1, -1):
            sigma = (1, sigma_2, sigma_3)
            u_vector = solve_dot_system(
                x,
                y,
                cross,
                radius_1 - sigma_2 * radius_2,
                radius_1 - sigma_3 * radius_3,
            )
            quadratic = dot_norm(u_vector, x, y, cross) - 1
            linear = (
                2 * dot_of_vectors(b_vector, u_vector, x, y, cross)
                - 2 * radius_1
            )
            discriminant = linear**2 - 4 * quadratic * constant
            if quadratic.contains_zero():
                raise AssertionError(f"{sigma} 的二次项可能为零")
            if not discriminant.strictly_positive():
                raise AssertionError(f"{sigma} 的判别式未证严格为正")
            square_root = square_root_bounds(discriminant)
            roots = (
                (-linear - square_root) / (2 * quadratic),
                (-linear + square_root) / (2 * quadratic),
            )
            for root in roots:
                for direction, radius in zip(
                    sigma,
                    (radius_1, radius_2, radius_3),
                    strict=True,
                ):
                    product = root * (root + direction * radius)
                    if not product.strictly_positive():
                        raise AssertionError(
                            f"{sigma} 的根没有满足物理半径符号条件"
                        )
            reports[sigma] = {
                "A": quadratic.approximate(),
                "Delta": discriminant.approximate(),
                "roots": tuple(root.approximate() for root in roots),
            }
    return gram_determinant, reports


def verify_merge_denominators(box: tuple[Interval, Interval, Interval]):
    a, b, c = box

    def projection(distance, first_to_third, second_to_third):
        return (
            distance**2
            + first_to_third**2
            - second_to_third**2
        ) / (2 * distance)

    projection_12 = projection(a, b, c)
    projection_13 = projection(b, a, c)
    projection_23 = projection(c, a, b)
    projection_21 = projection(a, c, b)
    projection_31 = projection(b, c, a)
    projection_32 = projection(c, b, a)
    denominators = {
        "12->3 P0 a2=a2_prime": a - 2 - projection_12 + 1,
        "21->3 P0 a2=a2_prime": a - 3 - projection_21 + 1,
        "13->2 P1 alpha2=alpha2_prime": b + 1 - projection_13 + 2,
        "31->2 P1 alpha2=alpha2_prime": b + 3 - projection_31 + 2,
        "23->1 P1 alpha2=alpha2_prime": c + 1 - projection_23 + 3,
        "32->1 P1 alpha2=alpha2_prime": c + 2 - projection_32 + 3,
    }
    for name, denominator in denominators.items():
        if denominator.contains_zero():
            raise AssertionError(f"{name} 的有理参数分母可能为零")
    return denominators


def verify_other_roles_distinct(box: tuple[Interval, Interval, Interval]):
    a, b, c = box

    def class_roles(
        distance,
        first_to_third,
        second_to_third,
        first_radius: int,
        second_radius: int,
        third_radius: int,
        class_id: str,
    ):
        projection = (
            distance**2
            + first_to_third**2
            - second_to_third**2
        ) / (2 * distance)
        height_squared = first_to_third**2 - projection**2
        height = square_root_bounds(height_squared)

        def point_a(source):
            return -height / (source - projection + third_radius)

        def point_b(source):
            return (source - projection - third_radius) / height

        if class_id == "P0":
            return {
                "a2": point_b(first_radius),
                "alpha2": point_a(-first_radius),
                "a2_prime": point_a(distance - second_radius),
                "alpha2_prime": point_b(distance + second_radius),
            }
        if class_id == "P1":
            return {
                "a2": point_a(first_radius),
                "alpha2": point_b(-first_radius),
                "a2_prime": point_b(distance - second_radius),
                "alpha2_prime": point_a(distance + second_radius),
            }
        raise ValueError(f"未实现方向类 {class_id}")

    cases = {
        "12->3": (
            class_roles(a, b, c, 3, 2, 1, "P0"),
            frozenset(("a2", "a2_prime")),
        ),
        "21->3": (
            class_roles(a, c, b, 2, 3, 1, "P0"),
            frozenset(("a2", "a2_prime")),
        ),
        "13->2": (
            class_roles(b, a, c, 3, 1, 2, "P1"),
            frozenset(("alpha2", "alpha2_prime")),
        ),
        "31->2": (
            class_roles(b, c, a, 1, 3, 2, "P1"),
            frozenset(("alpha2", "alpha2_prime")),
        ),
        "23->1": (
            class_roles(c, a, b, 2, 1, 3, "P1"),
            frozenset(("alpha2", "alpha2_prime")),
        ),
        "32->1": (
            class_roles(c, b, a, 1, 2, 3, "P1"),
            frozenset(("alpha2", "alpha2_prime")),
        ),
    }
    for orientation, (roles, merged_pair) in cases.items():
        names = tuple(roles)
        for first_index, first_name in enumerate(names):
            for second_name in names[first_index + 1 :]:
                pair = frozenset((first_name, second_name))
                difference = roles[first_name] - roles[second_name]
                if pair == merged_pair:
                    if not difference.contains_zero():
                        raise AssertionError(
                            f"{orientation} 的指定合并角色区间不相交"
                        )
                elif difference.contains_zero():
                    raise AssertionError(
                        f"{orientation} 的其它角色未证为不同点"
                    )
    return {
        orientation: {
            name: value.approximate() for name, value in roles.items()
        }
        for orientation, (roles, _) in cases.items()
    }


def main() -> None:
    scale = 10**6
    box = (
        Interval(F(7_796_641, scale), F(7_796_642, scale)),
        Interval(F(5_020_202, scale), F(5_020_203, scale)),
        Interval(F(4_149_601, scale), F(4_149_602, scale)),
    )
    image = verify_unique_root(box)
    gram_determinant, d8_reports = verify_d8(box)
    denominators = verify_merge_denominators(box)
    roles = verify_other_roles_distinct(box)
    print("root_box", tuple(interval.approximate() for interval in box))
    print("krawczyk_image", tuple(interval.approximate() for interval in image))
    print("gram_determinant", gram_determinant.approximate())
    print(
        "merge_denominators",
        {name: value.approximate() for name, value in denominators.items()},
    )
    print("role_parameters", roles)
    for sigma in sorted(d8_reports, reverse=True):
        print("D8", sigma, d8_reports[sigma])
    print(
        "three_direction_failures",
        {
            "12/21->3": "P0 a2=a2_prime, hence K=Kp",
            "13/31->2": "P1 alpha2=alpha2_prime, hence K=Kp",
            "23/32->1": "P1 alpha2=alpha2_prime, hence K=Kp",
        },
    )


if __name__ == "__main__":
    main()
