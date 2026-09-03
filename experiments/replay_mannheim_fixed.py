"""精确重放校准实例上的 Mannheim 65 E 候选。

本脚本只验证固定实例

    Gamma1=((0, 0), 4)
    Gamma2=((13, 0), 2)
    Gamma3=((4, 15), 1)

它用有理数和每个解对各自的实二次扩域做精确运算，不把浮点残差当作
正确性证据。小数只用于人类可读输出；切触等式和切向符号都在二次域中
精确判定。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt, sqrt


F = Fraction
RationalPoint = tuple[Fraction, Fraction]
RationalLine = tuple[Fraction, Fraction, Fraction]
RationalCircle = tuple[RationalPoint, Fraction]


@dataclass(frozen=True, slots=True)
class Quadratic:
    """表示 ``a + b sqrt(D)``，其中 ``D`` 是正有理数。"""

    a: Fraction
    b: Fraction
    discriminant: Fraction

    def __init__(
        self,
        a: int | Fraction = 0,
        b: int | Fraction = 0,
        discriminant: int | Fraction = 0,
    ) -> None:
        object.__setattr__(self, "a", F(a))
        object.__setattr__(self, "b", F(b))
        object.__setattr__(self, "discriminant", F(discriminant))

    def _coerce(self, other: object) -> Quadratic:
        if isinstance(other, Quadratic):
            if other.discriminant != self.discriminant:
                raise ValueError("不能混合不同的二次扩域")
            return other
        if isinstance(other, (int, Fraction)):
            return Quadratic(other, 0, self.discriminant)
        return NotImplemented  # type: ignore[return-value]

    def __add__(self, other: object) -> Quadratic:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return NotImplemented
        return Quadratic(
            self.a + other_value.a,
            self.b + other_value.b,
            self.discriminant,
        )

    __radd__ = __add__

    def __neg__(self) -> Quadratic:
        return Quadratic(-self.a, -self.b, self.discriminant)

    def __sub__(self, other: object) -> Quadratic:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return NotImplemented
        return self + (-other_value)

    def __rsub__(self, other: object) -> Quadratic:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return NotImplemented
        return other_value - self

    def __mul__(self, other: object) -> Quadratic:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return NotImplemented
        return Quadratic(
            self.a * other_value.a
            + self.b * other_value.b * self.discriminant,
            self.a * other_value.b + self.b * other_value.a,
            self.discriminant,
        )

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> Quadratic:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return NotImplemented
        denominator = (
            other_value.a**2
            - other_value.b**2 * self.discriminant
        )
        if denominator == 0:
            raise ZeroDivisionError("二次域除数为零")
        return Quadratic(
            (
                self.a * other_value.a
                - self.b * other_value.b * self.discriminant
            )
            / denominator,
            (self.b * other_value.a - self.a * other_value.b)
            / denominator,
            self.discriminant,
        )

    def __rtruediv__(self, other: object) -> Quadratic:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return NotImplemented
        return other_value / self

    def __eq__(self, other: object) -> bool:
        other_value = self._coerce(other)
        if other_value is NotImplemented:
            return False
        return self.a == other_value.a and self.b == other_value.b

    def sign(self) -> int:
        """精确返回该实二次域元素的符号。"""

        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if self.b == 0 or (self.a > 0) == (self.b > 0):
            return (self.a > 0) - (self.a < 0)
        squared_difference = self.a**2 - self.b**2 * self.discriminant
        if squared_difference == 0:
            return 0
        rational_sign = (squared_difference > 0) - (squared_difference < 0)
        return rational_sign if self.a > 0 else -rational_sign

    def approximate(self) -> float:
        return float(self.a) + float(self.b) * sqrt(float(self.discriminant))


QuadraticPoint = tuple[Quadratic, Quadratic]
QuadraticLine = tuple[Quadratic, Quadratic, Quadratic]


def add(first, second):
    return (first[0] + second[0], first[1] + second[1])


def subtract(first, second):
    return (first[0] - second[0], first[1] - second[1])


def multiply(scalar, point):
    return (scalar * point[0], scalar * point[1])


def dot(first, second):
    return first[0] * second[0] + first[1] * second[1]


def line_through(first, second):
    line = (
        first[1] - second[1],
        second[0] - first[0],
        first[0] * second[1] - second[0] * first[1],
    )
    if line[0] == 0 and line[1] == 0:
        raise ValueError("不能经过同一点画线")
    return line


def line_intersection(first, second):
    a, b, c = first
    d, e, f = second
    denominator = a * e - b * d
    if denominator == 0:
        raise ValueError("两条直线没有唯一有限交点")
    return (
        (b * f - c * e) / denominator,
        (c * d - a * f) / denominator,
    )


def same_line(first: RationalLine, second: RationalLine) -> bool:
    """判断两组三元系数是否表示同一条射影直线。"""

    a, b, c = first
    d, e, f = second
    return a * e == b * d and a * f == c * d and b * f == c * e


def rational_circle(center: RationalPoint, through: RationalPoint) -> RationalCircle:
    delta = subtract(center, through)
    radius_squared = dot(delta, delta)
    if radius_squared == 0:
        raise ValueError("基础圆的圆心和圆上一点重合")
    return (center, radius_squared)


def fraction_square_root(value: Fraction) -> Fraction:
    if value < 0:
        raise ValueError("有理数不是非负平方")
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        raise ValueError(f"分子 {value.numerator} 不是平方数")
    if denominator * denominator != value.denominator:
        raise ValueError(f"分母 {value.denominator} 不是平方数")
    return F(numerator, denominator)


def rational_line_circle_intersections(
    line: RationalLine,
    circle: RationalCircle,
) -> tuple[RationalPoint, RationalPoint]:
    a, b, c = line
    center, radius_squared = circle
    if b != 0:
        base = (F(0), -c / b)
    else:
        base = (-c / a, F(0))
    direction = (b, -a)
    offset = subtract(base, center)
    quadratic_a = dot(direction, direction)
    quadratic_b = 2 * dot(offset, direction)
    quadratic_c = dot(offset, offset) - radius_squared
    root_discriminant = (
        quadratic_b**2 - 4 * quadratic_a * quadratic_c
    )
    root_square = fraction_square_root(root_discriminant)
    roots = (
        (-quadratic_b - root_square) / (2 * quadratic_a),
        (-quadratic_b + root_square) / (2 * quadratic_a),
    )
    points = tuple(add(base, multiply(root, direction)) for root in roots)
    return tuple(sorted(points))  # type: ignore[return-value]


def other_point(
    points: tuple[RationalPoint, RationalPoint],
    known: RationalPoint,
) -> RationalPoint:
    candidates = tuple(point for point in points if point != known)
    if len(candidates) != 1:
        raise ValueError("不能唯一绑定第二交点")
    return candidates[0]


def lift_point(point: RationalPoint, discriminant: Fraction) -> QuadraticPoint:
    return (
        Quadratic(point[0], 0, discriminant),
        Quadratic(point[1], 0, discriminant),
    )


def approximate_point(point: QuadraticPoint) -> tuple[float, float]:
    return (point[0].approximate(), point[1].approximate())


@dataclass(slots=True)
class DependencyGraph:
    dependencies: dict[str, tuple[str, ...]]
    paid_order: list[str]
    paid_kinds: dict[str, str]

    @classmethod
    def initial(cls) -> DependencyGraph:
        return cls(
            dependencies={
                "O1": (),
                "O2": (),
                "O3": (),
                "Gamma1": ("O1",),
                "Gamma2": ("O2",),
                "Gamma3": ("O3",),
            },
            paid_order=[],
            paid_kinds={},
        )

    def _add(self, node_id: str, dependencies: tuple[str, ...]) -> None:
        if node_id in self.dependencies:
            raise ValueError(f"重复节点 {node_id}")
        missing = tuple(
            item for item in dependencies if item not in self.dependencies
        )
        if missing:
            raise ValueError(f"节点 {node_id} 引用未知依赖 {missing}")
        self.dependencies[node_id] = dependencies

    def point(self, node_id: str, *dependencies: str) -> None:
        self._add(node_id, tuple(dependencies))

    def draw(self, node_id: str, kind: str, *dependencies: str) -> int:
        if kind not in {"line", "circle"}:
            raise ValueError(f"未知计费对象类型 {kind}")
        self._add(node_id, tuple(dependencies))
        self.paid_order.append(node_id)
        self.paid_kinds[node_id] = kind
        return len(self.paid_order)

    def paid_ancestors(self, root: str) -> frozenset[str]:
        live: set[str] = set()
        pending = [root]
        while pending:
            node_id = pending.pop()
            if node_id in live:
                continue
            live.add(node_id)
            pending.extend(self.dependencies[node_id])
        return frozenset(item for item in live if item in self.paid_kinds)


@dataclass(frozen=True, slots=True)
class TargetResult:
    sign: str
    center: QuadraticPoint
    radius_squared: Quadratic
    output_id: str
    draw_index: int


class MannheimReplay:
    def __init__(self) -> None:
        self.graph = DependencyGraph.initial()
        self.o1: RationalPoint = (F(0), F(0))
        self.o2: RationalPoint = (F(13), F(0))
        self.o3: RationalPoint = (F(4), F(15))
        self.gamma1 = (self.o1, F(16))
        self.gamma2 = (self.o2, F(4))
        self.gamma3 = (self.o3, F(1))
        self.points: dict[str, RationalPoint] = {}
        self.batch_lines: dict[str, RationalLine] = {}
        self.batch_points: dict[str, RationalPoint] = {}
        self.similarity_centers: dict[str, RationalPoint] = {}
        self.targets: dict[str, TargetResult] = {}

    def build_global_prefix(self) -> None:
        ell = line_through(self.o1, self.o2)
        self.graph.draw("ell", "line", "O1", "O2")
        alpha, a = rational_line_circle_intersections(ell, self.gamma1)
        a1, alpha1 = rational_line_circle_intersections(ell, self.gamma2)
        self.points.update(
            {"alpha": alpha, "a": a, "a1": a1, "alpha1": alpha1}
        )
        self.graph.point("alpha", "ell", "Gamma1")
        self.graph.point("a", "ell", "Gamma1")
        self.graph.point("a1", "ell", "Gamma2")
        self.graph.point("alpha1", "ell", "Gamma2")

        prefix_discriminant = F(241)
        x = Quadratic(0, 1, prefix_discriminant)
        x_point = (x, Quadratic(0, 0, prefix_discriminant))
        o1 = lift_point(self.o1, prefix_discriminant)
        o3 = lift_point(self.o3, prefix_discriminant)
        c0_radius_squared = dot(subtract(o1, o3), subtract(o1, o3))
        self.graph.draw("parallel_c0", "circle", "O1", "O3")
        self.graph.point("parallel_X", "parallel_c0", "ell")
        c_x_radius_squared = dot(
            subtract(x_point, o3),
            subtract(x_point, o3),
        )
        self.graph.draw("parallel_cX", "circle", "parallel_X", "O3")
        q = lift_point((F(4), F(-15)), prefix_discriminant)
        if dot(subtract(q, o1), subtract(q, o1)) != c0_radius_squared:
            raise AssertionError("Q 不在 parallel_c0 上")
        if dot(subtract(q, x_point), subtract(q, x_point)) != c_x_radius_squared:
            raise AssertionError("Q 不在 parallel_cX 上")
        self.graph.point("parallel_Q", "parallel_c0", "parallel_cX")
        self.graph.draw("parallel_diameter", "line", "parallel_X", "parallel_Q")
        r = (
            2 * x - q[0],
            2 * x_point[1] - q[1],
        )
        if dot(subtract(r, x_point), subtract(r, x_point)) != c_x_radius_squared:
            raise AssertionError("R 不在 parallel_cX 上")
        self.graph.point("parallel_R", "parallel_diameter", "parallel_cX")
        parallel = line_through(o3, r)
        if parallel[0] != 0 or parallel[1] == 0:
            raise AssertionError("4 E 前缀没有产生水平线")
        if parallel[2] / parallel[1] != -15:
            raise AssertionError("4 E 前缀的水平线不经过 O3")
        self.graph.draw("ell3", "line", "O3", "parallel_R")

        ell3: RationalLine = (F(0), F(1), F(-15))
        a3, b3 = rational_line_circle_intersections(ell3, self.gamma3)
        self.points.update({"A": a3, "B": b3})
        self.graph.point("A", "ell3", "Gamma3")
        self.graph.point("B", "ell3", "Gamma3")

    def draw_batch(self, key: str) -> None:
        if key in self.batch_lines:
            return
        if key.startswith("alpha1"):
            source_key = "alpha1"
        elif key.startswith("alpha"):
            source_key = "alpha"
        elif key.startswith("a1"):
            source_key = "a1"
        elif key.startswith("a"):
            source_key = "a"
        else:
            raise ValueError(f"未知批量线标签 {key}")
        end_key = key[-1]
        source = self.points[source_key]
        end = self.points[end_key]
        batch_line = line_through(source, end)
        line_id = f"batch_{key}"
        self.graph.draw(line_id, "line", source_key, end_key)
        second = other_point(
            rational_line_circle_intersections(batch_line, self.gamma3),
            end,
        )
        expected_points = {
            "alphaA": (F(460, 137), F(2160, 137)),
            "alphaB": (F(76, 17), F(240, 17)),
            "aA": (F(340, 113), F(1680, 113)),
            "aB": (F(564, 113), F(1680, 113)),
            "a1A": (F(995, 289), F(4095, 289)),
            "a1B": (F(137, 29), F(455, 29)),
            "alpha1A": (F(155, 41), F(575, 41)),
            "alpha1B": (F(57, 13), F(207, 13)),
        }
        if second != expected_points[key]:
            raise AssertionError(f"批量线 {key} 的第二交点与校准值不符")
        point_id = f"second_{key}"
        self.graph.point(point_id, line_id, "Gamma3")
        self.batch_lines[key] = batch_line
        self.batch_points[key] = second

    def bind_similarity_center(self, kind: str) -> RationalPoint:
        node_id = f"H23_{kind}"
        if kind in self.similarity_centers:
            return self.similarity_centers[kind]
        if kind == "ext":
            parents = ("a1A", "alpha1B")
            expected = (F(-5), F(30))
        else:
            parents = ("a1B", "alpha1A")
            expected = (F(7), F(10))
        center = line_intersection(
            self.batch_lines[parents[0]],
            self.batch_lines[parents[1]],
        )
        if center != expected:
            raise AssertionError(f"{kind} 相似中心与校准值不符")
        self.graph.point(
            node_id,
            f"batch_{parents[0]}",
            f"batch_{parents[1]}",
        )
        self.similarity_centers[kind] = center
        return center

    def build_pair(
        self,
        class_id: str,
        roles: tuple[str, str, str, str],
        sign_roots: tuple[tuple[str, int], tuple[str, int]],
    ) -> None:
        alpha2, a2, a2_prime, alpha2_prime = (
            self.batch_points[key] for key in roles
        )
        first = line_through(a2, alpha2_prime)
        second = line_through(alpha2, a2_prime)
        self.graph.draw(
            f"{class_id}_K_line_1",
            "line",
            f"second_{roles[1]}",
            f"second_{roles[3]}",
        )
        self.graph.draw(
            f"{class_id}_K_line_2",
            "line",
            f"second_{roles[0]}",
            f"second_{roles[2]}",
        )
        k = line_intersection(first, second)
        self.graph.point(
            f"{class_id}_K",
            f"{class_id}_K_line_1",
            f"{class_id}_K_line_2",
        )

        third = line_through(a2, alpha2)
        fourth = line_through(a2_prime, alpha2_prime)
        self.graph.draw(
            f"{class_id}_Kp_line_1",
            "line",
            f"second_{roles[1]}",
            f"second_{roles[0]}",
        )
        self.graph.draw(
            f"{class_id}_Kp_line_2",
            "line",
            f"second_{roles[2]}",
            f"second_{roles[3]}",
        )
        k_prime = line_intersection(third, fourth)
        self.graph.point(
            f"{class_id}_Kp",
            f"{class_id}_Kp_line_1",
            f"{class_id}_Kp_line_2",
        )
        tau = line_through(k, k_prime)
        expected_tau = {
            "P0": (F(3308), F(1095), F(-30340)),
            "P1": (F(116), F(45), F(-1108)),
            "P2": (F(3164), F(1335), F(-30940)),
            "P3": (F(1092), F(325), F(-9772)),
        }
        if not same_line(tau, expected_tau[class_id]):
            raise AssertionError(f"{class_id} 的 tau 与校准值不符")
        self.graph.draw(
            f"{class_id}_tau",
            "line",
            f"{class_id}_K",
            f"{class_id}_Kp",
        )

        a, b, c = tau
        if b == 0:
            raise AssertionError("校准实例中的 tau 不应为竖线")
        linear_y_constant = -c - 15 * b
        quadratic_a = b**2 + a**2
        quadratic_b = -8 * b**2 - 2 * a * linear_y_constant
        quadratic_c = 15 * b**2 + linear_y_constant**2
        discriminant = (
            quadratic_b**2 - 4 * quadratic_a * quadratic_c
        )
        if discriminant <= 0:
            raise AssertionError("tau 与 Gamma3 没有两个不同实交点")
        for sign, root_sign in sign_roots:
            x = Quadratic(
                -quadratic_b / (2 * quadratic_a),
                F(root_sign, 1) / (2 * quadratic_a),
                discriminant,
            )
            y = (-a * x - c) / b
            point_id = f"{class_id}_M3_{sign}"
            self.graph.point(point_id, f"{class_id}_tau", "Gamma3")
            self._build_target(class_id, sign, (x, y), point_id)

    def _build_target(
        self,
        class_id: str,
        sign: str,
        contact_3: QuadraticPoint,
        contact_3_id: str,
    ) -> None:
        discriminant = contact_3[0].discriminant
        same_direction = sign[1] == sign[2]
        similarity_kind = "ext" if same_direction else "int"
        h_rational = self.bind_similarity_center(similarity_kind)
        h = lift_point(h_rational, discriminant)
        o1 = lift_point(self.o1, discriminant)
        o2 = lift_point(self.o2, discriminant)
        o3 = lift_point(self.o3, discriminant)

        contact_line_id = f"{class_id}_{sign}_contact_line"
        self.graph.draw(
            contact_line_id,
            "line",
            f"H23_{similarity_kind}",
            contact_3_id,
        )
        direction = subtract(contact_3, h)
        offset = subtract(h, o2)
        quadratic_a = dot(direction, direction)
        quadratic_b = 2 * dot(offset, direction)
        homologue_root = F(2 if same_direction else -2)
        anti_root = -quadratic_b / quadratic_a - homologue_root
        contact_2 = add(h, multiply(anti_root, direction))
        if dot(subtract(contact_2, o2), subtract(contact_2, o2)) != 4:
            raise AssertionError(f"{sign} 的 Gamma2 接触点不在输入圆上")
        contact_2_id = f"{class_id}_M2_{sign}"
        self.graph.point(contact_2_id, contact_line_id, "Gamma2")

        radius_3_id = f"{class_id}_{sign}_radius_3"
        radius_2_id = f"{class_id}_{sign}_radius_2"
        self.graph.draw(radius_3_id, "line", "O3", contact_3_id)
        self.graph.draw(radius_2_id, "line", "O2", contact_2_id)
        center = line_intersection(
            line_through(o3, contact_3),
            line_through(o2, contact_2),
        )
        center_id = f"{class_id}_center_{sign}"
        self.graph.point(center_id, radius_3_id, radius_2_id)
        radius_squared = dot(
            subtract(center, contact_3),
            subtract(center, contact_3),
        )
        output_id = f"target_{sign}"
        draw_index = self.graph.draw(
            output_id,
            "circle",
            center_id,
            contact_3_id,
        )

        for input_center, input_radius, direction_sign in zip(
            (o1, o2, o3),
            (4, 2, 1),
            sign,
            strict=True,
        ):
            center_distance_squared = dot(
                subtract(center, input_center),
                subtract(center, input_center),
            )
            signed_term = (
                center_distance_squared
                - radius_squared
                - input_radius**2
            )
            residual = (
                signed_term * signed_term
                - 4 * radius_squared * input_radius**2
            )
            if residual != 0:
                raise AssertionError(f"{sign} 的切触多项式不为零")
            expected_positive = direction_sign == "+"
            if (signed_term.sign() > 0) != expected_positive:
                raise AssertionError(f"{sign} 的切向符号错误")
            if (
                not expected_positive
                and (radius_squared - input_radius**2).sign() <= 0
            ):
                raise AssertionError(f"{sign} 的内含型半径条件错误")

        self.targets[sign] = TargetResult(
            sign=sign,
            center=center,
            radius_squared=radius_squared,
            output_id=output_id,
            draw_index=draw_index,
        )

    def run(self) -> None:
        self.build_global_prefix()
        for key in ("alphaA", "aB", "a1A", "alpha1B"):
            self.draw_batch(key)
        self.build_pair(
            "P0",
            ("alphaA", "aB", "a1A", "alpha1B"),
            (("+++", 1), ("---", -1)),
        )

        for key in ("alphaB", "aA"):
            self.draw_batch(key)
        self.build_pair(
            "P2",
            ("alphaB", "aA", "alpha1B", "a1A"),
            (("-++", 1), ("+--", -1)),
        )

        for key in ("a1B", "alpha1A"):
            self.draw_batch(key)
        self.build_pair(
            "P1",
            ("alphaB", "aA", "a1B", "alpha1A"),
            (("--+", 1), ("++-", -1)),
        )
        self.build_pair(
            "P3",
            ("alphaA", "aB", "alpha1A", "a1B"),
            (("-+-", -1), ("+-+", 1)),
        )
        self.audit()

    def audit(self) -> None:
        expected_signs = {
            "---",
            "+++",
            "--+",
            "++-",
            "-++",
            "+--",
            "-+-",
            "+-+",
        }
        if set(self.targets) != expected_signs:
            raise AssertionError("八个物理切向符号没有全部出现")
        line_count = sum(
            kind == "line" for kind in self.graph.paid_kinds.values()
        )
        circle_count = sum(
            kind == "circle" for kind in self.graph.paid_kinds.values()
        )
        if (line_count, circle_count) != (55, 10):
            raise AssertionError("65 E 的直线/圆分解错误")
        if self.targets["+++"].draw_index != 18:
            raise AssertionError("三重外切圆没有在第 18 E 首次出现")
        if len(self.graph.paid_order) != 65:
            raise AssertionError("八解没有在第 65 E 完成")

        ancestor_sets = {
            sign: self.graph.paid_ancestors(result.output_id)
            for sign, result in self.targets.items()
        }
        if {len(items) for items in ancestor_sets.values()} != {18}:
            raise AssertionError("并非每个目标都有 18 个计费祖先")
        union = frozenset().union(*ancestor_sets.values())
        reuse = sum(len(items) for items in ancestor_sets.values()) - len(union)
        if len(union) != 65 or reuse != 79:
            raise AssertionError("八目标联合祖先或复用量错误")

        approximate_targets = {
            sign: (
                *approximate_point(result.center),
                sqrt(result.radius_squared.approximate()),
            )
            for sign, result in self.targets.items()
        }
        rounded_targets = {
            tuple(round(value, 8) for value in row)
            for row in approximate_targets.values()
        }
        if len(rounded_targets) != 8:
            raise AssertionError("八个目标圆不是两两不同")

        print(
            "score",
            {
                "lines": line_count,
                "circles": circle_count,
                "total": len(self.graph.paid_order),
                "first_ext": self.targets["+++"].draw_index,
                "all_targets": len(self.graph.paid_order),
            },
        )
        print(
            "dependencies",
            {
                "per_target": {
                    sign: len(ancestor_sets[sign])
                    for sign in sorted(ancestor_sets)
                },
                "union": len(union),
                "reuse": reuse,
            },
        )
        for sign in sorted(approximate_targets):
            center_x, center_y, radius = approximate_targets[sign]
            print(
                "target",
                sign,
                {
                    "center": (round(center_x, 12), round(center_y, 12)),
                    "radius": round(radius, 12),
                    "draw_index": self.targets[sign].draw_index,
                    "exact_tangency": True,
                },
            )


def main() -> None:
    replay = MannheimReplay()
    replay.run()


if __name__ == "__main__":
    main()
