"""完整重放 Mannheim 三简单合并夹具的有序依赖图。

固定严格 ``D8`` 夹具为

    Gamma1=((0, 0), 8)
    Gamma2=((253, 0), 5)
    Gamma3=((155, 120), 3)

``P0``、``P2``、``P3`` 各发生一次简单对向角色合并，``P1`` 正规。
脚本逐项建立合法 4 E 平行前缀、八条批量线、四个方向类及八个目标
后缀，并按几何对象本身去重。结果区分 66 E 的保守结构上界、57 E 的
去重历史轨迹和目标联合祖先裁剪后的 49 E 程序。
"""

from __future__ import annotations

from fractions import Fraction

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
    subtract,
)
from replay_mannheim_centered_parallel_repair import (
    build_roles,
    verify_target_pair,
)
from replay_mannheim_fixed import DependencyGraph, fraction_square_root
from scan_mannheim_degeneracies import (
    ROLE_KEYS,
    ROLE_NAMES,
    canonical_line,
    is_d8,
)


F = Fraction
Point = tuple[Fraction, Fraction]
Line = tuple[Fraction, Fraction, Fraction]
Circle = tuple[Point, Fraction]
Algebraic = Fraction | Quadratic
AlgebraicPoint = tuple[Algebraic, Algebraic]
AlgebraicLine = tuple[Algebraic, Algebraic, Algebraic]
AlgebraicCircle = tuple[AlgebraicPoint, Algebraic]

CENTERS: tuple[Point, Point, Point] = (
    (F(0), F(0)),
    (F(253), F(0)),
    (F(155), F(120)),
)
RADII = (F(8), F(5), F(3))

SIGMAS = {
    "P0": (1, 1, 1),
    "P1": (1, 1, -1),
    "P2": (1, -1, -1),
    "P3": (1, -1, 1),
}


def scalar_key(value) -> tuple:
    if isinstance(value, Quadratic):
        return ("quadratic", value.discriminant, value.a, value.b)
    return ("rational", F(value))


def point_key(point) -> tuple:
    return tuple(scalar_key(value) for value in point)


def collapse_scalar(value) -> Algebraic:
    if not isinstance(value, Quadratic):
        return F(value)
    try:
        square_root = fraction_square_root(value.discriminant)
    except ValueError:
        return value
    return value.a + value.b * square_root


def collapse_point(point) -> AlgebraicPoint:
    return tuple(collapse_scalar(value) for value in point)  # type: ignore[return-value]


def collapse_line(line) -> AlgebraicLine:
    return tuple(collapse_scalar(value) for value in line)  # type: ignore[return-value]


def collapse_circle(circle) -> AlgebraicCircle:
    center, radius_squared = circle
    return (collapse_point(center), collapse_scalar(radius_squared))


def algebraic_line_key(line) -> tuple:
    simplified = tuple(collapse_scalar(value) for value in line)
    discriminants = {
        value.discriminant
        for value in simplified
        if isinstance(value, Quadratic)
    }
    if not discriminants:
        return ("rational", canonical_line(simplified))
    if len(discriminants) != 1:
        raise ValueError("一条直线不能混合不同的二次扩域")
    discriminant = next(iter(discriminants))
    lifted = tuple(
        value
        if isinstance(value, Quadratic)
        else Quadratic(value, 0, discriminant)
        for value in simplified
    )
    pivot = next((value for value in lifted if value != 0), None)
    if pivot is None:
        raise ValueError("齐次直线系数不能全为零")
    normalized = tuple(value / pivot for value in lifted)
    if all(value.b == 0 for value in normalized):
        return (
            "rational",
            canonical_line(tuple(value.a for value in normalized)),
        )
    return (
        "quadratic",
        discriminant,
        tuple((value.a, value.b) for value in normalized),
    )


def circle(center: Point, through: Point) -> Circle:
    radius_squared = dot(subtract(center, through), subtract(center, through))
    if radius_squared == 0:
        raise ValueError("圆心与圆上一点重合")
    return (center, radius_squared)


def on_circle(point: Point, value: Circle) -> bool:
    center, radius_squared = value
    return dot(subtract(point, center), subtract(point, center)) == radius_squared


def role_batch_keys(profile: str) -> dict[str, str]:
    named = dict(zip(ROLE_NAMES, ROLE_KEYS[profile], strict=True))
    return {
        "x": named["a2"],
        "y": named["alpha2"],
        "z": named["a2_prime"],
        "w": named["alpha2_prime"],
    }


class ExactObjectGraph:
    """带几何对象去重和逻辑别名的依赖图。"""

    def __init__(
        self,
        centers: tuple[Point, Point, Point],
        radii: tuple[Fraction, Fraction, Fraction],
    ) -> None:
        self.graph = DependencyGraph.initial()
        self.aliases: dict[str, str] = {}
        self.paid_aliases: dict[str, str] = {}
        self.point_registry = {
            point_key(centers[index]): f"O{index + 1}" for index in range(3)
        }
        self.object_registry: dict[tuple, str] = {}
        for index, (center, radius) in enumerate(
            zip(centers, radii, strict=True),
            start=1,
        ):
            self.object_registry[
                ("circle", point_key(center), scalar_key(radius**2))
            ] = f"Gamma{index}"

    def resolve(self, node_id: str) -> str:
        while node_id in self.aliases:
            node_id = self.aliases[node_id]
        return node_id

    def point(self, node_id: str, value, *dependencies: str) -> str:
        key = point_key(value)
        if key in self.point_registry:
            existing = self.point_registry[key]
            self.aliases[node_id] = existing
            return existing
        resolved = tuple(self.resolve(item) for item in dependencies)
        self.graph.point(node_id, *resolved)
        self.point_registry[key] = node_id
        return node_id

    def _draw(
        self,
        node_id: str,
        kind: str,
        key: tuple,
        dependencies: tuple[str, ...],
    ) -> str:
        if key in self.object_registry:
            existing = self.object_registry[key]
            self.aliases[node_id] = existing
            self.paid_aliases[node_id] = existing
            return existing
        resolved = tuple(self.resolve(item) for item in dependencies)
        self.graph.draw(node_id, kind, *resolved)
        self.object_registry[key] = node_id
        return node_id

    def line(
        self,
        node_id: str,
        value: AlgebraicLine,
        *dependencies: str,
        special_key: str | None = None,
    ) -> str:
        key = (
            ("special-line", special_key)
            if special_key is not None
            else ("line", algebraic_line_key(value))
        )
        return self._draw(node_id, "line", key, tuple(dependencies))

    def circle(
        self,
        node_id: str,
        value: AlgebraicCircle | None,
        *dependencies: str,
        special_key: str | None = None,
    ) -> str:
        key = (
            ("special-circle", special_key)
            if special_key is not None
            else (
                "circle",
                point_key(value[0]),
                scalar_key(value[1]),
            )
        )
        return self._draw(node_id, "circle", key, tuple(dependencies))


class ThreeBlockReplay:
    def __init__(
        self,
        centers: tuple[Point, Point, Point] = CENTERS,
        radii: tuple[Fraction, Fraction, Fraction] = RADII,
    ) -> None:
        self.centers = centers
        self.radii = radii
        self.objects = ExactObjectGraph(centers, radii)
        self.o1, self.o2, self.o3 = centers
        self.r1, self.r2, self.r3 = radii
        self.gamma3 = (self.o3, self.r3**2)
        self.named_points: dict[str, Point] = {}
        self.named_point_ids: dict[str, str] = {}
        self.batch_lines: dict[str, Line] = {}
        self.batch_line_ids: dict[str, str] = {}
        self.batch_points: dict[str, Point] = {}
        self.batch_point_ids: dict[str, str] = {}
        self.similarity_centers: dict[str, Point] = {}
        self.similarity_center_ids: dict[str, str] = {}
        self.targets: dict[str, dict] = {}

    def build_prefix(self) -> None:
        ell = line_through(self.o1, self.o2)
        self.objects.line("ell", ell, "O1", "O2")

        direction = (F(1), F(0))
        endpoints = {
            "alpha": add(self.o1, multiply(-self.r1, direction)),
            "a": add(self.o1, multiply(self.r1, direction)),
            "a1": add(self.o2, multiply(-self.r2, direction)),
            "alpha1": add(self.o2, multiply(self.r2, direction)),
        }
        for name, point in endpoints.items():
            circle_id = "Gamma1" if name in {"alpha", "a"} else "Gamma2"
            point_id = self.objects.point(name, point, "ell", circle_id)
            self.named_points[name] = point
            self.named_point_ids[name] = point_id

        distance_squared = dot(
            subtract(self.o3, self.o1),
            subtract(self.o3, self.o1),
        )
        try:
            fraction_square_root(distance_squared)
        except ValueError:
            pass
        else:
            raise AssertionError("夹具的平行前缀应进入非有理二次域")

        root = Quadratic(0, 1, distance_squared)
        x = (root, Quadratic(0, 0, distance_squared))
        self.objects.circle(
            "parallel_c0",
            (self.o1, distance_squared),
            "O1",
            "O3",
        )
        x_id = self.objects.point("parallel_X", x, "parallel_c0", "ell")
        self.objects.circle(
            "parallel_cX",
            None,
            x_id,
            "O3",
            special_key="center-X-through-O3",
        )
        reflected_o3 = (self.o3[0], -self.o3[1])
        q_id = self.objects.point(
            "parallel_Q",
            reflected_o3,
            "parallel_c0",
            "parallel_cX",
        )
        self.objects.line(
            "parallel_diameter",
            (F(1), F(0), F(0)),
            x_id,
            q_id,
            special_key="irrational-XQ",
        )
        r = (
            2 * root - reflected_o3[0],
            Quadratic(self.o3[1], 0, distance_squared),
        )
        r_id = self.objects.point(
            "parallel_R",
            r,
            "parallel_diameter",
            "parallel_cX",
        )
        ell3 = (F(0), F(1), -self.o3[1])
        self.objects.line("ell3", ell3, "O3", r_id)

        for name, sign in (("A", -1), ("B", 1)):
            point = add(self.o3, multiply(sign * self.r3, direction))
            point_id = self.objects.point(name, point, "ell3", "Gamma3")
            self.named_points[name] = point
            self.named_point_ids[name] = point_id

    def draw_batch(self, key: str) -> None:
        if key in self.batch_lines:
            return
        if key.startswith("alpha1"):
            source_name = "alpha1"
        elif key.startswith("alpha"):
            source_name = "alpha"
        elif key.startswith("a1"):
            source_name = "a1"
        elif key.startswith("a"):
            source_name = "a"
        else:
            raise ValueError(f"未知批量线 {key}")
        endpoint_name = key[-1]
        source = self.named_points[source_name]
        endpoint = self.named_points[endpoint_name]
        value = line_through(source, endpoint)
        line_id = self.objects.line(
            f"batch_{key}",
            value,
            self.named_point_ids[source_name],
            self.named_point_ids[endpoint_name],
        )
        second = second_circle_point(
            source,
            endpoint,
            (self.o3, self.r3),
            allow_tangent_alias=True,
        )
        point_id = self.objects.point(
            f"second_{key}",
            second,
            line_id,
            "Gamma3",
        )
        self.batch_lines[key] = value
        self.batch_line_ids[key] = line_id
        self.batch_points[key] = second
        self.batch_point_ids[key] = point_id

    def bind_similarity_center(self, kind: str) -> tuple[Point, str]:
        if kind in self.similarity_centers:
            return (
                self.similarity_centers[kind],
                self.similarity_center_ids[kind],
            )
        parents = (
            ("a1A", "alpha1B") if kind == "ext" else ("a1B", "alpha1A")
        )
        value = line_intersection(
            self.batch_lines[parents[0]],
            self.batch_lines[parents[1]],
        )
        if kind == "ext":
            denominator = self.r2 - self.r3
            expected = (
                (
                    self.r2 * self.o3[0] - self.r3 * self.o2[0]
                )
                / denominator,
                (
                    self.r2 * self.o3[1] - self.r3 * self.o2[1]
                )
                / denominator,
            )
        else:
            denominator = self.r2 + self.r3
            expected = (
                (
                    self.r2 * self.o3[0] + self.r3 * self.o2[0]
                )
                / denominator,
                (
                    self.r2 * self.o3[1] + self.r3 * self.o2[1]
                )
                / denominator,
            )
        if value != expected:
            raise AssertionError(f"{kind} 相似中心错误")
        node_id = f"H23_{kind}"
        point_id = self.objects.point(
            node_id,
            value,
            self.batch_line_ids[parents[0]],
            self.batch_line_ids[parents[1]],
        )
        self.similarity_centers[kind] = value
        self.similarity_center_ids[kind] = point_id
        return value, point_id

    def build_target(
        self,
        profile: str,
        sign: str,
        target: dict,
        contact_3_id: str,
    ) -> None:
        paid_lines = tuple(collapse_line(value) for value in target["paid_lines"])
        if len(paid_lines) != 3:
            raise AssertionError("三块合并夹具的目标后缀应有三条候选线")
        contact_line, radius_2, radius_3 = paid_lines
        contact_3 = collapse_point(target["contact_3"])
        target_center = collapse_point(target["center"])

        physical_sign = target.get("physical_sign", sign)
        similarity_kind = (
            "ext" if physical_sign[1] == physical_sign[2] else "int"
        )
        _, similarity_id = self.bind_similarity_center(similarity_kind)
        contact_line_id = self.objects.line(
            f"{profile}_{sign}_contact_line",
            contact_line,
            similarity_id,
            contact_3_id,
        )
        contact_2 = line_intersection(contact_line, radius_2)
        if not on_circle(contact_2, (self.o2, self.r2**2)):
            raise AssertionError(f"{profile} {sign} 的第二圆接触点错误")
        contact_2_id = self.objects.point(
            f"{profile}_{sign}_M2",
            contact_2,
            contact_line_id,
            "Gamma2",
        )
        radius_3_id = self.objects.line(
            f"{profile}_{sign}_radius_3",
            radius_3,
            "O3",
            contact_3_id,
        )
        radius_2_id = self.objects.line(
            f"{profile}_{sign}_radius_2",
            radius_2,
            "O2",
            contact_2_id,
        )
        recovered_center = line_intersection(radius_3, radius_2)
        if recovered_center != target_center:
            raise AssertionError(f"{profile} {sign} 没有恢复目标圆心")
        center_id = self.objects.point(
            f"{profile}_{sign}_center",
            target_center,
            radius_3_id,
            radius_2_id,
        )
        output_circle = collapse_circle(target["output_circle"])
        output_id = self.objects.circle(
            f"target_{sign}",
            output_circle,
            center_id,
            contact_3_id,
        )
        self.targets[sign] = {
            "profile": profile,
            "output_id": output_id,
            "circle": output_circle,
            "draw_index": self.objects.graph.paid_order.index(output_id) + 1,
        }

    def verify_pair(
        self,
        profile: str,
        tau: Line,
        *,
        allow_repeated_physical_signs: bool = False,
    ) -> dict[str, dict]:
        targets, _ = verify_target_pair(
            self.centers,
            self.radii,
            tau,
            sigma=SIGMAS[profile],
            allow_repeated_physical_signs=allow_repeated_physical_signs,
        )
        return targets

    def build_regular_pair(
        self,
        profile: str,
        target_order: tuple[str, str] | None,
        *,
        allow_repeated_physical_signs: bool = False,
    ) -> None:
        keys = role_batch_keys(profile)
        roles = build_roles(self.centers, self.radii, profile)
        x, y, z, w = (roles[name] for name in "xyzw")
        role_ids = {
            name: self.batch_point_ids[key] for name, key in keys.items()
        }

        k_left = line_through(x, w)
        k_right = line_through(y, z)
        k_left_id = self.objects.line(
            f"{profile}_K_line_1",
            k_left,
            role_ids["x"],
            role_ids["w"],
        )
        k_right_id = self.objects.line(
            f"{profile}_K_line_2",
            k_right,
            role_ids["y"],
            role_ids["z"],
        )
        k = line_intersection(k_left, k_right)
        k_id = self.objects.point(f"{profile}_K", k, k_left_id, k_right_id)

        kp_left = line_through(x, y)
        kp_right = line_through(z, w)
        kp_left_id = self.objects.line(
            f"{profile}_Kp_line_1",
            kp_left,
            role_ids["x"],
            role_ids["y"],
        )
        kp_right_id = self.objects.line(
            f"{profile}_Kp_line_2",
            kp_right,
            role_ids["z"],
            role_ids["w"],
        )
        k_prime = line_intersection(kp_left, kp_right)
        kp_id = self.objects.point(
            f"{profile}_Kp",
            k_prime,
            kp_left_id,
            kp_right_id,
        )
        tau = line_through(k, k_prime)
        tau_id = self.objects.line(
            f"{profile}_tau",
            tau,
            k_id,
            kp_id,
        )

        targets = self.verify_pair(
            profile,
            tau,
            allow_repeated_physical_signs=allow_repeated_physical_signs,
        )
        contact_ids = {}
        for sign, target in targets.items():
            contact = collapse_point(target["contact_3"])
            contact_ids[sign] = self.objects.point(
                f"{profile}_{sign}_M3",
                contact,
                tau_id,
                "Gamma3",
            )
        ordered_targets = tuple(targets) if target_order is None else target_order
        for sign in ordered_targets:
            logical_sign = (
                f"{profile}:{sign}" if allow_repeated_physical_signs else sign
            )
            self.build_target(
                profile,
                logical_sign,
                targets[sign],
                contact_ids[sign],
            )

    def build_parallel_pair(
        self,
        profile: str,
        infinite_diagonal: str,
        target_order: tuple[str, str] | None,
        *,
        allow_repeated_physical_signs: bool = False,
    ) -> None:
        """重放第 8.8 节的一般两圆一线平行修复。"""

        if infinite_diagonal not in {"K", "Kp"}:
            raise ValueError("无穷远对角点只能是 K 或 Kp")
        keys = role_batch_keys(profile)
        roles = build_roles(self.centers, self.radii, profile)
        role_ids = {
            name: self.batch_point_ids[key] for name, key in keys.items()
        }
        x, y, z, w = (roles[name] for name in "xyzw")

        if infinite_diagonal == "Kp":
            finite_pairs = (("x", "w"), ("y", "z"))
            infinite_pairs = (("x", "y"), ("z", "w"))
            first_name, second_name = "x", "y"
        else:
            finite_pairs = (("x", "y"), ("z", "w"))
            infinite_pairs = (("x", "w"), ("y", "z"))
            first_name, second_name = "x", "w"

        finite_lines = tuple(
            line_through(*(roles[name] for name in pair))
            for pair in finite_pairs
        )
        finite_line_ids = tuple(
            self.objects.line(
                f"{profile}_{infinite_diagonal}_finite_line_{index}",
                value,
                *(role_ids[name] for name in pair),
            )
            for index, (value, pair) in enumerate(
                zip(finite_lines, finite_pairs, strict=True),
                start=1,
            )
        )
        finite_point = line_intersection(*finite_lines)
        finite_point_id = self.objects.point(
            f"{profile}_{infinite_diagonal}_finite_point",
            finite_point,
            *finite_line_ids,
        )

        infinite_lines = tuple(
            line_through(*(roles[name] for name in pair))
            for pair in infinite_pairs
        )
        if same_line(*infinite_lines):
            raise AssertionError("平行分支的两条定义弦意外重合")
        if determinant(infinite_lines[0][:2], infinite_lines[1][:2]) != 0:
            raise AssertionError("声明的对角点没有位于无穷远")

        first_point = roles[first_name]
        first_point_id = role_ids[first_name]
        second_point = roles[second_name]
        second_point_id = role_ids[second_name]
        first_circle = circle(first_point, finite_point)
        first_circle_id = self.objects.circle(
            f"{profile}_{infinite_diagonal}_circle_first",
            first_circle,
            first_point_id,
            finite_point_id,
        )
        reflected = add(first_point, subtract(first_point, finite_point))
        reflected_id = self.objects.point(
            f"{profile}_{infinite_diagonal}_reflected",
            reflected,
            first_circle_id,
            finite_line_ids[0],
        )
        second_circle = circle(second_point, reflected)
        second_circle_id = self.objects.circle(
            f"{profile}_{infinite_diagonal}_circle_second",
            second_circle,
            second_point_id,
            reflected_id,
        )
        q = add(finite_point, subtract(first_point, second_point))
        if not on_circle(q, first_circle) or not on_circle(q, second_circle):
            raise AssertionError("两圆一线修复的第二公共点错误")
        q_id = self.objects.point(
            f"{profile}_{infinite_diagonal}_Q",
            q,
            first_circle_id,
            second_circle_id,
        )
        tau = line_through(finite_point, q)
        if determinant(tau[:2], infinite_lines[0][:2]) != 0:
            raise AssertionError("修复所得 tau 没有通过无穷远对角点")
        tau_id = self.objects.line(
            f"{profile}_tau",
            tau,
            finite_point_id,
            q_id,
        )

        targets = self.verify_pair(
            profile,
            tau,
            allow_repeated_physical_signs=allow_repeated_physical_signs,
        )
        contact_ids = {}
        for sign, target in targets.items():
            contact = collapse_point(target["contact_3"])
            contact_ids[sign] = self.objects.point(
                f"{profile}_{sign}_M3",
                contact,
                tau_id,
                "Gamma3",
            )
        ordered_targets = tuple(targets) if target_order is None else target_order
        for sign in ordered_targets:
            logical_sign = (
                f"{profile}:{sign}" if allow_repeated_physical_signs else sign
            )
            self.build_target(
                profile,
                logical_sign,
                targets[sign],
                contact_ids[sign],
            )

    def build_merge_pair(
        self,
        profile: str,
        target_order: tuple[str, str],
        preferred_merged_key: str | None = None,
    ) -> None:
        keys = role_batch_keys(profile)
        roles = build_roles(self.centers, self.radii, profile)
        equal_pairs = [
            (first, second)
            for index, first in enumerate("xyzw")
            for second in "xyzw"[index + 1 :]
            if roles[first] == roles[second]
        ]
        if len(equal_pairs) != 1 or equal_pairs[0] not in {
            ("x", "z"),
            ("y", "w"),
        }:
            raise AssertionError(f"{profile} 不是简单对向合并")
        first, second = equal_pairs[0]
        merged = roles[first]
        merged_candidates = (keys[first], keys[second])
        if preferred_merged_key is None:
            preferred_merged_key = merged_candidates[0]
        if preferred_merged_key not in merged_candidates:
            raise AssertionError("指定的合并点来源不属于合并角色")
        merged_id = self.batch_point_ids[preferred_merged_key]

        radial = line_through(self.o3, merged)
        radial_id = self.objects.line(
            f"{profile}_merged_radius",
            radial,
            "O3",
            merged_id,
        )
        circle_p = circle(merged, self.o3)
        circle_p_id = self.objects.circle(
            f"{profile}_circle_P",
            circle_p,
            merged_id,
            "O3",
        )
        reflected = add(merged, subtract(merged, self.o3))
        reflected_id = self.objects.point(
            f"{profile}_reflected",
            reflected,
            circle_p_id,
            radial_id,
        )
        shared_circle = circle(self.o3, reflected)
        shared_id = self.objects.circle(
            f"{profile}_shared_2r3",
            shared_circle,
            "O3",
            reflected_id,
        )
        circle_e = circle(reflected, self.o3)
        circle_e_id = self.objects.circle(
            f"{profile}_circle_E",
            circle_e,
            reflected_id,
            "O3",
        )

        radial_vector = subtract(merged, self.o3)
        tangent_points = (
            (
                Quadratic(merged[0], -radial_vector[1], 3),
                Quadratic(merged[1], radial_vector[0], 3),
            ),
            (
                Quadratic(merged[0], radial_vector[1], 3),
                Quadratic(merged[1], -radial_vector[0], 3),
            ),
        )
        tangent_point_ids = tuple(
            self.objects.point(
                f"{profile}_tangent_point_{index}",
                point,
                shared_id,
                circle_e_id,
            )
            for index, point in enumerate(tangent_points, start=1)
        )
        tangent = (
            merged[0] - self.o3[0],
            merged[1] - self.o3[1],
            -(
                (merged[0] - self.o3[0]) * merged[0]
                + (merged[1] - self.o3[1]) * merged[1]
            ),
        )
        tangent_id = self.objects.line(
            f"{profile}_tangent",
            tangent,
            *tangent_point_ids,
        )

        if (first, second) == ("x", "z"):
            other_names = ("y", "w")
        else:
            other_names = ("x", "z")
        other = line_through(*(roles[name] for name in other_names))
        other_id = self.objects.line(
            f"{profile}_other_chord",
            other,
            *(self.batch_point_ids[keys[name]] for name in other_names),
        )
        tangent_center = line_intersection(tangent, other)
        tangent_center_id = self.objects.point(
            f"{profile}_tangent_center",
            tangent_center,
            tangent_id,
            other_id,
        )
        tangent_circle = circle(tangent_center, merged)
        tangent_circle_id = self.objects.circle(
            f"{profile}_tangent_circle",
            tangent_circle,
            tangent_center_id,
            merged_id,
        )

        radius_squared = tangent_circle[1]
        tau = (
            2 * (tangent_center[0] - self.o3[0]),
            2 * (tangent_center[1] - self.o3[1]),
            self.o3[0] ** 2
            + self.o3[1] ** 2
            - self.r3**2
            - tangent_center[0] ** 2
            - tangent_center[1] ** 2
            + radius_squared,
        )
        targets = self.verify_pair(profile, tau)
        contacts = {
            sign: collapse_point(target["contact_3"])
            for sign, target in targets.items()
        }
        if merged not in contacts.values():
            raise AssertionError(f"{profile} 的合并点不是目标接触点")
        other_contacts = set(contacts.values()) - {merged}
        if len(other_contacts) != 1:
            raise AssertionError(f"{profile} 没有唯一的第二接触点")
        second_contact = next(iter(other_contacts))
        if not on_circle(second_contact, self.gamma3):
            raise AssertionError("第二接触点不在第三圆上")
        if not on_circle(second_contact, tangent_circle):
            raise AssertionError("第二接触点不在修复圆上")
        if not same_line(line_through(merged, second_contact), tau):
            raise AssertionError("修复圆公共弦不是目标接触弦")
        second_contact_id = self.objects.point(
            f"{profile}_second_contact",
            second_contact,
            tangent_circle_id,
            "Gamma3",
        )

        contact_ids = {
            sign: merged_id if contact == merged else second_contact_id
            for sign, contact in contacts.items()
        }
        for sign in target_order:
            self.build_target(profile, sign, targets[sign], contact_ids[sign])

    def run(self) -> None:
        if not is_d8(self.centers, self.radii):
            raise AssertionError("三块合并夹具不属于 D8")
        self.build_prefix()

        # P0 的 aB 与 a1A 合并。先只画 a1A，使三重外切目标尽早完成；
        # aB 稍后仍要服务 P3，所以不会从八解联合程序中删除。
        for key in ("alphaA", "a1A", "alpha1B"):
            self.draw_batch(key)
        self.build_merge_pair(
            "P0",
            ("+++", "---"),
            preferred_merged_key="a1A",
        )

        for key in ("aB", "alphaB", "aA", "a1B", "alpha1A"):
            self.draw_batch(key)
        self.build_regular_pair("P1", ("++-", "--+"))
        self.build_merge_pair("P2", ("+--", "-++"))
        self.build_merge_pair("P3", ("+-+", "-+-"))
        self.audit()

    def audit(self) -> None:
        expected_signs = {
            "+++",
            "---",
            "++-",
            "--+",
            "+--",
            "-++",
            "+-+",
            "-+-",
        }
        if set(self.targets) != expected_signs:
            raise AssertionError("八个物理切向符号没有全部出现")
        target_circles = {target["circle"] for target in self.targets.values()}
        if len(target_circles) != 8:
            raise AssertionError("八个目标圆不是两两不同")

        graph = self.objects.graph
        line_count = sum(kind == "line" for kind in graph.paid_kinds.values())
        circle_count = sum(
            kind == "circle" for kind in graph.paid_kinds.values()
        )
        if (line_count, circle_count) != (37, 20):
            raise AssertionError("三块合并夹具的精确对象分解错误")
        if len(graph.paid_order) != 57:
            raise AssertionError("三块合并夹具没有在 57 E 完成")
        if self.targets["+++"]["draw_index"] != 16:
            raise AssertionError("三重外切圆没有在第 16 E 完成")

        ancestor_sets = {
            sign: graph.paid_ancestors(target["output_id"])
            for sign, target in self.targets.items()
        }
        union = frozenset().union(*ancestor_sets.values())
        missing = set(graph.paid_order) - set(union)
        expected_missing = {
            "P1_K_line_1",
            "P1_K_line_2",
            "P1_Kp_line_1",
            "P1_Kp_line_2",
            "P1_tau",
            "P3_+-+_contact_line",
            "batch_aA",
            "batch_aB",
        }
        if len(union) != 49 or missing != expected_missing:
            raise AssertionError(
                f"八目标联合祖先为 {len(union)}，遗漏 {sorted(missing)}"
            )
        reuse = sum(len(items) for items in ancestor_sets.values()) - len(union)
        pruned_line_count = sum(
            graph.paid_kinds[node_id] == "line" for node_id in union
        )
        pruned_circle_count = sum(
            graph.paid_kinds[node_id] == "circle" for node_id in union
        )
        if (pruned_line_count, pruned_circle_count) != (29, 20):
            raise AssertionError("裁剪程序的直线/圆分解错误")

        structural_upper = 5 + 8 + 13 + 3 * 14 - 2
        if structural_upper != 66:
            raise AssertionError("保守结构上界计算错误")
        additional_fixture_reuse = structural_upper - len(graph.paid_order)
        if additional_fixture_reuse != 9:
            raise AssertionError("夹具的额外对象复用量错误")
        ancestor_pruning = len(graph.paid_order) - len(union)
        if ancestor_pruning != 8:
            raise AssertionError("目标联合祖先裁剪量错误")
        if len(self.objects.paid_aliases) != 14:
            raise AssertionError("计费对象别名数量错误")

        print(
            "score",
            {
                "lines": line_count,
                "circles": circle_count,
                "deduplicated_trace": len(graph.paid_order),
                "first_ext": self.targets["+++"]["draw_index"],
                "pruned_ext": len(ancestor_sets["+++"]),
                "pruned_all_targets": len(union),
                "pruned_lines": pruned_line_count,
                "pruned_circles": pruned_circle_count,
                "structural_upper": structural_upper,
                "additional_fixture_reuse": additional_fixture_reuse,
                "ancestor_pruning": ancestor_pruning,
            },
        )
        print(
            "dependencies",
            {
                "per_target": {
                    sign: len(ancestor_sets[sign]) for sign in sorted(ancestor_sets)
                },
                "union": len(union),
                "reuse": reuse,
                "non_ancestors": sorted(missing),
            },
        )
        print("paid_aliases", dict(sorted(self.objects.paid_aliases.items())))


def main() -> None:
    replay = ThreeBlockReplay()
    replay.run()


if __name__ == "__main__":
    main()
