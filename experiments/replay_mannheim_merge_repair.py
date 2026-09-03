"""精确重放 Mannheim 对向角色合并时的 1 E 增量修复。

固定夹具为

    Gamma1=((0, 0), 3)
    Gamma2=((7, 0), 2)
    Gamma3=((3, 5), 1)

``O1O2`` 方向的 ``P0`` 块满足 ``alpha2 = alpha2_prime``。脚本从该
合并点构造第三圆切线，再用一个以切线交点为圆心的圆同时恢复这一对解
在第三圆上的两个接触点。所有正确性断言均为 Fraction 或实二次域中的
精确等式。
"""

from __future__ import annotations

from fractions import Fraction

from check_mannheim_degenerate_fixture import (
    add,
    determinant,
    dot,
    line_intersection,
    line_through,
    multiply,
    second_circle_point,
    subtract,
)
from replay_mannheim_fixed import Quadratic
from scan_mannheim_degeneracies import canonical_line, is_d8


F = Fraction


def lift_point(point, discriminant):
    return (
        Quadratic(point[0], 0, discriminant),
        Quadratic(point[1], 0, discriminant),
    )


def lift_line(line, discriminant):
    return tuple(Quadratic(value, 0, discriminant) for value in line)


def on_line(point, line) -> bool:
    return line[0] * point[0] + line[1] * point[1] + line[2] == 0


def same_line(first, second) -> bool:
    return (
        first[0] * second[1] == first[1] * second[0]
        and first[0] * second[2] == first[2] * second[0]
        and first[1] * second[2] == first[2] * second[1]
    )


def contact_point(input_center, input_radius, target_center, target_radius, sign):
    delta = subtract(target_center, input_center)
    center_distance = target_radius + sign * input_radius
    return add(
        input_center,
        multiply(sign * input_radius / center_distance, delta),
    )


def verify_target(
    sign_name,
    target_center,
    target_radius,
    expected_third_contact,
    similarity_center,
    centers,
    radii,
):
    physical_sign = 1 if sign_name == "+++" else -1
    contacts = tuple(
        contact_point(
            input_center,
            input_radius,
            target_center,
            target_radius,
            physical_sign,
        )
        for input_center, input_radius in zip(centers, radii, strict=True)
    )
    if contacts[2] != expected_third_contact:
        raise AssertionError(f"{sign_name} 的第三圆接触点错误")
    if determinant(
        subtract(contacts[1], similarity_center),
        subtract(contacts[2], similarity_center),
    ) != 0:
        raise AssertionError(f"{sign_name} 的反对应接触点不经过相似中心")

    recovery_line = line_through(similarity_center, contacts[2])
    radial_2 = line_through(centers[1], contacts[1])
    radial_3 = line_through(centers[2], contacts[2])
    recovered_center = line_intersection(radial_2, radial_3)
    if recovered_center != target_center:
        raise AssertionError(f"{sign_name} 没有恢复目标圆心")
    if not on_line(contacts[1], recovery_line):
        raise AssertionError(f"{sign_name} 没有恢复第二圆接触点")
    radius_squared = dot(
        subtract(contacts[2], recovered_center),
        subtract(contacts[2], recovered_center),
    )
    if radius_squared != target_radius**2:
        raise AssertionError(f"{sign_name} 的最终圆半径错误")
    for input_center, input_radius in zip(centers, radii, strict=True):
        delta = subtract(target_center, input_center)
        expected_distance = target_radius + physical_sign * input_radius
        if dot(delta, delta) != expected_distance**2:
            raise AssertionError(f"{sign_name} 的切触等式错误")
    return {
        "center": target_center,
        "radius": target_radius,
        "contact_2": contacts[1],
        "contact_3": contacts[2],
        "paid_lines": (recovery_line, radial_2, radial_3),
        "output_circle": (target_center, radius_squared),
    }


def main() -> None:
    centers = ((F(0), F(0)), (F(7), F(0)), (F(3), F(5)))
    radii = (F(3), F(2), F(1))
    o1, o2, o3 = centers
    r1, r2, r3 = radii
    if not is_d8(centers, radii):
        raise AssertionError("合并修复夹具不属于 D8")
    direction = (F(1), F(0))
    named_points = {
        "alpha": add(o1, multiply(-r1, direction)),
        "a": add(o1, multiply(r1, direction)),
        "a1": add(o2, multiply(-r2, direction)),
        "alpha1": add(o2, multiply(r2, direction)),
        "A": add(o3, multiply(-r3, direction)),
        "B": add(o3, multiply(r3, direction)),
    }
    batch_points = {
        source_name + endpoint_name: second_circle_point(
            named_points[source_name],
            named_points[endpoint_name],
            (o3, r3),
            allow_tangent_alias=True,
        )
        for source_name in ("alpha", "a", "a1", "alpha1")
        for endpoint_name in ("A", "B")
    }
    roles = {
        "alpha2": batch_points["alphaA"],
        "a2": batch_points["aB"],
        "a2_prime": batch_points["a1A"],
        "alpha2_prime": batch_points["alpha1B"],
    }
    merged = (F(3), F(6))
    if roles["alpha2"] != merged or roles["alpha2_prime"] != merged:
        raise AssertionError("P0 没有按夹具声明发生对向角色合并")
    if len(set(roles.values())) != 3:
        raise AssertionError("P0 除指定合并外还有其它角色合并")

    radial_merged = line_through(o3, merged)
    if not same_line(radial_merged, (F(1), F(0), F(-3))):
        raise AssertionError("合并点半径线错误")
    reflected_center = add(merged, subtract(merged, o3))
    if reflected_center != (F(3), F(7)):
        raise AssertionError("Circle(P,O3) 没有给出 O3 关于 P 的对称点")

    discriminant = F(3)
    upper = (Quadratic(3, 1, discriminant), Quadratic(6, 0, discriminant))
    lower = (Quadratic(3, -1, discriminant), Quadratic(6, 0, discriminant))
    lifted_o3 = lift_point(o3, discriminant)
    lifted_reflected = lift_point(reflected_center, discriminant)
    for intersection in (upper, lower):
        if (
            dot(
                subtract(intersection, lifted_o3),
                subtract(intersection, lifted_o3),
            )
            != 4
        ):
            raise AssertionError("切线宏交点不在 Circle(O3,E) 上")
        if (
            dot(
                subtract(intersection, lifted_reflected),
                subtract(intersection, lifted_reflected),
            )
            != 4
        ):
            raise AssertionError("切线宏交点不在 Circle(E,O3) 上")
    tangent = line_through(upper, lower)
    if not same_line(tangent, lift_line((F(0), F(1), F(-6)), discriminant)):
        raise AssertionError("五步切线宏没有得到 Gamma3 在 P 处的切线")

    other_chord = line_through(roles["a2"], roles["a2_prime"])
    tangent_center = line_intersection(lift_line(other_chord, discriminant), tangent)
    expected_tangent_center = lift_point((F(39, 5), F(6)), discriminant)
    if tangent_center != expected_tangent_center:
        raise AssertionError("切线与另一角色弦的交点错误")

    second_contact = (F(2043, 601), F(2454, 601))
    tangent_center_rational = (tangent_center[0].a, tangent_center[1].a)
    tangent_radius_squared = dot(
        subtract(tangent_center_rational, merged),
        subtract(tangent_center_rational, merged),
    )
    for contact in (merged, second_contact):
        if dot(subtract(contact, o3), subtract(contact, o3)) != r3**2:
            raise AssertionError("接触点不在 Gamma3 上")
        if (
            dot(
                subtract(contact, tangent_center_rational),
                subtract(contact, tangent_center_rational),
            )
            != tangent_radius_squared
        ):
            raise AssertionError("接触点不在 Circle(J,P) 上")
    if merged == second_contact:
        raise AssertionError("Circle(J,P) 没有给出不同的第二接触点")

    similarity_center = (F(-1), F(10))
    targets = {
        "---": verify_target(
            "---",
            (F(3), F(0)),
            F(6),
            merged,
            similarity_center,
            centers,
            radii,
        ),
        "+++": verify_target(
            "+++",
            (F(4437, 1079), F(2640, 1079)),
            F(1926, 1079),
            second_contact,
            similarity_center,
            centers,
            radii,
        ),
    }

    repair_lines = (
        (F(0), F(1), F(-6)),
        other_chord,
    )
    repair_circles = (
        (merged, dot(subtract(merged, o3), subtract(merged, o3))),
        (o3, dot(subtract(o3, reflected_center), subtract(o3, reflected_center))),
        (
            reflected_center,
            dot(subtract(reflected_center, o3), subtract(reflected_center, o3)),
        ),
        (tangent_center_rational, tangent_radius_squared),
    )
    suffix_lines = tuple(
        line
        for target in targets.values()
        for line in target["paid_lines"]
    )
    suffix_circles = tuple(
        target["output_circle"] for target in targets.values()
    )
    paid_lines = repair_lines + suffix_lines
    paid_circles = repair_circles + suffix_circles
    if len({canonical_line(line) for line in paid_lines}) != len(paid_lines):
        raise AssertionError("合并修复模块含有重复直线")
    if len(set(paid_circles)) != len(paid_circles):
        raise AssertionError("合并修复模块含有重复圆")
    repaired_pair_cost = len(paid_lines) + len(paid_circles)
    if repaired_pair_cost != 14:
        raise AssertionError("合并修复模块不是 14 个不同计费对象")
    regular_pair_cost = 5 + 2 * 4
    if repaired_pair_cost != regular_pair_cost + 1:
        raise AssertionError("合并修复没有保持 1 E 增量")

    print(
        "merge",
        {
            "alpha2": roles["alpha2"],
            "alpha2_prime": roles["alpha2_prime"],
            "other_roles": (roles["a2"], roles["a2_prime"]),
        },
    )
    print(
        "repair",
        {
            "tangent_center": tangent_center_rational,
            "contacts": (merged, second_contact),
            "pair_cost": repaired_pair_cost,
            "delta_vs_regular": 1,
        },
    )
    print(
        "targets",
        {
            sign: {
                key: value
                for key, value in target.items()
                if key not in {"paid_lines", "output_circle"}
            }
            for sign, target in targets.items()
        },
    )


if __name__ == "__main__":
    main()
