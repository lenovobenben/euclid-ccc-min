"""精确重放 Mannheim 居中平行退化块的 3 E 修复。

主夹具为

    Gamma1=((0, 0), 3)
    Gamma2=((9, 0), 2)
    Gamma3=((6, 4), 1)

其 ``P0`` 块满足 ``K'`` 位于无穷远点且有限对角点 ``K = O3``。
脚本用两个以角色点为圆心、经过 ``O3`` 的圆及其公共弦恢复目标直径，
再在实二次域中验证两个目标圆和 9 E 的联合模块计数。

脚本还给出一个严格 ``D8`` 的有理反例，证明 ``K'`` 位于无穷远点本身
不推出 ``K = O3``，因而这个便宜程序只能用于居中子分支。
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

from check_mannheim_degenerate_fixture import (
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
from replay_mannheim_fixed import Quadratic
from scan_mannheim_degeneracies import ROLE_KEYS, ROLE_NAMES, is_d8


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


def parallel_through(point, reference):
    return (
        reference[0],
        reference[1],
        -reference[0] * point[0] - reference[1] * point[1],
    )


def build_roles(centers, radii, profile="P0"):
    o1, o2, o3 = centers
    r1, r2, r3 = radii
    direction = (F(1), F(0))
    named = {
        "alpha": add(o1, multiply(-r1, direction)),
        "a": add(o1, multiply(r1, direction)),
        "a1": add(o2, multiply(-r2, direction)),
        "alpha1": add(o2, multiply(r2, direction)),
        "A": add(o3, multiply(-r3, direction)),
        "B": add(o3, multiply(r3, direction)),
    }
    batch_points = {
        source_name + endpoint_name: second_circle_point(
            named[source_name],
            named[endpoint_name],
            (o3, r3),
            allow_tangent_alias=True,
        )
        for source_name in ("alpha", "a", "a1", "alpha1")
        for endpoint_name in ("A", "B")
    }
    named_roles = dict(
        zip(
            ROLE_NAMES,
            (batch_points[key] for key in ROLE_KEYS[profile]),
            strict=True,
        )
    )
    return {
        "x": named_roles["a2"],
        "y": named_roles["alpha2"],
        "z": named_roles["a2_prime"],
        "w": named_roles["alpha2_prime"],
    }


def diagonal_data(roles):
    x, y, z, w = (roles[name] for name in ("x", "y", "z", "w"))
    k_left = line_through(x, w)
    k_right = line_through(y, z)
    kp_left = line_through(x, y)
    kp_right = line_through(z, w)
    if determinant(kp_left[:2], kp_right[:2]) != 0:
        raise AssertionError("K' 的两条定义弦不平行")
    if same_line(kp_left, kp_right):
        raise AssertionError("K' 的两条定义弦意外重合")
    return (
        line_intersection(k_left, k_right),
        k_left,
        k_right,
        kp_left,
        kp_right,
    )


def verify_target_pair(
    centers,
    radii,
    tau,
    sigma=(1, 1, 1),
    expected_discriminant=None,
    allow_repeated_physical_signs=False,
):
    e2 = subtract(centers[1], centers[0])
    e3 = subtract(centers[2], centers[0])
    b_vector = solve_dot_system(
        e2,
        e3,
        (dot(e2, e2) - radii[1] ** 2 + radii[0] ** 2) / 2,
        (dot(e3, e3) - radii[2] ** 2 + radii[0] ** 2) / 2,
    )
    u_vector = solve_dot_system(
        e2,
        e3,
        radii[0] - sigma[1] * radii[1],
        radii[0] - sigma[2] * radii[2],
    )
    quadratic = dot(u_vector, u_vector) - 1
    linear = 2 * dot(b_vector, u_vector) - 2 * radii[0]
    constant = dot(b_vector, b_vector) - radii[0] ** 2
    discriminant = linear**2 - 4 * quadratic * constant
    if (
        expected_discriminant is not None
        and discriminant != expected_discriminant
    ):
        raise AssertionError("P0 目标二次式的判别式错误")

    lifted_centers = tuple(
        lift_point(center, discriminant) for center in centers
    )
    lifted_b = lift_point(b_vector, discriminant)
    lifted_u = lift_point(u_vector, discriminant)
    lifted_tau = lift_line(tau, discriminant)
    if sigma[1] == sigma[2]:
        similarity_denominator = radii[1] - radii[2]
        similarity_center = (
            (
                radii[1] * centers[2][0]
                - radii[2] * centers[1][0]
            )
            / similarity_denominator,
            (
                radii[1] * centers[2][1]
                - radii[2] * centers[1][1]
            )
            / similarity_denominator,
        )
    else:
        similarity_denominator = radii[1] + radii[2]
        similarity_center = (
            (
                radii[1] * centers[2][0]
                + radii[2] * centers[1][0]
            )
            / similarity_denominator,
            (
                radii[1] * centers[2][1]
                + radii[2] * centers[1][1]
            )
            / similarity_denominator,
        )
    lifted_similarity = lift_point(similarity_center, discriminant)

    targets = {}
    for root_sign in (-1, 1):
        signed_radius = Quadratic(
            -linear / (2 * quadratic),
            F(root_sign, 2) / quadratic,
            discriminant,
        )
        signed_radius_sign = signed_radius.sign()
        radius = (
            signed_radius
            if signed_radius_sign > 0
            else -signed_radius
        )
        physical = tuple(signed_radius_sign * item for item in sigma)
        sign_name = "".join("+" if item > 0 else "-" for item in physical)
        center = add(lifted_b, multiply(signed_radius, lifted_u))
        contacts = []
        for input_center, input_radius, physical_sign in zip(
            lifted_centers,
            radii,
            physical,
            strict=True,
        ):
            delta = subtract(center, input_center)
            center_distance = radius + physical_sign * input_radius
            contact = add(
                input_center,
                multiply(
                    physical_sign * input_radius / center_distance,
                    delta,
                ),
            )
            if dot(delta, delta) != center_distance * center_distance:
                raise AssertionError(f"{sign_name} 的切触等式错误")
            if dot(
                subtract(contact, input_center),
                subtract(contact, input_center),
            ) != input_radius**2:
                raise AssertionError(f"{sign_name} 的接触点不在输入圆上")
            contacts.append(contact)

        contact_2, contact_3 = contacts[1], contacts[2]
        if not on_line(contact_3, lifted_tau):
            raise AssertionError(f"{sign_name} 的第三圆接触点不在 tau 上")
        if determinant(
            subtract(contact_2, lifted_similarity),
            subtract(contact_3, lifted_similarity),
        ) != 0:
            raise AssertionError(f"{sign_name} 的反对应接触点错误")

        recovery_line = line_through(lifted_similarity, contact_3)
        radial_2 = line_through(lifted_centers[1], contact_2)
        radial_3 = line_through(lifted_centers[2], contact_3)
        recovered_center = line_intersection(radial_2, radial_3)
        if recovered_center != center:
            raise AssertionError(f"{sign_name} 没有恢复目标圆心")
        radius_squared = dot(
            subtract(center, contact_3),
            subtract(center, contact_3),
        )
        if radius_squared != radius * radius:
            raise AssertionError(f"{sign_name} 的输出圆半径错误")
        paid_lines = (recovery_line, radial_2)
        if not same_line(radial_3, lifted_tau):
            paid_lines += (radial_3,)
        target_key = (
            f"{sign_name}@{root_sign:+d}"
            if allow_repeated_physical_signs
            else sign_name
        )
        targets[target_key] = {
            "physical_sign": sign_name,
            "center": center,
            "radius": radius,
            "contact_3": contact_3,
            "paid_lines": paid_lines,
            "output_circle": (center, radius_squared),
        }

    if allow_repeated_physical_signs:
        if len(targets) != 2:
            raise AssertionError("有向二次式没有恢复两个目标根")
    else:
        expected_signs = {
            "".join("+" if item > 0 else "-" for item in sigma),
            "".join("+" if item < 0 else "-" for item in sigma),
        }
        if set(targets) != expected_signs:
            raise AssertionError("没有恢复预期的一对物理切向符号")
    tau_through_center = on_line(lifted_centers[2], lifted_tau)
    contacts_are_antipodal = add(
        *(target["contact_3"] for target in targets.values()),
    ) == multiply(2, lifted_centers[2])
    if tau_through_center and not contacts_are_antipodal:
        raise AssertionError("两个第三圆接触点不是 tau 的对径点")
    return targets, discriminant


def assert_distinct_lines(lines) -> None:
    for first, second in combinations(lines, 2):
        if same_line(first, second):
            raise AssertionError("模块含有重复计费直线")


def main() -> None:
    centers = ((F(0), F(0)), (F(9), F(0)), (F(6), F(4)))
    radii = (F(3), F(2), F(1))
    if not is_d8(centers, radii):
        raise AssertionError("居中平行夹具不属于 D8")
    roles = build_roles(centers, radii)
    expected_roles = {
        "x": (F(6), F(3)),
        "y": (F(33, 5), F(24, 5)),
        "z": (F(27, 5), F(16, 5)),
        "w": (F(6), F(5)),
    }
    if roles != expected_roles:
        raise AssertionError("P0 角色点与预期值不符")
    k, _, _, kp_left, _ = diagonal_data(roles)
    o3 = centers[2]
    if k != o3:
        raise AssertionError("主夹具的有限对角点 K 不等于 O3")
    if add(roles["x"], roles["w"]) != multiply(2, o3):
        raise AssertionError("x,w 不是第三圆对径点")
    if add(roles["y"], roles["z"]) != multiply(2, o3):
        raise AssertionError("y,z 不是第三圆对径点")

    circle_x = (roles["x"], radii[2] ** 2)
    circle_z = (roles["z"], radii[2] ** 2)
    q = add(add(roles["x"], roles["z"]), multiply(-1, o3))
    for circle in (circle_x, circle_z):
        center, radius_squared = circle
        for point in (o3, q):
            if dot(subtract(point, center), subtract(point, center)) != radius_squared:
                raise AssertionError("辅助圆没有通过声明的公共交点")
    if q == o3:
        raise AssertionError("两个辅助圆没有不同于 O3 的第二交点")
    tau = line_through(o3, q)
    if not same_line(tau, (F(3), F(-1), F(-14))):
        raise AssertionError("三步修复没有得到预期直径")
    if determinant(tau[:2], kp_left[:2]) != 0:
        raise AssertionError("三步修复所得直径不平行于 K' 的定义弦")

    targets, discriminant = verify_target_pair(
        centers,
        radii,
        tau,
        expected_discriminant=F(640, 9),
    )
    lifted_tau = lift_line(tau, discriminant)
    paid_lines = (lifted_tau,) + tuple(
        line
        for target in targets.values()
        for line in target["paid_lines"]
    )
    paid_circles = (
        (
            lift_point(circle_x[0], discriminant),
            Quadratic(circle_x[1], 0, discriminant),
        ),
        (
            lift_point(circle_z[0], discriminant),
            Quadratic(circle_z[1], 0, discriminant),
        ),
    ) + tuple(target["output_circle"] for target in targets.values())
    assert_distinct_lines(paid_lines)
    if len(set(paid_circles)) != len(paid_circles):
        raise AssertionError("模块含有重复计费圆")
    repaired_pair_cost = len(paid_lines) + len(paid_circles)
    if repaired_pair_cost != 9:
        raise AssertionError("居中平行修复与双后缀不是 9 E")
    if repaired_pair_cost != 5 + 2 * 4 - 4:
        raise AssertionError("居中平行分支没有相对正规模块节省 4 E")

    counter_centers = (
        (F(0), F(0)),
        (F(28), F(0)),
        (F(13), F(27, 2)),
    )
    counter_radii = (F(5), F(4), F(1))
    if not is_d8(counter_centers, counter_radii):
        raise AssertionError("非居中平行反例不属于 D8")
    counter_roles = build_roles(counter_centers, counter_radii)
    counter_k, _, _, counter_kp_left, _ = diagonal_data(counter_roles)
    if counter_k == counter_centers[2]:
        raise AssertionError("非居中平行反例意外满足 K=O3")
    counter_q = add(
        add(counter_roles["x"], counter_roles["z"]),
        multiply(-1, counter_centers[2]),
    )
    cheap_line = line_through(counter_centers[2], counter_q)
    true_tau = parallel_through(counter_k, counter_kp_left)
    if same_line(cheap_line, true_tau):
        raise AssertionError("三步程序意外覆盖非居中平行反例")

    print(
        "centered_parallel",
        {
            "roles": roles,
            "K": k,
            "Q": q,
            "tau": tau,
            "pair_cost": repaired_pair_cost,
            "delta_vs_regular_pair": -4,
        },
    )
    print(
        "targets",
        {
            sign: {
                "center": target["center"],
                "radius": target["radius"],
                "contact_3": target["contact_3"],
            }
            for sign, target in targets.items()
        },
    )
    print(
        "scope_counterexample",
        {
            "centers": counter_centers,
            "radii": counter_radii,
            "K": counter_k,
            "O3": counter_centers[2],
        },
    )


if __name__ == "__main__":
    main()
