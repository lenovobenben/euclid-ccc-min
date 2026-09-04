"""精确重放 Mannheim 的一般 ``K`` 平行 5 E 修复。

固定夹具为

    Gamma1=((0, 0), 3)
    Gamma2=((7, 0), 2)
    Gamma3=((8, 4), 1)

其 ``P2`` 块的 ``K`` 位于无穷远点，而 ``K'`` 有限。脚本先用两条弦
取得 ``K'``，再以两圆一线作出经过 ``K'`` 且平行于 ``xw`` 的目标弦；
最后在实二次域中验证该弦对应的两个目标圆和 13 E 联合模块计数。
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
    same_line,
    subtract,
)
from replay_mannheim_centered_parallel_repair import (
    assert_distinct_lines,
    build_roles,
    lift_line,
    lift_point,
    verify_target_pair,
)
from replay_mannheim_fixed import Quadratic
from scan_mannheim_degeneracies import canonical_line, is_d8


F = Fraction


def on_circle(point, circle) -> bool:
    center, radius_squared = circle
    delta = subtract(point, center)
    return dot(delta, delta) == radius_squared


def main() -> None:
    centers = ((F(0), F(0)), (F(7), F(0)), (F(8), F(4)))
    radii = (F(3), F(2), F(1))
    if not is_d8(centers, radii):
        raise AssertionError("K 平行夹具不属于 D8")

    roles = build_roles(centers, radii, profile="P2")
    expected_roles = {
        "x": (F(8), F(5)),
        "y": (F(36, 5), F(17, 5)),
        "z": (F(9), F(4)),
        "w": (F(37, 5), F(24, 5)),
    }
    if roles != expected_roles:
        raise AssertionError("P2 角色点与预期值不符")
    x, y, z, w = (roles[name] for name in ("x", "y", "z", "w"))
    k_left = line_through(x, w)
    k_right = line_through(y, z)
    if determinant(k_left[:2], k_right[:2]) != 0:
        raise AssertionError("K 的两条定义弦不平行")
    if same_line(k_left, k_right):
        raise AssertionError("K 的两条定义弦意外重合")

    kp_left = line_through(x, y)
    kp_right = line_through(z, w)
    kp = line_intersection(kp_left, kp_right)
    expected_kp = (F(39, 5), F(23, 5))
    if kp != expected_kp:
        raise AssertionError("有限对角点 K' 与预期值不符")

    kpx_squared = dot(subtract(x, kp), subtract(x, kp))
    kpw_squared = dot(subtract(w, kp), subtract(w, kp))
    if kpx_squared != kpw_squared:
        raise AssertionError("圆内相交弦没有给出 K'x=K'w")

    circle_x = (x, kpx_squared)
    reflected = add(add(x, x), multiply(-1, kp))
    if reflected != (F(41, 5), F(27, 5)):
        raise AssertionError("Circle(x,K') 的对径交点错误")
    if not on_circle(kp, circle_x) or not on_circle(reflected, circle_x):
        raise AssertionError("第一个辅助圆的交点错误")
    if determinant(
        subtract(reflected, kp),
        subtract(x, kp),
    ) != 0:
        raise AssertionError("反射点不在已有弦 K'x 上")

    circle_w = (
        w,
        dot(subtract(reflected, w), subtract(reflected, w)),
    )
    q = add(add(kp, x), multiply(-1, w))
    if q != (F(42, 5), F(24, 5)):
        raise AssertionError("两个辅助圆的另一交点错误")
    for point in (reflected, q):
        if not on_circle(point, circle_x) or not on_circle(point, circle_w):
            raise AssertionError("辅助圆没有通过两个声明交点")
    if reflected == q:
        raise AssertionError("两个辅助圆没有不同的第二交点")

    tau = line_through(kp, q)
    if canonical_line(tau) != (1, -3, 6):
        raise AssertionError("两圆一线宏没有得到预期 tau")
    if determinant(tau[:2], k_left[:2]) != 0:
        raise AssertionError("所得 tau 不平行于 K 的定义弦")

    targets, discriminant = verify_target_pair(
        centers,
        radii,
        tau,
        sigma=(1, -1, -1),
    )
    paid_lines = tuple(
        lift_line(line, discriminant)
        for line in (kp_left, kp_right, tau)
    ) + tuple(
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
            lift_point(circle_w[0], discriminant),
            Quadratic(circle_w[1], 0, discriminant),
        ),
    ) + tuple(target["output_circle"] for target in targets.values())
    assert_distinct_lines(paid_lines)
    if len(set(paid_circles)) != len(paid_circles):
        raise AssertionError("模块含有重复计费圆")
    repaired_pair_cost = len(paid_lines) + len(paid_circles)
    if repaired_pair_cost != 13:
        raise AssertionError("K 平行修复与双后缀不是 13 E")
    if repaired_pair_cost != 5 + 2 * 4:
        raise AssertionError("K 平行分支没有保持正规模块成本")

    print(
        "parallel_k",
        {
            "K_lines": (
                canonical_line(k_left),
                canonical_line(k_right),
            ),
            "Kp": kp,
            "reflected": reflected,
            "Q": q,
            "tau": canonical_line(tau),
            "pair_cost": repaired_pair_cost,
            "delta_vs_regular_pair": 0,
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


if __name__ == "__main__":
    main()
