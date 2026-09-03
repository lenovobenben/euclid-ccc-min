"""精确重放 Mannheim 的一般 ``K'`` 平行 5 E 修复。

固定夹具为

    Gamma1=((0, 0), 5)
    Gamma2=((28, 0), 4)
    Gamma3=((13, 27/2), 1)

其 ``P0`` 块的 ``K'`` 位于无穷远点，而有限对角点 ``K != O3``。
脚本先用两条弦取得 ``K``，再以两圆一线作出经过 ``K`` 且平行于
``xy`` 的目标弦；最后在实二次域中验证该弦对应的两个目标圆和 13 E
联合模块计数。
"""

from __future__ import annotations

from fractions import Fraction

from check_mannheim_degenerate_fixture import (
    add,
    determinant,
    dot,
    line_through,
    multiply,
    same_line,
    subtract,
)
from replay_mannheim_centered_parallel_repair import (
    assert_distinct_lines,
    build_roles,
    diagonal_data,
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
    centers = (
        (F(0), F(0)),
        (F(28), F(0)),
        (F(13), F(27, 2)),
    )
    radii = (F(5), F(4), F(1))
    if not is_d8(centers, radii):
        raise AssertionError("K' 平行夹具不属于 D8")

    roles = build_roles(centers, radii)
    expected_roles = {
        "x": (F(174, 13), F(327, 26)),
        "y": (F(24932, 1885), F(54567, 3770)),
        "z": (F(1868, 145), F(3627, 290)),
        "w": (F(318, 25), F(723, 50)),
    }
    if roles != expected_roles:
        raise AssertionError("P0 角色点与预期值不符")
    k, k_left, k_right, kp_left, kp_right = diagonal_data(roles)
    expected_k = (F(457, 35), F(1418, 105))
    if k != expected_k or k == centers[2]:
        raise AssertionError("有限对角点 K 与预期值不符")
    if canonical_line(kp_left) != (312, 26, -4503):
        raise AssertionError("第一条 K' 定义弦错误")
    if canonical_line(kp_right) != (120, 10, -1671):
        raise AssertionError("第二条 K' 定义弦错误")

    x, y = roles["x"], roles["y"]
    kx_squared = dot(subtract(x, k), subtract(x, k))
    ky_squared = dot(subtract(y, k), subtract(y, k))
    if kx_squared != ky_squared:
        raise AssertionError("圆内相交弦没有给出 Kx=Ky")

    circle_x = (x, kx_squared)
    reflected = add(add(x, x), multiply(-1, k))
    if reflected != (F(6239, 455), F(15901, 1365)):
        raise AssertionError("Circle(x,K) 的对径交点错误")
    if not on_circle(k, circle_x) or not on_circle(reflected, circle_x):
        raise AssertionError("第一个辅助圆的交点错误")
    if determinant(
        subtract(reflected, k),
        subtract(x, k),
    ) != 0:
        raise AssertionError("反射点不在已有弦 Kx 上")

    circle_y = (
        y,
        dot(subtract(reflected, y), subtract(reflected, y)),
    )
    q = add(add(k, x), multiply(-1, y))
    if q != (F(34875, 2639), F(91898, 7917)):
        raise AssertionError("两个辅助圆的另一交点错误")
    for point in (reflected, q):
        if not on_circle(point, circle_x) or not on_circle(point, circle_y):
            raise AssertionError("辅助圆没有通过两个声明交点")
    if reflected == q:
        raise AssertionError("两个辅助圆没有不同的第二交点")

    tau = line_through(k, q)
    if canonical_line(tau) != (252, 21, -3574):
        raise AssertionError("两圆一线宏没有得到预期 tau")
    if determinant(tau[:2], kp_left[:2]) != 0:
        raise AssertionError("所得 tau 不平行于 K' 的定义弦")

    targets, discriminant = verify_target_pair(centers, radii, tau)
    paid_lines = tuple(
        lift_line(line, discriminant)
        for line in (k_left, k_right, tau)
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
            lift_point(circle_y[0], discriminant),
            Quadratic(circle_y[1], 0, discriminant),
        ),
    ) + tuple(target["output_circle"] for target in targets.values())
    assert_distinct_lines(paid_lines)
    if len(set(paid_circles)) != len(paid_circles):
        raise AssertionError("模块含有重复计费圆")
    repaired_pair_cost = len(paid_lines) + len(paid_circles)
    if repaired_pair_cost != 13:
        raise AssertionError("K' 平行修复与双后缀不是 13 E")
    if repaired_pair_cost != 5 + 2 * 4:
        raise AssertionError("K' 平行分支没有保持正规模块成本")

    print(
        "parallel_kp",
        {
            "K": k,
            "Kp_lines": (
                canonical_line(kp_left),
                canonical_line(kp_right),
            ),
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
