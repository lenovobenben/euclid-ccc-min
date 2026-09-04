"""精确重放 Mannheim 双对角色合并时的 1 E 修复。

固定夹具为

    Gamma1=((0, 0), 5)
    Gamma2=((16, 0), 3)
    Gamma3=((12, 6), 1)

其 ``P2`` 块满足 ``x = z`` 且 ``y = w``，四条对角定义弦全部重合。
退化极线就是连接两个不同角色点的弦；脚本在实二次域中验证两个目标圆。
一般联合模块至多 9 E，而这个固定夹具还复用一条定心半径，实际为 8 E。
"""

from __future__ import annotations

from fractions import Fraction

from check_mannheim_degenerate_fixture import same_line
from replay_mannheim_centered_parallel_repair import (
    assert_distinct_lines,
    build_roles,
    lift_line,
    verify_target_pair,
)
from scan_mannheim_degeneracies import canonical_line, is_d8


F = Fraction


def main() -> None:
    centers = ((F(0), F(0)), (F(16), F(0)), (F(12), F(6)))
    radii = (F(5), F(3), F(1))
    if not is_d8(centers, radii):
        raise AssertionError("双对合并夹具不属于 D8")

    roles = build_roles(centers, radii, profile="P2")
    expected_roles = {
        "x": (F(12), F(7)),
        "y": (F(56, 5), F(27, 5)),
        "z": (F(12), F(7)),
        "w": (F(56, 5), F(27, 5)),
    }
    if roles != expected_roles:
        raise AssertionError("P2 双对角色点与预期值不符")
    if roles["x"] != roles["z"] or roles["y"] != roles["w"]:
        raise AssertionError("夹具没有发生双对角色合并")
    if roles["x"] == roles["y"]:
        raise AssertionError("双对合并只剩一个角色点")

    tau = (
        roles["x"][1] - roles["y"][1],
        roles["y"][0] - roles["x"][0],
        roles["x"][0] * roles["y"][1]
        - roles["y"][0] * roles["x"][1],
    )
    if canonical_line(tau) != (2, -1, -17):
        raise AssertionError("双对合并弦与预期 tau 不符")

    targets, discriminant = verify_target_pair(
        centers,
        radii,
        tau,
        sigma=(1, -1, -1),
        expected_discriminant=F(144),
    )
    lifted_tau = lift_line(tau, discriminant)
    candidate_lines = (lifted_tau,) + tuple(
        line
        for target in targets.values()
        for line in target["paid_lines"]
    )
    paid_lines = []
    for line in candidate_lines:
        if any(same_line(line, existing) for existing in paid_lines):
            continue
        paid_lines.append(line)
    paid_circles = tuple(
        target["output_circle"] for target in targets.values()
    )
    assert_distinct_lines(paid_lines)
    if len(set(paid_circles)) != len(paid_circles):
        raise AssertionError("两个输出圆意外重合")
    repaired_pair_cost = len(paid_lines) + len(paid_circles)
    if repaired_pair_cost != 8:
        raise AssertionError("固定双对合并夹具不是 8 E")
    if repaired_pair_cost > 1 + 2 * 4:
        raise AssertionError("双对合并分支超过 9 E 上界")

    print(
        "double_merge",
        {
            "x=z": roles["x"],
            "y=w": roles["y"],
            "tau": canonical_line(tau),
            "pair_cost": repaired_pair_cost,
            "delta_vs_regular_pair": -5,
            "generic_pair_upper_bound": 9,
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
