"""验证同一 Mannheim 方向中三个简单合并块可以同时出现。

固定严格 ``D8`` 夹具为

    Gamma1=((0, 0), 8)
    Gamma2=((253, 0), 5)
    Gamma3=((155, 120), 3)

其 ``P0``、``P2``、``P3`` 分别发生一个简单对向角色合并，``P1``
保持四点不同。三个修复共享同一个以 ``O3`` 为圆心、半径为 ``2*r3``
的圆，所以联合保守成本是 66 E。全部判定均使用 Fraction 和精确
``D8`` 根符号检查。
"""

from __future__ import annotations

from fractions import Fraction

from replay_mannheim_centered_parallel_repair import build_roles
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
PROFILES = ("P0", "P1", "P2", "P3")


def equal_pairs(roles) -> set[tuple[str, str]]:
    names = ("x", "y", "z", "w")
    return {
        (first, second)
        for first_index, first in enumerate(names)
        for second in names[first_index + 1 :]
        if roles[first] == roles[second]
    }


def main() -> None:
    centers = (
        (F(0), F(0)),
        (F(253), F(0)),
        (F(155), F(120)),
    )
    radii = (F(8), F(5), F(3))
    if not is_d8(centers, radii):
        raise AssertionError("三块合并夹具不属于 D8")

    expected = {
        "P0": {("x", "z")},
        "P1": set(),
        "P2": {("x", "z")},
        "P3": {("y", "w")},
    }
    actual = {}
    merged_points = {}
    for profile in PROFILES:
        roles = build_roles(centers, radii, profile=profile)
        actual[profile] = equal_pairs(roles)
        if actual[profile] != expected[profile]:
            raise AssertionError(f"{profile} 的角色合并型错误")
        if actual[profile]:
            first_name, _ = next(iter(actual[profile]))
            merged_points[profile] = roles[first_name]

    events = analyze_fixture(centers, radii)
    expected_events = {
        "P0:merge:a2=a2_prime",
        "P0:collapse:K=Kp",
        "P2:merge:a2=a2_prime",
        "P2:collapse:K=Kp",
        "P3:merge:alpha2=alpha2_prime",
        "P3:collapse:K=Kp",
    }
    structural_events = {
        event
        for event in events
        if event.startswith(PROFILES)
        and any(
            marker in event
            for marker in (
                ":merge:",
                ":parallel:",
                ":coincident:",
                ":collapse:",
                ":undefined:",
            )
        )
    }
    if structural_events != expected_events:
        raise AssertionError("三块合并夹具含有未声明的结构事件")

    normalized = {
        "p": radii[0] / radii[2],
        "q": radii[1] / radii[2],
        "c": -centers[2][0] / radii[2],
        "d": centers[1][0] / radii[2],
        "h": -centers[2][1] / radii[2],
    }
    expected_normalized = {
        "p": F(8, 3),
        "q": F(5, 3),
        "c": F(-155, 3),
        "d": F(253, 3),
        "h": F(-40),
    }
    if normalized != expected_normalized:
        raise AssertionError("归一化参数错误")

    simple_blocks = sum(bool(pairs) for pairs in actual.values())
    shared_circles = set()
    for merged in merged_points.values():
        reflected = (
            2 * merged[0] - centers[2][0],
            2 * merged[1] - centers[2][1],
        )
        radius_squared = (
            (reflected[0] - centers[2][0]) ** 2
            + (reflected[1] - centers[2][1]) ** 2
        )
        shared_circles.add((centers[2], radius_squared))
    expected_shared_circle = {(centers[2], F(36))}
    if shared_circles != expected_shared_circle:
        raise AssertionError("三个简单合并修复没有复用半径 2*r3 的圆")

    per_block_upper = 65 + simple_blocks
    shared_branch_upper = 65 + int(bool(simple_blocks))
    if (
        simple_blocks != 3
        or per_block_upper != 68
        or shared_branch_upper != 66
    ):
        raise AssertionError("三块合并的保守模块成本错误")

    print(
        "three_block_merge",
        {
            "centers": centers,
            "radii": radii,
            "normalized": normalized,
            "equal_pairs": actual,
            "simple_blocks": simple_blocks,
            "per_block_upper": per_block_upper,
            "shared_circle": next(iter(shared_circles)),
            "shared_branch_upper": shared_branch_upper,
        },
    )


if __name__ == "__main__":
    main()
