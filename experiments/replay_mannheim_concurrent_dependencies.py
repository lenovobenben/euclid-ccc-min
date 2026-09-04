"""精确重放 Mannheim 共点优化后的正规与双 ``K'`` 八解程序。

四个方向类的接触弦 ``tau`` 恒过同一个 Mannheim 共点 ``S``。正规夹具
先构造 P0、P2，再让 P1、P3 复用公共 ``K'`` 弦；双 ``K'`` 平行夹具
先构造 P0、P1，再令 P2、P3 各画两条有限 K 弦并连接 S 与 K。两种程序
都把前两条 ``tau`` 的交点作为免费点 S。
"""

from __future__ import annotations

from fractions import Fraction

from check_mannheim_degenerate_fixture import line_intersection, line_through
from replay_mannheim_centered_parallel_repair import build_roles
from replay_mannheim_three_block_dependencies import (
    ThreeBlockReplay,
    collapse_point,
    role_batch_keys,
)
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
RADII = (F(4), F(2), F(1))
FIXTURES = {
    "regular": (
        ((F(0), F(0)), (F(13), F(0)), (F(4), F(15))),
        set(),
        ("P0", "P2"),
        ("P1", "P3"),
        "Kp",
        False,
        (47, 10, 57),
    ),
    "double_Kp": (
        ((F(0), F(0)), (F(45), F(0)), (F(30), F(21))),
        {"P2:parallel:Kp", "P3:parallel:Kp"},
        ("P0", "P1"),
        ("P2", "P3"),
        "K",
        True,
        (51, 10, 61),
    ),
}

TARGET_ORDERS = {
    "P0": ("+++", "---"),
    "P1": ("++-", "--+"),
    "P2": ("+--", "-++"),
    "P3": ("+-+", "-+-"),
}


def regular_tau(replay: ThreeBlockReplay, profile: str):
    roles = build_roles(replay.centers, replay.radii, profile)
    x, y, z, w = (roles[name] for name in "xyzw")
    k = line_intersection(line_through(x, w), line_through(y, z))
    k_prime = line_intersection(line_through(x, y), line_through(z, w))
    return line_through(k, k_prime)


def build_pair_from_mannheim_point(
    replay: ThreeBlockReplay,
    profile: str,
    mannheim_point,
    mannheim_point_id: str,
    finite_diagonal: str,
    *,
    allow_repeated_physical_signs: bool,
) -> None:
    roles = build_roles(replay.centers, replay.radii, profile)
    role_ids = {
        name: replay.batch_point_ids[key]
        for name, key in role_batch_keys(profile).items()
    }
    if finite_diagonal == "Kp":
        pairs = (("x", "y"), ("z", "w"))
    elif finite_diagonal == "K":
        pairs = (("x", "w"), ("y", "z"))
    else:
        raise ValueError("finite_diagonal 必须是 K 或 Kp")

    diagonal_lines = tuple(
        line_through(*(roles[name] for name in pair)) for pair in pairs
    )
    diagonal_line_ids = tuple(
        replay.objects.line(
            f"{profile}_{finite_diagonal}_S_line_{index}",
            value,
            *(role_ids[name] for name in pair),
        )
        for index, (value, pair) in enumerate(
            zip(diagonal_lines, pairs, strict=True),
            start=1,
        )
    )
    diagonal_point = line_intersection(*diagonal_lines)
    diagonal_point_id = replay.objects.point(
        f"{profile}_{finite_diagonal}_S_point",
        diagonal_point,
        *diagonal_line_ids,
    )
    tau = line_through(mannheim_point, diagonal_point)
    tau_id = replay.objects.line(
        f"{profile}_tau_from_S",
        tau,
        mannheim_point_id,
        diagonal_point_id,
    )

    targets = replay.verify_pair(
        profile,
        tau,
        allow_repeated_physical_signs=allow_repeated_physical_signs,
    )
    contact_ids = {
        sign: replay.objects.point(
            f"{profile}_{sign}_M3_from_B",
            collapse_point(target["contact_3"]),
            tau_id,
            "Gamma3",
        )
        for sign, target in targets.items()
    }
    target_order = (
        tuple(targets)
        if allow_repeated_physical_signs
        else TARGET_ORDERS[profile]
    )
    for sign in target_order:
        logical_sign = (
            f"{profile}:{sign}" if allow_repeated_physical_signs else sign
        )
        replay.build_target(
            profile,
            logical_sign,
            targets[sign],
            contact_ids[sign],
        )


def run_fixture(fixture_name: str) -> dict:
    (
        centers,
        expected_events,
        seed_profiles,
        remaining_profiles,
        finite_diagonal,
        repeated_signs,
        expected_count,
    ) = FIXTURES[fixture_name]
    if not is_d8(centers, RADII):
        raise AssertionError(f"{fixture_name} 不属于严格 D8")
    if analyze_fixture(centers, RADII) != expected_events:
        raise AssertionError(f"{fixture_name} 的退化型错误")

    replay = ThreeBlockReplay(centers, RADII)
    replay.build_prefix()
    for key in (
        "alphaA",
        "aB",
        "a1A",
        "alpha1B",
    ):
        replay.draw_batch(key)
    first_seed, second_seed = seed_profiles
    if first_seed != "P0":
        raise AssertionError("首个方向类必须保留 P0 以达到 18 E 外切解")
    replay.build_regular_pair(
        first_seed,
        None if repeated_signs else TARGET_ORDERS[first_seed],
        allow_repeated_physical_signs=repeated_signs,
    )
    for key in (
        "alphaB",
        "aA",
        "a1B",
        "alpha1A",
    ):
        replay.draw_batch(key)
    replay.build_regular_pair(
        second_seed,
        None if repeated_signs else TARGET_ORDERS[second_seed],
        allow_repeated_physical_signs=repeated_signs,
    )

    first_tau = regular_tau(replay, first_seed)
    second_tau = regular_tau(replay, second_seed)
    mannheim_point = line_intersection(first_tau, second_tau)
    mannheim_point_id = replay.objects.point(
        "Mannheim_S",
        mannheim_point,
        replay.objects.resolve(f"{first_seed}_tau"),
        replay.objects.resolve(f"{second_seed}_tau"),
    )
    for profile in remaining_profiles:
        build_pair_from_mannheim_point(
            replay,
            profile,
            mannheim_point,
            mannheim_point_id,
            finite_diagonal,
            allow_repeated_physical_signs=repeated_signs,
        )

    graph = replay.objects.graph
    if len(replay.targets) != 8:
        raise AssertionError("没有恢复八个目标圆")
    if len({target["circle"] for target in replay.targets.values()}) != 8:
        raise AssertionError("八个目标圆不是两两不同")
    ancestors = {
        sign: graph.paid_ancestors(target["output_id"])
        for sign, target in replay.targets.items()
    }
    union = frozenset().union(*ancestors.values())
    if set(union) != set(graph.paid_order):
        raise AssertionError("存在不属于八目标联合祖先的计费对象")
    line_count = sum(graph.paid_kinds[node] == "line" for node in union)
    circle_count = sum(graph.paid_kinds[node] == "circle" for node in union)
    if (line_count, circle_count, len(union)) != expected_count:
        raise AssertionError(f"{fixture_name} 的共点程序成本错误")

    reuse = sum(len(items) for items in ancestors.values()) - len(union)
    report = {
        "fixture": fixture_name,
        "lines": line_count,
        "circles": circle_count,
        "all_targets": len(union),
        "first_ext": replay.targets[
            "+++" if not repeated_signs else "P0:+++@-1"
        ]["draw_index"],
        "per_target": {
            sign: len(items) for sign, items in sorted(ancestors.items())
        },
        "reuse": reuse,
        "paid_aliases": dict(sorted(replay.objects.paid_aliases.items())),
    }
    print("concurrent_replay", report)
    return report


def main() -> None:
    reports = {
        fixture_name: run_fixture(fixture_name) for fixture_name in FIXTURES
    }
    if reports["regular"]["all_targets"] != 57:
        raise AssertionError("正规共点程序没有达到 57 E")
    if reports["double_Kp"]["all_targets"] != 61:
        raise AssertionError("双 K' 共点程序没有达到 61 E")


if __name__ == "__main__":
    main()
