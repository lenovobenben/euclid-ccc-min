"""扫描 Mannheim 四个方向类同时退化时的局部成本分布。

脚本复用 ``scan_mannheim_degeneracies`` 的严格 ``D8`` 判定和 Fraction
中间对象。共点程序先完整构造 P0、P2，并保证两条 ``tau`` 都已画出；
每块保守记 13 E，双对合并仍只需 9 E。二者免费给出根心 S。后行的
P1、P3 若四点不同，只需两条有限对角弦和 ``Line(S,D)``，连两个目标
后缀至多 11 E；简单或双对合并直接用合并点和 S 画 ``tau``，至多
9 E。偶然对象重合、公共弦和居中平行的额外节省均忽略。

这是有界扫描，不是连续参数域上“至多两个简单合并”的证明。
"""

from __future__ import annotations

import argparse
from collections import Counter

from scan_mannheim_degeneracies import analyze_fixture, is_d8, make_fixture


PROFILES = ("P0", "P1", "P2", "P3")


def classify(events: set[str]) -> tuple[tuple[str, ...], ...]:
    simple = []
    double = []
    parallel = []
    for profile in PROFILES:
        merge_xz = f"{profile}:merge:a2=a2_prime" in events
        merge_yw = f"{profile}:merge:alpha2=alpha2_prime" in events
        if merge_xz and merge_yw:
            double.append(profile)
        elif merge_xz or merge_yw:
            simple.append(profile)
        if (
            f"{profile}:parallel:K" in events
            or f"{profile}:parallel:Kp" in events
        ):
            parallel.append(profile)
    return tuple(simple), tuple(double), tuple(parallel)


def profile_cost(profile: str, events: set[str]) -> int:
    merge_xz = f"{profile}:merge:a2=a2_prime" in events
    merge_yw = f"{profile}:merge:alpha2=alpha2_prime" in events
    double_merge = merge_xz and merge_yw
    simple_merge = merge_xz != merge_yw
    if profile in {"P0", "P2"}:
        if double_merge:
            return 9
        return 13
    if double_merge or simple_merge:
        return 9
    return 11


def scan(max_radius: int, max_coordinate: int) -> None:
    tested = 0
    in_domain = 0
    distribution: Counter[tuple[int, int, int, int]] = Counter()
    maximum_cost = -1
    maximum_rows = []
    for r1 in range(3, max_radius + 1):
        for r2 in range(2, r1):
            for r3 in range(1, r2):
                for distance in range(
                    r1 + r2 + 1,
                    max_coordinate + 1,
                ):
                    for u in range(-3, distance + 4):
                        for v in range(1, max_coordinate + 1):
                            tested += 1
                            fixture = make_fixture(
                                r1,
                                r2,
                                r3,
                                distance,
                                u,
                                v,
                            )
                            if not is_d8(*fixture):
                                continue
                            in_domain += 1
                            events = analyze_fixture(*fixture)
                            simple, double, parallel = classify(events)
                            branch_cost = 13 + sum(
                                profile_cost(profile, events)
                                for profile in PROFILES
                            )
                            distribution[
                                (
                                    len(simple),
                                    len(double),
                                    len(parallel),
                                    branch_cost,
                                )
                            ] += 1
                            row = {
                                "centers": fixture[0],
                                "radii": fixture[1],
                                "simple": simple,
                                "double": double,
                                "parallel": parallel,
                                "branch_cost_upper": branch_cost,
                            }
                            if branch_cost > maximum_cost:
                                maximum_cost = branch_cost
                                maximum_rows = [row]
                            elif (
                                branch_cost == maximum_cost
                                and len(maximum_rows) < 10
                            ):
                                maximum_rows.append(row)

    print(
        "scan",
        {
            "tested": tested,
            "in_D8": in_domain,
            "maximum_branch_cost_upper": maximum_cost,
        },
    )
    print("distribution")
    for key, count in sorted(distribution.items()):
        simple, double, parallel, branch_cost = key
        print(
            {
                "simple": simple,
                "double": double,
                "parallel": parallel,
                "branch_cost_upper": branch_cost,
                "count": count,
            }
        )
    print("first_maximum_fixtures")
    for row in maximum_rows:
        print(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-radius", type=int, default=6)
    parser.add_argument("--max-coordinate", type=int, default=18)
    args = parser.parse_args()
    scan(args.max_radius, args.max_coordinate)


if __name__ == "__main__":
    main()
