"""扫描 Mannheim 圆心线程序同时退化时的局部成本分布。

脚本复用 ``scan_mannheim_degeneracies`` 的严格 ``D8`` 判定和 Fraction
中间对象。无有限简单合并时，P0、P2 各用至多 5 E 构造 ``tau`` 核心，
P1、P3 各用至多 3 E，再让四个方向类各用 7 E 圆心线后缀。有限简单
合并先行块保守记完整 13 E，但同时扣除一条强制冗余批量线；若两块都
简单合并，再扣除它们共享的 ``2r3`` 圆。尾部合并块按 1 E ``tau`` 加
7 E 后缀计。居中分支和其它偶然对象重合的额外节省均忽略。

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
            return 8
        if simple_merge:
            return 13
        return 12
    if double_merge or simple_merge:
        return 8
    return 10


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
                            seed_simple_count = sum(
                                (
                                    f"{profile}:merge:a2=a2_prime" in events
                                )
                                != (
                                    f"{profile}:merge:alpha2=alpha2_prime"
                                    in events
                                )
                                for profile in ("P0", "P2")
                            )
                            branch_cost = 13 + sum(
                                profile_cost(profile, events)
                                for profile in PROFILES
                            )
                            branch_cost -= seed_simple_count
                            branch_cost -= max(0, seed_simple_count - 1)
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
