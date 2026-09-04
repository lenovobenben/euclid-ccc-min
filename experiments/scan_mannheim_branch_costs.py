"""扫描 Mannheim 有限 ``K'`` 圆心线程序的分支成本分布。

脚本复用 ``scan_mannheim_degeneracies`` 的严格 ``D8`` 判定和 Fraction
中间对象。有限 ``K'`` 且根心不在 ``O1O2`` 上时，每个双目标后缀由
7 E 降为 6 E。若 P0 或 P2 有有限简单合并，则先独立完成该解对，再由
两个目标圆心的连线与接触弦取得根心；其余方向类全部按后行块计数。
若根心位于 ``O1O2``，四个非退化 ``K'`` 会合到根心，改用一个 5 E
接触弦核心和三个至多 3 E 后行核心。其它偶然对象重合均忽略。

这是有界扫描，不代替连续参数域上的分支证明。
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


def merge_kind(profile: str, events: set[str]) -> str:
    merge_xz = f"{profile}:merge:a2=a2_prime" in events
    merge_yw = f"{profile}:merge:alpha2=alpha2_prime" in events
    if merge_xz and merge_yw:
        return "double"
    if merge_xz or merge_yw:
        return "simple"
    return "none"


def root_center_on_ell(fixture) -> bool:
    centers, radii = fixture
    o1, o2, o3 = centers
    r1, r2, r3 = radii
    distance = o2[0] - o1[0]
    root_x = (
        o2[0] ** 2
        + o2[1] ** 2
        - r2**2
        - o1[0] ** 2
        - o1[1] ** 2
        + r1**2
    ) / (2 * distance)
    root_y = (
        o3[0] ** 2
        + o3[1] ** 2
        - r3**2
        - o1[0] ** 2
        - o1[1] ** 2
        + r1**2
        - 2 * (o3[0] - o1[0]) * root_x
    ) / (2 * (o3[1] - o1[1]))
    return root_y == o1[1]


def optimized_cost(fixture, events: set[str]) -> int:
    seed_simple = tuple(
        profile
        for profile in ("P0", "P2")
        if merge_kind(profile, events) == "simple"
    )
    if seed_simple:
        selected = seed_simple[0]
        # 13 E 公共前缀，13 E 完整合并块，1 E 圆心线；合并又强制
        # 删除一条批量线。其余方向类作为后行块。
        cost = 13 + 13 + 1 - 1
        for profile in PROFILES:
            if profile == selected:
                continue
            if merge_kind(profile, events) != "none":
                cost += 8
            elif f"{profile}:parallel:Kp" in events:
                cost += 10
            else:
                cost += 9
        return cost

    base = 13 + sum(profile_cost(profile, events) for profile in PROFILES)
    if root_center_on_ell(fixture):
        return min(base, 55)
    finite_kp_savings = sum(
        merge_kind(profile, events) == "none"
        and f"{profile}:parallel:Kp" not in events
        for profile in PROFILES
    )
    return base - finite_kp_savings


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
                            branch_cost = optimized_cost(fixture, events)
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
