"""在整数盒中逐夹具执行 Mannheim 61 E 有序分支证书。

与只统计局部类别成本的 ``scan_mannheim_branch_costs.py`` 不同，本脚本
对每个 ``D8`` 输入真正建立精确对象依赖图，并构造八个目标圆。扫描仍是
有限回归，不承担连续参数覆盖证明。
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict

from replay_mannheim_ordered_branches import OrderedBranchReplay
from scan_mannheim_degeneracies import is_d8, make_fixture


def scan(max_radius: int, max_coordinate: int) -> None:
    tested = 0
    in_domain = 0
    trace_distribution: Counter[int] = Counter()
    union_distribution: Counter[int] = Counter()
    branch_distribution: Counter[tuple[str, str, str, str]] = Counter()
    maximum_trace = -1
    maximum_union = -1
    rows_by_cost = defaultdict(list)

    for r1 in range(3, max_radius + 1):
        for r2 in range(2, r1):
            for r3 in range(1, r2):
                for distance in range(r1 + r2 + 1, max_coordinate + 1):
                    for u in range(-3, distance + 4):
                        for v in range(1, max_coordinate + 1):
                            tested += 1
                            centers, radii = make_fixture(
                                r1, r2, r3, distance, u, v
                            )
                            if not is_d8(centers, radii):
                                continue
                            in_domain += 1
                            report = OrderedBranchReplay(
                                f"scan_{in_domain}",
                                centers,
                                radii,
                                emit=False,
                            ).run()
                            trace = report["trace"]
                            union = report["all_targets"]
                            branch_key = tuple(
                                report["branches"][profile]
                                for profile in ("P0", "P2", "P1", "P3")
                            )
                            trace_distribution[trace] += 1
                            union_distribution[union] += 1
                            branch_distribution[branch_key] += 1
                            maximum_trace = max(maximum_trace, trace)
                            maximum_union = max(maximum_union, union)
                            if len(rows_by_cost[(trace, union)]) < 10:
                                rows_by_cost[(trace, union)].append(
                                    {
                                        "centers": centers,
                                        "radii": radii,
                                        "branches": branch_key,
                                        "trace": trace,
                                        "all_targets": union,
                                    }
                                )

    print(
        "ordered_branch_scan",
        {
            "tested": tested,
            "in_D8": in_domain,
            "maximum_trace": maximum_trace,
            "maximum_all_targets": maximum_union,
        },
    )
    print("trace_distribution", dict(sorted(trace_distribution.items())))
    print("union_distribution", dict(sorted(union_distribution.items())))
    print("branch_distribution")
    for branches, count in sorted(branch_distribution.items()):
        print({"branches": branches, "count": count})
    print("first_joint_maximum_fixtures")
    for row in rows_by_cost[(maximum_trace, maximum_union)]:
        print(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-radius", type=int, default=6)
    parser.add_argument("--max-coordinate", type=int, default=18)
    args = parser.parse_args()
    scan(args.max_radius, args.max_coordinate)


if __name__ == "__main__":
    main()
