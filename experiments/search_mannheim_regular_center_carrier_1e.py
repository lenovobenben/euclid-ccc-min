"""搜索一步可画且同时经过一对目标圆心的载体对象。

正规 6 E 双目标后缀含两个输出圆和四条辅助线。若在接触弦核心完成后，
可以一步画出一条直线或一个圆并使其同时经过该方向类的两个目标圆心，
再画两条已知第三圆接触半径即可分别取得圆心，最后画两个输出圆：

    1 个双圆心载体 + 2 条接触半径 + 2 个输出圆 = 5 E。

脚本在三个严格正规 ``D8`` 夹具上同步枚举由两个已有确定点定义的全部
直线和圆，不允许任意点。命中必须另作精确重放和符号证明；零结果只是
当前多夹具点绑定模型中的一步载体筛查，不是 5 E 不可能性证明。
"""

from __future__ import annotations

import argparse
from itertools import combinations
from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    candidates,
    on_drawable,
)
from search_mannheim_regular_sequential_locus_reuse import (
    PROFILES,
    build_state,
)


def carries_centers(drawable, centers) -> bool:
    center_rows = tuple(centers.values())
    if len(center_rows) != 2:
        raise AssertionError("每个正规方向类必须恰有两个目标圆心")
    return all(
        on_drawable(first, drawable.kind, value)
        and on_drawable(second, drawable.kind, value)
        for first, second, value in zip(
            center_rows[0],
            center_rows[1],
            drawable.values,
            strict=True,
        )
    )


def search_state(selected_tuple, *, progress_interval: int):
    (
        state,
        _,
        target_centers,
        _,
        _,
        _,
        _,
    ) = build_state(frozenset(selected_tuple))
    old_carriers = {
        profile: tuple(
            drawable.drawable_id
            for drawable in state.drawables
            if carries_centers(drawable, centers)
        )
        for profile, centers in target_centers.items()
    }
    start = monotonic()
    checked = 0
    hits = []
    for candidate in candidates(state):
        checked += 1
        profiles = tuple(
            profile
            for profile, centers in target_centers.items()
            if carries_centers(candidate.drawable, centers)
        )
        if profiles:
            hits.append(
                {
                    "candidate": candidate.describe(),
                    "kind": candidate.kind,
                    "profiles": profiles,
                }
            )
        if progress_interval and checked % progress_interval == 0:
            print(
                "progress",
                {
                    "checked": checked,
                    "completed": selected_tuple,
                    "hits": len(hits),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    report = {
        "completed": selected_tuple,
        "remaining": tuple(target_centers),
        "initial_points": len(state.points),
        "initial_drawables": len(state.drawables),
        "candidates": checked,
        "old_carriers": old_carriers,
        "hits": len(hits),
        "elapsed_seconds": round(monotonic() - start, 3),
    }
    print("center_carrier_state", report, flush=True)
    for hit in hits:
        print("candidate", {"completed": selected_tuple, **hit}, flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--completed-count",
        type=int,
        choices=(0, 1, 2, 3),
        default=0,
    )
    parser.add_argument("--completed-profiles")
    parser.add_argument("--progress-interval", type=int, default=250_000)
    args = parser.parse_args()

    selected_sets = tuple(combinations(PROFILES, args.completed_count))
    if args.completed_profiles:
        selected = tuple(args.completed_profiles.split(","))
        if (
            len(selected) not in {1, 2, 3}
            or len(set(selected)) != len(selected)
            or any(profile not in PROFILES for profile in selected)
        ):
            raise SystemExit(
                "--completed-profiles 必须给出一至三个不同方向类"
            )
        selected_sets = (selected,)

    reports = [
        search_state(selected, progress_interval=args.progress_interval)
        for selected in selected_sets
    ]
    print(
        "regular_center_carrier_1e_search",
        {
            "samples": 3,
            "completed_count": len(selected_sets[0]),
            "states": len(reports),
            "candidates": sum(report["candidates"] for report in reports),
            "old_carriers": sum(
                sum(bool(rows) for rows in report["old_carriers"].values())
                for report in reports
            ),
            "hits": sum(report["hits"] for report in reports),
            "elapsed_seconds": round(
                sum(report["elapsed_seconds"] for report in reports), 3
            ),
        },
    )


if __name__ == "__main__":
    main()
