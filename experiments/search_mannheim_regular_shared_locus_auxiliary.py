"""搜索正规 Mannheim 目标对之间共享的圆心线辅助对象。

先用当前 6 E 后缀完成一或两个方向类，再枚举当时由两个已有确定点可
画出的全部直线和圆。若同一个候选对象分别与已有对象相交，在两个未完成
方向类的目标圆心线上各产生一个非根心点，则再画两条圆心线只需 2 E，
两套原 ``辅助对象 + 圆心线`` 的 4 E 合计可以压成 3 E，给出 48 E
正规程序候选。

脚本同步使用三个严格正规 ``D8`` 夹具，不允许任意点。它是多夹具浮点
筛查，不是下界；命中仍须精确重放和符号证明。零结果只覆盖当前状态、
点绑定方式和“一条共享辅助对象 + 两条圆心线”的程序形状。
"""

from __future__ import annotations

import argparse
from itertools import combinations
from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    build_witnesses,
    candidate_has_witness,
    candidates,
    same_line,
)
from search_mannheim_regular_sequential_locus_reuse import (
    PROFILES,
    build_state,
)


def search_state(selected_tuple, *, progress_interval: int):
    selected = frozenset(selected_tuple)
    (
        state,
        target_lines,
        _,
        _,
        root_centers,
        _,
    ) = build_state(selected)
    if len(target_lines) < 2:
        raise AssertionError("共享辅助对象至少需要两个未完成方向类")

    witnesses = build_witnesses(state, target_lines, root_centers)
    start = monotonic()
    checked = 0
    profile_hits = {profile: 0 for profile in target_lines}
    multi_hits = []
    for candidate in candidates(state):
        checked += 1
        assisted = {}
        for profile, lines in target_lines.items():
            direct = candidate.kind == "line" and all(
                same_line(value, line)
                for value, line in zip(
                    candidate.drawable.values,
                    lines,
                    strict=True,
                )
            )
            witness = candidate_has_witness(candidate, witnesses[profile])
            if direct or witness is not None:
                profile_hits[profile] += 1
                assisted[profile] = {
                    "direct": direct,
                    "witness_object": witness,
                }
        if len(assisted) >= 2:
            multi_hits.append(
                {
                    "candidate": candidate.describe(),
                    "assisted": assisted,
                }
            )
        if progress_interval and checked % progress_interval == 0:
            print(
                "progress",
                {
                    "completed": selected_tuple,
                    "checked": checked,
                    "single_profile_hits": sum(profile_hits.values()),
                    "multi_hits": len(multi_hits),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    report = {
        "completed": selected_tuple,
        "remaining": tuple(target_lines),
        "initial_points": len(state.points),
        "initial_drawables": len(state.drawables),
        "witness_objects": {
            profile: len(rows) for profile, rows in witnesses.items()
        },
        "candidates": checked,
        "profile_hits": profile_hits,
        "multi_hits": len(multi_hits),
        "elapsed_seconds": round(monotonic() - start, 3),
    }
    print("shared_locus_state", report, flush=True)
    for hit in multi_hits:
        print("candidate", {"completed": selected_tuple, **hit}, flush=True)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--completed-count",
        type=int,
        choices=(1, 2),
        default=1,
        help="先完成的方向类数量；默认 1。",
    )
    parser.add_argument(
        "--completed-profiles",
        help="只检查一个逗号分隔的方向类组合，例如 P0,P2。",
    )
    parser.add_argument("--progress-interval", type=int, default=250_000)
    args = parser.parse_args()

    selected_sets = tuple(combinations(PROFILES, args.completed_count))
    if args.completed_profiles:
        selected = tuple(args.completed_profiles.split(","))
        if (
            len(selected) not in {1, 2}
            or len(set(selected)) != len(selected)
            or any(profile not in PROFILES for profile in selected)
        ):
            raise SystemExit("--completed-profiles 必须给出一或两个不同方向类")
        selected_sets = (selected,)

    reports = [
        search_state(selected, progress_interval=args.progress_interval)
        for selected in selected_sets
    ]
    print(
        "regular_shared_locus_auxiliary_search",
        {
            "samples": 3,
            "completed_count": len(selected_sets[0]),
            "states": len(reports),
            "candidates": sum(report["candidates"] for report in reports),
            "single_profile_hits": sum(
                sum(report["profile_hits"].values()) for report in reports
            ),
            "multi_hits": sum(report["multi_hits"] for report in reports),
            "elapsed_seconds": round(
                sum(report["elapsed_seconds"] for report in reports), 3
            ),
        },
    )


if __name__ == "__main__":
    main()
