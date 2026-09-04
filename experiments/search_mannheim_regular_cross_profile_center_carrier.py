"""搜索正规核心上连接任意两个目标圆心的短载体程序。

此前的载体搜索只要求对象同时经过同一方向类的两个圆心。本脚本把目标
扩展为八个目标圆心中的任意一对，以寻找跨方向共享定心基础设施。

初态一步候选只有八条必需第三圆接触半径能经过目标圆心。脚本先在画出
其中一条半径后枚举第二步，检查它能否直接成为跨圆心载体，并核对所有
经过其它圆心的第二步是否仍只是另一条必需半径。随后对八条半径的 28 个
无序组合逐一应用确定交点闭包，枚举第三步双圆心载体。只同步使用三个
严格正规 ``D8`` 夹具；命中须另作精确重放和符号证明，零结果不是下界。
"""

from __future__ import annotations

from itertools import combinations
from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    candidates,
    on_drawable,
    same_drawable,
)
from search_mannheim_regular_radius_first_center_carrier import add_radius
from search_mannheim_regular_sequential_locus_reuse import build_state


def carries_point(drawable, points) -> bool:
    return all(
        on_drawable(point, drawable.kind, value)
        for point, value in zip(points, drawable.values, strict=True)
    )


def carries_pair(drawable, first, second) -> bool:
    return carries_point(drawable, first["center"]) and carries_point(
        drawable,
        second["center"],
    )


def is_radius_for(drawable, target) -> bool:
    return drawable.kind == "line" and all(
        same_drawable("line", value, expected)
        for value, expected in zip(
            drawable.values,
            target["radius_values"],
            strict=True,
        )
    )


def target_name(target) -> str:
    return f'{target["profile"]}:{target["key"]}'


def main() -> None:
    (
        initial_state,
        _,
        target_centers,
        required_radii,
        _,
        _,
        _,
    ) = build_state(frozenset())
    targets = tuple(
        {
            "profile": profile,
            "key": key,
            "center": target_centers[profile][key],
            "contact_values": required_radii[profile][key][0],
            "radius_values": required_radii[profile][key][1],
        }
        for profile in target_centers
        for key in target_centers[profile]
    )
    if len(targets) != 8:
        raise AssertionError("正规核心必须恰有八个目标圆心")

    start = monotonic()
    one_radius_states = {}
    second_candidates = 0
    second_center_hits = []
    second_carrier_hits = []
    unexpected_second_hits = []
    for first in targets:
        first_name = target_name(first)
        state_one = add_radius(
            initial_state,
            first["profile"],
            first["key"],
            first["contact_values"],
            first["radius_values"],
        )
        one_radius_states[first_name] = state_one
        for second in candidates(state_one):
            second_candidates += 1
            other_hits = tuple(
                other
                for other in targets
                if other is not first
                and carries_point(second.drawable, other["center"])
            )
            if not other_hits:
                continue
            second_center_hits.append(
                {
                    "first": first_name,
                    "candidate": second.describe(),
                    "centers": tuple(target_name(other) for other in other_hits),
                }
            )
            if carries_point(second.drawable, first["center"]):
                second_carrier_hits.append(
                    {
                        "first": first_name,
                        "candidate": second.describe(),
                        "centers": tuple(target_name(other) for other in other_hits),
                    }
                )
            for other in other_hits:
                if not is_radius_for(second.drawable, other):
                    unexpected_second_hits.append(
                        {
                            "first": first_name,
                            "candidate": second.describe(),
                            "center": target_name(other),
                        }
                    )

    third_candidates = 0
    old_third_carriers = []
    new_third_carriers = []
    pair_reports = []
    for pair_index, (first, second) in enumerate(combinations(targets, 2), start=1):
        state_one = one_radius_states[target_name(first)]
        state_two = add_radius(
            state_one,
            second["profile"],
            second["key"],
            second["contact_values"],
            second["radius_values"],
        )
        old_hits = tuple(
            drawable.drawable_id
            for drawable in state_two.drawables
            if carries_pair(drawable, first, second)
        )
        pair_candidates = 0
        pair_hits = []
        for third in candidates(state_two):
            pair_candidates += 1
            third_candidates += 1
            if carries_pair(third.drawable, first, second):
                pair_hits.append(third.describe())
                new_third_carriers.append(
                    {
                        "centers": (target_name(first), target_name(second)),
                        "carrier": third.describe(),
                    }
                )
        if old_hits:
            old_third_carriers.append(
                {
                    "centers": (target_name(first), target_name(second)),
                    "carriers": old_hits,
                }
            )
        pair_reports.append(
            {
                "centers": (target_name(first), target_name(second)),
                "same_profile": first["profile"] == second["profile"],
                "candidates": pair_candidates,
                "old_hits": len(old_hits),
                "new_hits": len(pair_hits),
            }
        )
        print(
            "progress",
            {
                "pair": pair_index,
                "pair_total": 28,
                **pair_reports[-1],
                "elapsed_seconds": round(monotonic() - start, 3),
            },
            flush=True,
        )

    print(
        "regular_cross_profile_center_carrier_search",
        {
            "samples": 3,
            "targets": len(targets),
            "one_radius_states": len(one_radius_states),
            "second_candidates": second_candidates,
            "second_center_hits": len(second_center_hits),
            "unexpected_second_hits": len(unexpected_second_hits),
            "second_carrier_hits": len(second_carrier_hits),
            "two_radius_states": len(pair_reports),
            "same_profile_states": sum(
                report["same_profile"] for report in pair_reports
            ),
            "cross_profile_states": sum(
                not report["same_profile"] for report in pair_reports
            ),
            "third_candidates": third_candidates,
            "old_third_carriers": len(old_third_carriers),
            "new_third_carriers": len(new_third_carriers),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for hit in unexpected_second_hits:
        print("unexpected_second_center_object", hit)
    for hit in second_carrier_hits:
        print("second_step_carrier", hit)
    for hit in old_third_carriers:
        print("old_third_step_carrier", hit)
    for hit in new_third_carriers:
        print("third_step_carrier", hit)


if __name__ == "__main__":
    main()
