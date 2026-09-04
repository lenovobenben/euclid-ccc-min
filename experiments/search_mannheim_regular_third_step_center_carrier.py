"""搜索正规 5 E 双目标后缀中最后画出的双圆心载体。

若两个输出圆都放在最后，5 E 后缀只剩三个辅助对象。初态没有对象经过
任一目标圆心，所以每个圆心都必须是至少两个新对象的交点。三个辅助
对象中必有一个同时经过两个圆心。

一步载体与第二步载体已有独立筛查。本脚本覆盖剩余顺序：第一步经过
一个目标圆心，第二步经过同方向类的另一个圆心，第三步才成为双圆心
载体。第一步候选由正规核心的已有点定义；以后每一步只使用新旧对象的
确定有限实交点。命中必须另作精确重放和符号证明；零结果只覆盖三个
严格正规 ``D8`` 夹具和当前多夹具点绑定，不是 5 E 下界。
"""

from __future__ import annotations

from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    candidates,
    on_drawable,
)
from search_mannheim_regular_center_carrier_1e import carries_centers
from search_mannheim_regular_sequential_locus_reuse import build_state
from search_parallel_3e import apply_candidate


def carries_point(drawable, points) -> bool:
    return all(
        on_drawable(point, drawable.kind, value)
        for point, value in zip(points, drawable.values, strict=True)
    )


def main() -> None:
    (
        initial_state,
        _,
        target_centers,
        _,
        _,
        _,
        _,
    ) = build_state(frozenset())

    start = monotonic()
    first_candidates = tuple(candidates(initial_state))
    first_hits = []
    for first in first_candidates:
        centers_hit = tuple(
            (profile, key)
            for profile, centers in target_centers.items()
            for key, points in centers.items()
            if carries_point(first.drawable, points)
        )
        if centers_hit:
            first_hits.append((first, centers_hit))

    second_candidates = 0
    second_hits = []
    third_candidates = 0
    old_carrier_hits = []
    carrier_hits = []
    for first, centers_hit in first_hits:
        state_one = apply_candidate(
            initial_state,
            first,
            len(initial_state.program) + 1,
        )
        for second in candidates(state_one):
            second_candidates += 1
            opposite_hits = []
            for profile, first_key in centers_hit:
                for second_key, points in target_centers[profile].items():
                    if second_key != first_key and carries_point(second.drawable, points):
                        opposite_hits.append((profile, first_key, second_key))
            if not opposite_hits:
                continue
            second_hits.append(
                {
                    "first": first.describe(),
                    "second": second.describe(),
                    "pairs": tuple(opposite_hits),
                }
            )
            state_two = apply_candidate(
                state_one,
                second,
                len(state_one.program) + 1,
            )
            relevant_profiles = tuple(
                dict.fromkeys(profile for profile, _, _ in opposite_hits)
            )
            for drawable in state_two.drawables:
                profiles = tuple(
                    profile
                    for profile in relevant_profiles
                    if carries_centers(drawable, target_centers[profile])
                )
                if profiles:
                    old_carrier_hits.append(
                        {
                            "first": first.describe(),
                            "second": second.describe(),
                            "carrier": drawable.drawable_id,
                            "profiles": profiles,
                        }
                    )
            for third in candidates(state_two):
                third_candidates += 1
                profiles = tuple(
                    profile
                    for profile in relevant_profiles
                    if carries_centers(third.drawable, target_centers[profile])
                )
                if profiles:
                    carrier_hits.append(
                        {
                            "first": first.describe(),
                            "second": second.describe(),
                            "carrier": third.describe(),
                            "profiles": profiles,
                        }
                    )

    print(
        "regular_third_step_center_carrier_search",
        {
            "samples": 3,
            "initial_points": len(initial_state.points),
            "initial_drawables": len(initial_state.drawables),
            "first_candidates": len(first_candidates),
            "first_center_hits": len(first_hits),
            "second_candidates": second_candidates,
            "second_opposite_center_hits": len(second_hits),
            "third_candidates": third_candidates,
            "old_carrier_hits": len(old_carrier_hits),
            "new_carrier_hits": len(carrier_hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for hit in second_hits:
        print("opposite_center_pair", hit)
    for hit in old_carrier_hits:
        print("old_carrier", hit)
    for hit in carrier_hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
