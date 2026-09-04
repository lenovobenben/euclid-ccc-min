"""搜索正规 5 E 交错后缀中的两步双圆心载体。

考虑以下必要形状：前两步相交得到第一个目标圆心，第三步画出第一个
目标圆，第四步由双圆心载体与另一对象相交得到第二个圆心，第五步画出
第二个目标圆。由于目标圆彼此不经过对方圆心，前两步之一必须同时经过
两个圆心。

一步载体已由 ``search_mannheim_regular_center_carrier_1e.py`` 排除。
本脚本枚举所有一步可画且经过任一目标圆心的首对象；应用它的确定交点
闭包后，再枚举所有一步可画的双圆心载体。只使用三个严格正规 ``D8``
夹具中的已有点和有限实交点。零结果不是 5 E 下界。
"""

from __future__ import annotations

from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    candidates,
    on_drawable,
)
from search_mannheim_regular_center_carrier_1e import carries_centers
from search_mannheim_regular_sequential_locus_reuse import (
    PROFILES,
    build_state,
)
from search_parallel_3e import apply_candidate


def carries_point(drawable, points) -> bool:
    return all(
        on_drawable(point, drawable.kind, value)
        for point, value in zip(points, drawable.values, strict=True)
    )


def main() -> None:
    (
        state,
        _,
        target_centers,
        _,
        _,
        _,
        _,
    ) = build_state(frozenset())
    old_center_objects = {
        f"{profile}:{key}": tuple(
            drawable.drawable_id
            for drawable in state.drawables
            if carries_point(drawable, points)
        )
        for profile, centers in target_centers.items()
        for key, points in centers.items()
    }

    start = monotonic()
    first_candidates = tuple(candidates(state))
    first_hits = []
    for candidate in first_candidates:
        centers_hit = tuple(
            (profile, key)
            for profile, centers in target_centers.items()
            for key, points in centers.items()
            if carries_point(candidate.drawable, points)
        )
        if centers_hit:
            first_hits.append((candidate, centers_hit))

    second_candidates = 0
    carrier_hits = []
    for first, centers_hit in first_hits:
        state_one = apply_candidate(state, first, len(state.program) + 1)
        for second in candidates(state_one):
            second_candidates += 1
            profiles = tuple(
                profile
                for profile, centers in target_centers.items()
                if carries_centers(second.drawable, centers)
                and any(hit_profile == profile for hit_profile, _ in centers_hit)
            )
            if profiles:
                carrier_hits.append(
                    {
                        "first": first.describe(),
                        "first_centers": centers_hit,
                        "carrier": second.describe(),
                        "profiles": profiles,
                    }
                )

    report = {
        "samples": 3,
        "initial_points": len(state.points),
        "initial_drawables": len(state.drawables),
        "old_center_objects": old_center_objects,
        "first_candidates": len(first_candidates),
        "first_center_hits": len(first_hits),
        "second_candidates": second_candidates,
        "carrier_hits": len(carrier_hits),
        "elapsed_seconds": round(monotonic() - start, 3),
    }
    print("regular_two_step_center_carrier_search", report)
    for first, centers_hit in first_hits:
        print(
            "first_center_object",
            {
                "candidate": first.describe(),
                "centers": centers_hit,
            },
        )
    for hit in carrier_hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
