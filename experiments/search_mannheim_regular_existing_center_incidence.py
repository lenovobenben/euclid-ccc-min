"""筛查已完成目标对的对象是否经过尚未完成的目标圆心。

正规 49 E 程序允许任意调整四个方向类的完成顺序。若先完成的目标圆、
接触半径或圆心线经过后续目标圆心，它就能充当后续定心所需的一个已有
载体，并改变 5 E 后缀的组合分类。

本脚本同步重放三个严格正规 ``D8`` 夹具的全部 15 个真子集状态，检查
尚未完成的目标圆心是否已经作为确定点出现，以及每个已有直线或圆是否
经过它。只做多夹具筛查；零结果不是一般参数下的不入射证明。
"""

from __future__ import annotations

from itertools import combinations

from search_mannheim_kp_parallel_center_locus_2e import on_drawable
from search_mannheim_regular_sequential_locus_reuse import (
    PROFILES,
    build_state,
    same_point,
)


def main() -> None:
    reports = []
    total_states = 0
    total_centers = 0
    point_hits = 0
    drawable_hits = 0
    for selected_count in range(len(PROFILES)):
        for selected_tuple in combinations(PROFILES, selected_count):
            state, _, target_centers, _, _, _, _ = build_state(
                frozenset(selected_tuple)
            )
            state_point_hits = []
            state_drawable_hits = []
            for profile, centers in target_centers.items():
                for key, values in centers.items():
                    total_centers += 1
                    matching_points = tuple(
                        point.point_id
                        for point in state.points
                        if all(
                            same_point(value, expected)
                            for value, expected in zip(
                                point.values,
                                values,
                                strict=True,
                            )
                        )
                    )
                    matching_drawables = tuple(
                        drawable.drawable_id
                        for drawable in state.drawables
                        if all(
                            on_drawable(value, drawable.kind, drawable_value)
                            for value, drawable_value in zip(
                                values,
                                drawable.values,
                                strict=True,
                            )
                        )
                    )
                    if matching_points:
                        state_point_hits.append(
                            (profile, key, matching_points)
                        )
                    if matching_drawables:
                        state_drawable_hits.append(
                            (profile, key, matching_drawables)
                        )
            total_states += 1
            point_hits += len(state_point_hits)
            drawable_hits += len(state_drawable_hits)
            if state_point_hits or state_drawable_hits:
                reports.append(
                    {
                        "completed": selected_tuple,
                        "point_hits": tuple(state_point_hits),
                        "drawable_hits": tuple(state_drawable_hits),
                    }
                )

    print(
        "regular_existing_center_incidence_search",
        {
            "samples": 3,
            "states": total_states,
            "remaining_center_checks": total_centers,
            "point_hits": point_hits,
            "drawable_hits": drawable_hits,
        },
    )
    for report in reports:
        print("candidate", report)


if __name__ == "__main__":
    main()
