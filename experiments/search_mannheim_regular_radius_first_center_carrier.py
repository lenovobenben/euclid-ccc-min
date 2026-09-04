"""搜索接触半径先行的 5 E 正规双目标后缀。

对一个正规方向类，先画一个或两个本来就必须画出的第三圆接触半径，
把它们与已有对象的有限实交点加入状态。随后枚举一步可画对象；若它同时
经过该方向类的两个目标圆心，就能以这三个非输出对象取得两个圆心，再画
两个输出圆，总成本为 5 E。

脚本同步使用三个严格正规 ``D8`` 夹具，只使用确定交点。零结果只覆盖
当前公共核心状态、点绑定方式和“接触半径先行 + 一步双圆心载体”的
程序形状，不是 5 E 下界。
"""

from __future__ import annotations

from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import candidates
from search_mannheim_regular_center_carrier_1e import carries_centers
from search_mannheim_regular_sequential_locus_reuse import (
    PROFILES,
    build_state,
    find_point_id,
)
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    State,
    apply_candidate,
)


def add_radius(
    state: State,
    profile: str,
    key: str,
    contact_values,
    radius_values,
) -> State:
    contact_id = find_point_id(state, contact_values)
    candidate = Candidate(
        "line",
        "O3",
        contact_id,
        DrawableBundle(
            f"required_radius_{profile}_{key}",
            "line",
            radius_values,
        ),
    )
    return apply_candidate(state, candidate, len(state.program) + 1)


def carrier_hits(state, centers):
    checked = 0
    hits = []
    for candidate in candidates(state):
        checked += 1
        if carries_centers(candidate.drawable, centers):
            hits.append(
                {
                    "candidate": candidate.describe(),
                    "kind": candidate.kind,
                }
            )
    return checked, hits


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
    start = monotonic()
    reports = []
    total_candidates = 0
    total_hits = 0

    for profile in PROFILES:
        radius_keys = tuple(required_radii[profile])
        if len(radius_keys) != 2:
            raise AssertionError(f"{profile} 必须恰有两条接触半径")
        one_radius_states = {}
        for key in radius_keys:
            state = add_radius(
                initial_state,
                profile,
                key,
                *required_radii[profile][key],
            )
            one_radius_states[key] = state
            checked, hits = carrier_hits(state, target_centers[profile])
            total_candidates += checked
            total_hits += len(hits)
            reports.append(
                {
                    "profile": profile,
                    "radii_drawn": (key,),
                    "points": len(state.points),
                    "candidates": checked,
                    "hits": hits,
                }
            )

        first_key, second_key = radius_keys
        state = add_radius(
            one_radius_states[first_key],
            profile,
            second_key,
            *required_radii[profile][second_key],
        )
        checked, hits = carrier_hits(state, target_centers[profile])
        total_candidates += checked
        total_hits += len(hits)
        reports.append(
            {
                "profile": profile,
                "radii_drawn": radius_keys,
                "points": len(state.points),
                "candidates": checked,
                "hits": hits,
            }
        )

    print(
        "regular_radius_first_center_carrier_search",
        {
            "samples": 3,
            "states": len(reports),
            "candidates": total_candidates,
            "hits": total_hits,
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for report in reports:
        print("radius_first_state", report)


if __name__ == "__main__":
    main()
