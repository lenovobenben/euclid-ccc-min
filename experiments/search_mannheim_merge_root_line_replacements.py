"""搜索简单合并饱和分支中定根线的等成本替换。

``P0`` 两目标和接触弦完成后，当前程序画两目标圆心线，与接触弦相交
取得根心 ``S``。三个不同半径夹具的同步筛查还能找到其它由已有两点
直接画出、并恒过 ``S`` 的直线。每条候选都可用同样 1 E 替换圆心线，
但会带来一组不同的免费交点。

本脚本逐一加入这些候选定根线，再检查 ``P2/P1/P3`` 的后行接触弦能否
由“一条辅助对象 + 目标线”在 2 E 内完成。命中必须另作精确重放；零
命中只排除这组等成本替换与当前两步程序形状。
"""

from __future__ import annotations

from time import monotonic

from search_mannheim_kp_parallel_center_locus_2e import (
    build_witnesses,
    candidate_has_witness,
    candidates,
    on_drawable,
    same_line,
    target_points,
)
from search_mannheim_merge_tail_tau_2e import build_state
from search_parallel_3e import apply_candidate


def main() -> None:
    start = monotonic()
    state, target_lines, root_centers, paid_prefix = build_state(
        include_root_line=False
    )
    root_candidates = tuple(
        candidate
        for candidate in candidates(state)
        if all(
            on_drawable(root, candidate.kind, value)
            for root, value in zip(
                root_centers,
                candidate.drawable.values,
                strict=True,
            )
        )
    )
    if len(root_candidates) != 4:
        raise AssertionError(
            f"预期四条候选定根线，实际得到 {len(root_candidates)}"
        )

    reports = []
    for replacement_index, replacement in enumerate(root_candidates, start=1):
        replaced = apply_candidate(
            state,
            replacement,
            len(paid_prefix) + 1,
        )
        initial_target_points = target_points(
            replaced,
            target_lines,
            root_centers,
        )
        witnesses = build_witnesses(
            replaced,
            target_lines,
            root_centers,
        )
        hits = []
        checked = 0
        for candidate in candidates(replaced):
            checked += 1
            direct_profiles = tuple(
                profile
                for profile, lines in target_lines.items()
                if candidate.kind == "line"
                and all(
                    same_line(value, line)
                    for value, line in zip(
                        candidate.drawable.values,
                        lines,
                        strict=True,
                    )
                )
            )
            witness_profiles = {
                profile: candidate_has_witness(candidate, witnesses[profile])
                for profile in target_lines
            }
            if direct_profiles or any(witness_profiles.values()):
                hits.append(
                    {
                        "candidate": candidate.describe(),
                        "direct_profiles": direct_profiles,
                        "witnesses": witness_profiles,
                    }
                )
        report = {
            "replacement": replacement.describe(),
            "initial_points": len(replaced.points),
            "initial_drawables": len(replaced.drawables),
            "initial_target_points": initial_target_points,
            "candidates": checked,
            "hits": hits,
        }
        reports.append(report)
        print(
            "replacement_progress",
            {
                "index": replacement_index,
                "replacement": replacement.describe(),
                "candidates": checked,
                "hits": len(hits),
                "elapsed_seconds": round(monotonic() - start, 3),
            },
            flush=True,
        )

    print(
        "mannheim_merge_root_line_replacements",
        {
            "samples": len(root_centers),
            "paid_prefix_objects": len(paid_prefix),
            "replacement_candidates": len(root_candidates),
            "total_candidates": sum(report["candidates"] for report in reports),
            "hits": sum(len(report["hits"]) for report in reports),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for report in reports:
        for hit in report["hits"]:
            print(
                "candidate",
                {"replacement": report["replacement"], **hit},
            )


if __name__ == "__main__":
    main()
