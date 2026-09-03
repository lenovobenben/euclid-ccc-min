"""筛查先行简单合并后，后行接触弦能否由 3 E 降为 2 E。

在三个不同半径的严格 55 E 饱和夹具上，先精确重放 ``P0`` 简单合并
方向类的两个目标、其接触弦、两目标圆心线与根心 ``S``。此时正式程序
对 ``P2/P1/P3`` 各画两条对角弦取得有限对角点，再画经过 ``S`` 的接触
弦，共 3 E。

本脚本补全先行状态的全部免费有限实交点，再枚举任意已有两点可画的
一条线或圆。若候选本身就是某条后行接触弦，或它与旧对象给出该接触
弦上不同于 ``S`` 的第二点，便可在至多 2 E 内完成该核心并把完整上界
降低 1 E。

搜索使用三个夹具上的同步浮点筛查。命中必须精确重放；零命中只排除
当前对象闭包与“一条辅助对象加目标线”的程序形状。
"""

from __future__ import annotations

from time import monotonic

from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from replay_mannheim_ordered_branches import branch_data
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_double_kp_global_sequential_2e import float_line
from search_mannheim_kp_parallel_center_locus_2e import (
    MERGE_FIXTURES,
    build_witnesses,
    candidate_has_witness,
    candidates,
    same_line,
    target_points,
)
from search_mannheim_root_center_2e import add_known_drawable
from search_parallel_3e import DrawableBundle, PointBundle, State


PROFILES = ("P2", "P1", "P3")


def build_state(*, include_root_line: bool = True):
    replays = []
    prefixes = []
    for index, (centers, radii) in enumerate(MERGE_FIXTURES):
        replay = KpCenterLocusReplay(
            f"merge_tail_tau_{index}",
            centers=centers,
            radii=radii,
            emit=False,
        )
        report = replay.run()
        if report["all_targets"] != 55:
            raise AssertionError("搜索夹具没有饱和当前 55 E 程序")
        paid = tuple(replay.objects.graph.paid_order)
        cutoff_name = (
            "P0_completed_pair_center_locus"
            if include_root_line
            else "P0_tau_for_S"
        )
        cutoff = paid.index(cutoff_name) + 1
        prefixes.append(paid[:cutoff])
        replays.append(replay)
    if len({len(prefix) for prefix in prefixes}) != 1:
        raise AssertionError("三个合并夹具的先行计费对象数不一致")

    points = tuple(
        PointBundle(
            f"O{index + 1}",
            tuple(
                tuple(float(coordinate) for coordinate in replay.centers[index])
                for replay in replays
            ),
        )
        for index in range(3)
    )
    input_drawables = []
    for node_id in ("Gamma1", "Gamma2", "Gamma3"):
        rows = tuple(object_value(replay.objects, node_id) for replay in replays)
        input_drawables.append(
            DrawableBundle(
                node_id,
                "circle",
                tuple(value for _, value in rows),
            )
        )
    state = State(points, tuple(input_drawables))
    for move, node_ids in enumerate(zip(*prefixes, strict=True), start=1):
        rows = tuple(
            object_value(replay.objects, node_id)
            for replay, node_id in zip(replays, node_ids, strict=True)
        )
        kinds = {kind for kind, _ in rows}
        if len(kinds) != 1:
            raise AssertionError(f"第 {move} 个对象的种类不一致：{node_ids}")
        state = add_known_drawable(
            state,
            node_ids[0],
            rows[0][0],
            tuple(value for _, value in rows),
            move,
        )

    root_centers = tuple(
        point_value(replay.objects, "Mannheim_S_from_completed_merge_pair")
        for replay in replays
    )
    target_lines = {
        profile: tuple(
            float_line(branch_data(replay, profile)["tau"])
            for replay in replays
        )
        for profile in PROFILES
    }
    return state, target_lines, root_centers, prefixes[0]


def main() -> None:
    start = monotonic()
    state, target_lines, root_centers, paid_prefix = build_state()
    initial_target_points = target_points(state, target_lines, root_centers)
    witnesses = build_witnesses(state, target_lines, root_centers)
    hits = []
    checked = 0
    for candidate in candidates(state):
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
        if checked % 100_000 == 0:
            print(
                "progress",
                {
                    "checked": checked,
                    "hits": len(hits),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    print(
        "mannheim_merge_tail_tau_2e_search",
        {
            "samples": len(MERGE_FIXTURES),
            "paid_prefix_objects": len(paid_prefix),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "witness_objects": {
                profile: len(rows) for profile, rows in witnesses.items()
            },
            "candidates": checked,
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    print("initial_target_points", initial_target_points)
    for hit in hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
