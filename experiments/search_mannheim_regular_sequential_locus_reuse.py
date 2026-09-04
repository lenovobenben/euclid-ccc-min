"""筛查正规 49 E 程序能否从已完成目标对免费取得后续圆心线上的点。

正规 Mannheim 程序的每个方向类使用 6 E 后缀：先用
``Line(O3,K')`` 和 ``Line(S,J)`` 画出目标圆心线，再画两条第三圆
接触半径和两个输出圆。若先完成若干方向类后，新增对象的自动交点闭包
已经在另一方向类的圆心线上给出一个不同于根心 ``S`` 的点，那么后者
可以直接画 ``Line(S,J)``，省去 ``Line(O3,K')``，把总上界从 49 E
降到 48 E。

本脚本在三个严格正规 ``D8`` 夹具上同步检查四个方向类的全部真子集，
其中包括尚未完成任何目标对的初始状态。除免费点复用外，还检查本来就
必须画出的第三圆接触半径能否帮助任一圆心线，以及标准 ``O3-K'``
辅助线能否跨方向复用。所有状态均精确重放当前合法程序；搜索点只来自
已有对象的有限实交点，不引入任意点。数值多夹具零结果不是下界，命中
必须另作符号证明和精确重放。
"""

from __future__ import annotations

from itertools import combinations
from math import hypot

from check_mannheim_degenerate_fixture import line_through as exact_line_through
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from replay_mannheim_ordered_branches import branch_data
from replay_mannheim_three_block_dependencies import collapse_point
from search_mannheim_center_locus_2e import (
    object_value,
    point_value,
)
from search_mannheim_double_kp_global_sequential_2e import (
    float_line,
    float_scalar,
)
from search_mannheim_kp_parallel_center_locus_2e import CoreOnlyReplay
from search_mannheim_root_center_2e import FIXTURES, add_known_drawable
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
)


PROFILES = ("P0", "P2", "P1", "P3")
TOLERANCE = 1e-7


class SelectedPairsReplay(CoreOnlyReplay):
    """完成指定方向类，其余方向类停在接触弦核心之后。"""

    selected_profiles: frozenset[str] = frozenset()

    def build_pair(self, profile, *args, **kwargs) -> None:
        if profile in self.selected_profiles:
            KpCenterLocusReplay.build_pair(self, profile, *args, **kwargs)


def same_point(first, second) -> bool:
    return hypot(first[0] - second[0], first[1] - second[1]) <= TOLERANCE


def on_line(point, line) -> bool:
    residual = abs(line[0] * point[0] + line[1] * point[1] + line[2])
    return residual <= TOLERANCE * max(1.0, hypot(*point))


def build_state(selected_profiles: frozenset[str]):
    replays = []
    for index, (centers, radii) in enumerate(FIXTURES):
        replay = SelectedPairsReplay(
            f"regular_sequential_locus_{index}",
            centers=centers,
            radii=radii,
            emit=False,
        )
        replay.selected_profiles = selected_profiles
        replay.run()
        branches = {
            profile: branch_data(replay, profile)["kind"]
            for profile in PROFILES
        }
        if set(branches.values()) != {"regular"}:
            raise AssertionError(f"顺序复用夹具必须严格正规：{branches}")
        if len(replay.targets) != 2 * len(selected_profiles):
            raise AssertionError("已完成目标圆数与所选方向类不符")
        replays.append(replay)

    paid_sequences = tuple(
        tuple(replay.objects.graph.paid_order) for replay in replays
    )
    if any(sequence != paid_sequences[0] for sequence in paid_sequences[1:]):
        raise AssertionError("三个夹具的计费对象顺序不一致")
    paid_sequence = paid_sequences[0]
    expected_paid = 25 + 6 * len(selected_profiles)
    if len(paid_sequence) != expected_paid:
        raise AssertionError(
            f"所选方向类应产生 {expected_paid} 个计费对象，"
            f"实际为 {len(paid_sequence)}"
        )

    points = tuple(
        PointBundle(
            f"O{center_index + 1}",
            tuple(
                tuple(float(coordinate) for coordinate in centers[center_index])
                for centers, _ in FIXTURES
            ),
        )
        for center_index in range(3)
    )
    input_drawables = []
    for node_id in ("Gamma1", "Gamma2", "Gamma3"):
        rows = tuple(object_value(replay.objects, node_id) for replay in replays)
        if {kind for kind, _ in rows} != {"circle"}:
            raise AssertionError("输入对象不是圆")
        input_drawables.append(
            DrawableBundle(
                node_id,
                "circle",
                tuple(value for _, value in rows),
            )
        )
    state = State(points, tuple(input_drawables))
    for move, node_ids in enumerate(zip(*paid_sequences, strict=True), start=1):
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
        point_value(replay.objects, "Mannheim_S_center_locus")
        for replay in replays
    )
    target_lines = {}
    required_radii = {}
    kp_auxiliaries = {}
    for profile in PROFILES:
        if profile in selected_profiles:
            continue
        rows = []
        sample_targets = []
        sample_data = []
        for replay in replays:
            data = replay.core_data[profile]
            sample_data.append(data)
            targets = replay.verify_pair(
                profile,
                data["tau"],
                allow_repeated_physical_signs=True,
            )
            centers = tuple(
                collapse_point(target["center"])
                for target in targets.values()
            )
            rows.append(float_line(exact_line_through(*centers)))
            sample_targets.append(targets)
        target_lines[profile] = tuple(rows)
        target_keys = tuple(sorted(sample_targets[0]))
        if any(tuple(sorted(targets)) != target_keys for targets in sample_targets):
            raise AssertionError(f"{profile} 的目标键在夹具间不一致")
        required_radii[profile] = {}
        for key in target_keys:
            contacts = tuple(
                collapse_point(targets[key]["contact_3"])
                for targets in sample_targets
            )
            contact_values = tuple(
                tuple(float_scalar(coordinate) for coordinate in contact)
                for contact in contacts
            )
            radius_values = tuple(
                float_line(exact_line_through(replay.o3, contact))
                for replay, contact in zip(replays, contacts, strict=True)
            )
            required_radii[profile][key] = (contact_values, radius_values)
        kps = tuple(data["Kp"][2] for data in sample_data)
        if any(kp is None for kp in kps):
            raise AssertionError(f"正规 {profile} 的 K' 必须有限")
        kp_values = tuple(
            tuple(float_scalar(coordinate) for coordinate in kp)
            for kp in kps
        )
        kp_radial_values = tuple(
            float_line(exact_line_through(replay.o3, kp))
            for replay, kp in zip(replays, kps, strict=True)
        )
        kp_auxiliaries[profile] = (kp_values, kp_radial_values)
    aligned_contacts = []
    for profile, keyed_radii in required_radii.items():
        for key, (contact_values, _) in keyed_radii.items():
            matches = tuple(
                point.point_id
                for point in state.points
                if all(
                    same_point(value, expected)
                    for value, expected in zip(
                        point.values,
                        contact_values,
                        strict=True,
                    )
                )
            )
            if len(matches) > 1:
                raise AssertionError(
                    f"{profile}:{key} 接触点有多个跨夹具绑定：{matches}"
                )
            if not matches:
                aligned_contacts.append(
                    PointBundle(
                        f"{profile}_{key}_contact_3",
                        contact_values,
                    )
                )
    for profile, (kp_values, _) in kp_auxiliaries.items():
        matches = tuple(
            point.point_id
            for point in state.points
            if all(
                same_point(value, expected)
                for value, expected in zip(
                    point.values,
                    kp_values,
                    strict=True,
                )
            )
        )
        if len(matches) > 1:
            raise AssertionError(
                f"{profile}:K' 有多个跨夹具绑定：{matches}"
            )
        if not matches:
            aligned_contacts.append(
                PointBundle(f"{profile}_Kp", kp_values)
            )
    if aligned_contacts:
        state = State(
            state.points + tuple(aligned_contacts),
            state.drawables,
            state.program,
        )
    return (
        state,
        target_lines,
        required_radii,
        kp_auxiliaries,
        root_centers,
        paid_sequence,
    )


def source_names(point_id: str, paid_sequence: tuple[str, ...]):
    source_ids = point_id.split("[", 1)[0].split("&")
    result = []
    for source_id in source_ids:
        if source_id.startswith("object_"):
            move = int(source_id.removeprefix("object_"))
            result.append(paid_sequence[move - 1])
        else:
            result.append(source_id)
    return tuple(result)


def find_hits(state, target_lines, root_centers, paid_sequence):
    hits = {profile: [] for profile in target_lines}
    for point in state.points:
        if all(
            same_point(value, root)
            for value, root in zip(point.values, root_centers, strict=True)
        ):
            continue
        for profile, lines in target_lines.items():
            if all(
                on_line(value, line)
                for value, line in zip(point.values, lines, strict=True)
            ):
                hits[profile].append(
                    {
                        "point": point.point_id,
                        "sources": source_names(point.point_id, paid_sequence),
                    }
                )
    return {profile: rows for profile, rows in hits.items() if rows}


def find_point_id(state, values):
    matches = tuple(
        point.point_id
        for point in state.points
        if all(
            same_point(value, expected)
            for value, expected in zip(point.values, values, strict=True)
        )
    )
    if len(matches) != 1:
        raise AssertionError(f"接触点应有唯一状态绑定，实际为 {matches}")
    return matches[0]


def find_required_radius_hits(
    state,
    target_lines,
    required_radii,
    root_centers,
    paid_sequence,
):
    hits = []
    move = len(paid_sequence) + 1
    for profile, keyed_radii in required_radii.items():
        for key, (contact_values, radius_values) in keyed_radii.items():
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
            state_one = apply_candidate(state, candidate, move)
            new_points = state_one.points[len(state.points) :]
            for target_profile, lines in target_lines.items():
                witnesses = tuple(
                    {
                        "point": point.point_id,
                        "sources": source_names(
                            point.point_id,
                            paid_sequence
                            + (f"required_radius_{profile}_{key}",),
                        ),
                    }
                    for point in new_points
                    if not all(
                        same_point(value, root)
                        for value, root in zip(
                            point.values,
                            root_centers,
                            strict=True,
                        )
                    )
                    and all(
                        on_line(value, line)
                        for value, line in zip(
                            point.values,
                            lines,
                            strict=True,
                        )
                    )
                )
                if witnesses:
                    hits.append(
                        {
                            "radius_profile": profile,
                            "target": key,
                            "assisted_profile": target_profile,
                            "radius": candidate.describe(),
                            "witnesses": witnesses,
                        }
                    )
    return hits


def find_kp_auxiliary_cross_hits(
    state,
    target_lines,
    kp_auxiliaries,
    root_centers,
    paid_sequence,
):
    hits = []
    move = len(paid_sequence) + 1
    for profile, (kp_values, radial_values) in kp_auxiliaries.items():
        kp_id = find_point_id(state, kp_values)
        candidate = Candidate(
            "line",
            "O3",
            kp_id,
            DrawableBundle(
                f"kp_auxiliary_{profile}",
                "line",
                radial_values,
            ),
        )
        state_one = apply_candidate(state, candidate, move)
        new_points = state_one.points[len(state.points) :]
        assisted = {}
        for target_profile, lines in target_lines.items():
            witnesses = tuple(
                {
                    "point": point.point_id,
                    "sources": source_names(
                        point.point_id,
                        paid_sequence + (f"kp_auxiliary_{profile}",),
                    ),
                }
                for point in new_points
                if not all(
                    same_point(value, root)
                    for value, root in zip(
                        point.values,
                        root_centers,
                        strict=True,
                    )
                )
                and all(
                    on_line(value, line)
                    for value, line in zip(
                        point.values,
                        lines,
                        strict=True,
                    )
                )
            )
            if witnesses:
                assisted[target_profile] = witnesses
        if profile not in assisted:
            raise AssertionError(
                f"{profile} 的标准 O3-K' 辅助线没有恢复自身圆心线截点"
            )
        cross = {
            target_profile: witnesses
            for target_profile, witnesses in assisted.items()
            if target_profile != profile
        }
        if cross:
            hits.append(
                {
                    "auxiliary_profile": profile,
                    "line": candidate.describe(),
                    "assisted_profiles": cross,
                }
            )
    return hits


def main() -> None:
    reports = []
    total_states = 0
    total_points = 0
    total_hits = 0
    required_radius_candidates = 0
    required_radius_hits = 0
    kp_auxiliary_candidates = 0
    kp_auxiliary_cross_hits = 0
    for selected_count in range(0, len(PROFILES)):
        for selected_tuple in combinations(PROFILES, selected_count):
            selected = frozenset(selected_tuple)
            (
                state,
                target_lines,
                required_radii,
                kp_auxiliaries,
                root_centers,
                paid_sequence,
            ) = build_state(selected)
            hits = find_hits(
                state,
                target_lines,
                root_centers,
                paid_sequence,
            )
            kp_hits = find_kp_auxiliary_cross_hits(
                state,
                target_lines,
                kp_auxiliaries,
                root_centers,
                paid_sequence,
            )
            radius_hits = find_required_radius_hits(
                state,
                target_lines,
                required_radii,
                root_centers,
                paid_sequence,
            )
            total_states += 1
            total_points += len(state.points)
            total_hits += sum(len(rows) for rows in hits.values())
            required_radius_candidates += sum(
                len(rows) for rows in required_radii.values()
            )
            required_radius_hits += len(radius_hits)
            kp_auxiliary_candidates += len(kp_auxiliaries)
            kp_auxiliary_cross_hits += len(kp_hits)
            if hits or radius_hits or kp_hits:
                reports.append(
                    {
                        "completed": selected_tuple,
                        "paid_objects": len(paid_sequence),
                        "points": len(state.points),
                        "free_point_hits": hits,
                        "required_radius_hits": radius_hits,
                        "kp_auxiliary_cross_hits": kp_hits,
                    }
                )

    print(
        "regular_sequential_locus_reuse_search",
        {
            "samples": len(FIXTURES),
            "states": total_states,
            "bundled_points_checked": total_points,
            "free_point_hits": total_hits,
            "required_radius_candidates": required_radius_candidates,
            "required_radius_hits": required_radius_hits,
            "kp_auxiliary_candidates": kp_auxiliary_candidates,
            "kp_auxiliary_cross_hits": kp_auxiliary_cross_hits,
        },
    )
    for report in reports:
        print("candidate", report)


if __name__ == "__main__":
    main()
