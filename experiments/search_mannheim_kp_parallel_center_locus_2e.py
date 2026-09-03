"""筛查 ``K'`` 无穷远块能否用 2 E 提前画出目标圆心线。

第 8.19 节在有限 ``K'`` 时用 ``Line(O3,K')`` 与 ``O1O2`` 的交点
取得圆心线上的第二点。这里针对唯一可能同时出现的
``P2/P3:Kp_parallel`` 分支，精确重放当前 55 E 程序在所有目标后缀
之前已经画出的对象，并补全它们之间的免费有限实交点。

随后枚举任意两个已有点可画的一条直线或圆。若新对象本身就是目标
圆心线，或它与某个旧对象的免费交点给出圆心线上严格不同于根心 ``S``
的第二点，再画 ``Line(S,J)`` 即可在 2 E 内取得圆心线，使该方向类的
双目标后缀由 7 E 降为 6 E。

这是三个严格双平行样本上的浮点筛查，不是 2 E 下界。命中必须另作
精确重放；零命中只排除当前对象闭包与上述“一条辅助对象加圆心线”
程序形状。
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
from math import hypot
from time import monotonic

from check_mannheim_degenerate_fixture import line_through as exact_line_through
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from replay_mannheim_ordered_branches import branch_data
from replay_mannheim_three_block_dependencies import collapse_point
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_double_kp_global_sequential_2e import float_line
from search_mannheim_root_center_2e import add_known_drawable
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    PointBundle,
    State,
    circle_through,
    drawable_signature,
    intersections,
    line_through,
)


PROFILES = ("P2", "P3")
TOLERANCE = 1e-7
F = Fraction
FIXTURES = (
    (
        ((F(0), F(0)), (F(13, 2), F(0)), (F(9, 2), F(2))),
        (F(3), F(4, 3), F(1)),
    ),
    (
        ((F(0), F(0)), (F(25, 4), F(0)), (F(4), F(2))),
        (F(8, 3), F(3, 2), F(1)),
    ),
    (
        (
            (F(0), F(0)),
            (F(143, 12), F(0)),
            (F(33, 4), F(9, 2)),
        ),
        (F(9, 2), F(2), F(1)),
    ),
)


def merge_fixture(p, q, distance, height):
    return (
        (
            (F(0), F(0)),
            (distance, F(0)),
            (distance * p / (p + q), height),
        ),
        (p, q, F(1)),
    )


MERGE_FIXTURES = (
    merge_fixture(F(13, 9), F(17, 15), F(21112, 2115), F(220, 47)),
    merge_fixture(F(13, 7), F(25, 19), F(2743, 399), F(17, 6)),
    merge_fixture(F(37, 18), F(13, 8), F(1325, 216), F(20, 9)),
)


class CoreOnlyReplay(KpCenterLocusReplay):
    """沿用正式程序，但在目标后缀开始处停止。"""

    def build_pair(self, *args, **kwargs) -> None:
        return None

    def audit_center_locus(
        self,
        data,
        input_cost,
        seed_costs,
        tail_core_costs,
        suffix_costs,
    ):
        self.core_data = data
        return {
            "input_cost": input_cost,
            "seed_costs": seed_costs,
            "tail_core_costs": tail_core_costs,
            "suffix_costs": suffix_costs,
        }


class FinitePairsFirstReplay(CoreOnlyReplay):
    """在平行方向类之前完成已有 6 E 程序的两个有限方向类。"""

    def build_pair(self, profile, *args, **kwargs) -> None:
        if profile in {"P0", "P1"}:
            KpCenterLocusReplay.build_pair(self, profile, *args, **kwargs)


def same_point(first, second) -> bool:
    return hypot(first[0] - second[0], first[1] - second[1]) <= TOLERANCE


def on_line(point, line) -> bool:
    residual = abs(line[0] * point[0] + line[1] * point[1] + line[2])
    return residual <= TOLERANCE * max(1.0, hypot(*point))


def same_line(first, second) -> bool:
    return max(abs(a - b) for a, b in zip(first, second, strict=True)) <= 1e-7


def same_drawable(kind, first, second) -> bool:
    if kind == "line":
        return same_line(first, second)
    return max(abs(a - b) for a, b in zip(first, second, strict=True)) <= 1e-6


def on_drawable(point, kind, value) -> bool:
    if kind == "line":
        return on_line(point, value)
    center_x, center_y, radius_squared = value
    residual = abs(
        (point[0] - center_x) ** 2
        + (point[1] - center_y) ** 2
        - radius_squared
    )
    return residual <= TOLERANCE * max(1.0, radius_squared)


def build_state(
    *,
    fixtures=FIXTURES,
    include_finite_targets: bool = False,
    merge_branch: bool = False,
):
    replay_class = FinitePairsFirstReplay if include_finite_targets else CoreOnlyReplay
    replays = []
    for index, (centers, radii) in enumerate(fixtures):
        replay = replay_class(
            f"kp_parallel_center_locus_{index}",
            centers=centers,
            radii=radii,
            emit=False,
        )
        replay.run()
        branches = {
            profile: branch_data(replay, profile)["kind"]
            for profile in ("P0", "P1", "P2", "P3")
        }
        if branches != {
            "P0": "simple_merge" if merge_branch else "regular",
            "P1": "regular",
            "P2": "Kp_parallel",
            "P3": "Kp_parallel",
        }:
            raise AssertionError(f"搜索夹具分支错误：{branches}")
        replays.append(replay)

    paid_sequences = tuple(
        tuple(replay.objects.graph.paid_order) for replay in replays
    )
    if len({len(sequence) for sequence in paid_sequences}) != 1:
        raise AssertionError("三个双平行夹具的计费对象数不一致")

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
        point_value(
            replay.objects,
            (
                "Mannheim_S_from_completed_merge_pair"
                if merge_branch
                else "Mannheim_S_center_locus"
            ),
        )
        for replay in replays
    )
    target_lines = {}
    for profile in PROFILES:
        rows = []
        for replay in replays:
            data = replay.core_data[profile]
            targets = replay.verify_pair(
                profile,
                data["tau"],
                allow_repeated_physical_signs=True,
            )
            centers = tuple(
                collapse_point(target["center"]) for target in targets.values()
            )
            rows.append(float_line(exact_line_through(*centers)))
        target_lines[profile] = tuple(rows)
    return state, target_lines, root_centers, paid_sequences[0]


def target_points(state, target_lines, root_centers):
    result = {}
    for profile, lines in target_lines.items():
        result[profile] = tuple(
            point.point_id
            for point in state.points
            if all(
                on_line(value, line)
                for value, line in zip(point.values, lines, strict=True)
            )
        )
        non_root = tuple(
            point_id
            for point_id in result[profile]
            if not all(
                same_point(value, root)
                for value, root in zip(
                    next(
                        point.values
                        for point in state.points
                        if point.point_id == point_id
                    ),
                    root_centers,
                    strict=True,
                )
            )
        )
        if non_root:
            raise AssertionError(
                f"{profile} 在搜索前已有圆心线第二点：{non_root}"
            )
    return result


def build_witnesses(state, target_lines, root_centers):
    """预计算每条目标线与每个旧对象的非根心交点。"""

    witnesses = {profile: [] for profile in target_lines}
    for profile in target_lines:
        lines = target_lines[profile]
        for drawable in state.drawables:
            sample_rows = []
            for line, kind, value, root in zip(
                lines,
                (drawable.kind for _ in lines),
                drawable.values,
                root_centers,
                strict=True,
            ):
                roots = intersections("line", line, kind, value)
                roots = tuple(point for point in roots if not same_point(point, root))
                if not roots:
                    break
                sample_rows.append(roots)
            else:
                witnesses[profile].append(
                    (
                        drawable.drawable_id,
                        drawable.kind,
                        drawable.values,
                        tuple(sample_rows),
                    )
                )
    return witnesses


def candidate_has_witness(candidate, profile_witnesses) -> str | None:
    for drawable_id, kind, drawable_values, sample_rows in profile_witnesses:
        if candidate.kind == kind and all(
            same_drawable(kind, candidate_value, existing_value)
            for candidate_value, existing_value in zip(
                candidate.drawable.values,
                drawable_values,
                strict=True,
            )
        ):
            continue
        if all(
            any(on_drawable(point, candidate.kind, value) for point in points)
            for value, points in zip(
                candidate.drawable.values,
                sample_rows,
                strict=True,
            )
        ):
            return drawable_id
    return None


def candidates(state):
    existing = {drawable_signature(item) for item in state.drawables}
    seen = set()
    for first, second in combinations(state.points, 2):
        if any(
            old.kind == "line"
            and all(
                on_line(first_value, old_value)
                and on_line(second_value, old_value)
                for first_value, second_value, old_value in zip(
                    first.values,
                    second.values,
                    old.values,
                    strict=True,
                )
            )
            for old in state.drawables
        ):
            continue
        values = tuple(
            line_through(a, b)
            for a, b in zip(first.values, second.values, strict=True)
        )
        if any(value is None for value in values):
            continue
        drawable = DrawableBundle("move_1", "line", values)
        signature = drawable_signature(drawable)
        if (
            signature in existing
            or signature in seen
            or any(
                old.kind == "line"
                and all(
                    same_line(value, old_value)
                    for value, old_value in zip(
                        drawable.values,
                        old.values,
                        strict=True,
                    )
                )
                for old in state.drawables
            )
        ):
            continue
        seen.add(signature)
        yield Candidate("line", first.point_id, second.point_id, drawable)
    for center in state.points:
        for through in state.points:
            if center is through:
                continue
            values = tuple(
                circle_through(a, b)
                for a, b in zip(center.values, through.values, strict=True)
            )
            if any(value is None for value in values):
                continue
            drawable = DrawableBundle("move_1", "circle", values)
            signature = drawable_signature(drawable)
            if (
                signature in existing
                or signature in seen
                or any(
                    old.kind == "circle"
                    and all(
                        same_drawable("circle", value, old_value)
                        for value, old_value in zip(
                            drawable.values,
                            old.values,
                            strict=True,
                        )
                    )
                    for old in state.drawables
                )
            ):
                continue
            seen.add(signature)
            yield Candidate("circle", center.point_id, through.point_id, drawable)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-finite-targets", action="store_true")
    parser.add_argument("--merge-fixtures", action="store_true")
    args = parser.parse_args()
    start = monotonic()
    fixtures = MERGE_FIXTURES if args.merge_fixtures else FIXTURES
    state, target_lines, root_centers, paid_sequence = build_state(
        fixtures=fixtures,
        include_finite_targets=args.include_finite_targets,
        merge_branch=args.merge_fixtures,
    )
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
        "kp_parallel_center_locus_2e_search",
        {
            "samples": len(fixtures),
            "include_finite_targets": args.include_finite_targets,
            "merge_fixtures": args.merge_fixtures,
            "paid_core_objects": len(paid_sequence),
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
