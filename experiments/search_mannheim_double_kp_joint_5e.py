"""搜索双 ``K'`` 平行修复能否共享第一个对象。

在三个不同的有理双平行夹具上同步枚举。初始状态包含合法平行前缀、
八条批量线、两个有限对角点、两块各自的有限弦与方向弦，以及这些
对象直接给出的命名点。候选形状固定为：

1. 两个方向类共享第一个新对象；
2. 每个方向类各画一个独立辅助对象；
3. 各用两个已经严格分离的点画出目标平行线。

命中即把原来两套 ``2 圆 + 1 线`` 的 6 E 联合修复降为 5 E。无命中只
排除上述依赖形状，不构成 5 E 不可能性证明。搜索使用浮点数发现候选；
任何命中仍须另作精确重放。点对必须在每个样本中相距超过阈值，防止
同一精确交点因浮点误差被登记两次而产生伪命中。
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
from math import hypot, sqrt
from time import monotonic

from check_mannheim_degenerate_fixture import (
    line_intersection as exact_line_intersection,
    line_through as exact_line_through,
)
from replay_mannheim_centered_parallel_repair import build_roles
from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8
from search_parallel_3e import (
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
    circle_through,
    drawable_signature,
    generate_candidates,
    intersections,
    line_through,
)
from search_parallel_block_3e import (
    is_target_chord_circle,
    is_target_line,
    normalized_parallel_through,
)


F = Fraction
FIXTURES = (
    (
        (
            (F(0), F(0)),
            (F(225, 17), F(0)),
            (F(150, 17), F(93, 17)),
        ),
        (F(4), F(2), F(1)),
    ),
    (
        (
            (F(0), F(0)),
            (F(117, 7), F(0)),
            (F(78, 7), F(51, 7)),
        ),
        (F(4), F(2), F(1)),
    ),
    (
        ((F(0), F(0)), (F(45), F(0)), (F(30), F(21))),
        (F(4), F(2), F(1)),
    ),
)

EXPECTED_EVENTS = {"P2:parallel:Kp", "P3:parallel:Kp"}

BATCH_KEYS = (
    "alphaA",
    "aB",
    "a1A",
    "alpha1B",
    "alphaB",
    "aA",
    "a1B",
    "alpha1A",
)

POINT_KEYS = (
    "O1",
    "O2",
    "O3",
    "alpha",
    "a",
    "a1",
    "alpha1",
    "A",
    "B",
    "parallel_X",
    "parallel_Q",
    "parallel_R",
    "H23_ext",
    "H23_int",
    *BATCH_KEYS,
    "P2_K",
    "P3_K",
)


def float_point(point) -> tuple[float, float]:
    return (float(point[0]), float(point[1]))


def float_exact_line(value) -> tuple[float, float, float]:
    a, b, c = (float(coordinate) for coordinate in value)
    norm = hypot(a, b)
    a, b, c = a / norm, b / norm, c / norm
    if a < -1e-9 or (abs(a) <= 1e-9 and b < 0):
        a, b, c = -a, -b, -c
    return (a, b, c)


def batch_source_name(key: str) -> str:
    if key.startswith("alpha1"):
        return "alpha1"
    if key.startswith("alpha"):
        return "alpha"
    if key.startswith("a1"):
        return "a1"
    return "a"


def build_state() -> tuple[State, dict[str, tuple[tuple[float, float, float], ...]]]:
    point_rows = {key: [] for key in POINT_KEYS}
    drawable_rows: dict[str, list] = {}
    drawable_kinds: dict[str, str] = {}
    target_rows = {"P2": [], "P3": []}

    def add_drawable(key: str, kind: str, value) -> None:
        drawable_rows.setdefault(key, []).append(value)
        drawable_kinds[key] = kind

    for centers, radii in FIXTURES:
        if not is_d8(centers, radii):
            raise AssertionError("搜索夹具不属于严格 D8")
        if analyze_fixture(centers, radii) != EXPECTED_EVENTS:
            raise AssertionError("搜索夹具含有额外的 Mannheim 退化")
        replay = ThreeBlockReplay(centers, radii)
        replay.build_prefix()
        for key in BATCH_KEYS:
            replay.draw_batch(key)
        replay.bind_similarity_center("ext")
        replay.bind_similarity_center("int")

        o1, o2, o3 = centers
        distance_13 = sqrt(
            float((o3[0] - o1[0]) ** 2 + (o3[1] - o1[1]) ** 2)
        )
        parallel_x = (distance_13, 0.0)
        parallel_q = (float(o3[0]), -float(o3[1]))
        parallel_r = (
            2 * distance_13 - parallel_q[0],
            -parallel_q[1],
        )

        points = {
            "O1": float_point(o1),
            "O2": float_point(o2),
            "O3": float_point(o3),
            **{
                key: float_point(value)
                for key, value in replay.named_points.items()
            },
            "parallel_X": parallel_x,
            "parallel_Q": parallel_q,
            "parallel_R": parallel_r,
            "H23_ext": float_point(replay.similarity_centers["ext"]),
            "H23_int": float_point(replay.similarity_centers["int"]),
            **{
                key: float_point(value)
                for key, value in replay.batch_points.items()
            },
        }

        add_drawable(
            "Gamma1",
            "circle",
            (float(o1[0]), float(o1[1]), float(radii[0] ** 2)),
        )
        add_drawable(
            "Gamma2",
            "circle",
            (float(o2[0]), float(o2[1]), float(radii[1] ** 2)),
        )
        add_drawable(
            "Gamma",
            "circle",
            (float(o3[0]), float(o3[1]), float(radii[2] ** 2)),
        )
        add_drawable("ell", "line", line_through(points["O1"], points["O2"]))
        add_drawable(
            "parallel_c0",
            "circle",
            (float(o1[0]), float(o1[1]), distance_13**2),
        )
        add_drawable(
            "parallel_cX",
            "circle",
            circle_through(parallel_x, points["O3"]),
        )
        add_drawable(
            "parallel_diameter",
            "line",
            line_through(parallel_x, parallel_q),
        )
        add_drawable(
            "ell3",
            "line",
            line_through(points["O3"], parallel_r),
        )

        for key in BATCH_KEYS:
            source = points[batch_source_name(key)]
            endpoint = points[key[-1]]
            add_drawable(f"batch_{key}", "line", line_through(source, endpoint))

        for profile in ("P2", "P3"):
            roles = build_roles(centers, radii, profile)
            x, y, z, w = (roles[name] for name in "xyzw")
            finite_lines = (
                exact_line_through(x, w),
                exact_line_through(y, z),
            )
            direction_lines = (
                exact_line_through(x, y),
                exact_line_through(z, w),
            )
            finite_point = exact_line_intersection(*finite_lines)
            points[f"{profile}_K"] = float_point(finite_point)
            for index, value in enumerate(finite_lines, start=1):
                add_drawable(
                    f"{profile}_finite_{index}",
                    "line",
                    float_exact_line(value),
                )
            for index, value in enumerate(direction_lines, start=1):
                add_drawable(
                    f"{profile}_direction_{index}",
                    "line",
                    float_exact_line(value),
                )
            target_rows[profile].append(
                normalized_parallel_through(
                    points[f"{profile}_K"],
                    float_exact_line(direction_lines[0]),
                )
            )

        for key in POINT_KEYS:
            point_rows[key].append(points[key])

    point_bundles = tuple(
        PointBundle(key, tuple(values)) for key, values in point_rows.items()
    )
    raw_drawables = tuple(
        DrawableBundle(
            key,
            drawable_kinds[key],
            tuple(values),
        )
        for key, values in drawable_rows.items()
    )
    seen = set()
    drawables = []
    for drawable in raw_drawables:
        signature = drawable_signature(drawable)
        if signature in seen:
            continue
        seen.add(signature)
        drawables.append(drawable)
    targets = {key: tuple(values) for key, values in target_rows.items()}
    return State(point_bundles, tuple(drawables)), targets


def target_point_pair(state: State, targets) -> tuple[str, str] | None:
    points_on_target = []
    for point in state.points:
        if all(
            abs(a * x + b * y + c) <= 1e-7
            for (x, y), (a, b, c) in zip(
                point.values,
                targets,
                strict=True,
            )
        ):
            points_on_target.append(point)
    for first, second in combinations(points_on_target, 2):
        if all(
            hypot(
                first_value[0] - second_value[0],
                first_value[1] - second_value[1],
            )
            > 1e-6
            for first_value, second_value in zip(
                first.values,
                second.values,
                strict=True,
            )
        ):
            return first.point_id, second.point_id
    return None


def candidate_target_witness(
    candidate,
    state: State,
    targets,
    finite_point: PointBundle,
) -> tuple[str, str] | None:
    if is_target_line(candidate.drawable, targets):
        return (candidate.describe(), "target-line")
    if is_target_chord_circle(candidate.drawable, state, targets):
        return (candidate.describe(), "target-chord-circle")

    for existing in state.drawables:
        rows = tuple(
            intersections(
                candidate.drawable.kind,
                candidate_value,
                existing.kind,
                existing_value,
            )
            for candidate_value, existing_value in zip(
                candidate.drawable.values,
                existing.values,
                strict=True,
            )
        )
        root_counts = {len(row) for row in rows}
        if len(root_counts) != 1 or not rows[0]:
            continue
        for root_index in range(len(rows[0])):
            values = tuple(row[root_index] for row in rows)
            if not all(
                abs(a * x + b * y + c) <= 1e-7
                for (x, y), (a, b, c) in zip(
                    values,
                    targets,
                    strict=True,
                )
            ):
                continue
            if not all(
                hypot(
                    value[0] - known[0],
                    value[1] - known[1],
                )
                > 1e-6
                for value, known in zip(
                    values,
                    finite_point.values,
                    strict=True,
                )
            ):
                continue
            return (
                candidate.describe(),
                f"move_2&{existing.drawable_id}[{root_index}]",
            )
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--first-kind",
        choices=("line", "circle", "both"),
        default="circle",
    )
    parser.add_argument(
        "--second-kind",
        choices=("line", "circle", "both"),
        default="circle",
    )
    args = parser.parse_args()
    state, targets = build_state()
    first_candidates = generate_candidates(state, 1)
    if args.first_kind != "both":
        first_candidates = tuple(
            candidate
            for candidate in first_candidates
            if candidate.kind == args.first_kind
        )

    start = monotonic()
    second_candidates = 0
    hits = []
    finite_points = {
        profile: next(
            point for point in state.points if point.point_id == f"{profile}_K"
        )
        for profile in targets
    }
    for first_index, first in enumerate(first_candidates, start=1):
        state_one = apply_candidate(state, first, 1)
        found: dict[str, tuple[tuple[str, ...], tuple[str, str]] | None] = {
            profile: None for profile in targets
        }
        for profile, target in targets.items():
            point_pair = target_point_pair(state_one, target)
            if point_pair is not None:
                found[profile] = (state_one.program, point_pair)
        current_second_candidates = generate_candidates(state_one, 2)
        if args.second_kind != "both":
            current_second_candidates = tuple(
                candidate
                for candidate in current_second_candidates
                if candidate.kind == args.second_kind
            )
        for second in current_second_candidates:
            second_candidates += 1
            for profile, target in targets.items():
                if found[profile] is not None:
                    continue
                witness = candidate_target_witness(
                    second,
                    state_one,
                    target,
                    finite_points[profile],
                )
                if witness is not None:
                    found[profile] = (
                        state_one.program + (second.describe(),),
                        witness,
                    )
            if all(value is not None for value in found.values()):
                break
        if all(value is not None for value in found.values()):
            hits.append((first.describe(), found))
            break
        if first_index % 50 == 0:
            print(
                "progress",
                {
                    "first": first_index,
                    "first_total": len(first_candidates),
                    "second_candidates": second_candidates,
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    print(
        "summary",
        {
            "samples": len(FIXTURES),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "first_kind": args.first_kind,
            "second_kind": args.second_kind,
            "first_candidates": len(first_candidates),
            "second_candidates": second_candidates,
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    for shared, found in hits:
        print("candidate", {"shared": shared, "branches": found})


if __name__ == "__main__":
    main()
