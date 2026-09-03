"""筛查 Mannheim 圆心线能否在双目标后缀前用 2 E 画出。

第 8.18 节先完成一个目标圆，再由该圆心与根心 ``S`` 画出方向类的
圆心线。若能在接触弦核心完成后用“一个辅助对象 + 圆心线”提前得到
同一条线，就可只用两条第三圆半径线恢复该方向类的两个目标圆心，令
双目标后缀由 7 E 降为 6 E。

搜索在三个严格正规 ``D8`` 夹具上同步重放 25 个已画对象：合法平行
前缀、八条批量线、``P0/P2`` 接触弦核心、根心，以及 ``P1/P3`` 接触
弦核心。随后枚举所有由现有点可画的直线和圆，并检查它与任一旧对象
的交点是否在指定方向类的圆心线上。若是，第二步即可连接该交点与
``S``。

这是多夹具浮点筛查，不是 2 E 下界。命中必须另作精确重放；零命中只
排除当前夹具套件、点绑定和“第一步与一个旧对象求交”的程序形状。
"""

from __future__ import annotations

import math
from time import monotonic

from replay_mannheim_center_locus_dependencies import CenterLocusReplay
from search_mannheim_root_center_2e import (
    FIXTURES,
    add_known_drawable,
)
from search_parallel_3e import (
    DrawableBundle,
    PointBundle,
    State,
    circle_through,
    generate_candidates,
    intersections,
    line_through,
)


PROFILES = ("P0", "P1", "P2", "P3")
PRE_SUFFIX_OBJECTS = 25
TOLERANCE = 1e-7


def scalar_value(key) -> float:
    if key[0] == "rational":
        return float(key[1])
    _, discriminant, rational, radical = key
    return float(rational) + float(radical) * math.sqrt(float(discriminant))


def point_value(objects, node_id: str) -> tuple[float, float]:
    resolved = objects.resolve(node_id)
    key = next(
        key for key, value in objects.point_registry.items() if value == resolved
    )
    return tuple(scalar_value(coordinate) for coordinate in key)


def normalize_line(value) -> tuple[float, float, float]:
    norm = math.hypot(value[0], value[1])
    result = tuple(coordinate / norm for coordinate in value)
    if result[0] < -1e-9 or (
        abs(result[0]) <= 1e-9 and result[1] < 0
    ):
        result = tuple(-coordinate for coordinate in result)
    return result


def object_value(objects, node_id: str):
    resolved = objects.resolve(node_id)
    key = next(
        key for key, value in objects.object_registry.items() if value == resolved
    )
    if key[0] == "special-circle":
        return "circle", circle_through(
            point_value(objects, "parallel_X"),
            point_value(objects, "O3"),
        )
    if key[0] == "special-line":
        return "line", line_through(
            point_value(objects, "parallel_X"),
            point_value(objects, "parallel_Q"),
        )
    if key[0] == "line":
        line_key = key[1]
        if line_key[0] == "rational":
            coefficients = tuple(float(item) for item in line_key[1])
        else:
            _, discriminant, pairs = line_key
            coefficients = tuple(
                float(rational)
                + float(radical) * math.sqrt(float(discriminant))
                for rational, radical in pairs
            )
        return "line", normalize_line(coefficients)
    if key[0] != "circle":
        raise AssertionError(f"无法恢复对象 {node_id}: {key[0]}")
    center = tuple(scalar_value(coordinate) for coordinate in key[1])
    return "circle", (center[0], center[1], scalar_value(key[2]))


def on_line(point, line) -> bool:
    residual = abs(line[0] * point[0] + line[1] * point[1] + line[2])
    return residual <= TOLERANCE * max(1.0, math.hypot(*point))


def same_point(first, second) -> bool:
    return math.hypot(first[0] - second[0], first[1] - second[1]) <= TOLERANCE


def make_state():
    replays = []
    for index, (centers, radii) in enumerate(FIXTURES):
        replay = CenterLocusReplay(
            f"center_locus_search_{index}",
            centers=centers,
            radii=radii,
            emit=False,
        )
        report = replay.run()
        if set(report["branches"].values()) != {"regular"}:
            raise AssertionError("圆心线搜索夹具必须全部严格正规")
        replays.append(replay)

    paid_sequences = tuple(
        tuple(replay.objects.graph.paid_order[:PRE_SUFFIX_OBJECTS])
        for replay in replays
    )
    if len(set(paid_sequences)) != 1:
        raise AssertionError("三个夹具的接触弦核心对象顺序不一致")
    paid_sequence = paid_sequences[0]

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
    initial_drawables = []
    for node_id in ("Gamma1", "Gamma2", "Gamma3"):
        values = []
        for replay in replays:
            kind, value = object_value(replay.objects, node_id)
            if kind != "circle":
                raise AssertionError("输入圆对象种类错误")
            values.append(value)
        initial_drawables.append(
            DrawableBundle(node_id, "circle", tuple(values))
        )
    state = State(points, tuple(initial_drawables))

    for move, node_id in enumerate(paid_sequence, start=1):
        rows = tuple(object_value(replay.objects, node_id) for replay in replays)
        kinds = {kind for kind, _ in rows}
        if len(kinds) != 1:
            raise AssertionError(f"{node_id} 的对象种类不一致")
        kind = rows[0][0]
        state = add_known_drawable(
            state,
            node_id,
            kind,
            tuple(value for _, value in rows),
            move,
        )

    target_lines = {
        profile: tuple(
            object_value(replay.objects, f"{profile}_center_locus")[1]
            for replay in replays
        )
        for profile in PROFILES
    }
    root_centers = tuple(
        point_value(replay.objects, "Mannheim_S_center_locus")
        for replay in replays
    )
    return state, target_lines, root_centers


def point_on_target(point_values, target_lines) -> bool:
    return all(
        on_line(point, line)
        for point, line in zip(point_values, target_lines, strict=True)
    )


def point_is_root_center(point_values, root_centers) -> bool:
    return all(
        same_point(point, center)
        for point, center in zip(point_values, root_centers, strict=True)
    )


def main() -> None:
    start = monotonic()
    state, target_lines, root_centers = make_state()
    initial_target_points = {
        profile: tuple(
            point.point_id
            for point in state.points
            if point_on_target(point.values, target_lines[profile])
        )
        for profile in PROFILES
    }
    if any(len(points) != 1 for points in initial_target_points.values()):
        raise AssertionError(
            f"每条圆心线应只有根心一个已有点: {initial_target_points}"
        )

    candidates = generate_candidates(state, PRE_SUFFIX_OBJECTS + 1)
    hits = []
    for candidate_index, candidate in enumerate(candidates, start=1):
        for existing in state.drawables:
            rows = tuple(
                intersections(
                    candidate.kind,
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
                point_values = tuple(row[root_index] for row in rows)
                for profile in PROFILES:
                    if point_is_root_center(point_values, root_centers):
                        continue
                    if point_on_target(point_values, target_lines[profile]):
                        hits.append(
                            (
                                profile,
                                candidate.describe(),
                                existing.drawable_id,
                                root_index,
                            )
                        )
        if candidate_index % 20000 == 0:
            print(
                "progress",
                {
                    "checked": candidate_index,
                    "total": len(candidates),
                    "hits": len(hits),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    print(
        "center_locus_2e_search",
        {
            "samples": len(FIXTURES),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "first_candidates": len(candidates),
            "hits": len(hits),
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )
    print("initial_target_points", initial_target_points)
    for hit in hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
