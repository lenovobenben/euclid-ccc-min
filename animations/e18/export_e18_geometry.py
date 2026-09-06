"""从精确 Mannheim 依赖图分别导出全外切、全内切的 18 E 单目标轨迹。"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for path in (ROOT, ROOT / "experiments"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from animations.e49.export_geometry import (
    ANIMATION_CENTERS as E49_CENTERS,
    ANIMATION_RADII as E49_RADII,
)
from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8
from search_mannheim_center_locus_2e import object_value, point_value


OUTPUT_DIRECTORY = Path(__file__).parent
MODES = {"external": "+++", "internal": "---"}
# Keep the existing layout, with a larger second circle so H23 remains nearby.
ANIMATION_CENTERS = E49_CENTERS
ANIMATION_RADII = (E49_RADII[0], Fraction(4), E49_RADII[2])


def build_replay() -> ThreeBlockReplay:
    if not is_d8(ANIMATION_CENTERS, ANIMATION_RADII):
        raise AssertionError("动画输入不属于 D8")
    if analyze_fixture(ANIMATION_CENTERS, ANIMATION_RADII):
        raise AssertionError("动画输入不是严格正规夹具")
    replay = ThreeBlockReplay(ANIMATION_CENTERS, ANIMATION_RADII)
    replay.build_prefix()
    for key in ("alphaA", "aB", "a1A", "alpha1B"):
        replay.draw_batch(key)
    # This method verifies both physical tangency signs in exact quadratic fields.
    replay.build_regular_pair("P0", ("+++", "---"))
    return replay


def approximate(value) -> float:
    return value.approximate() if hasattr(value, "approximate") else float(value)


def exact_contacts(replay: ThreeBlockReplay, sign: str) -> list[dict]:
    """核对相切、内含方向和最终三个免费切点，不用浮点残差作证明。"""
    center, radius_squared = replay.targets[sign]["circle"]
    polarity = 1 if sign == "+++" else -1
    contacts = []
    for index, (origin, input_radius) in enumerate(zip(replay.centers, replay.radii, strict=True), 1):
        delta = tuple(center[i] - origin[i] for i in range(2))
        distance_squared = sum(value * value for value in delta)
        signed_term = distance_squared - radius_squared - input_radius**2
        radius = signed_term / (2 * polarity * input_radius)
        if radius.sign() <= 0 or radius * radius != radius_squared:
            raise AssertionError("目标圆的精确相切方向错误")
        if polarity < 0 and (radius - input_radius).sign() <= 0:
            raise AssertionError("全内切目标必须包住输入圆")
        factor = polarity * input_radius / (radius + polarity * input_radius)
        point = tuple(origin[i] + factor * delta[i] for i in range(2))
        to_input = tuple(point[i] - origin[i] for i in range(2))
        to_target = tuple(point[i] - center[i] for i in range(2))
        if sum(value * value for value in to_input) != input_radius**2:
            raise AssertionError("切点不在输入圆上")
        if sum(value * value for value in to_target) != radius_squared:
            raise AssertionError("切点不在目标圆上")
        contacts.append({
            "input_id": f"Gamma{index}",
            "at": [approximate(value) for value in point],
            "available_after": 18,
            "dependencies": [replay.targets[sign]["output_id"], f"Gamma{index}"],
        })
    return contacts


def build_export(mode: str) -> dict:
    sign = MODES[mode]
    replay = build_replay()
    graph = replay.objects.graph
    target = replay.targets[sign]
    paid = graph.paid_ancestors(target["output_id"])
    order = [node for node in graph.paid_order if node in paid]
    kinds = [graph.paid_kinds[node] for node in order]
    if len(order) != 18 or kinds.count("line") != 15 or kinds.count("circle") != 3:
        raise AssertionError("单目标轨迹不再是 15 线 + 3 圆 = 18 E")
    paid_index = {node: index for index, node in enumerate(order, 1)}
    available = {}

    def available_after(node: str) -> int:
        if node not in available:
            if node in graph.paid_kinds:
                if node not in paid_index:
                    raise AssertionError(f"轨迹依赖被裁剪掉的付费对象：{node}")
                available[node] = paid_index[node]
            else:
                available[node] = max((available_after(parent) for parent in graph.dependencies[node]), default=0)
        return available[node]

    # Include the two free intersections of tau and Gamma3. The unused target's
    # paid suffix is excluded; each movie starts independently from the inputs.
    point_ids = sorted({
        node for node in replay.objects.point_registry.values()
        if graph.paid_ancestors(node).issubset(paid)
    })
    points = {
        node: {
            "at": list(point_value(replay.objects, node)),
            "available_after": available_after(node),
            "dependencies": list(graph.dependencies[node]),
        }
        for node in point_ids
    }
    events = []
    for step, node in enumerate(order, 1):
        references = list(graph.dependencies[node])
        if len(references) != 2 or any(point not in points for point in references):
            raise AssertionError("每一步必须由两个已有点定位")
        if any(points[point]["available_after"] >= step for point in references):
            raise AssertionError("定位点尚未可用")
        kind, value = object_value(replay.objects, node)
        geometry = (
            dict(zip(("a", "b", "c"), value, strict=True))
            if kind == "line" else
            {"center": list(value[:2]), "radius": math.sqrt(value[2])}
        )
        event = {"e_move": step, "id": node, "op": kind, "references": references, "geometry": geometry}
        if node == target["output_id"]:
            event["target"] = {"key": sign, "profile": "P0"}
        events.append(event)
    if events[-1]["id"] != target["output_id"]:
        raise AssertionError("目标圆不是最后一步")
    _, geometry = object_value(replay.objects, target["output_id"])
    targets = [{
        "display_index": 1, "key": sign, "profile": "P0", "output_id": target["output_id"],
        "draw_index": 18, "center": list(geometry[:2]), "radius": math.sqrt(geometry[2]),
    }]
    contacts = exact_contacts(replay, sign)
    payload = {"points": points, "events": events, "targets": targets, "contacts": contacts}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
    return {
        "schema": "euclid-ccc-min-manim-e18/v1",
        "mode": mode,
        "source": {
            "replay": "experiments/replay_mannheim_three_block_dependencies.py",
            "replay_class": "ThreeBlockReplay",
            "selection": "paid ancestors of the selected P0 target",
            "trace_sha256": digest,
            "centers_exact": [[str(value) for value in center] for center in replay.centers],
            "radii_exact": [str(radius) for radius in replay.radii],
        },
        "verified_result": {
            "coverage": "gen", "physical_sign": sign, "e_move": 18,
            "line_draws": 15, "circle_draws": 3, "target_draws": 1,
        },
        "initial": {"circles": [
            {"id": f"Gamma{i}", "center_id": f"O{i}", "center": list(map(float, center)), "radius": float(radius)}
            for i, (center, radius) in enumerate(zip(replay.centers, replay.radii, strict=True), 1)
        ]},
        **payload,
    }


def main() -> None:
    for mode in MODES:
        path = OUTPUT_DIRECTORY / f"geometry_{mode}.json"
        path.write_text(json.dumps(build_export(mode), ensure_ascii=False, indent=2) + "\n")
        print(path)


if __name__ == "__main__":
    main()
