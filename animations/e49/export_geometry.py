"""从严格正规 Mannheim 重放导出 49 E 动画几何数据。"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTS_PATH = REPOSITORY_ROOT / "experiments"
if str(EXPERIMENTS_PATH) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS_PATH))

from replay_mannheim_kp_center_locus_dependencies import (  # noqa: E402
    KpCenterLocusReplay,
)
from search_mannheim_center_locus_2e import (  # noqa: E402
    object_value,
    point_value,
)


OUTPUT_PATH = Path(__file__).with_name("geometry.json")
FIXTURE_NAME = "regular"
TARGET_STEPS = (29, 31, 35, 37, 41, 43, 47, 49)


def _fraction_text(value) -> str:
    return str(value)


def _phase(e_move: int) -> str:
    if e_move <= 5:
        return "parallel_prefix"
    if e_move <= 13:
        return "batch_lines"
    if e_move <= 25:
        return "contact_chords"
    return "target_pairs"


def _trace_sha256(payload: dict) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_export() -> dict:
    replay = KpCenterLocusReplay(FIXTURE_NAME, emit=False)
    report = replay.run()
    graph = replay.objects.graph

    if report["branches"] != {
        "P0": "regular",
        "P2": "regular",
        "P1": "regular",
        "P3": "regular",
    }:
        raise RuntimeError("动画夹具不再位于 Mannheim 严格正规分支")
    if report["all_targets"] != 49 or report["trace"] != 49:
        raise RuntimeError("权威重放不再给出完整 49 E 轨迹")
    if report["union_lines"] != 39 or report["union_circles"] != 10:
        raise RuntimeError("49 E 的直线/圆台账发生变化")
    if report["non_ancestors"]:
        raise RuntimeError("动画轨迹含有不属于八目标联合祖先的对象")

    paid_index = {
        node_id: index
        for index, node_id in enumerate(graph.paid_order, start=1)
    }
    availability_cache: dict[str, int] = {}

    def available_after(node_id: str) -> int:
        if node_id in availability_cache:
            return availability_cache[node_id]
        if node_id in paid_index:
            result = paid_index[node_id]
        else:
            dependencies = graph.dependencies[node_id]
            result = max(
                (available_after(parent) for parent in dependencies),
                default=0,
            )
        availability_cache[node_id] = result
        return result

    point_ids = sorted(set(replay.objects.point_registry.values()))
    points = {
        point_id: {
            "at": list(point_value(replay.objects, point_id)),
            "available_after": available_after(point_id),
            "dependencies": list(graph.dependencies[point_id]),
        }
        for point_id in point_ids
    }

    target_by_output = {
        target["output_id"]: (target_key, target)
        for target_key, target in replay.targets.items()
    }
    events: list[dict] = []
    for e_move, node_id in enumerate(graph.paid_order, start=1):
        kind = graph.paid_kinds[node_id]
        references = list(graph.dependencies[node_id])
        if len(references) != 2 or any(
            reference not in points for reference in references
        ):
            raise RuntimeError(f"{node_id} 不是由两个已有点定位的基础操作")
        if any(points[reference]["available_after"] >= e_move for reference in references):
            raise RuntimeError(f"{node_id} 使用了尚未出现的定位点")

        recovered_kind, value = object_value(replay.objects, node_id)
        if recovered_kind != kind:
            raise RuntimeError(f"{node_id} 的几何类型与依赖图不一致")
        if kind == "line":
            geometry = {"a": value[0], "b": value[1], "c": value[2]}
        else:
            center_x, center_y, radius_squared = value
            geometry = {
                "center": [center_x, center_y],
                "radius": math.sqrt(radius_squared),
            }

        event = {
            "e_move": e_move,
            "id": node_id,
            "op": kind,
            "phase": _phase(e_move),
            "references": references,
            "geometry": geometry,
        }
        if node_id in target_by_output:
            target_key, target = target_by_output[node_id]
            event["target"] = {
                "key": target_key,
                "profile": target["profile"],
            }
        events.append(event)

    targets = []
    for display_index, (target_key, target) in enumerate(
        sorted(replay.targets.items(), key=lambda item: item[1]["draw_index"]),
        start=1,
    ):
        _, value = object_value(replay.objects, target["output_id"])
        center_x, center_y, radius_squared = value
        targets.append(
            {
                "display_index": display_index,
                "key": target_key,
                "profile": target["profile"],
                "output_id": target["output_id"],
                "draw_index": target["draw_index"],
                "center": [center_x, center_y],
                "radius": math.sqrt(radius_squared),
            }
        )
    if tuple(target["draw_index"] for target in targets) != TARGET_STEPS:
        raise RuntimeError("八个目标圆的画出时刻发生变化")

    fixture = {
        "name": FIXTURE_NAME,
        "centers_exact": [
            [_fraction_text(coordinate) for coordinate in center]
            for center in replay.centers
        ],
        "radii_exact": [_fraction_text(radius) for radius in replay.radii],
    }
    initial = {
        "circles": [
            {
                "id": f"Gamma{index}",
                "center_id": f"O{index}",
                "center": [float(center[0]), float(center[1])],
                "radius": float(radius),
            }
            for index, (center, radius) in enumerate(
                zip(replay.centers, replay.radii, strict=True),
                start=1,
            )
        ]
    }
    trace_payload = {"points": points, "events": events, "targets": targets}
    return {
        "schema": "euclid-ccc-min-manim-e49/v1",
        "source": {
            "replay": "experiments/replay_mannheim_kp_center_locus_dependencies.py",
            "replay_class": "KpCenterLocusReplay",
            "fixture": fixture,
            "trace_sha256": _trace_sha256(trace_payload),
        },
        "verified_result": {
            "profile": "CCC-ALL-8",
            "coverage": "gen",
            "e_move": report["all_targets"],
            "line_draws": report["union_lines"],
            "circle_draws": report["union_circles"],
            "target_draws": len(targets),
            "target_steps": list(TARGET_STEPS),
        },
        "initial": initial,
        "points": points,
        "events": events,
        "targets": targets,
    }


def main() -> int:
    OUTPUT_PATH.write_text(
        json.dumps(build_export(), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT_PATH)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
