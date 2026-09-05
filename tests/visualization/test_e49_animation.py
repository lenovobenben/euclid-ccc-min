from __future__ import annotations

import configparser
import json
import math
import sys
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
ANIMATION_PATH = REPOSITORY_ROOT / "animations" / "e49"
GEOMETRY_PATH = ANIMATION_PATH / "geometry.json"
MANIM_CONFIG_PATH = ANIMATION_PATH / "manim.cfg"
sys.path.insert(0, str(ANIMATION_PATH))

from export_geometry import build_export  # noqa: E402


class E49AnimationDataTests(unittest.TestCase):
    def test_release_render_profile_is_4k_30fps(self) -> None:
        parser = configparser.ConfigParser()
        parser.read(MANIM_CONFIG_PATH, encoding="utf-8")
        cli = parser["CLI"]

        self.assertEqual(cli.getint("pixel_width"), 3840)
        self.assertEqual(cli.getint("pixel_height"), 2160)
        self.assertEqual(cli.getint("frame_rate"), 30)

    def test_snapshot_is_the_current_exact_replay_export(self) -> None:
        stored = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, build_export())

    def test_trace_has_49_legal_draws_and_eight_targets(self) -> None:
        exported = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))
        result = exported["verified_result"]
        events = exported["events"]
        points = exported["points"]

        self.assertEqual(exported["schema"], "euclid-ccc-min-manim-e49/v1")
        self.assertEqual(result["e_move"], 49)
        self.assertEqual(result["line_draws"], 39)
        self.assertEqual(result["circle_draws"], 10)
        self.assertEqual(result["target_draws"], 8)
        self.assertEqual(result["target_steps"], [29, 31, 35, 37, 41, 43, 47, 49])
        self.assertEqual([event["e_move"] for event in events], list(range(1, 50)))

        for event in events:
            self.assertEqual(len(event["references"]), 2, event["id"])
            for reference in event["references"]:
                self.assertIn(reference, points, event["id"])
                self.assertLess(
                    points[reference]["available_after"],
                    event["e_move"],
                    event["id"],
                )

    def test_all_exported_targets_are_distinct_and_tangent_to_three_inputs(self) -> None:
        exported = json.loads(GEOMETRY_PATH.read_text(encoding="utf-8"))
        targets = exported["targets"]
        inputs = exported["initial"]["circles"]

        identities = {
            (
                round(target["center"][0], 12),
                round(target["center"][1], 12),
                round(target["radius"], 12),
            )
            for target in targets
        }
        self.assertEqual(len(identities), 8)

        for target in targets:
            x, y = target["center"]
            radius = target["radius"]
            self.assertGreater(radius, 0)
            for given in inputs:
                gx, gy = given["center"]
                distance = math.hypot(x - gx, y - gy)
                residual = min(
                    abs(distance - (radius + given["radius"])),
                    abs(distance - abs(radius - given["radius"])),
                )
                self.assertLess(residual, 1e-10, (target["key"], given["id"]))


if __name__ == "__main__":
    unittest.main()
