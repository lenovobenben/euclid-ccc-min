"""两条独立 18 E 轨迹及其镜头范围的回归检查。"""

from __future__ import annotations

import importlib.util
import json
import math
import tempfile
import unittest

from animations.e18.export_e18_geometry import MODES, OUTPUT_DIRECTORY, build_export

HAS_MANIM = importlib.util.find_spec("manim") is not None
if HAS_MANIM:
    from manim import tempconfig
    from animations.e18.e18_progress import E18External, E18Internal
    from animations.e49.e49_progress import GEOMETRY_SCALE, logical_to_scene


class E18AnimationDataTests(unittest.TestCase):
    def test_snapshots_match_exact_replay(self) -> None:
        for mode in MODES:
            with self.subTest(mode=mode):
                stored = json.loads((OUTPUT_DIRECTORY / f"geometry_{mode}.json").read_text())
                self.assertEqual(stored, build_export(mode))

    def test_each_movie_has_18_independently_legal_draws(self) -> None:
        for mode in MODES:
            with self.subTest(mode=mode):
                data = build_export(mode)
                available = {"O1", "O2", "O3", "Gamma1", "Gamma2", "Gamma3"}
                remaining = dict(data["points"])
                kinds = []
                target_events = []
                for step, event in enumerate(data["events"], 1):
                    while True:
                        free = [node for node, point in remaining.items() if set(point["dependencies"]) <= available]
                        if not free:
                            break
                        for node in free:
                            available.add(node)
                            del remaining[node]
                    self.assertEqual(event["e_move"], step)
                    self.assertTrue(set(event["references"]) <= available)
                    self.assertEqual(len(event["references"]), 2)
                    self.assertNotIn(event["id"], available)
                    available.add(event["id"])
                    kinds.append(event["op"])
                    if "target" in event:
                        target_events.append(step)
                self.assertEqual((kinds.count("line"), kinds.count("circle")), (15, 3))
                self.assertEqual(target_events, [18])
                self.assertEqual(data["events"][-1]["target"]["key"], MODES[mode])

    def test_first_14_draws_are_shared_and_suffixes_are_distinct(self) -> None:
        external, internal = (build_export(mode) for mode in MODES)
        self.assertEqual(external["initial"], internal["initial"])
        self.assertEqual(external["events"][:14], internal["events"][:14])
        self.assertNotEqual(external["events"][14:], internal["events"][14:])
        external_ids = {event["id"] for event in external["events"]}
        internal_ids = {event["id"] for event in internal["events"]}
        self.assertEqual(len(external_ids & internal_ids), 14)
        self.assertNotIn("target_+++", internal_ids)
        self.assertNotIn("target_---", external_ids)

    def test_exported_tangency_signs_and_contact_marks(self) -> None:
        for mode in MODES:
            with self.subTest(mode=mode):
                data = build_export(mode)
                target = data["targets"][0]
                for initial, contact in zip(data["initial"]["circles"], data["contacts"], strict=True):
                    distance = math.dist(target["center"], initial["center"])
                    if mode == "external":
                        self.assertAlmostEqual(distance, target["radius"] + initial["radius"], places=8)
                    else:
                        self.assertGreater(target["radius"], initial["radius"])
                        self.assertAlmostEqual(distance + initial["radius"], target["radius"], places=8)
                    self.assertAlmostEqual(math.dist(contact["at"], initial["center"]), initial["radius"], places=8)
                    self.assertAlmostEqual(math.dist(contact["at"], target["center"]), target["radius"], places=8)
                    self.assertEqual(contact["available_after"], 18)


@unittest.skipUnless(HAS_MANIM, "requires the Manim rendering environment")
class E18ViewportTests(unittest.TestCase):
    def test_both_movies_keep_references_and_complete_final_circles_visible(self) -> None:
        with tempfile.TemporaryDirectory() as directory, tempconfig({
            "dry_run": True, "progress_bar": "none", "media_dir": directory,
            "frame_rate": 60,
        }):
            for scene_class in (E18External, E18Internal):
                with self.subTest(scene=scene_class.__name__):
                    scene = scene_class(skip_animations=True)
                    scene.data = json.loads(scene.data_path.read_text())
                    scene.setup_construction()
                    for event in scene.data["events"]:
                        scene.adjust_camera(event["e_move"])
                        scene.assert_references_visible(event)
                        if event["op"] == "line":
                            full = scene.make_drawable(event)
                            visible = scene.visible_line(event, full)
                            self.assertLess(visible.get_length(), full.get_length())
                    frame = scene.camera.frame
                    for circle in [*scene.data["initial"]["circles"], scene.data["targets"][0]]:
                        center = logical_to_scene(circle["center"]) - frame.get_center()
                        radius = circle["radius"] * GEOMETRY_SCALE
                        self.assertLess(abs(center[0]) + radius, frame.width / 2)
                        self.assertLess(abs(center[1]) + radius, frame.height / 2)


if __name__ == "__main__":
    unittest.main()
