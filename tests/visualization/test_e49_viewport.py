"""Run these viewport checks inside the animation's Manim image."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ANIMATION_PATH = Path(__file__).resolve().parents[2] / "animations" / "e49"
sys.path.insert(0, str(ANIMATION_PATH))
HAS_MANIM = importlib.util.find_spec("manim") is not None

if HAS_MANIM:
    from manim import tempconfig
    from e49_progress import CAMERA_CUES, E49Progress, GOLD, INPUT_FRAME_WIDTH, logical_to_scene


@unittest.skipUnless(HAS_MANIM, "requires the Manim rendering environment")
class E49ViewportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.media = tempfile.TemporaryDirectory()
        self.addCleanup(self.media.cleanup)
        self.configuration = tempconfig(
            {"dry_run": True, "progress_bar": "none", "media_dir": self.media.name}
        )
        self.configuration.__enter__()
        self.addCleanup(self.configuration.__exit__, None, None, None)
        self.scene = E49Progress(skip_animations=True)
        self.scene.data = json.loads((ANIMATION_PATH / "geometry.json").read_text())
        self.scene.setup_construction()

    def test_reference_marks_stay_small_across_all_camera_scales(self) -> None:
        widths = (INPUT_FRAME_WIDTH, 16.2, 3.2, 2.7, 15.8)
        sizes = []
        for width in widths:
            self.scene.camera.frame.set(width=width)
            marker = self.scene.reference_marker("O3", GOLD)
            ring, dot = marker
            sizes.append((ring.width / width, dot.width / width))
            # A solid dot must not obscure a dense intersection, even zoomed in.
            self.assertLess(dot.width / width * 1920, 4)
            self.assertLess(ring.width / width * 1920, 14)
        for size in sizes[1:]:
            self.assertAlmostEqual(size[0], sizes[0][0])
            self.assertAlmostEqual(size[1], sizes[0][1])

    def test_zoom_preserves_screen_styles_without_per_object_updaters(self) -> None:
        drawable = self.scene.make_drawable(self.scene.data["events"][0])
        center_group = self.scene.input_centers[2]
        dot = center_group[0]
        original_center = dot.get_center().copy()
        original_label_ratio = center_group[1].width / self.scene.camera.frame.width
        original_stroke_ratio = drawable.get_stroke_width() / self.scene.camera.frame.width
        for width in (3.2, 16.2, 2.7, INPUT_FRAME_WIDTH):
            self.scene.camera.frame.set(width=width)
            self.scene.sync_viewport()
            self.assertAlmostEqual(drawable.get_stroke_width() / width, original_stroke_ratio)
            self.assertAlmostEqual(center_group[1].width / width, original_label_ratio)
            for actual, expected in zip(dot.get_center(), original_center):
                self.assertAlmostEqual(actual, expected)
            for item in (drawable, center_group, self.scene.counter):
                self.assertFalse(any(member.get_updaters() for member in item.get_family()))

    def test_zoomed_line_animation_uses_the_visible_segment(self) -> None:
        event = self.scene.data["events"][13]
        center, width = CAMERA_CUES[14]
        self.scene.camera.frame.set(width=width).move_to(logical_to_scene(center))
        full_line = self.scene.make_drawable(event)
        visible = self.scene.visible_line(event, full_line)
        self.assertLess(visible.get_length(), full_line.get_length() / 2)
        frame = self.scene.camera.frame
        for endpoint in (visible.get_start(), visible.get_end()):
            relative = endpoint - frame.get_center()
            self.assertLessEqual(abs(relative[0]), frame.width * 0.511)
            self.assertLessEqual(abs(relative[1]), frame.height * 0.511)
        # The replacement must still be the same infinite geometric line.
        direction = full_line.get_end() - full_line.get_start()
        for endpoint in (visible.get_start(), visible.get_end()):
            delta = endpoint - full_line.get_start()
            self.assertAlmostEqual(direction[0] * delta[1] - direction[1] * delta[0], 0)


if __name__ == "__main__":
    unittest.main()
