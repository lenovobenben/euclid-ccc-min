"""分别演示全外切与包围型全内切圆的正规 18 E 构造。"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from manim import Circle, Create, Dot, FadeIn, FadeOut, LaggedStart, Rectangle, VGroup, Write, linear
from manim import DOWN, LEFT, RIGHT, UP
from animations.e49.e49_progress import (
    AUX_STROKE, BACKGROUND, CENTER_DOT_RADIUS, CIRCLE_BLUE,
    E49Progress, FOREGROUND, GOLD, INPUT_FRAME_CENTER,
    INPUT_FRAME_WIDTH, LINE_BLUE, MARKER_STROKE, MUTED, REFERENCE_RADIUS,
    TARGET, logical_to_scene,
)


DATA_DIRECTORY = Path(__file__).parent
CORE_FRAME_CENTER = (6.0, 13.3)
CORE_FRAME_WIDTH = 4.7
RECOVERY_FRAME_CENTER = (5.5, 5.0)
RECOVERY_FRAME_WIDTH = 13.2


class E18SingleTarget(E49Progress):
    total_moves = 18
    mode = "external"
    title = "全外切圆"
    meaning = "与三个输入圆分别外切"
    # Covers the overview, both auxiliary circles, and the external similarity
    # center. Existing lines must extend far enough when the camera pulls back.
    logical_bounds = (-50.0, 50.0, -35.0, 45.0)
    data_path = DATA_DIRECTORY / "geometry_external.json"

    def play_intro(self) -> None:
        title = self.make_text(f"18 E · {self.title}", font_size=56, weight="MEDIUM")
        title.move_to(UP * 2.3)
        meaning = self.make_text(self.meaning, font_size=28, color=TARGET)
        meaning.next_to(title, DOWN, buff=0.35)
        rules = self.make_text(
            "输入的三个圆与三个圆心免费\n"
            "每画一条直线或一个圆计 1 E，已有对象的交点免费\n"
            "不任选点，不搬运圆规开度",
            font_size=24, line_spacing=1.15,
        ).move_to(DOWN * 0.05)
        count = self.make_text("15 条直线 + 3 个圆 = 18 E", font_size=32, color=GOLD)
        count.move_to(DOWN * 1.6)
        scope = self.make_text("严格一般位置的构造演示 · 未声称最优", font_size=20, color=MUTED)
        scope.move_to(DOWN * 2.45)
        page = VGroup(title, meaning, rules, count, scope)
        self.play(Write(title), FadeIn(meaning), run_time=0.9)
        self.play(FadeIn(rules), run_time=0.6)
        self.play(FadeIn(count), FadeIn(scope), run_time=0.5)
        self.wait(2.8)
        self.play(FadeOut(page), run_time=0.5)

    def setup_construction(self) -> None:
        super().setup_construction()
        self.caption = None
        self.caption_base_width = 0.0

    def sync_viewport(self) -> None:
        before = self.viewport_key
        super().sync_viewport()
        if before != self.viewport_key and getattr(self, "caption", None) is not None:
            self.pin_caption()

    def pin_caption(self) -> None:
        scale = self.camera.frame.width / self.base_frame_width
        self.caption.set(width=self.caption_base_width * scale)
        self.caption.move_to(
            self.camera.frame.get_corner(UP + LEFT)
            + RIGHT * (self.caption.width / 2 + 0.28 * scale)
            + DOWN * (self.caption.height / 2 + 0.24 * scale)
        )

    def set_caption(self, content: str) -> None:
        if self.caption is not None:
            self.play(FadeOut(self.caption), run_time=0.2)
        label = self.make_text(content, font_size=22, color=FOREGROUND)
        backdrop = Rectangle(
            width=label.width + 0.25, height=label.height + 0.2,
            stroke_width=0, fill_color=BACKGROUND, fill_opacity=0.9,
        ).move_to(label)
        self.caption = VGroup(backdrop, label).set_z_index(18)
        self.caption_base_width = self.caption.width
        self.pin_caption()
        self.play(FadeIn(self.caption), run_time=0.2)

    def adjust_camera(self, e_move: int) -> None:
        cues = {
            2: ((6.7, 0.0), 16.2),
            5: (INPUT_FRAME_CENTER, INPUT_FRAME_WIDTH),
            10: (CORE_FRAME_CENTER, CORE_FRAME_WIDTH),
            15: ((5.5, 11.2), 18.0),
            16: (RECOVERY_FRAME_CENTER, RECOVERY_FRAME_WIDTH),
            18: self.final_frame(),
        }
        if e_move in cues:
            self.move_camera(*cues[e_move])

    def final_frame(self):
        if self.mode == "internal":
            return (self.data["targets"][0]["center"], 14.3)
        return (INPUT_FRAME_CENTER, 12.1)

    def label_point(self, point_id: str, label: str, direction=RIGHT, color=FOREGROUND) -> None:
        if point_id in self.key_point_groups:
            return
        point = self.point_at(point_id)
        dot = Dot(point, radius=self.screen_stroke_width(0.012), color=color)
        text = self.make_text(label, font_size=15, color=color)
        text.scale(self.camera.frame.width / self.base_frame_width)
        text.next_to(dot, direction, buff=self.screen_stroke_width(0.07))
        group = VGroup(dot, text).set_z_index(11)
        self.keep_screen_size(group, point)
        self.key_point_groups[point_id] = group
        self.play(FadeIn(group), run_time=0.2)

    def hide_point_labels(self) -> None:
        if self.key_point_groups:
            self.play(FadeOut(VGroup(*self.key_point_groups.values())), run_time=0.2)
            self.key_point_groups.clear()

    def reveal_key_points(self, e_move: int) -> None:
        if e_move == 1:
            for node, label in (("alpha", "α"), ("a", "a"), ("a1", "a₁"), ("alpha1", "α₁")):
                self.label_point(node, label, DOWN)
        elif e_move == 5:
            self.label_point("A", "A", LEFT)
            self.label_point("B", "B", RIGHT)
        elif e_move == 11:
            self.label_point("P0_K", "K", UP + RIGHT)
        elif e_move == 13:
            self.label_point("P0_Kp", "K′", RIGHT)
        elif e_move == 14:
            self.show_contact_choice()
        elif e_move == 15:
            self.label_point(self.data["events"][16]["references"][1], "M₂", RIGHT)
        elif e_move == 17:
            self.label_point(self.data["events"][-1]["references"][0], "C", DOWN + LEFT, TARGET)

    def show_contact_choice(self) -> None:
        chosen = self.data["events"][-1]["references"][1]
        other = "P0_---_M3" if self.mode == "external" else "P0_+++_M3"
        caption = "14 E · 选取全外切接触点" if self.mode == "external" else "14 E · 选取全内切接触点"
        self.set_caption(caption)
        markers = VGroup(self.reference_marker(chosen, GOLD), self.reference_marker(other, MUTED))
        self.play(FadeIn(markers), run_time=0.25)
        self.label_point(chosen, "M₃", RIGHT, GOLD)
        self.wait(0.9)
        self.play(FadeOut(markers), run_time=0.25)

    def make_target_marker(self, event: dict) -> VGroup:
        existing = self.key_point_groups.pop(event["references"][0], None)
        if existing is not None:
            return existing
        point = self.point_at(event["references"][0])
        dot = Dot(point, radius=self.screen_stroke_width(CENTER_DOT_RADIUS), color=TARGET)
        label = self.make_text("C", font_size=18, color=TARGET)
        label.scale(self.camera.frame.width / self.base_frame_width)
        label.next_to(dot, DOWN + LEFT, buff=self.screen_stroke_width(0.08))
        group = VGroup(dot, label).set_z_index(12)
        self.keep_screen_size(group, point)
        return group

    def play_construction(self) -> None:
        phase_captions = {
            1: "1–5 E · 构造平行线",
            6: "6–9 E · 取得第三圆上的四个点",
            10: "10–14 E · 构造接触弦",
            15: "15–17 E · 恢复目标圆心",
            18: f"18 E · 画出{self.title}",
        }
        for event in self.data["events"]:
            step = event["e_move"]
            if step in (10, 15, 16):
                self.hide_point_labels()
            if step in phase_captions:
                self.set_caption(phase_captions[step])
            self.adjust_camera(step)
            self.assert_references_visible(event)
            if step == 15:
                self.label_point("H23_ext", "H", UP + LEFT)
            if step == 16:
                self.label_point(event["references"][1], "M₃", RIGHT, GOLD)
            if step == 17:
                self.label_point(event["references"][1], "M₂", RIGHT)
            references = self.show_references(event)
            drawable = self.make_drawable(event)
            visible = self.visible_line(event, drawable) if event["op"] == "line" else drawable
            self.update_counter(step)
            self.play(Create(visible, rate_func=linear), run_time=1.6 if step == 18 else (1.1 if event["op"] == "circle" else 0.8))
            if visible is not drawable:
                self.remove(visible)
                self.add(drawable)
            if step == 18:
                already_visible = event["references"][0] in self.key_point_groups
                marker = self.make_target_marker(event)
                self.target_drawables.add(drawable)
                self.target_markers.add(marker)
                finish = [FadeOut(references)]
                if not already_visible:
                    finish.append(FadeIn(marker))
                self.play(*finish, run_time=0.25)
            else:
                self.aux_drawables.add(drawable)
                drawable.screen_stroke_width = AUX_STROKE
                self.play(
                    drawable.animate.set_stroke(
                        color=CIRCLE_BLUE if event["op"] == "circle" else LINE_BLUE,
                        width=self.screen_stroke_width(AUX_STROKE), opacity=0.28,
                    ),
                    FadeOut(references), run_time=0.25,
                )
            self.reveal_key_points(step)
            self.wait(0.2)

    def finish_construction(self) -> None:
        self.hide_point_labels()
        for drawable in self.aux_drawables:
            drawable.screen_stroke_width = 0.3
        self.play(self.aux_drawables.animate.set_stroke(width=self.screen_stroke_width(0.3), opacity=0.045), run_time=0.6)
        self.set_caption(f"{self.title} · 18 E\n15 条直线 + 3 个圆")
        contacts = VGroup()
        for item in self.data["contacts"]:
            point = logical_to_scene(item["at"])
            mark = Circle(
                radius=self.screen_stroke_width(REFERENCE_RADIUS), color=GOLD,
                stroke_width=self.screen_stroke_width(MARKER_STROKE), fill_opacity=0,
            ).move_to(point)
            contacts.add(mark)
        # Tangency points are free intersections of the finished target and inputs.
        self.play(LaggedStart(*(FadeIn(mark) for mark in contacts), lag_ratio=0.3), run_time=0.9)
        self.wait(3.0)


class E18External(E18SingleTarget):
    """与三个输入圆都外切。"""


class E18Internal(E18SingleTarget):
    """包围三个输入圆，并与每个输入圆内切。"""

    mode = "internal"
    title = "全内切圆"
    meaning = "包住三个输入圆，并与它们分别内切"
    data_path = DATA_DIRECTORY / "geometry_internal.json"
