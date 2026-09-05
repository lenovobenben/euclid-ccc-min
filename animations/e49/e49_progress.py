"""Euclid-CCC-Min 的 Mannheim 正规 49 E 八解演示动画。"""

from __future__ import annotations

import json
import math
from pathlib import Path

from manim import (
    Arrow,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    GrowFromCenter,
    LaggedStart,
    Line,
    MovingCameraScene,
    Rectangle,
    Text,
    Transform,
    VGroup,
    Write,
    config,
    DOWN,
    LEFT,
    RIGHT,
    UP,
)


DATA_PATH = Path(__file__).with_name("geometry.json")
FONT = "Noto Sans CJK SC"

BACKGROUND = "#07101F"
FOREGROUND = "#E8F0F7"
MUTED = "#8397A8"
CIRCLE_BLUE = "#56B4D3"
LINE_BLUE = "#8AC6D9"
GOLD = "#FFC857"
TARGET = "#5DE2A5"
ALERT = "#FF7A90"
INPUT_CIRCLE = "#F45BFF"

LOGICAL_BOUNDS = (-42.0, 38.0, -23.0, 25.0)
GEOMETRY_SCALE = 0.32
GEOMETRY_ORIGIN = (5.5, 4.75)
INPUT_FRAME_CENTER = (5.0, 4.75)
INPUT_FRAME_WIDTH = 11.7
TARGET_FRAME_CENTER = (5.45, 4.65)
TARGET_FRAME_WIDTH = 15.8

SUBSCRIPTS = "₀₁₂₃₄₅₆₇₈₉"
TARGET_DIRECTIONS = (
    RIGHT,
    LEFT,
    RIGHT,
    LEFT,
    UP + LEFT,
    DOWN + LEFT,
    UP + RIGHT,
    DOWN + RIGHT,
)
KEY_POINTS = {
    "alpha": ("α", DOWN + LEFT),
    "a": ("a", DOWN + RIGHT),
    "a1": ("a₁", DOWN + LEFT),
    "alpha1": ("α₁", DOWN + RIGHT),
    "A": ("A", LEFT),
    "B": ("B", RIGHT),
    "Mannheim_S_center_locus": ("S", DOWN + RIGHT),
}

CAMERA_CUES = {
    2: ((6.7, 0.0), 16.2),
    5: (INPUT_FRAME_CENTER, INPUT_FRAME_WIDTH),
    14: ((6.0, 12.2), 3.2),
    18: ((6.45, 15.55), 4.3),
    19: ((5.8, 11.8), 3.2),
    22: ((5.3, 9.8), 2.7),
    23: ((6.75, 11.1), 2.7),
    24: ((6.0, 8.3), 4.6),
    26: ((6.2, 12.3), 3.0),
    27: ((2.2, 2.55), 5.5),
    28: (TARGET_FRAME_CENTER, TARGET_FRAME_WIDTH),
    32: ((5.5, 11.5), 3.0),
    33: ((-1.15, 2.55), 5.6),
    34: (TARGET_FRAME_CENTER, TARGET_FRAME_WIDTH),
    38: ((5.85, 11.1), 3.0),
    39: ((5.1, 2.55), 3.6),
    40: (TARGET_FRAME_CENTER, TARGET_FRAME_WIDTH),
    44: ((6.6, 11.75), 3.0),
    45: ((20.95, 2.55), 10.4),
    46: (TARGET_FRAME_CENTER, TARGET_FRAME_WIDTH),
}


def logical_to_scene(point: list[float] | tuple[float, float]):
    return (
        RIGHT * ((point[0] - GEOMETRY_ORIGIN[0]) * GEOMETRY_SCALE)
        + UP * ((point[1] - GEOMETRY_ORIGIN[1]) * GEOMETRY_SCALE)
    )


def clipped_line(geometry: dict) -> tuple[list[float], list[float]]:
    a, b, c = geometry["a"], geometry["b"], geometry["c"]
    x_min, x_max, y_min, y_max = LOGICAL_BOUNDS
    candidates: list[tuple[float, float]] = []
    if abs(b) > 1e-10:
        for x in (x_min, x_max):
            y = (-a * x - c) / b
            if y_min - 1e-8 <= y <= y_max + 1e-8:
                candidates.append((x, y))
    if abs(a) > 1e-10:
        for y in (y_min, y_max):
            x = (-b * y - c) / a
            if x_min - 1e-8 <= x <= x_max + 1e-8:
                candidates.append((x, y))
    unique: list[tuple[float, float]] = []
    for point in candidates:
        if not any(math.dist(point, other) < 1e-8 for other in unique):
            unique.append(point)
    if len(unique) < 2:
        raise ValueError(f"line does not cross viewport: {geometry}")
    first, second = max(
        (
            (first, second)
            for first in unique
            for second in unique
            if first != second
        ),
        key=lambda pair: math.dist(*pair),
    )
    return list(first), list(second)


class E49Progress(MovingCameraScene):
    def construct(self) -> None:
        config.background_color = BACKGROUND
        self.data = json.loads(DATA_PATH.read_text(encoding="utf-8"))
        self.aux_drawables = VGroup()
        self.target_drawables = VGroup()
        self.target_markers = VGroup()
        self.key_point_groups: dict[str, VGroup] = {}

        self.play_intro()
        self.setup_construction()
        self.play_construction()
        self.finish_construction()

    def make_text(self, content: str, **kwargs) -> Text:
        return Text(
            content,
            font=FONT,
            color=kwargs.pop("color", FOREGROUND),
            **kwargs,
        )

    def play_intro(self) -> None:
        title = self.make_text(
            "三个圆，八个切圆",
            font_size=58,
            weight="MEDIUM",
        ).move_to(UP * 2.45)
        subtitle = self.make_text(
            "阿波罗尼乌斯 CCC 问题",
            font_size=31,
            color=GOLD,
            weight="MEDIUM",
        ).next_to(title, DOWN, buff=0.3)
        description = self.make_text(
            "给定三个两两外离、大小不一、圆心不共线的圆\n"
            "画出全部八个与三圆分别内切或外切的目标圆",
            font_size=26,
            line_spacing=1.18,
            color=FOREGROUND,
        ).move_to(DOWN * 0.15)
        scope = self.make_text(
            "演示 Mannheim 构造的正规一般位置",
            font_size=22,
            color=MUTED,
        ).move_to(DOWN * 1.85)

        page_one = VGroup(title, subtitle, description, scope)
        self.play(Write(title), FadeIn(subtitle, shift=UP * 0.08), run_time=0.9)
        self.play(FadeIn(description, shift=UP * 0.08), run_time=0.7)
        self.play(FadeIn(scope), run_time=0.4)
        self.wait(1.5)
        self.play(FadeOut(page_one), run_time=0.5)

        metrics_title = self.make_text(
            "怎样计算一步",
            font_size=52,
            weight="MEDIUM",
        ).move_to(UP * 2.55)
        rules = self.make_text(
            "输入的三个圆与三个圆心：0 E\n"
            "每画一条直线或一个圆：1 E\n"
            "已有对象的确定有限实交点：0 E\n"
            "不允许任取点，也不把高级作图当作一步",
            font_size=25,
            line_spacing=1.05,
        ).move_to(UP * 0.55)
        domain = self.make_text(
            "49 E 覆盖 D₈ 中开稠密、满测度的 Mannheim 一般位置",
            font_size=21,
            color=MUTED,
        ).move_to(DOWN * 1.22)
        old_score = self.make_text("合法朴素展开 65 E", font_size=36, color=MUTED)
        new_score = self.make_text("当前上界 49 E", font_size=46, color=GOLD, weight="MEDIUM")
        scores = VGroup(old_score, new_score).arrange(RIGHT, buff=1.7).move_to(DOWN * 2.25)
        arrow = Arrow(
            old_score.get_right() + RIGHT * 0.12,
            new_score.get_left() + LEFT * 0.12,
            color=CIRCLE_BLUE,
            stroke_width=3.0,
            buff=0,
        )
        saving = self.make_text("减少 16 E", font_size=24, color=TARGET, weight="MEDIUM")
        saving.next_to(scores, DOWN, buff=0.34)
        page_two = VGroup(metrics_title, rules, domain, old_score, arrow, new_score, saving)

        self.play(Write(metrics_title), FadeIn(rules, shift=UP * 0.08), run_time=0.9)
        self.play(FadeIn(domain), run_time=0.35)
        self.play(FadeIn(old_score), Create(arrow), GrowFromCenter(new_score), run_time=0.75)
        self.play(FadeIn(saving, shift=UP * 0.06), run_time=0.35)
        self.wait(1.8)
        self.play(FadeOut(page_two), run_time=0.55)

    def setup_construction(self) -> None:
        self.camera.frame.move_to(logical_to_scene(INPUT_FRAME_CENTER)).set(
            width=INPUT_FRAME_WIDTH
        )
        self.base_frame_width = INPUT_FRAME_WIDTH

        self.counter = self.make_text(
            "00 / 49",
            font_size=42,
            color=GOLD,
            weight="MEDIUM",
        )
        self.counter_base_width = self.counter.width
        self.counter.set_stroke(BACKGROUND, width=7, background=True).set_z_index(20)
        self.counter_backdrop = Rectangle(
            width=self.counter.width + 0.42,
            height=self.counter.height + 0.24,
            stroke_width=0,
            fill_color=BACKGROUND,
            fill_opacity=0.88,
        ).set_z_index(19)
        backdrop_base_width = self.counter_backdrop.width
        backdrop_base_height = self.counter_backdrop.height

        def pin_counter(counter: Text) -> None:
            scale = self.camera.frame.width / self.base_frame_width
            counter.set(width=self.counter_base_width * scale)
            counter.move_to(
                self.camera.frame.get_corner(UP + RIGHT)
                + LEFT * (counter.width / 2 + 0.30 * scale)
                + DOWN * (counter.height / 2 + 0.22 * scale)
            )

        self.counter.add_updater(pin_counter)

        def pin_counter_backdrop(backdrop: Rectangle) -> None:
            scale = self.camera.frame.width / self.base_frame_width
            backdrop.set(
                width=backdrop_base_width * scale,
                height=backdrop_base_height * scale,
            )
            backdrop.move_to(self.counter)

        self.counter_backdrop.add_updater(pin_counter_backdrop)
        input_circles = VGroup()
        input_centers = VGroup()
        for index, record in enumerate(self.data["initial"]["circles"], start=1):
            circle = Circle(
                radius=record["radius"] * GEOMETRY_SCALE,
                color=INPUT_CIRCLE,
                stroke_width=self.screen_stroke_width(2.4),
                fill_opacity=0.025,
                fill_color=INPUT_CIRCLE,
            ).move_to(logical_to_scene(record["center"]))
            self.keep_screen_stroke(circle, 2.4)
            circle.set_z_index(3)
            center = Dot(logical_to_scene(record["center"]), radius=0.052, color=FOREGROUND)
            direction = DOWN if index < 3 else UP
            label = self.make_text(f"O{SUBSCRIPTS[index]}", font_size=17)
            label.next_to(center, direction, buff=0.08)
            gamma = self.make_text(f"Γ{SUBSCRIPTS[index]}", font_size=17, color=INPUT_CIRCLE)
            gamma.next_to(circle, UP if index < 3 else LEFT, buff=0.08)
            input_circles.add(circle, gamma)
            input_centers.add(center, label)
        self.input_circles = input_circles
        self.input_centers = input_centers
        self.play(
            FadeIn(VGroup(self.counter_backdrop, self.counter)),
            LaggedStart(*(Create(item) for item in input_circles[::2]), lag_ratio=0.15),
            FadeIn(VGroup(*input_circles[1::2], *input_centers)),
            run_time=1.0,
        )
        self.wait(0.3)

    def point_at(self, point_id: str):
        return logical_to_scene(self.data["points"][point_id]["at"])

    def screen_stroke_width(self, width: float) -> float:
        """把屏幕线宽换算成当前移动镜头中的场景线宽。"""

        return width * self.camera.frame.width / self.base_frame_width

    def keep_screen_stroke(self, drawable, width: float) -> None:
        """缩放镜头时保持线在屏幕上的粗细稳定。"""

        drawable.screen_stroke_width = width
        drawable.add_updater(
            lambda item: item.set_stroke(
                width=self.screen_stroke_width(item.screen_stroke_width)
            )
        )

    def make_drawable(self, event: dict):
        is_target = "target" in event
        color = TARGET if is_target else GOLD
        width = 3.0 if is_target else 1.8
        if event["op"] == "line":
            start, end = clipped_line(event["geometry"])
            drawable = Line(
                logical_to_scene(start),
                logical_to_scene(end),
                color=color,
                stroke_width=self.screen_stroke_width(width),
            )
        else:
            geometry = event["geometry"]
            drawable = Circle(
                radius=geometry["radius"] * GEOMETRY_SCALE,
                color=color,
                stroke_width=self.screen_stroke_width(width),
                fill_opacity=0,
            ).move_to(logical_to_scene(geometry["center"]))
        self.keep_screen_stroke(drawable, width)
        return drawable

    def reference_marker(self, point_id: str, color: str) -> VGroup:
        position = self.point_at(point_id)
        return VGroup(
            Circle(
                radius=0.105,
                color=color,
                stroke_width=self.screen_stroke_width(1.8),
            ).move_to(position),
            Dot(position, radius=0.044, color=color),
        ).set_z_index(15)

    def show_references(self, event: dict) -> VGroup:
        first_id, second_id = event["references"]
        if event["op"] == "line":
            overlay = VGroup(
                self.reference_marker(first_id, GOLD),
                self.reference_marker(second_id, GOLD),
            )
        else:
            center = self.point_at(first_id)
            through = self.point_at(second_id)
            radius = Line(
                center,
                through,
                color=ALERT,
                stroke_width=self.screen_stroke_width(1.4),
            ).set_z_index(14)
            overlay = VGroup(
                radius,
                self.reference_marker(first_id, ALERT),
                self.reference_marker(second_id, GOLD),
            )
        self.play(FadeIn(overlay, scale=0.75), run_time=0.13)
        return overlay

    def update_counter(self, e_move: int):
        new_counter = self.make_text(
            f"{e_move:02d} / 49",
            font_size=42,
            color=GOLD,
            weight="MEDIUM",
        )
        new_counter.scale(self.camera.frame.width / self.base_frame_width)
        new_counter.move_to(self.counter)
        new_counter.set_stroke(BACKGROUND, width=7, background=True).set_z_index(20)
        return Transform(self.counter, new_counter)

    def adjust_camera(self, e_move: int) -> None:
        cue = CAMERA_CUES.get(e_move)
        if cue is None:
            return
        center, width = cue
        self.play(
            self.camera.frame.animate.move_to(logical_to_scene(center)).set(width=width),
            run_time=0.42,
        )

    def assert_references_visible(self, event: dict) -> None:
        """渲染时验证本步的两个尺规定位点没有被局部镜头裁掉。"""

        frame = self.camera.frame
        center = frame.get_center()
        half_width = frame.width * 0.46
        half_height = frame.height * 0.43
        for point_id in event["references"]:
            point = self.point_at(point_id)
            if (
                abs(point[0] - center[0]) > half_width
                or abs(point[1] - center[1]) > half_height
            ):
                raise ValueError(
                    f"E{event['e_move']:02d} 的定位点 {point_id} 位于镜头安全区外"
                )

    def reveal_key_points(self, e_move: int) -> None:
        for point_id, (label_text, direction) in KEY_POINTS.items():
            record = self.data["points"][point_id]
            if record["available_after"] != e_move or point_id in self.key_point_groups:
                continue
            dot = Dot(self.point_at(point_id), radius=0.035, color=FOREGROUND)
            label = self.make_text(label_text, font_size=14, color=FOREGROUND)
            label.next_to(dot, direction, buff=0.055)
            group = VGroup(dot, label).set_z_index(7)
            self.key_point_groups[point_id] = group
            self.play(FadeIn(group, scale=0.7), run_time=0.12)

    def make_target_marker(self, event: dict) -> VGroup:
        display_index = next(
            item["display_index"]
            for item in self.data["targets"]
            if item["output_id"] == event["id"]
        )
        center_id = event["references"][0]
        dot = Dot(self.point_at(center_id), radius=0.052, color=TARGET)
        label = self.make_text(
            f"K{SUBSCRIPTS[display_index]}",
            font_size=15,
            color=TARGET,
            weight="MEDIUM",
        )
        label.next_to(dot, TARGET_DIRECTIONS[display_index - 1], buff=0.07)
        return VGroup(dot, label).set_z_index(12)

    def play_construction(self) -> None:
        for event in self.data["events"]:
            e_move = event["e_move"]
            if e_move == 14 and self.key_point_groups:
                endpoint_groups = VGroup(
                    *(
                        group
                        for point_id, group in self.key_point_groups.items()
                        if point_id != "Mannheim_S_center_locus"
                    )
                )
                if len(endpoint_groups):
                    self.play(FadeOut(endpoint_groups), run_time=0.2)

            self.adjust_camera(e_move)
            self.assert_references_visible(event)
            references = self.show_references(event)
            drawable = self.make_drawable(event)
            is_target = "target" in event
            animations = [Create(drawable), self.update_counter(e_move)]
            target_marker = None
            if is_target:
                target_marker = self.make_target_marker(event)
                animations.append(FadeIn(target_marker, scale=0.7))
            self.play(
                *animations,
                run_time=0.56 if is_target else 0.27,
            )
            if is_target:
                self.target_drawables.add(drawable)
                self.target_markers.add(target_marker)
                self.play(FadeOut(references), run_time=0.1)
                self.wait(0.12)
            else:
                self.aux_drawables.add(drawable)
                faded_color = CIRCLE_BLUE if event["op"] == "circle" else LINE_BLUE
                drawable.screen_stroke_width = 0.8
                self.play(
                    drawable.animate.set_stroke(
                        color=faded_color,
                        opacity=0.22,
                    ),
                    FadeOut(references),
                    run_time=0.09,
                )
            self.reveal_key_points(e_move)

    def finish_construction(self) -> None:
        self.play(
            self.camera.frame.animate.move_to(logical_to_scene(TARGET_FRAME_CENTER)).set(
                width=TARGET_FRAME_WIDTH
            ),
            run_time=0.65,
        )
        remaining_points = VGroup()
        root_center_id = "Mannheim_S_center_locus"
        if root_center_id in self.key_point_groups:
            remaining_points.add(self.key_point_groups[root_center_id])
        for drawable in self.aux_drawables:
            drawable.screen_stroke_width = 0.55
        for drawable in self.target_drawables:
            drawable.screen_stroke_width = 2.8
        finish_animations = [
            self.aux_drawables.animate.set_stroke(opacity=0.045),
            self.target_drawables.animate.set_stroke(color=TARGET, opacity=0.95),
        ]
        if len(remaining_points):
            finish_animations.append(FadeOut(remaining_points))
        self.play(*finish_animations, run_time=0.65)

        pulses = []
        for target in self.data["targets"]:
            pulse = Circle(
                radius=0.12,
                color=TARGET,
                stroke_width=self.screen_stroke_width(2.0),
                fill_opacity=0,
            ).move_to(logical_to_scene(target["center"]))
            pulses.append(pulse.animate.scale(3.0).set_stroke(opacity=0))
        self.play(LaggedStart(*pulses, lag_ratio=0.09), run_time=1.25)

        profile = self.make_text(
            "CCC-ALL-8 · 一般位置",
            font_size=24,
            color=FOREGROUND,
            weight="MEDIUM",
        )
        total = self.make_text(
            "39 条直线 + 10 个圆 = 49 E",
            font_size=31,
            color=TARGET,
            weight="MEDIUM",
        )
        conclusion = VGroup(profile, total).arrange(DOWN, buff=0.16)
        scale = self.camera.frame.width / self.base_frame_width
        conclusion.scale(scale)
        conclusion.move_to(
            self.camera.frame.get_corner(DOWN + LEFT)
            + RIGHT * (conclusion.width / 2 + 0.35 * scale)
            + UP * (conclusion.height / 2 + 0.28 * scale)
        )
        conclusion.set_stroke(BACKGROUND, width=8, background=True).set_z_index(20)
        self.play(FadeIn(conclusion, shift=UP * 0.08), run_time=0.55)
        self.wait(3.0)
