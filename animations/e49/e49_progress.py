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
    UpdateFromAlphaFunc,
    VGroup,
    Write,
    config,
    linear,
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

# Widths and marker radii are specified at the overview camera scale. They stay
# the same size on screen when the camera moves into the dense construction.
INPUT_STROKE = 0.95
ACTIVE_STROKE = 0.85
TARGET_STROKE = 1.05
AUX_STROKE = 0.40
MARKER_STROKE = 0.60
REFERENCE_RADIUS = 0.038
REFERENCE_DOT_RADIUS = 0.010
CENTER_DOT_RADIUS = 0.016

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


def clipped_line(
    geometry: dict,
    bounds: tuple[float, float, float, float] = LOGICAL_BOUNDS,
) -> tuple[list[float], list[float]]:
    a, b, c = geometry["a"], geometry["b"], geometry["c"]
    x_min, x_max, y_min, y_max = bounds
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
        self.screen_strokes = []
        self.screen_anchors = []
        self.viewport_key = None

        self.counter = self.make_text(
            "00 / 49",
            font_size=42,
            color=GOLD,
            weight="MEDIUM",
        )
        self.counter_base_width = self.counter.width
        self.counter.set_z_index(20)
        self.counter_backdrop = Rectangle(
            width=self.counter.width + 0.42,
            height=self.counter.height + 0.24,
            stroke_width=0,
            fill_color=BACKGROUND,
            fill_opacity=0.88,
        ).set_z_index(19)
        self.backdrop_base_width = self.counter_backdrop.width
        self.backdrop_base_height = self.counter_backdrop.height
        self.sync_viewport()
        input_circles = VGroup()
        input_centers = VGroup()
        for index, record in enumerate(self.data["initial"]["circles"], start=1):
            circle = Circle(
                radius=record["radius"] * GEOMETRY_SCALE,
                color=INPUT_CIRCLE,
                stroke_width=self.screen_stroke_width(INPUT_STROKE),
                fill_opacity=0,
                fill_color=INPUT_CIRCLE,
            ).move_to(logical_to_scene(record["center"]))
            self.keep_screen_stroke(circle, INPUT_STROKE)
            circle.set_z_index(3)
            position = logical_to_scene(record["center"])
            center = Dot(position, radius=CENTER_DOT_RADIUS, color=FOREGROUND)
            direction = DOWN if index < 3 else UP
            label = self.make_text(f"O{SUBSCRIPTS[index]}", font_size=17)
            label.next_to(center, direction, buff=0.08)
            gamma = self.make_text(f"Γ{SUBSCRIPTS[index]}", font_size=17, color=INPUT_CIRCLE)
            gamma.next_to(circle, UP if index < 3 else LEFT, buff=0.08)
            self.keep_screen_size(gamma, gamma.get_center())
            center_group = VGroup(center, label)
            self.keep_screen_size(center_group, position)
            input_circles.add(circle, gamma)
            input_centers.add(center_group)
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
        """记录屏幕线宽；仅在镜头变化时更新，静止背景可以缓存。"""

        drawable.screen_stroke_width = width
        self.screen_strokes.append(drawable)

    def keep_screen_size(self, marker, anchor) -> None:
        self.screen_anchors.append(
            [marker, anchor.copy(), self.camera.frame.width / self.base_frame_width]
        )

    def update_mobjects(self, dt: float) -> None:
        super().update_mobjects(dt)
        if hasattr(self, "base_frame_width"):
            self.sync_viewport()

    def sync_viewport(self) -> None:
        """镜头运动时同步线宽、点标记和计数器，避免逐对象常驻 updater。"""

        frame = self.camera.frame
        key = (frame.width, *frame.get_center())
        if key == self.viewport_key:
            return
        self.viewport_key = key
        scale = frame.width / self.base_frame_width
        for drawable in self.screen_strokes:
            drawable.set_stroke(width=drawable.screen_stroke_width * scale)
        for record in self.screen_anchors:
            marker, anchor, previous_scale = record
            marker.scale(scale / previous_scale, about_point=anchor)
            record[2] = scale
        self.pin_counter()

    def pin_counter(self) -> None:
        scale = self.camera.frame.width / self.base_frame_width
        self.counter.set(width=self.counter_base_width * scale)
        self.counter.move_to(
            self.camera.frame.get_corner(UP + RIGHT)
            + LEFT * (self.counter.width / 2 + 0.30 * scale)
            + DOWN * (self.counter.height / 2 + 0.22 * scale)
        )
        self.counter_backdrop.set(
            width=self.backdrop_base_width * scale,
            height=self.backdrop_base_height * scale,
        ).move_to(self.counter)

    def make_drawable(self, event: dict):
        is_target = "target" in event
        color = TARGET if is_target else GOLD
        width = TARGET_STROKE if is_target else ACTIVE_STROKE
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
                radius=self.screen_stroke_width(REFERENCE_RADIUS),
                color=color,
                stroke_width=self.screen_stroke_width(MARKER_STROKE),
                fill_opacity=0,
            ).move_to(position),
            Dot(position, radius=self.screen_stroke_width(REFERENCE_DOT_RADIUS), color=color),
        ).set_z_index(15)

    def visible_line(self, event: dict, drawable: Line) -> Line:
        """把绘制时间用在可见线段上，避免局部镜头只看到几帧扫过。"""

        frame = self.camera.frame
        center = frame.get_center()
        x = center[0] / GEOMETRY_SCALE + GEOMETRY_ORIGIN[0]
        y = center[1] / GEOMETRY_SCALE + GEOMETRY_ORIGIN[1]
        # Extend just beyond the viewport so the replacement is invisible.
        half_width = frame.width * 0.51 / GEOMETRY_SCALE
        half_height = frame.height * 0.51 / GEOMETRY_SCALE
        start, end = clipped_line(
            event["geometry"],
            (x - half_width, x + half_width, y - half_height, y + half_height),
        )
        return Line(
            logical_to_scene(start),
            logical_to_scene(end),
            color=drawable.get_color(),
            stroke_width=drawable.get_stroke_width(),
        )

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
                stroke_width=self.screen_stroke_width(MARKER_STROKE),
            ).set_z_index(14)
            overlay = VGroup(
                radius,
                self.reference_marker(first_id, ALERT),
                self.reference_marker(second_id, GOLD),
            )
        self.play(FadeIn(overlay), run_time=0.2)
        self.wait(0.1)
        return overlay

    def update_counter(self, e_move: int) -> None:
        new_counter = self.make_text(
            f"{e_move:02d} / 49",
            font_size=42,
            color=GOLD,
            weight="MEDIUM",
        )
        new_counter.scale(self.camera.frame.width / self.base_frame_width)
        new_counter.move_to(self.counter)
        new_counter.set_z_index(20)
        # E is an integer. Morphing Text glyphs made intermediate numbers tear
        # and kept the HUD moving for every drawing animation.
        self.counter.become(new_counter)
        self.pin_counter()

    def move_camera(self, center, width: float) -> None:
        frame = self.camera.frame
        start_center = frame.get_center().copy()
        end_center = logical_to_scene(center)
        start_width = frame.width
        zoom = math.log(width / start_width)
        travel = math.dist(start_center, end_center) / min(start_width, width)
        if abs(zoom) < 1e-10 and travel < 1e-10:
            return
        seconds = min(1.6, 0.8 + 0.25 * abs(zoom) + 0.15 * travel)
        seconds = round(seconds * config.frame_rate) / config.frame_rate

        def interpolate_camera(camera_frame, alpha):
            camera_frame.set(width=start_width * math.exp(zoom * alpha))
            camera_frame.move_to((1 - alpha) * start_center + alpha * end_center)

        self.play(UpdateFromAlphaFunc(frame, interpolate_camera), run_time=seconds)
        self.sync_viewport()
        self.wait(0.1)

    def adjust_camera(self, e_move: int) -> None:
        cue = CAMERA_CUES.get(e_move)
        if cue is None:
            return
        center, width = cue
        self.move_camera(center, width)

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
            position = self.point_at(point_id)
            dot = Dot(position, radius=self.screen_stroke_width(0.012), color=FOREGROUND)
            label = self.make_text(label_text, font_size=14, color=FOREGROUND)
            label.scale(self.camera.frame.width / self.base_frame_width)
            label.next_to(dot, direction, buff=self.screen_stroke_width(0.055))
            group = VGroup(dot, label).set_z_index(7)
            self.keep_screen_size(group, position)
            self.key_point_groups[point_id] = group
            self.play(FadeIn(group), run_time=0.2)

    def make_target_marker(self, event: dict) -> VGroup:
        display_index = next(
            item["display_index"]
            for item in self.data["targets"]
            if item["output_id"] == event["id"]
        )
        center_id = event["references"][0]
        position = self.point_at(center_id)
        dot = Dot(position, radius=self.screen_stroke_width(CENTER_DOT_RADIUS), color=TARGET)
        label = self.make_text(
            f"K{SUBSCRIPTS[display_index]}",
            font_size=15,
            color=TARGET,
            weight="MEDIUM",
        )
        label.scale(self.camera.frame.width / self.base_frame_width)
        label.next_to(dot, TARGET_DIRECTIONS[display_index - 1], buff=self.screen_stroke_width(0.07))
        marker = VGroup(dot, label).set_z_index(12)
        self.keep_screen_size(marker, position)
        return marker

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
            self.update_counter(e_move)
            visible = self.visible_line(event, drawable) if event["op"] == "line" else drawable
            animations = [Create(visible, rate_func=linear)]
            target_marker = None
            if is_target:
                target_marker = self.make_target_marker(event)
                animations.append(FadeIn(target_marker))
            self.play(
                *animations,
                run_time=1.0 if is_target else 0.6,
            )
            if visible is not drawable:
                self.remove(visible)
                self.add(drawable)
            if is_target:
                self.target_drawables.add(drawable)
                self.target_markers.add(target_marker)
                self.play(FadeOut(references), run_time=0.2)
                self.wait(0.2)
            else:
                self.aux_drawables.add(drawable)
                faded_color = CIRCLE_BLUE if event["op"] == "circle" else LINE_BLUE
                drawable.screen_stroke_width = AUX_STROKE
                self.play(
                    drawable.animate.set_stroke(
                        color=faded_color,
                        width=self.screen_stroke_width(AUX_STROKE),
                        opacity=0.28,
                    ),
                    FadeOut(references),
                    run_time=0.2,
                )
            self.reveal_key_points(e_move)

    def finish_construction(self) -> None:
        self.move_camera(TARGET_FRAME_CENTER, TARGET_FRAME_WIDTH)
        remaining_points = VGroup()
        root_center_id = "Mannheim_S_center_locus"
        if root_center_id in self.key_point_groups:
            remaining_points.add(self.key_point_groups[root_center_id])
        for drawable in self.aux_drawables:
            drawable.screen_stroke_width = 0.3
        for drawable in self.target_drawables:
            drawable.screen_stroke_width = 0.95
        finish_animations = [
            self.aux_drawables.animate.set_stroke(width=self.screen_stroke_width(0.3), opacity=0.045),
            self.target_drawables.animate.set_stroke(width=self.screen_stroke_width(0.95), color=TARGET, opacity=0.95),
        ]
        if len(remaining_points):
            finish_animations.append(FadeOut(remaining_points))
        self.play(*finish_animations, run_time=0.65)

        pulses = []
        for target in self.data["targets"]:
            pulse = Circle(
                radius=self.screen_stroke_width(REFERENCE_RADIUS),
                color=TARGET,
                stroke_width=self.screen_stroke_width(MARKER_STROKE),
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
        conclusion.set_stroke(
            BACKGROUND, width=self.screen_stroke_width(2.5), background=True
        ).set_z_index(20)
        self.play(FadeIn(conclusion, shift=UP * 0.08), run_time=0.55)
        self.wait(3.0)
