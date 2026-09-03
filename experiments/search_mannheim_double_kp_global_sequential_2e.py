"""搜索已完成 P0/P1 正规模块后的双 ``K'`` 全局复用。

对三个仅含 ``P2:parallel:Kp`` 与 ``P3:parallel:Kp`` 的严格 D8 有理夹具，
先同步重放公共前缀、八条批量线以及 P0/P1 的四个目标圆。状态加入这些
正规模块的全部计费对象、显式构造点，以及新对象与此前对象产生的有限
实交点。程序先检查 P2、P3 的目标线是否已经各有两个已知点；若没有，
再先做一块的标准 3 E 修复并枚举一个辅助对象。

命中必须精确重放。当前直接命中的共同点由 P0、P1 两条 ``tau`` 相交
产生，其一般恒等和精确 E 证书另见共点验证与重放脚本。
"""

from __future__ import annotations

from math import hypot
from time import monotonic

from replay_mannheim_fixed import Quadratic
from replay_mannheim_three_block_dependencies import (
    ExactObjectGraph,
    ThreeBlockReplay,
)
from search_mannheim_double_kp_joint_5e import (
    BATCH_KEYS,
    FIXTURES,
    build_state,
    candidate_target_witness,
    target_point_pair,
)
from search_mannheim_double_kp_sequential_2e import (
    build_after_repair,
    point,
)
from search_parallel_3e import (
    Candidate,
    DrawableBundle,
    PointBundle,
    State,
    apply_candidate,
    drawable_signature,
    generate_candidates,
    point_signature,
)


def float_scalar(value) -> float:
    if isinstance(value, Quadratic):
        return value.approximate()
    return float(value)


def float_line(value) -> tuple[float, float, float]:
    a, b, c = (float_scalar(coordinate) for coordinate in value)
    norm = hypot(a, b)
    a, b, c = a / norm, b / norm, c / norm
    if a < -1e-9 or (abs(a) <= 1e-9 and b < 0):
        a, b, c = -a, -b, -c
    return (a, b, c)


class RecordingObjectGraph(ExactObjectGraph):
    """记录精确重放器交给对象注册器的实际几何值。"""

    def __init__(self, centers, radii) -> None:
        super().__init__(centers, radii)
        self.point_values: dict[str, tuple[float, float]] = {}
        self.drawable_values: dict[str, tuple[str, tuple[float, ...]]] = {}

    def point(self, node_id: str, value, *dependencies: str) -> str:
        physical_id = super().point(node_id, value, *dependencies)
        self.point_values[physical_id] = tuple(
            float_scalar(coordinate) for coordinate in value
        )
        return physical_id

    def line(
        self,
        node_id: str,
        value,
        *dependencies: str,
        special_key: str | None = None,
    ) -> str:
        physical_id = super().line(
            node_id,
            value,
            *dependencies,
            special_key=special_key,
        )
        if value is not None:
            self.drawable_values[physical_id] = ("line", float_line(value))
        return physical_id

    def circle(
        self,
        node_id: str,
        value,
        *dependencies: str,
        special_key: str | None = None,
    ) -> str:
        physical_id = super().circle(
            node_id,
            value,
            *dependencies,
            special_key=special_key,
        )
        if value is not None:
            center, radius_squared = value
            self.drawable_values[physical_id] = (
                "circle",
                (
                    float_scalar(center[0]),
                    float_scalar(center[1]),
                    float_scalar(radius_squared),
                ),
            )
        return physical_id


def replay_regular_blocks(centers, radii):
    replay = ThreeBlockReplay(centers, radii)
    recorder = RecordingObjectGraph(centers, radii)
    replay.objects = recorder
    replay.build_prefix()
    for key in BATCH_KEYS:
        replay.draw_batch(key)
    replay.build_regular_pair(
        "P0",
        None,
        allow_repeated_physical_signs=True,
    )
    replay.build_regular_pair(
        "P1",
        None,
        allow_repeated_physical_signs=True,
    )
    return recorder


def build_global_state():
    state, targets = build_state()
    recorders = tuple(
        replay_regular_blocks(centers, radii) for centers, radii in FIXTURES
    )
    paid_order = recorders[0].graph.paid_order
    if any(recorder.graph.paid_order != paid_order for recorder in recorders[1:]):
        raise AssertionError("三个夹具的 P0/P1 计费轨迹不一致")
    if any(
        set(recorder.point_values) != set(recorders[0].point_values)
        for recorder in recorders[1:]
    ):
        raise AssertionError("三个夹具的 P0/P1 显式点轨迹不一致")

    existing_drawables = {
        drawable_signature(drawable) for drawable in state.drawables
    }
    next_move = 100
    added_drawables = 0
    for node_id in paid_order:
        rows = tuple(
            recorder.drawable_values.get(node_id) for recorder in recorders
        )
        if any(row is None for row in rows):
            if node_id != "parallel_cX":
                raise AssertionError(f"缺少计费对象的几何值：{node_id}")
            continue
        kinds = {row[0] for row in rows if row is not None}
        if len(kinds) != 1:
            raise AssertionError(f"计费对象类型不一致：{node_id}")
        kind = kinds.pop()
        values = tuple(row[1] for row in rows if row is not None)
        drawable = DrawableBundle(node_id, kind, values)
        signature = drawable_signature(drawable)
        if signature in existing_drawables:
            continue
        candidate = Candidate(kind, node_id, node_id, drawable)
        state = apply_candidate(state, candidate, next_move)
        existing_drawables.add(signature)
        next_move += 1
        added_drawables += 1

    known_points = {point_signature(item) for item in state.points}
    explicit_points_added = 0
    for point_id in recorders[0].point_values:
        point_bundle = PointBundle(
            point_id,
            tuple(
                recorder.point_values[point_id] for recorder in recorders
            ),
        )
        signature = point_signature(point_bundle)
        if signature in known_points:
            continue
        known_points.add(signature)
        state = State(
            state.points + (point_bundle,),
            state.drawables,
            state.program,
        )
        explicit_points_added += 1
    return state, targets, added_drawables, explicit_points_added


def main() -> None:
    initial_state, targets, added_drawables, explicit_points_added = (
        build_global_state()
    )
    print(
        "global_state",
        {
            "points": len(initial_state.points),
            "drawables": len(initial_state.drawables),
            "regular_drawables_added": added_drawables,
            "explicit_points_added": explicit_points_added,
        },
        flush=True,
    )

    direct_pairs = {
        profile: target_point_pair(initial_state, targets[profile])
        for profile in ("P2", "P3")
    }
    print("direct_target_pairs", direct_pairs, flush=True)
    if all(pair is not None for pair in direct_pairs.values()):
        return

    for first_profile, remaining_profile in (("P2", "P3"), ("P3", "P2")):
        state, _ = build_after_repair(
            first_profile,
            (initial_state, targets),
        )
        finite = point(state, f"{remaining_profile}_K")
        if target_point_pair(state, targets[remaining_profile]) is not None:
            raise AssertionError(
                f"{remaining_profile} 在枚举辅助对象前已经完成"
            )

        start = monotonic()
        candidates = generate_candidates(state, 4)
        hits = []
        for candidate_index, candidate in enumerate(candidates, start=1):
            witness = candidate_target_witness(
                candidate,
                state,
                targets[remaining_profile],
                finite,
            )
            if witness is not None:
                hits.append((candidate.describe(), witness))
            if candidate_index % 100_000 == 0:
                print(
                    "progress",
                    {
                        "first_profile": first_profile,
                        "checked": candidate_index,
                        "total": len(candidates),
                        "hits": len(hits),
                        "elapsed_seconds": round(monotonic() - start, 3),
                    },
                    flush=True,
                )
        print(
            "summary",
            {
                "samples": len(FIXTURES),
                "first_profile": first_profile,
                "remaining_profile": remaining_profile,
                "points": len(state.points),
                "drawables": len(state.drawables),
                "candidates": len(candidates),
                "hits": len(hits),
                "elapsed_seconds": round(monotonic() - start, 3),
            },
        )
        for hit in hits:
            print("candidate", hit)


if __name__ == "__main__":
    main()
