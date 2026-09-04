"""筛查正规 Mannheim 单根的动态 3 E 直接接触点程序。

初态是合法 5 E 平行前缀与 ``P0`` 四条批量线组成的 9 E 状态。第一步
任取由初态可画的对象 ``A``；第二步对象 ``B`` 必须使用 ``A`` 产生的
新点。若第三步才首次画出经过 ``+++`` 接触点 ``M3`` 的轨迹，它的两个
定义点中至少有一个必须由 ``B`` 新产生，否则该轨迹在第二步前已经
可画。

本轮覆盖其中另一个定义点 ``P`` 在画 ``B`` 前已经存在的情形。对每个
``P`` 预先建立 ``Line(P,M3)``、``Circle(P,M3)`` 和 ``P M3`` 的垂直
平分线，并与画 ``B`` 前的对象求交。第二步只须检查它与旧对象的新
交点是否命中这些有限目标签名，命中即可用第三步画目标直线或圆。

两个定义点都由 ``B`` 同时产生的目标轨迹，以及前两步互相独立但
第三步使用两个派生点的程序不在本轮范围内。这是多实例浮点筛查，
不是 3 E 下界；任何命中仍须用精确算术独立重放。
"""

from __future__ import annotations

from dataclasses import dataclass
from time import monotonic

from search_mannheim_ext_contact_static_3e import (
    MODES,
    contains_targets,
    final_drawable,
    support_locus,
    target_points,
)
from search_mannheim_ext_tau_4e import (
    INITIAL_COST,
    build_state,
    generate_incremental_candidates,
    point_signature,
    tolerant_drawable_signature,
)
from search_parallel_3e import (
    Candidate,
    apply_candidate,
    generate_candidates,
    intersections,
)


PROGRESS_INTERVAL = 100


@dataclass(frozen=True, slots=True)
class TargetRequirement:
    mode: str
    support_id: str
    final_signature: tuple | None


def target_requirements(state, targets):
    known_points = {
        point_signature(point.values)
        for point in state.points
    }
    existing_drawables = {
        tolerant_drawable_signature(drawable)
        for drawable in state.drawables
    }
    by_drawable = [dict() for _ in state.drawables]
    unstable_final_signatures = 0

    for mode in MODES:
        for support in state.points:
            locus_rows = tuple(
                support_locus(mode, support_value, target)
                for support_value, target in zip(
                    support.values,
                    targets,
                    strict=True,
                )
            )
            if any(value is None for _, value in locus_rows):
                continue
            kinds = {kind for kind, _ in locus_rows}
            if len(kinds) != 1:
                continue
            locus_kind = locus_rows[0][0]
            locus_values = tuple(value for _, value in locus_rows)

            for drawable_index, drawable in enumerate(state.drawables):
                rows = tuple(
                    intersections(
                        locus_kind,
                        locus_value,
                        drawable.kind,
                        drawable_value,
                    )
                    for locus_value, drawable_value in zip(
                        locus_values,
                        drawable.values,
                        strict=True,
                    )
                )
                root_counts = {len(row) for row in rows}
                if len(root_counts) != 1 or not rows[0]:
                    continue
                for root_index in range(len(rows[0])):
                    values = tuple(row[root_index] for row in rows)
                    signature = point_signature(values)
                    if signature in known_points:
                        continue
                    final = final_drawable(mode, support, values)
                    if final is None:
                        continue
                    if contains_targets(final, targets):
                        final_signature = tolerant_drawable_signature(final)
                        if final_signature in existing_drawables:
                            continue
                    else:
                        # 近乎平行或相切的浮点交点会放大重建误差。保留
                        # 必要点签名；若它实际命中，再用精确算术判定。
                        final_signature = None
                        unstable_final_signatures += 1
                    by_drawable[drawable_index].setdefault(
                        signature,
                        TargetRequirement(
                            mode,
                            support.point_id,
                            final_signature,
                        ),
                    )
    return tuple(by_drawable), unstable_final_signatures


def candidate_hits_requirements(
    candidate: Candidate,
    state,
    requirements,
):
    for drawable, requirement_map in zip(
        state.drawables,
        requirements,
        strict=True,
    ):
        if not requirement_map:
            continue
        rows = tuple(
            intersections(
                candidate.kind,
                candidate_value,
                drawable.kind,
                drawable_value,
            )
            for candidate_value, drawable_value in zip(
                candidate.drawable.values,
                drawable.values,
                strict=True,
            )
        )
        root_counts = {len(row) for row in rows}
        if len(root_counts) != 1 or not rows[0]:
            continue
        for root_index in range(len(rows[0])):
            values = tuple(row[root_index] for row in rows)
            signature = point_signature(values)
            requirement = requirement_map.get(signature)
            if requirement is not None:
                return requirement
    return None


def main() -> None:
    start = monotonic()
    state, _ = build_state()
    targets = target_points()
    first_candidates = generate_candidates(state, INITIAL_COST + 1)

    incremental_second_candidates = 0
    requirement_signatures = 0
    unstable_final_signatures = 0
    hits = []

    for first_index, first in enumerate(first_candidates, start=1):
        state_one = apply_candidate(state, first, INITIAL_COST + 1)
        requirements, unstable_count = target_requirements(
            state_one,
            targets,
        )
        requirement_signatures += sum(len(row) for row in requirements)
        unstable_final_signatures += unstable_count
        existing_drawables = {
            tolerant_drawable_signature(drawable)
            for drawable in state_one.drawables
        }

        for second in generate_incremental_candidates(
            state,
            state_one,
            INITIAL_COST + 2,
        ):
            incremental_second_candidates += 1
            second_signature = tolerant_drawable_signature(second.drawable)
            if second_signature in existing_drawables:
                continue
            requirement = candidate_hits_requirements(
                second,
                state_one,
                requirements,
            )
            if requirement is None:
                continue
            if (
                requirement.final_signature is not None
                and requirement.final_signature == second_signature
            ):
                continue
            hits.append(
                (
                    first.describe(),
                    second.describe(),
                    requirement.mode,
                    requirement.support_id,
                )
            )
            if len(hits) >= 3:
                break
        if hits:
            break

        if first_index % PROGRESS_INTERVAL == 0:
            print(
                "progress",
                {
                    "first": first_index,
                    "first_total": len(first_candidates),
                    "requirement_signatures": requirement_signatures,
                    "unstable_final_signatures": (
                        unstable_final_signatures
                    ),
                    "incremental_second_candidates": (
                        incremental_second_candidates
                    ),
                    "elapsed_seconds": round(monotonic() - start, 3),
                },
                flush=True,
            )

    if hits:
        raise AssertionError(
            "发现动态 3 E 直接接触轨迹，需要精确重放："
            f"{hits}"
        )

    print(
        "mannheim_ext_contact_dynamic_3e_search",
        {
            "samples": len(targets),
            "initial_points": len(state.points),
            "initial_drawables": len(state.drawables),
            "first_candidates": len(first_candidates),
            "requirement_signatures": requirement_signatures,
            "unstable_final_signatures": unstable_final_signatures,
            "incremental_second_candidates": incremental_second_candidates,
            "three_e_contact_locus_hits": 0,
            "elapsed_seconds": round(monotonic() - start, 3),
        },
    )


if __name__ == "__main__":
    main()
