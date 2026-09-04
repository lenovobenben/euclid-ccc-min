"""筛查 Mannheim 目标圆是否经过其它目标圆心。

若先画出的目标圆经过另一个目标圆心，它可能与另一已有对象相交并免费
产生该圆心，从而支持 5 E 交错后缀。本脚本在三个严格正规 ``D8``
夹具上检查全部 56 个有向“目标圆—其它目标圆心”关系。

同一方向类的两个根位于同一个实二次域，本脚本对这 8 个有向关系作精确
判定；不同方向类通常位于不同二次扩域，当前轻量精确数类不能直接混合
它们，仍只做多夹具浮点筛查。跨方向零结果不是一般参数证明。
"""

from __future__ import annotations

from replay_mannheim_three_block_dependencies import (
    collapse_circle,
    collapse_point,
)
from search_mannheim_double_kp_global_sequential_2e import float_scalar
from search_mannheim_kp_parallel_center_locus_2e import CoreOnlyReplay
from search_mannheim_regular_sequential_locus_reuse import PROFILES
from search_mannheim_root_center_2e import FIXTURES


TOLERANCE = 1e-7


def float_point(point):
    return tuple(float_scalar(coordinate) for coordinate in point)


def float_circle(value):
    center, radius_squared = value
    center_x, center_y = float_point(center)
    return center_x, center_y, float_scalar(radius_squared)


def on_circle(point, value) -> bool:
    center_x, center_y, radius_squared = value
    residual = abs(
        (point[0] - center_x) ** 2
        + (point[1] - center_y) ** 2
        - radius_squared
    )
    return residual <= TOLERANCE * max(1.0, radius_squared)


def on_circle_exact(point, value) -> bool:
    center, radius_squared = value
    delta_x = point[0] - center[0]
    delta_y = point[1] - center[1]
    return delta_x * delta_x + delta_y * delta_y == radius_squared


def targets_for_fixture(index, centers, radii):
    replay = CoreOnlyReplay(
        f"target_circle_center_incidence_{index}",
        centers=centers,
        radii=radii,
        emit=False,
    )
    replay.run()
    targets = {}
    for profile in PROFILES:
        data = replay.core_data[profile]
        for key, target in replay.verify_pair(
            profile,
            data["tau"],
            allow_repeated_physical_signs=True,
        ).items():
            targets[f"{profile}:{key}"] = {
                "profile": profile,
                "exact_center": collapse_point(target["center"]),
                "exact_circle": collapse_circle(target["output_circle"]),
            }
    if len(targets) != 8:
        raise AssertionError("严格正规夹具没有恰好八个目标")
    return {
        key: {
            **target,
            "center": float_point(target["exact_center"]),
            "circle": float_circle(target["exact_circle"]),
        }
        for key, target in targets.items()
    }


def main() -> None:
    samples = tuple(
        targets_for_fixture(index, centers, radii)
        for index, (centers, radii) in enumerate(FIXTURES)
    )
    keys = tuple(samples[0])
    if any(tuple(sample) != keys for sample in samples[1:]):
        raise AssertionError("三个夹具的目标键不一致")

    hits = []
    fixture_incidences = []
    exact_same_profile_incidences = []
    for sample in samples:
        incidences = tuple(
            (source, target)
            for source in keys
            for target in keys
            if source != target
            and on_circle(sample[target]["center"], sample[source]["circle"])
        )
        fixture_incidences.append(incidences)
        exact_same_profile_incidences.append(
            tuple(
                (source, target)
                for source in keys
                for target in keys
                if source != target
                and sample[source]["profile"] == sample[target]["profile"]
                and on_circle_exact(
                    sample[target]["exact_center"],
                    sample[source]["exact_circle"],
                )
            )
        )

    for source in keys:
        for target in keys:
            if source == target:
                continue
            if all(
                on_circle(sample[target]["center"], sample[source]["circle"])
                for sample in samples
            ):
                hits.append((source, target))

    print(
        "mannheim_target_circle_center_incidence_search",
        {
            "samples": len(samples),
            "directed_pairs": len(keys) * (len(keys) - 1),
            "exact_same_profile_pairs": sum(
                1
                for source in keys
                for target in keys
                if source != target
                and samples[0][source]["profile"]
                == samples[0][target]["profile"]
            ),
            "exact_same_profile_incidences": tuple(
                len(rows) for rows in exact_same_profile_incidences
            ),
            "fixture_incidences": tuple(
                len(rows) for rows in fixture_incidences
            ),
            "common_hits": len(hits),
        },
    )
    for hit in hits:
        print("candidate", hit)


if __name__ == "__main__":
    main()
