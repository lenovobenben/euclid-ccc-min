"""逐类验证简单合并目标的三条强制直线复用。

八个严格 ``D8`` 整数夹具分别覆盖四个方向类的 ``x=z`` 与 ``y=w``。
对以合并点为第三圆接触点的目标，脚本断言：

* 接触恢复线复用产生该点的 Gamma2 侧批量线；
* 第三圆半径复用合并修复开头的半径线；
* 第二圆半径复用公共圆心线 ``ell``。
"""

from __future__ import annotations

from fractions import Fraction

from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
FIXTURES = {
    ("P0", "xz"): (
        ((F(0), F(0)), (F(7), F(0)), (F(4), F(2))),
        (F(3), F(2), F(1)),
        "P0:+++@-1",
        "batch_a1A",
    ),
    ("P0", "yw"): (
        ((F(0), F(0)), (F(7), F(0)), (F(3), F(5))),
        (F(3), F(2), F(1)),
        "P0:---@+1",
        "batch_alpha1B",
    ),
    ("P1", "xz"): (
        ((F(0), F(0)), (F(9), F(0)), (F(5), F(1))),
        (F(3), F(2), F(1)),
        "P1:++-@-1",
        "batch_a1B",
    ),
    ("P1", "yw"): (
        ((F(0), F(0)), (F(6), F(0)), (F(5), F(6))),
        (F(3), F(2), F(1)),
        "P1:--+@+1",
        "batch_alpha1A",
    ),
    ("P2", "xz"): (
        ((F(0), F(0)), (F(8), F(0)), (F(5), F(2))),
        (F(3), F(2), F(1)),
        "P2:+--@-1",
        "batch_alpha1B",
    ),
    ("P2", "yw"): (
        ((F(0), F(0)), (F(7), F(0)), (F(-3), F(3))),
        (F(3), F(2), F(1)),
        "P2:-++@+1",
        "batch_a1A",
    ),
    ("P3", "xz"): (
        ((F(0), F(0)), (F(7), F(0)), (F(6), F(4))),
        (F(3), F(2), F(1)),
        "P3:+-+@-1",
        "batch_alpha1A",
    ),
    ("P3", "yw"): (
        ((F(0), F(0)), (F(9), F(0)), (F(2), F(4))),
        (F(3), F(2), F(1)),
        "P3:-+-@+1",
        "batch_a1B",
    ),
}

BATCH_KEYS = (
    "alphaA",
    "aB",
    "a1A",
    "alpha1B",
    "alphaB",
    "aA",
    "a1B",
    "alpha1A",
)


def verify_fixture(profile: str, merge_kind: str) -> dict[str, str]:
    centers, radii, target_key, batch_line = FIXTURES[(profile, merge_kind)]
    if not is_d8(centers, radii):
        raise AssertionError(f"{profile} {merge_kind} 夹具不属于 D8")
    events = analyze_fixture(centers, radii)
    merge_token = "a2=a2_prime" if merge_kind == "xz" else "alpha2=alpha2_prime"
    if f"{profile}:merge:{merge_token}" not in events:
        raise AssertionError(f"{profile} {merge_kind} 合并事件缺失")

    replay = ThreeBlockReplay(centers, radii)
    replay.build_prefix()
    for key in BATCH_KEYS:
        replay.draw_batch(key)
    replay.build_merge_pair(
        profile,
        None,
        allow_repeated_physical_signs=True,
    )

    expected = {
        f"{profile}_{target_key}_contact_line": batch_line,
        f"{profile}_{target_key}_radius_3": f"{profile}_merged_radius",
        f"{profile}_{target_key}_radius_2": "ell",
    }
    actual = replay.objects.paid_aliases
    for logical_id, paid_id in expected.items():
        if actual.get(logical_id) != paid_id:
            raise AssertionError(
                f"{profile} {merge_kind} 缺少复用 {logical_id} -> {paid_id}"
            )
    return expected


def main() -> None:
    reports = {
        f"{profile}_{merge_kind}": verify_fixture(profile, merge_kind)
        for profile, merge_kind in FIXTURES
    }
    print("mannheim_merge_target_aliases", reports)


if __name__ == "__main__":
    main()
