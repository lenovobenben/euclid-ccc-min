"""符号验证 Mannheim 四类公共 K' 弦与平行条件。

脚本在 ``Z[c,d,h,p,q]`` 中完成两项检查：

1. 四个方向类的八条 ``K'`` 定义弦恒等为四条几何直线；
2. 每类两条 ``K'`` 弦的方向行列式，除去 ``D8`` 上必非零的输入圆
   相切因子后，分别等于文档第 8.13 节的 ``G0,...,G3``。

实现复用接触弦恒等脚本的稀疏整数多项式，不需要外部计算机代数库。
"""

from __future__ import annotations

from verify_mannheim_contact_chord_identity import (
    PROFILES,
    c,
    conic_point,
    cross,
    d,
    h,
    p,
    parameter,
    q,
)


T = c**2 * d + c * d**2 + c * p**2 - c * q**2 + d * h**2 + d * p**2 + d
G = {
    "P0": T - 2 * c * p + 2 * c * q - 2 * d * p,
    "P1": T + 2 * c * p - 2 * c * q + 2 * d * p,
    "P2": T + 2 * c * p + 2 * c * q + 2 * d * p,
    "P3": T - 2 * c * p - 2 * c * q - 2 * d * p,
}

DISTANCE_13_SQUARED = c**2 + h**2
DISTANCE_23_SQUARED = (c + d) ** 2 + h**2
TANGENCY_FACTORS = {
    "P0": (
        (p - 1) ** 2 - DISTANCE_13_SQUARED,
        (q - 1) ** 2 - DISTANCE_23_SQUARED,
    ),
    "P1": (
        (p + 1) ** 2 - DISTANCE_13_SQUARED,
        (q + 1) ** 2 - DISTANCE_23_SQUARED,
    ),
    "P2": (
        (p + 1) ** 2 - DISTANCE_13_SQUARED,
        (q - 1) ** 2 - DISTANCE_23_SQUARED,
    ),
    "P3": (
        (p - 1) ** 2 - DISTANCE_13_SQUARED,
        (q + 1) ** 2 - DISTANCE_23_SQUARED,
    ),
}


def kp_lines(profile_id: str):
    roles = PROFILES[profile_id]["roles"]
    x, y, z, w = tuple(
        conic_point(parameter(kind, source)) for kind, source in roles
    )
    return cross(x, y), cross(z, w)


def same_projective_line(left, right) -> bool:
    return all(coordinate.is_zero for coordinate in cross(left, right))


def verify_shared_chords() -> None:
    lines = {profile: kp_lines(profile) for profile in PROFILES}
    identities = (
        (("P0", 0), ("P3", 0)),
        (("P0", 1), ("P2", 1)),
        (("P1", 0), ("P2", 0)),
        (("P1", 1), ("P3", 1)),
    )
    for left, right in identities:
        if not same_projective_line(lines[left[0]][left[1]], lines[right[0]][right[1]]):
            raise AssertionError(f"公共 K' 弦恒等失败：{left} != {right}")


def verify_parallel_factors() -> None:
    for profile_id in PROFILES:
        first, second = kp_lines(profile_id)
        direction_determinant = first[0] * second[1] - first[1] * second[0]
        first_factor, second_factor = TANGENCY_FACTORS[profile_id]
        expected = 8 * h * first_factor * second_factor * G[profile_id]
        if direction_determinant != expected:
            raise AssertionError(f"{profile_id} 的 K' 平行因式分解失败")


def verify_pair_relations() -> None:
    identities = (
        (G["P2"] - G["P0"], 4 * p * (c + d)),
        (G["P1"] - G["P3"], 4 * p * (c + d)),
        (G["P0"] - G["P3"], 4 * c * q),
        (G["P1"] - G["P2"], -4 * c * q),
        (G["P0"] + G["P1"], 2 * T),
        (G["P2"] + G["P3"], 2 * T),
    )
    if any(left != right for left, right in identities):
        raise AssertionError("K' 双平行组合的消元恒等式失败")


def main() -> None:
    verify_shared_chords()
    verify_parallel_factors()
    verify_pair_relations()
    print(
        "mannheim_kp_parallel_classification",
        {
            "shared_chords": 4,
            "profiles_factored": len(PROFILES),
            "pair_relations": 6,
        },
    )


if __name__ == "__main__":
    main()
