"""符号验证 Mannheim 四条接触弦恒过同一射影点。

沿用 ``verify_mannheim_contact_chord_identity`` 的归一化和稀疏整数
多项式环。四个方向类的 ``tau`` 已在该脚本中证明分别恒等于
``contact_line(sigma2, sigma3)``；这里证明四条线都通过同一个与符号
无关的射影点 ``B``，并独立检查四个三线行列式恒为零。
"""

from __future__ import annotations

from itertools import combinations

from verify_mannheim_contact_chord_identity import (
    PROFILES,
    c,
    contact_line,
    cross,
    d,
    h,
    p,
    q,
)


def dot(left, right):
    return sum(
        (first * second for first, second in zip(left, right, strict=True)),
        start=0,
    )


NUMERATOR_BX = 2 * c * d + d**2 + p**2 - q**2
NUMERATOR_BY = (
    c**2 * d
    + c * d**2
    + c * p**2
    - c * q**2
    - d * h**2
    + d * p**2
    - d
)

MANNHEIM_POINT_S = (
    h * NUMERATOR_BX,
    -NUMERATOR_BY,
    2 * d * h,
)

P0_P2_DISTINCT_FACTOR = (
    2
    * d**2
    * h**2
    * p
    * ((q - 1) ** 2 - ((c + d) ** 2 + h**2))
)


def main() -> None:
    lines = {
        profile_id: contact_line(profile["sigma"][1], profile["sigma"][2])
        for profile_id, profile in PROFILES.items()
    }
    for profile_id, line in lines.items():
        if not dot(MANNHEIM_POINT_S, line).is_zero:
            raise AssertionError(f"{profile_id} 的 tau 不通过 Mannheim 共点 S")

    triple_count = 0
    for profile_ids in combinations(lines, 3):
        first, second, third = (lines[profile_id] for profile_id in profile_ids)
        determinant = dot(cross(first, second), third)
        if not determinant.is_zero:
            raise AssertionError(f"三条 tau 不共点：{profile_ids}")
        triple_count += 1

    p0_p2_intersection = cross(lines["P0"], lines["P2"])
    expected_intersection = tuple(
        P0_P2_DISTINCT_FACTOR * coordinate
        for coordinate in MANNHEIM_POINT_S
    )
    if p0_p2_intersection != expected_intersection:
        raise AssertionError("P0、P2 的交点没有按相切因子化为 S")

    x_coordinate, y_coordinate, homogeneous = MANNHEIM_POINT_S
    power_3 = (
        x_coordinate**2 + y_coordinate**2 - homogeneous**2
    )
    power_1 = (
        (x_coordinate - c * homogeneous) ** 2
        + (y_coordinate - h * homogeneous) ** 2
        - (p * homogeneous) ** 2
    )
    power_2 = (
        (x_coordinate - (c + d) * homogeneous) ** 2
        + (y_coordinate - h * homogeneous) ** 2
        - (q * homogeneous) ** 2
    )
    if power_1 != power_3 or power_2 != power_3:
        raise AssertionError("Mannheim 共点 S 不是三个输入圆的根心")

    print(
        "mannheim_tau_concurrence",
        {
            "profiles": len(lines),
            "incidence_identities": len(lines),
            "triple_determinants": triple_count,
            "guaranteed_distinct_pair": "P0,P2",
            "radical_center_identities": 2,
            "point_maximum_degree": max(
                coordinate.total_degree for coordinate in MANNHEIM_POINT_S
            ),
        },
    )


if __name__ == "__main__":
    main()
