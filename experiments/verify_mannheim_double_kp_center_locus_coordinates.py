"""精确验证双 ``K'`` 平行族的两个圆心线截点公式。

沿用第三圆为原点单位圆、前两圆圆心为 ``(c,h)``、``(c+d,h)`` 的
归一化。唯一可能的双 ``K'`` 平行组合是 ``P2/P3``，其条件可写成

    c(p+q) + dp = 0,
    T = 0.

令 ``A,B`` 为 ``P0`` 的两条 ``K'`` 弦与 ``ell=O1O2`` 的交点，令
``J2,J3`` 为 ``P2,P3`` 目标圆心线与 ``ell`` 的交点。脚本清除所有
齐次分母，并在 ``Z[c,d,h,p,q]`` 中验证

    B-A  = (p+q)^2 / (d p q),
    J2-A = q(p+1)/(p+q) * (B-A),
    J3-A = -q(p-1)/(p+q) * (B-A).

验证器以双平行条件消去 ``c`` 与 ``h^2``；不使用浮点数或样本代入。
"""

from __future__ import annotations

from verify_mannheim_center_locus_kp_identity import (
    ELL_12,
    MANNHEIM_POINT_S,
    center_direction,
)
from verify_mannheim_contact_chord_identity import (
    PROFILES,
    Polynomial,
    c,
    conic_point,
    cross,
    d,
    h,
    p,
    parameter,
    q,
)


RADIUS_SUM = p + q
H2_NUMERATOR = p * q * (d**2 - RADIUS_SUM**2) - RADIUS_SUM**2


def vanishes_on_double_kp(value: Polynomial) -> bool:
    """代入 ``c=-dp/(p+q)`` 与双平行给出的 ``h^2``。"""

    grouped: dict[int, list[tuple[tuple[int, ...], int, int]]] = {
        0: [],
        1: [],
    }
    for exponent, coefficient in value.terms.items():
        c_power, d_power, h_power, p_power, q_power = exponent
        h_square_power = h_power // 2
        denominator_power = c_power + 2 * h_square_power
        grouped[h_power % 2].append(
            (exponent, coefficient, denominator_power)
        )

    for terms in grouped.values():
        if not terms:
            continue
        common_denominator_power = max(item[2] for item in terms)
        reduced = Polynomial()
        for exponent, coefficient, denominator_power in terms:
            c_power, d_power, h_power, p_power, q_power = exponent
            h_square_power = h_power // 2
            reduced += (
                coefficient
                * (-d * p) ** c_power
                * d**d_power
                * p**p_power
                * q**q_power
                * H2_NUMERATOR**h_square_power
                * RADIUS_SUM ** (
                    common_denominator_power - denominator_power
                )
            )
        if not reduced.is_zero:
            return False
    return True


def p0_chord_intercepts():
    roles = PROFILES["P0"]["roles"]
    x, y, z, w = tuple(
        conic_point(parameter(kind, source)) for kind, source in roles
    )
    first = cross(cross(x, y), ELL_12)
    second = cross(cross(z, w), ELL_12)
    return first, second


def center_locus_intercept(profile: str):
    center_line = cross(MANNHEIM_POINT_S, center_direction(profile))
    return cross(center_line, ELL_12)


def affine_difference_numerator(first, second):
    """返回 ``first_x-second_x`` 型仿射差的齐次分子。"""

    return first[0] * second[2] - second[0] * first[2]


def dot(first, second):
    return sum(
        (left * right for left, right in zip(first, second, strict=True)),
        start=Polynomial(),
    )


def main() -> None:
    a, b = p0_chord_intercepts()
    j2 = center_locus_intercept("P2")
    j3 = center_locus_intercept("P3")

    b_minus_a = affine_difference_numerator(b, a)
    distance_identity = (
        d * p * q * b_minus_a
        - RADIUS_SUM**2 * a[2] * b[2]
    )
    j2_identity = (
        RADIUS_SUM
        * affine_difference_numerator(j2, a)
        * b[2]
        - q
        * (p + 1)
        * b_minus_a
        * j2[2]
    )
    j3_identity = (
        RADIUS_SUM
        * affine_difference_numerator(j3, a)
        * b[2]
        + q
        * (p - 1)
        * b_minus_a
        * j3[2]
    )

    identities = {
        "B_minus_A": distance_identity,
        "J2_ratio": j2_identity,
        "J3_ratio": j3_identity,
    }
    failed = tuple(
        name
        for name, value in identities.items()
        if not vanishes_on_double_kp(value)
    )
    if failed:
        raise AssertionError(f"双 K' 平行截点恒等式失败：{failed}")

    # 三个同半径样本曾给出 Circle(A,B) 的伪通用候选。该圆与 ell 的
    # 另一交点是 2A-B；以下因式分解精确记录它只在额外子族成立。
    reflected = (
        2 * a[0] * b[2] - b[0] * a[2],
        2 * a[1] * b[2] - b[1] * a[2],
        a[2] * b[2],
    )
    p3_center_line = cross(MANNHEIM_POINT_S, center_direction("P3"))
    reflected_incidence = dot(p3_center_line, reflected)
    first_tangency_factor = -c**2 - h**2 + p**2 - 2 * p + 1
    second_tangency_factor = (
        -c**2 - 2 * c * d - d**2 - h**2 + q**2 - 2 * q + 1
    )
    special_factor = 2 * c * p + c * q - c + 2 * d * p - 2 * d
    root_on_ell = (
        c**2 * d
        + c * d**2
        + c * p**2
        - c * q**2
        + d * h**2
        + d * p**2
        - d
    )
    reflected_expected = (
        16
        * d
        * h**3
        * first_tangency_factor
        * special_factor
        * second_tangency_factor
        * root_on_ell
    )
    if reflected_incidence != reflected_expected:
        raise AssertionError("Circle(A,B) 候选的入射因式分解失败")
    if not vanishes_on_double_kp(
        RADIUS_SUM * special_factor - d * (p * q - p - 2 * q)
    ):
        raise AssertionError("Circle(A,B) 的额外半径条件化简失败")

    print(
        "mannheim_double_kp_center_locus_coordinates",
        {
            "identities": len(identities),
            "rejected_special_candidate_checks": 2,
            "maximum_degree": max(
                value.total_degree
                for value in (*identities.values(), reflected_incidence)
            ),
        },
    )


if __name__ == "__main__":
    main()
