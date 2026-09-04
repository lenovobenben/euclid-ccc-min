"""符号验证有限 ``K'`` 给出 Mannheim 方向类圆心线。

沿用 ``verify_mannheim_contact_chord_identity`` 的归一化：第三圆是原点
处单位圆，前两圆圆心为 ``(c,h)``、``(c+d,h)``，半径为 ``p,q``。
对每个方向类令 ``K'`` 为圆内接四边形的第二个对角点，再令

    J = Line(O3, K') ∩ Line(O1, O2).

脚本在 ``Z[c,d,h,p,q]`` 中验证 ``J`` 恒位于经过根心 ``S`` 的目标
圆心线上。于是只要 ``K'`` 有限、``K' != O3`` 且 ``J != S``，就能用
两条直线 ``O3 K'``、``S J`` 在完成任一目标前画出该圆心线。
"""

from __future__ import annotations

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
    determinant,
)
from verify_mannheim_tau_concurrence import MANNHEIM_POINT_S


O3 = (0 * c, 0 * c, 1 + 0 * c)
ELL_12 = (0 * c, 1 + 0 * c, -h)
ROOT_ON_ELL = (
    c**2 * d
    + c * d**2
    + c * p**2
    - c * q**2
    + d * h**2
    + d * p**2
    - d
)
ROOT_ON_ELL_CONSTANT = ROOT_ON_ELL - d * h**2


def dot(left, right):
    return sum(
        (first * second for first, second in zip(left, right, strict=True)),
        start=0,
    )


def kp_point(profile_id: str):
    roles = PROFILES[profile_id]["roles"]
    x, y, z, w = tuple(
        conic_point(parameter(kind, source)) for kind, source in roles
    )
    return cross(cross(x, y), cross(z, w))


def merge_parameters(profile_id: str):
    return tuple(
        parameter(kind, source)
        for kind, source in PROFILES[profile_id]["roles"]
    )


def center_direction(profile_id: str):
    _, sigma2, sigma3 = PROFILES[profile_id]["sigma"]
    numerator_uy = c * p - sigma2 * c * q + d * p - sigma3 * d
    return (
        h * (p - sigma2 * q),
        -numerator_uy,
        0 * c,
    )


def vanishes_when_root_on_ell(value: Polynomial) -> bool:
    """以 ``h^2=-R0/d`` 代入并清除分母，检查条件余式为零。"""

    parity_coefficients: dict[int, dict[int, Polynomial]] = {0: {}, 1: {}}
    for exponent, coefficient in value.terms.items():
        parity = exponent[2] % 2
        power = exponent[2] // 2
        reduced = (exponent[0], exponent[1], 0, exponent[3], exponent[4])
        term = Polynomial({reduced: coefficient})
        parity_coefficients[parity][power] = (
            parity_coefficients[parity].get(power, Polynomial()) + term
        )
    for coefficients in parity_coefficients.values():
        if not coefficients:
            continue
        degree = max(coefficients)
        remainder = sum(
            (
                coefficient
                * (-ROOT_ON_ELL_CONSTANT) ** power
                * d ** (degree - power)
                for power, coefficient in coefficients.items()
            ),
            start=Polynomial(),
        )
        if not remainder.is_zero:
            return False
    return True


def main() -> None:
    if dot(MANNHEIM_POINT_S, ELL_12) != -ROOT_ON_ELL:
        raise AssertionError("根心位于 O1O2 的条件多项式错误")
    maximum_degree = 0
    for profile_id in PROFILES:
        kp = kp_point(profile_id)
        radial = cross(O3, kp)
        j = cross(radial, ELL_12)
        center_line = cross(
            MANNHEIM_POINT_S,
            center_direction(profile_id),
        )
        incidence = dot(center_line, j)
        if not incidence.is_zero:
            raise AssertionError(f"{profile_id} 的 J 不在目标圆心线上")
        if any(
            not vanishes_when_root_on_ell(coordinate)
            for coordinate in cross(kp, MANNHEIM_POINT_S)
        ):
            raise AssertionError(
                f"{profile_id} 在 S 属于 O1O2 时没有满足 K'=S"
            )
        maximum_degree = max(
            maximum_degree,
            *(coordinate.total_degree for coordinate in j),
            *(coordinate.total_degree for coordinate in center_line),
        )

    seed_merge_checks = (
        (
            "P0",
            c * p - c * q + d * p - d,
            d - p - q,
            d + p + q,
        ),
        (
            "P2",
            c * p + c * q + d * p + d,
            -d + p - q,
            -d - p + q,
        ),
    )
    for profile_id, common, first_factor, second_factor in seed_merge_checks:
        x, y, z, w = merge_parameters(profile_id)
        merge_xz = determinant(x, z)
        merge_yw = determinant(y, w)
        for actual, expected in (
            (d * merge_xz, first_factor * common),
            (d * merge_yw, second_factor * common),
        ):
            if not vanishes_when_root_on_ell(actual - expected):
                raise AssertionError(
                    f"{profile_id} 的根心在线条件没有强制双对合并"
                )

    print(
        "mannheim_center_locus_kp_identity",
        {
            "profiles": len(PROFILES),
            "incidence_identities": len(PROFILES),
            "root_on_ell_kp_identities": len(PROFILES),
            "root_on_ell_seed_merge_checks": 4,
            "maximum_degree": maximum_degree,
        },
    )


if __name__ == "__main__":
    main()
