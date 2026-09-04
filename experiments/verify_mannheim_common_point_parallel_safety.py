"""验证后行平行块不能以有限对角点等于根心逃逸 3 E 上界。

共点程序先由 P0、P2 得到根心 S，再处理 P1、P3。四个角色两两不同时，
若一个对角点位于无穷远，就要用另一个有限对角点 D 画 ``Line(S,D)``。
本脚本在 ``Z[c,d,h,p,q]`` 中验证：同时要求 D=S 会落入输入圆相切、
非实参数或 P3 双对合并，因而不会产生遗漏的严格 D8 四点分支。

脚本只验证文档使用的多项式因式分解；D8 不等式排除各因子的论证见
``docs/MANNHEIM.md`` 第 8.15 节。
"""

from __future__ import annotations

from verify_mannheim_contact_chord_identity import (
    PROFILES,
    Polynomial,
    c,
    conic_point,
    cross,
    d,
    determinant,
    h,
    p,
    parameter,
    q,
)


def dot(left, right):
    return sum(
        (first * second for first, second in zip(left, right, strict=True)),
        start=0,
    )


def split_h_squared(value: Polynomial) -> tuple[Polynomial, Polynomial]:
    constant_terms = {}
    coefficient_terms = {}
    for exponent, coefficient in value.terms.items():
        if exponent[2] == 0:
            constant_terms[exponent] = coefficient
        elif exponent[2] == 2:
            reduced = (exponent[0], exponent[1], 0, exponent[3], exponent[4])
            coefficient_terms[reduced] = coefficient
        else:
            raise AssertionError("多项式不是 h^2 的一次式")
    return Polynomial(constant_terms), Polynomial(coefficient_terms)


def resultant_h_squared(left: Polynomial, right: Polynomial) -> Polynomial:
    left_constant, left_coefficient = split_h_squared(left)
    right_constant, right_coefficient = split_h_squared(right)
    return left_coefficient * right_constant - right_coefficient * left_constant


N_X = 2 * c * d + d**2 + p**2 - q**2
N_Y = (
    c**2 * d
    + c * d**2
    + c * p**2
    - c * q**2
    - d * h**2
    + d * p**2
    - d
)
S = (h * N_X, -N_Y, 2 * d * h)

T = c**2 * d + c * d**2 + c * p**2 - c * q**2 + d * h**2 + d * p**2 + d
G1 = T + 2 * c * p - 2 * c * q + 2 * d * p
G3 = T - 2 * c * p - 2 * c * q - 2 * d * p
H1 = (
    -c**2 * p
    - c**2 * q
    - 2 * c**2
    - 2 * c * d * p
    - 2 * c * d
    - d**2 * p
    - d**2
    - h**2 * p
    - h**2 * q
    - 2 * h**2
    + p**2 * q
    + p**2
    + p * q**2
    + 4 * p * q
    + 3 * p
    + q**2
    + 3 * q
    + 2
)
H3 = (
    c**2 * p
    - c**2 * q
    - 2 * c**2
    + 2 * c * d * p
    - 2 * c * d
    + d**2 * p
    - d**2
    + h**2 * p
    - h**2 * q
    - 2 * h**2
    + p**2 * q
    + p**2
    - p * q**2
    - 4 * p * q
    - 3 * p
    + q**2
    + 3 * q
    + 2
)

R = c**2 * d + c * d**2 + c * p**2 - c * q**2 + d * h**2 + d * p**2 - d
A1 = (
    -c**2 * p
    - c**2 * q
    - 2 * c * d * p
    - d**2 * p
    - d**2
    - h**2 * p
    - h**2 * q
    + p**2 * q
    + p**2
    + p * q**2
    + 2 * p * q
    + p
    + q**2
    + q
)
A3 = (
    c**2 * p
    - c**2 * q
    + 2 * c * d * p
    + d**2 * p
    - d**2
    + h**2 * p
    - h**2 * q
    + p**2 * q
    + p**2
    - p * q**2
    - 2 * p * q
    - p
    + q**2
    + q
)

C1 = c * p - c * q + d * p + d
D1 = -d**2 + (p + q) * (p + q + 2)
C3 = c * p + c * q + d * p - d
D3 = -d**2 + (p - q) * (p - q - 2)


def role_points(profile_id: str):
    return tuple(
        conic_point(parameter(kind, source))
        for kind, source in PROFILES[profile_id]["roles"]
    )


def line(first, second):
    return cross(first, second)


def direction_determinant(first, second):
    return first[0] * second[1] - first[1] * second[0]


def verify_profile_p1() -> None:
    x, y, z, w = role_points("P1")
    k_lines = (line(x, w), line(y, z))
    kp_lines = (line(x, y), line(z, w))
    tangency_13 = (p + 1) ** 2 - c**2 - h**2
    tangency_23 = (q + 1) ** 2 - (c + d) ** 2 - h**2

    expected = (
        (dot(k_lines[0], S), -2 * h**2 * (-d + p - q) ** 2 * A1),
        (dot(k_lines[1], S), -2 * h**2 * (d + p - q) ** 2 * A1),
        (dot(kp_lines[0], S), 2 * tangency_13**2 * R),
        (dot(kp_lines[1], S), -2 * tangency_23**2 * R),
        (
            direction_determinant(*k_lines),
            -8 * h**3 * (-d + p - q) * (d + p - q) * H1,
        ),
        (
            direction_determinant(*kp_lines),
            8 * h * tangency_13 * tangency_23 * G1,
        ),
        (resultant_h_squared(A1, G1), -C1 * D1),
        (resultant_h_squared(R, H1), C1 * D1),
    )
    if any(actual != wanted for actual, wanted in expected):
        raise AssertionError("P1 的有限对角点安全因式分解失败")


def verify_profile_p3() -> None:
    x, y, z, w = role_points("P3")
    k_lines = (line(x, w), line(y, z))
    kp_lines = (line(x, y), line(z, w))
    tangency_13 = (p - 1) ** 2 - c**2 - h**2
    tangency_23 = (q + 1) ** 2 - (c + d) ** 2 - h**2

    expected = (
        (dot(k_lines[0], S), -2 * h**2 * (-d + p + q) ** 2 * A3),
        (dot(k_lines[1], S), -2 * h**2 * (d + p + q) ** 2 * A3),
        (dot(kp_lines[0], S), -2 * tangency_13**2 * R),
        (dot(kp_lines[1], S), 2 * tangency_23**2 * R),
        (
            direction_determinant(*k_lines),
            8 * h**3 * (-d + p + q) * (d + p + q) * H3,
        ),
        (
            direction_determinant(*kp_lines),
            8 * h * tangency_13 * tangency_23 * G3,
        ),
        (resultant_h_squared(A3, G3), C3 * D3),
        (resultant_h_squared(R, H3), -C3 * D3),
    )
    if any(actual != wanted for actual, wanted in expected):
        raise AssertionError("P3 的有限对角点安全因式分解失败")

    merge_parameters = tuple(
        parameter(kind, source) for kind, source in PROFILES["P3"]["roles"]
    )
    merge_xz = determinant(merge_parameters[0], merge_parameters[2])
    merge_yw = determinant(merge_parameters[1], merge_parameters[3])
    merge_expected = (
        (resultant_h_squared(A3, merge_xz), (-d + p - q) * C3),
        (resultant_h_squared(A3, merge_yw), (d + p - q) * C3),
        (resultant_h_squared(R, merge_xz), -(-d + p - q) * C3),
        (resultant_h_squared(R, merge_yw), (d + p - q) * C3),
    )
    if any(actual != wanted for actual, wanted in merge_expected):
        raise AssertionError("P3 的残余因子没有强制双对合并")


def main() -> None:
    verify_profile_p1()
    verify_profile_p3()
    print(
        "mannheim_common_point_parallel_safety",
        {
            "profiles": ("P1", "P3"),
            "parallel_types": ("K", "Kp"),
            "incidence_factorizations": 8,
            "eliminants": 4,
            "double_merge_checks": 4,
        },
    )


if __name__ == "__main__":
    main()
