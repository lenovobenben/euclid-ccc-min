"""符号验证 Mannheim 极线恒等于对应 Apollonius 接触弦。

脚本只使用 Python 标准库。它在整数系数多项式环
``Z[c,d,h,p,q]`` 中展开四个方向类，验证文档第 8.12 节的恒等式

    d * tau = epsilon * 8 * D_xz * D_yw * lambda_sigma.

这里 ``tau`` 是由圆内接四边形的第三个对角点取极线得到的齐次直线，
``lambda_sigma`` 是相应有向 Apollonius 二次式两根的固定接触弦，
``D_xz``、``D_yw`` 是两个对向角色的齐次参数行列式。恒等式没有数值
代入、浮点容差或外部计算机代数依赖。
"""

from __future__ import annotations

from dataclasses import dataclass


Exponent = tuple[int, int, int, int, int]
ZERO_EXPONENT: Exponent = (0, 0, 0, 0, 0)


@dataclass(frozen=True, slots=True)
class Polynomial:
    """五变量稀疏整数系数多项式。"""

    terms: dict[Exponent, int]

    def __init__(self, value: int | dict[Exponent, int] = 0) -> None:
        if isinstance(value, int):
            terms = {} if value == 0 else {ZERO_EXPONENT: value}
        else:
            terms = {
                exponent: coefficient
                for exponent, coefficient in value.items()
                if coefficient != 0
            }
        object.__setattr__(self, "terms", terms)

    @staticmethod
    def coerce(value: int | Polynomial) -> Polynomial:
        return value if isinstance(value, Polynomial) else Polynomial(value)

    def __add__(self, other: int | Polynomial) -> Polynomial:
        result = dict(self.terms)
        for exponent, coefficient in self.coerce(other).terms.items():
            result[exponent] = result.get(exponent, 0) + coefficient
            if result[exponent] == 0:
                del result[exponent]
        return Polynomial(result)

    __radd__ = __add__

    def __neg__(self) -> Polynomial:
        return Polynomial(
            {exponent: -coefficient for exponent, coefficient in self.terms.items()}
        )

    def __sub__(self, other: int | Polynomial) -> Polynomial:
        return self + (-self.coerce(other))

    def __rsub__(self, other: int | Polynomial) -> Polynomial:
        return self.coerce(other) - self

    def __mul__(self, other: int | Polynomial) -> Polynomial:
        right = self.coerce(other)
        result: dict[Exponent, int] = {}
        for left_exponent, left_coefficient in self.terms.items():
            for right_exponent, right_coefficient in right.terms.items():
                exponent = tuple(
                    left_power + right_power
                    for left_power, right_power in zip(
                        left_exponent,
                        right_exponent,
                        strict=True,
                    )
                )
                result[exponent] = (
                    result.get(exponent, 0)
                    + left_coefficient * right_coefficient
                )
        return Polynomial(result)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> Polynomial:
        if exponent < 0:
            raise ValueError("多项式指数不能为负")
        result = Polynomial(1)
        factor = self
        power = exponent
        while power:
            if power % 2:
                result = result * factor
            factor = factor * factor
            power //= 2
        return result

    @property
    def is_zero(self) -> bool:
        return not self.terms

    @property
    def total_degree(self) -> int:
        return max((sum(exponent) for exponent in self.terms), default=-1)


Vector = tuple[Polynomial, Polynomial, Polynomial]
Parameter = tuple[Polynomial, Polynomial]


def variable(index: int) -> Polynomial:
    exponent = [0, 0, 0, 0, 0]
    exponent[index] = 1
    return Polynomial({tuple(exponent): 1})  # type: ignore[arg-type]


c, d, h, p, q = (variable(index) for index in range(5))


def cross(left: Vector, right: Vector) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def parameter(kind: str, source: Polynomial) -> Parameter:
    if kind == "A":
        return (source + 1, h)
    if kind == "B":
        return (h, 1 - source)
    raise ValueError(f"未知割线端点 {kind}")


def conic_point(value: Parameter) -> Vector:
    first, second = value
    return (
        first**2 - second**2,
        2 * first * second,
        first**2 + second**2,
    )


def determinant(left: Parameter, right: Parameter) -> Polynomial:
    return left[0] * right[1] - right[0] * left[1]


PROFILES = {
    "P0": {
        "sigma": (1, 1, 1),
        "epsilon": 1,
        "roles": (
            ("B", c + p),
            ("A", c - p),
            ("A", c + d - q),
            ("B", c + d + q),
        ),
    },
    "P1": {
        "sigma": (1, 1, -1),
        "epsilon": -1,
        "roles": (
            ("A", c + p),
            ("B", c - p),
            ("B", c + d - q),
            ("A", c + d + q),
        ),
    },
    "P2": {
        "sigma": (1, -1, -1),
        "epsilon": -1,
        "roles": (
            ("A", c + p),
            ("B", c - p),
            ("B", c + d + q),
            ("A", c + d - q),
        ),
    },
    "P3": {
        "sigma": (1, -1, 1),
        "epsilon": 1,
        "roles": (
            ("B", c + p),
            ("A", c - p),
            ("A", c + d + q),
            ("B", c + d - q),
        ),
    },
}


def contact_line(sigma2: int, sigma3: int) -> Vector:
    """返回清除分母后的 ``Line(B, sigma3*U)``。"""

    numerator_bx = 2 * c * d + d**2 + p**2 - q**2
    numerator_by = (
        c**2 * d
        + c * d**2
        + c * p**2
        - c * q**2
        - d * h**2
        + d * p**2
        - d
    )
    numerator_uy = c * p - sigma2 * c * q + d * p - sigma3 * d
    point_b = (
        h * numerator_bx,
        -numerator_by,
        2 * d * h,
    )
    point_sigma3_u = (
        sigma3 * h * (p - sigma2 * q),
        -sigma3 * numerator_uy,
        d * h,
    )
    return cross(point_b, point_sigma3_u)


def verify_profile(profile_id: str) -> dict[str, int | str]:
    profile = PROFILES[profile_id]
    sigma1, sigma2, sigma3 = profile["sigma"]
    if sigma1 != 1:
        raise AssertionError("代表元必须满足 sigma1=1")

    role_parameters = tuple(
        parameter(kind, source) for kind, source in profile["roles"]
    )
    x_parameter, y_parameter, z_parameter, w_parameter = role_parameters
    x, y, z, w = tuple(conic_point(value) for value in role_parameters)

    chord_xz = cross(x, z)
    chord_yw = cross(y, w)
    diagonal_point = cross(chord_xz, chord_yw)
    tau = (
        diagonal_point[0],
        diagonal_point[1],
        -diagonal_point[2],
    )
    lambda_sigma = contact_line(sigma2, sigma3)
    determinant_xz = determinant(x_parameter, z_parameter)
    determinant_yw = determinant(y_parameter, w_parameter)

    proportionality = (
        profile["epsilon"]
        * 8
        * determinant_xz
        * determinant_yw
    )
    for tau_coefficient, contact_coefficient in zip(
        tau,
        lambda_sigma,
        strict=True,
    ):
        if d * tau_coefficient != proportionality * contact_coefficient:
            raise AssertionError(f"{profile_id} 的比例恒等式失败")

    if not all(value.is_zero for value in cross(tau, lambda_sigma)):
        raise AssertionError(f"{profile_id} 的两条齐次直线不恒等")
    if determinant_xz.is_zero or determinant_yw.is_zero:
        raise AssertionError(f"{profile_id} 的合并因子不应恒为零")

    all_polynomials = (
        *tau,
        *lambda_sigma,
        determinant_xz,
        determinant_yw,
    )
    return {
        "profile": profile_id,
        "sigma": "".join("+" if sign > 0 else "-" for sign in profile["sigma"]),
        "epsilon": profile["epsilon"],
        "maximum_degree": max(value.total_degree for value in all_polynomials),
        "maximum_terms": max(len(value.terms) for value in all_polynomials),
    }


def main() -> None:
    reports = tuple(verify_profile(profile_id) for profile_id in PROFILES)
    print("mannheim_contact_chord_identity", reports)


if __name__ == "__main__":
    main()
