"""符号验证 Mannheim 三重外切根的两标号覆盖。

固定 ``Gamma1``，比较以下两个有序 Mannheim 程序：

* 以 ``Gamma1,Gamma2`` 为前两圆、``Gamma3`` 为第三圆；
* 以 ``Gamma1,Gamma3`` 为前两圆、``Gamma2`` 为第三圆。

当前 19 E 单根分支只会在 ``P0`` 的 ``y=w`` 简单合并中出现。本脚本
在整数系数稀疏多项式环中验证：上述两个方向不可能同时满足 ``y=w``。
证明只使用半径为正、前两条圆心边严格长于对应半径和，以及三个圆心
不共线；不使用浮点数或外部计算机代数系统。
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from replay_mannheim_ordered_branches import OrderedBranchReplay


@dataclass(frozen=True, slots=True)
class Polynomial:
    """任意固定变量数的稀疏整数系数多项式。"""

    terms: dict[tuple[int, ...], int]
    variable_count: int

    def __init__(
        self,
        terms: int | dict[tuple[int, ...], int],
        variable_count: int,
    ) -> None:
        if isinstance(terms, int):
            values = (
                {}
                if terms == 0
                else {(0,) * variable_count: terms}
            )
        else:
            values = {
                exponent: coefficient
                for exponent, coefficient in terms.items()
                if coefficient != 0
            }
            if any(len(exponent) != variable_count for exponent in values):
                raise ValueError("多项式指数维数错误")
        object.__setattr__(self, "terms", values)
        object.__setattr__(self, "variable_count", variable_count)

    @classmethod
    def constant(cls, value: int, variable_count: int) -> Polynomial:
        return cls(value, variable_count)

    @classmethod
    def variable(cls, index: int, variable_count: int) -> Polynomial:
        exponent = [0] * variable_count
        exponent[index] = 1
        return cls({tuple(exponent): 1}, variable_count)

    def coerce(self, value: int | Polynomial) -> Polynomial:
        if isinstance(value, Polynomial):
            if value.variable_count != self.variable_count:
                raise ValueError("不能混合不同变量数的多项式")
            return value
        return Polynomial.constant(value, self.variable_count)

    def __add__(self, other: int | Polynomial) -> Polynomial:
        result = dict(self.terms)
        for exponent, coefficient in self.coerce(other).terms.items():
            result[exponent] = result.get(exponent, 0) + coefficient
            if result[exponent] == 0:
                del result[exponent]
        return Polynomial(result, self.variable_count)

    __radd__ = __add__

    def __neg__(self) -> Polynomial:
        return Polynomial(
            {
                exponent: -coefficient
                for exponent, coefficient in self.terms.items()
            },
            self.variable_count,
        )

    def __sub__(self, other: int | Polynomial) -> Polynomial:
        return self + (-self.coerce(other))

    def __rsub__(self, other: int | Polynomial) -> Polynomial:
        return self.coerce(other) - self

    def __mul__(self, other: int | Polynomial) -> Polynomial:
        right = self.coerce(other)
        result: dict[tuple[int, ...], int] = {}
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
        return Polynomial(result, self.variable_count)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> Polynomial:
        if exponent < 0:
            raise ValueError("多项式指数不能为负")
        result = Polynomial.constant(1, self.variable_count)
        factor = self
        power = exponent
        while power:
            if power % 2:
                result *= factor
            factor *= factor
            power //= 2
        return result

    @property
    def is_zero(self) -> bool:
        return not self.terms


def main() -> None:
    # x=d12, y=d13, a=r1, b=r2, c=r3, k=cos(angle 213).
    x, y, a, b, c, k = (
        Polynomial.variable(index, 6) for index in range(6)
    )

    # P0:y=w 的圆参数行列式。把第三圆到公共方向的投影分别写成
    # y*k 和 x*k 后，两个方向的合并条件化为 residual_3=0 与
    # residual_2=0。二次 k 项在减去高的平方时完全抵消。
    raw_3 = (y * k + a - c) * (x - y * k + b - c) - y**2 * (
        1 - k**2
    )
    residual_3 = y * k * (x + b - a) + (a - c) * (x + b - c) - y**2
    raw_2 = (x * k + a - b) * (y - x * k + c - b) - x**2 * (
        1 - k**2
    )
    residual_2 = x * k * (y + c - a) + (a - b) * (y + c - b) - x**2
    if raw_3 != residual_3 or raw_2 != residual_2:
        raise AssertionError("P0:y=w 的线性余式展开错误")

    numerator_3 = y**2 - (a - c) * (x + b - c)
    numerator_2 = x**2 - (a - b) * (y + c - b)
    eliminated = (
        numerator_3 * x * (y + c - a)
        - numerator_2 * y * (x + b - a)
    )
    positive_factor = (
        -a**2 * x
        - a**2 * y
        + 2 * a * b * y
        + 2 * a * c * x
        - b**2 * y
        - c**2 * x
        + x**2 * y
        + x * y**2
    )
    linear_factor = x - y + b - c
    if eliminated != -linear_factor * positive_factor:
        raise AssertionError("两个合并条件的消元因式分解错误")

    # 外离条件写成 x=a+b+X、y=a+c+Y，X,Y>0。代入后
    # positive_factor 的每个系数都严格为正，故该因子不可能为零。
    X, Y, ra, rb, rc, unused = (
        Polynomial.variable(index, 6) for index in range(6)
    )
    side_x = ra + rb + X
    side_y = ra + rc + Y
    positive_on_domain = (
        -ra**2 * side_x
        - ra**2 * side_y
        + 2 * ra * rb * side_y
        + 2 * ra * rc * side_x
        - rb**2 * side_y
        - rc**2 * side_x
        + side_x**2 * side_y
        + side_x * side_y**2
    )
    if not positive_on_domain.terms or any(
        coefficient <= 0
        for coefficient in positive_on_domain.terms.values()
    ):
        raise AssertionError("外离域上的正因子证书含非正系数")
    if any(exponent[5] != 0 for exponent in positive_on_domain.terms):
        raise AssertionError("正因子证书意外使用了占位变量")

    # 消元迫使 y=x+b-c。代回任一合并余式，只剩
    # y*(x+b-a)*(k-1)。外离使前两个因子为正，而非共线使 k<1。
    forced_y = x + b - c
    forced_residual = (
        forced_y * k * (x + b - a)
        + (a - c) * (x + b - c)
        - forced_y**2
    )
    if forced_residual != forced_y * (x + b - a) * (k - 1):
        raise AssertionError("消元后的非共线矛盾因式分解错误")

    # 依赖图校准两个会决定最坏分数的 P0 分支。y=w 时 +++ 是另一根，
    # 恰有 19 个祖先；x=z 时 +++ 本身就是合并接触点根，只需 9 个祖先。
    F = Fraction
    merge_fixtures = {
        "xz": (
            ((F(0), F(0)), (F(7), F(0)), (F(4), F(2))),
            (F(3), F(2), F(1)),
        ),
        "yw": (
            ((F(0), F(0)), (F(7), F(0)), (F(3), F(5))),
            (F(3), F(2), F(1)),
        ),
    }
    replay_costs = {}
    for merge_kind, (centers, radii) in merge_fixtures.items():
        report = OrderedBranchReplay(
            f"ext_relabel_{merge_kind}",
            centers,
            radii,
            emit=False,
        ).run()
        replay_costs[merge_kind] = report["first_ext"]
    if replay_costs != {"xz": 9, "yw": 19}:
        raise AssertionError("P0 简单合并的 +++ 祖先成本发生变化")

    print(
        "mannheim_ext_relabel",
        {
            "candidate_labelings": 2,
            "merge_conditions": 2,
            "positive_factor_terms": len(positive_on_domain.terms),
            "simultaneous_bad_labelings": 0,
            "simple_merge_ext_costs": replay_costs,
            "full_domain_ext_upper": 18,
        },
    )


if __name__ == "__main__":
    main()
