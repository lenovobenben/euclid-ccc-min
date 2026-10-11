"""Gergonne 极点接口：符号等价检查与 19 E 条件极线模块的精确重放。

模块输入含已构造的四个相似中心、第三圆和两组直径端点。
此条件成本不是从 CCC 初态开始的成绩；坐标公式仅用于独立核对。
只需标准库；不采用数值容差或任意点。
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from check_mannheim_degenerate_fixture import (
    add, dot, line_intersection, line_through, multiply, same_line,
    second_circle_point, subtract,
)
from replay_mannheim_fixed import fraction_square_root
from verify_mannheim_contact_chord_identity import (
    PROFILES, c, contact_line, cross, d, h, p, q, verify_profile,
)
from verify_mannheim_center_locus_kp_identity import center_direction
from verify_mannheim_tau_concurrence import MANNHEIM_POINT_S


F = Fraction


def similarity_center(centers, radii, index, sign):
    # 圆 i 的有向半径为 r_i，第三圆为 sign*r_3。
    denominator = radii[index] - sign * radii[2]
    return multiply(1 / denominator, subtract(
        multiply(radii[index], centers[2]), multiply(sign * radii[2], centers[index])))


def polar(point, center, radius):
    dx, dy = subtract(point, center)
    return dx, dy, -dx * center[0] - dy * center[1] - radius**2


def verify_symbolic_interface():
    reports, poles = [], {}
    for profile, data in PROFILES.items():
        _, sigma2, sigma3 = data["sigma"]
        first = (-sigma3 * c, -sigma3 * h, p - sigma3)
        second = (-sigma3 * (c + d), -sigma3 * h, sigma2 * q - sigma3)
        axis = cross(first, second)
        pole = (-axis[0], -axis[1], axis[2])
        expected = (sigma3 * h * (p - sigma2 * q),
                    -sigma3 * (c * p - sigma2 * c * q + d * p - sigma3 * d),
                    d * h)
        if not all(value.is_zero for value in cross(pole, expected)):
            raise AssertionError(f"{profile} 的相似轴极点公式失败")
        if axis[2] != -d * h:
            raise AssertionError("相似轴不能经过第三圆心")
        chord = cross(MANNHEIM_POINT_S, pole)
        if not all(value.is_zero for value in cross(chord, contact_line(sigma2, sigma3))):
            raise AssertionError(f"{profile} 的 Gergonne 弦与 Mannheim 弦不等价")
        direction = center_direction(profile)
        if not (axis[0] * direction[1] - axis[1] * direction[0]).is_zero:
            raise AssertionError("相似轴法向量与目标圆心线方向不平行")
        verify_profile(profile)  # 同时运行原 Mannheim 弦恒等式。
        poles[(sigma2, sigma3)] = expected
        reports.append({"profile": profile, "pole_formula": True,
                        "same_contact_chord": True, "same_center_direction": True})
    for sigma2 in (-1, 1):
        one, two = poles[(sigma2, 1)], poles[(sigma2, -1)]
        if one[0] + two[0] != 0 * c or one[1] + two[1] != 2 * d:
            raise AssertionError("两对相似轴极点没有共同中点")
    return {"ring": "Z[c,d,h,p,q]", "profiles": reports,
            "pole_parallelogram_center": "(0, 1/h)", "score_certificate": False}


def replay_polar_gadget(centers, radii):
    """条件输入的合法 Line 重放；新增点逐个验证有限交点身份。"""
    origin, radius = centers[2], radii[2]
    points = {}
    lines, program = {}, []
    for index, prefix in ((0, "A"), (1, "C")):
        direction = subtract(centers[index], origin)
        length = fraction_square_root(dot(direction, direction))
        for sign, suffix in ((1, "+"), (-1, "-")):
            points[prefix + suffix] = add(origin, multiply(sign * radius / length, direction))
            points[f"H{index + 1}{suffix}"] = similarity_center(centers, radii, index, sign)
        lines[f"diameter_{prefix}"] = line_through(origin, centers[index])
    inputs = tuple(sorted(points))

    def draw_line(name, first, second):
        value = line_through(points[first], points[second])
        for previous, old in lines.items():
            if same_line(value, old):
                return previous
        lines[name] = value
        program.append({"op": "Line", "id": name, "points": [first, second]})
        return name

    def intersect_lines(name, first, second):
        value = line_intersection(lines[first], lines[second])
        if any(line[0] * value[0] + line[1] * value[1] + line[2] != 0
               for line in (lines[first], lines[second])):
            raise AssertionError("新点不在其两条父直线上")
        points[name] = value
        program.append({"op": "Intersect", "id": name, "objects": [first, second]})

    polars = {}
    for index, prefix, other in ((0, "A", "C+"), (1, "C", "A+")):
        # 两组共用 A+C+；四条基础交叉弦实际只绘制三条。
        shared_plus = draw_line(f"shared_{prefix}+", prefix + "+", other)
        shared_minus = draw_line(f"shared_{prefix}-", prefix + "-", other)
        for sign, suffix in ((1, "+"), (-1, "-")):
            key = f"H{index + 1}{suffix}"
            secant = draw_line(f"secant_{key}", key, other)
            root = second_circle_point(points[key], points[other], (origin, radius))
            line = lines[secant]
            if line[0] * root[0] + line[1] * root[1] + line[2] != 0:
                raise AssertionError("第二圆交点不在割线上")
            if dot(subtract(root, origin), subtract(root, origin)) != radius**2:
                raise AssertionError("第二圆交点不在输入圆上")
            if root == points[other]:
                raise AssertionError("条件模块排除割线相切")
            points[f"D_{key}"] = root
            program.append({"op": "IntersectOther", "id": f"D_{key}",
                            "objects": [secant, "Gamma3"], "known": other})
            cross_plus = draw_line(f"cross_plus_{key}", prefix + "+", f"D_{key}")
            cross_minus = draw_line(f"cross_minus_{key}", prefix + "-", f"D_{key}")
            intersect_lines(f"X_{key}", shared_plus, cross_minus)
            intersect_lines(f"Y_{key}", shared_minus, cross_plus)
            result = draw_line(f"polar_{key}", f"X_{key}", f"Y_{key}")
            if not same_line(lines[result], polar(points[key], origin, radius)):
                raise AssertionError("完全四边形程序没有给出目标极线")
            polars[(index, sign)] = result
    pole_results = []
    for profile, data in PROFILES.items():
        _, sigma2, sigma3 = data["sigma"]
        # H13 使用 sign=sigma3，H23 使用 sign=sigma3/sigma2。
        left, right = polars[(0, sigma3)], polars[(1, sigma3 * sigma2)]
        intersect_lines(f"pole_{profile}", left, right)
        axis = line_through(points[f"H1{'+' if sigma3 > 0 else '-'}"],
                            points[f"H2{'+' if sigma3 * sigma2 > 0 else '-'}"])
        pole = points[f"pole_{profile}"]
        # 检查点的极线就是轴；不把解析公式作为新点来源。
        if not same_line(polar(pole, origin, radius), axis):
            raise AssertionError("两条相似中心极线的交点不是相似轴的极点")
        pole_results.append({"profile": profile, "point": [str(value) for value in pole]})
    paid = sum(step["op"] == "Line" for step in program)
    if paid > 19:
        raise AssertionError("条件极线模块超过 19 E")
    return {"centers": [[str(value) for value in point] for point in centers],
            "radii": [str(value) for value in radii], "inputs": inputs,
            "input_objects": ["Gamma3", "diameter_A", "diameter_C"],
            "conditional_e": paid, "full_ccc_score": False,
            "poles": pole_results, "program": program}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    centers = ((F(0), F(0)), (F(15), F(0)), (F(0), F(20)))
    radii = (F(4), F(2), F(1))
    transforms = (lambda point: point,
                  lambda point: (7 - point[1], point[0] - 3),
                  lambda point: (2 * point[0] + 1, 2 * point[1] + 5))
    fixtures = [replay_polar_gadget(tuple(transform(point) for point in centers),
                                  tuple(scale * radius for radius in radii))
                for transform, scale in zip(transforms, (1, 1, 2), strict=True)]
    report = {"schema": "gergonne-polar-interface/v1",
              "symbolic": verify_symbolic_interface(), "gadget_fixtures": fixtures}
    Path(args.output).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"symbolic": report["symbolic"],
                      "conditional_e": [item["conditional_e"] for item in fixtures]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
