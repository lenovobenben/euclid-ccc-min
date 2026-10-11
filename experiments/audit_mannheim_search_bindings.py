"""校验旧排序绑定的漏检实例，以及现有圆心线关系的阳性对照。"""

from __future__ import annotations

import json
from dataclasses import asdict
from fractions import Fraction
from itertools import permutations
from math import sqrt

from mannheim_search_calibration import BindingAudit, root_bundles
from search_parallel_3e import PointBundle, intersections
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from search_mannheim_center_locus_2e import object_value, point_value
from search_mannheim_root_center_2e import FIXTURES
from search_parallel_3e import line_through


def known_root_control():
    # 同一个合法平行前缀，三个坐标表示：原样、旋转 180°、旋转 90°。
    transforms = (lambda p: p, lambda p: (-p[0], -p[1]), lambda p: (-p[1], p[0]))
    rows, known, expected = [], [], []
    for transform in transforms:
        o1, o3 = transform((0.0, 0.0)), transform((4.0, 15.0))
        x, q = transform((sqrt(241.0), 0.0)), transform((4.0, -15.0))
        radius_squared = (x[0] - o3[0]) ** 2 + (x[1] - o3[1]) ** 2
        rows.append(intersections("circle", (*o1, 241.0), "circle", (*x, radius_squared)))
        known.append(o3)
        expected.append(q)
    known_point = PointBundle("O3", tuple(known))
    old_bundles = tuple(tuple(row[i] for row in rows) for i in range(2))
    matches = lambda values: all(sum((a - b) ** 2 for a, b in zip(p, q)) < 1e-16
                                  for p, q in zip(values, expected))
    audit = BindingAudit()
    corrected = root_bundles(tuple(rows), (known_point,), audit=audit)
    old_hits = sum(matches(values) for values in old_bundles)
    new_hits = sum(matches(values) for _, values in corrected)
    assert old_hits == 0 and new_hits == 1
    return {"legacy_hits": old_hits, "known_first_hits": new_hits, "audit": asdict(audit)}


def independent_root_control():
    rows = (((-1.0, 0.0), (1.0, 0.0)), ((-1.0, 0.0), (1.0, 0.0)),
            ((-1.0, 0.0), (1.0, 0.0)))
    bundles = root_bundles(rows, policy="all")
    assert len(bundles) == 8
    assert (0, 1, 0) in {indices for indices, _ in bundles}
    audit = BindingAudit()
    assert not root_bundles((rows[0], (), rows[2]), policy="all", audit=audit)
    assert audit.root_count_mismatches == 1
    return {"diagnostic_bundles": len(bundles), "mismatched_real_roots_recorded": 1}


def center_locus_controls():
    rows = []
    for fixture_index, (centers, radii) in enumerate(FIXTURES):
        replay = KpCenterLocusReplay("binding_control", centers=centers, radii=radii, emit=False)
        report = replay.run()
        assert report["all_targets"] == 49
        for profile in ("P0", "P2", "P1", "P3"):
            kp = replay.objects.graph.dependencies[f"{profile}_O3_Kp"][1]
            first = line_through(point_value(replay.objects, "O3"), point_value(replay.objects, kp))
            _, ell = object_value(replay.objects, "ell")
            _, locus = object_value(replay.objects, f"{profile}_center_locus_from_Kp")
            j, = intersections("line", first, "line", ell)
            residual = abs(locus[0] * j[0] + locus[1] * j[1] + locus[2])
            assert residual < 1e-10
            rows.append({"fixture": fixture_index, "profile": profile, "residual": residual})
    return {"known_positive_controls": len(rows), "maximum_residual": max(row["residual"] for row in rows)}


def kp_concyclicity_controls():
    """以有理点的 4×4 行列式否定四个 K' 恒共圆这一候选关系。"""
    rows = []
    for fixture_index, (centers, radii) in enumerate(FIXTURES):
        replay = KpCenterLocusReplay("concyclicity_control", centers=centers, radii=radii, emit=False)
        replay.run()
        matrix = []
        for profile in ("P0", "P2", "P1", "P3"):
            node = replay.objects.resolve(replay.objects.graph.dependencies[f"{profile}_O3_Kp"][1])
            key = next(key for key, value in replay.objects.point_registry.items() if value == node)
            assert all(coordinate[0] == "rational" for coordinate in key)
            x, y = (coordinate[1] for coordinate in key)
            matrix.append((x * x + y * y, x, y, Fraction(1)))
        determinant = Fraction(0)
        for columns in permutations(range(4)):
            inversions = sum(columns[i] > columns[j] for i in range(4) for j in range(i + 1, 4))
            term = Fraction((-1) ** inversions)
            for i, j in enumerate(columns):
                term *= matrix[i][j]
            determinant += term
        assert determinant != 0
        rows.append({"fixture": fixture_index, "determinant": str(determinant)})
    return rows


def main():
    print(json.dumps({"known_root": known_root_control(), "independent_roots": independent_root_control(),
                      "center_locus": center_locus_controls(),
                      "kp_concyclicity": kp_concyclicity_controls()}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
