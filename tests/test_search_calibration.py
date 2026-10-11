"""已证明的几何实例和穷举阳性对照，验证新的搜索诊断。"""

from __future__ import annotations

import sys
import unittest
from time import monotonic
from math import sqrt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

from audit_mannheim_search_bindings import independent_root_control, known_root_control
from mannheim_search_calibration import apply_calibrated_candidate
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from search_mannheim_center_locus_2e import object_value
try:
    import numpy
except ImportError:
    numpy = None

if numpy is not None:
    from search_mannheim_joint_root_kp import build_joint_state, fast_target_candidates, third_target_candidates
    from search_mannheim_shared_circle_replacement import available, dual_circle_candidates, two_step_root_paths
    from mannheim_virtual_witnesses import witness_helpers
from search_mannheim_root_center_2e import FIXTURES, drawable_contains_targets
from search_parallel_3e import Candidate, DrawableBundle, PointBundle, State, drawable_signature, generate_candidates, line_through, intersections
from mannheim_search_calibration import BindingAudit, root_bundles


class SearchCalibrationTests(unittest.TestCase):
    def test_known_intersection_identity_survives_rotation(self):
        result = known_root_control()
        self.assertEqual(result["legacy_hits"], 0)
        self.assertEqual(result["known_first_hits"], 1)

    def test_diagnostic_root_pairings_and_missing_real_branches(self):
        result = independent_root_control()
        self.assertEqual(result["diagnostic_bundles"], 8)
        self.assertEqual(result["mismatched_real_roots_recorded"], 1)

    @unittest.skipIf(numpy is None, "联合搜索的 NumPy 剪枝需要 numpy")
    def test_goal_pruning_recovers_seed_tau_and_matches_unpruned_candidates(self):
        state, targets, kps, _ = build_joint_state()
        self.assertEqual(len(state.program), 17)
        self.assertEqual(len(kps), 4)
        replays = []
        for centers, radii in FIXTURES:
            replay = KpCenterLocusReplay("positive_control", centers=centers, radii=radii, emit=False)
            replay.run()
            replays.append(replay)
        count = len(state.points)
        for move, node in enumerate(("P0_K_core_line_1", "P0_K_core_line_2"), 18):
            values = tuple(object_value(replay.objects, node)[1] for replay in replays)
            candidate = Candidate("line", "positive_control_a", "positive_control_b",
                                  DrawableBundle(node, "line", values))
            state = apply_calibrated_candidate(state, candidate, move)
        required = frozenset(point.point_id for point in state.points[count:])
        fast = fast_target_candidates(state, targets, required, 20)
        exhaustive = tuple(candidate for candidate in generate_candidates(state, 20)
                           if (candidate.first_id in required or candidate.second_id in required)
                           and drawable_contains_targets(candidate.drawable, targets))
        fast_signatures = {drawable_signature(candidate.drawable) for candidate in fast}
        self.assertEqual(fast_signatures, {drawable_signature(candidate.drawable) for candidate in exhaustive})
        expected = DrawableBundle("P0_tau", "line", tuple(object_value(replay.objects, "P0_tau")[1]
                                  for replay in replays))
        self.assertIn(drawable_signature(expected), fast_signatures)

    @unittest.skipIf(numpy is None, "联合搜索的 NumPy 剪枝需要 numpy")
    def test_third_object_can_use_only_first_helper_points(self):
        # 第一辅助对象产生 B、D 后，可以先画 AB 再画 CD；第二步无需产生新点。
        bundle = lambda name, point: PointBundle(name, (point,) * 3)
        a, c = bundle("A", (-1.0, 0.0)), bundle("C", (0.0, -1.0))
        b, d = bundle("B", (1.0, 0.0)), bundle("D", (0.0, 1.0))
        base = State((a, c), ())
        horizontal = DrawableBundle("second", "line", ((0.0, 1.0, 0.0),) * 3)
        two = State((a, c, b, d), (horizontal,))
        candidates = third_target_candidates(base, two, ((0.0, 0.0),) * 3, 3)
        vertical = DrawableBundle("expected", "line", ((1.0, 0.0, 0.0),) * 3)
        self.assertIn(drawable_signature(vertical), {drawable_signature(item.drawable) for item in candidates})

    @unittest.skipIf(numpy is None, "联合搜索需要 numpy")
    def test_shared_circle_restores_two_geometric_intersections_under_rotation(self):
        # 合法 Circle(O,P)，与两条已知切线各有一个有限交点。
        transforms = (lambda p: p, lambda p: (-p[0], -p[1]), lambda p: (-p[1], p[0]))
        values = lambda point: tuple(transform(point) for transform in transforms)
        origin, through = PointBundle("O", values((0.0, 0.0))), PointBundle("P", values((0.0, 1.0)))
        supports = tuple(DrawableBundle(str(x), "line", tuple(line_through(transform((x, 0.0)), transform((x, 1.0)))
                         for transform in transforms)) for x in (-1.0, 1.0))
        state = State((origin, through), supports)
        left, right = values((-1.0, 0.0)), values((1.0, 0.0))
        self.assertFalse(available(state, left))
        self.assertFalse(available(state, right))
        circles = dual_circle_candidates(state, left, right, 1)
        self.assertEqual(len(circles), 1)
        constructed = apply_calibrated_candidate(state, circles[0], 1)
        self.assertTrue(available(constructed, left))
        self.assertTrue(available(constructed, right))

    @unittest.skipIf(numpy is None, "联合搜索需要 numpy")
    def test_two_step_root_completion_finds_intersection_without_injecting_it(self):
        points = tuple(PointBundle(name, (value,) * 3) for name, value in (
            ("A", (1.0, 0.0)), ("B", (1.0, 2.0)), ("C", (0.0, 1.0)), ("D", (2.0, 1.0))))
        state = State(points, ())
        targets = ((1.0, 1.0),) * 3
        self.assertFalse(available(state, targets))
        paths = two_step_root_paths(state, targets, 1, BindingAudit())
        self.assertIn(("Line(A, B)", "Line(C, D)"), paths)
        self.assertTrue(all(len(path) == 2 for path in paths))

    @unittest.skipIf(numpy is None, "目标见证剪枝需要 numpy")
    def test_virtual_witness_pruning_matches_brute_force_under_rotation(self):
        transforms = (lambda p: p, lambda p: (-p[0], -p[1]), lambda p: (-p[1], p[0]))
        values = lambda point: tuple(transform(point) for transform in transforms)
        points = tuple(PointBundle(name, values(point)) for name, point in (
            ("O", (0.0, 0.0)), ("R", (2.0, 0.0)), ("A", (0.0, 1.0)), ("B", (3.0, 1.0))))
        old = DrawableBundle("input_circle", "circle", ((0.0, 0.0, 4.0),) * 3)
        state = State(points, (old,))
        goal = DrawableBundle("virtual", "line", tuple(line_through(transform((sqrt(3), 0.0)),
                                 transform((sqrt(3), 1.0))) for transform in transforms))
        helpers, report = witness_helpers(state, [goal], 1, monotonic() + 10)
        expected = set()
        for candidate in generate_candidates(state, 1):
            rows = tuple(intersections(candidate.kind, p, old.kind, q)
                         for p, q in zip(candidate.drawable.values, old.values, strict=True))
            for _, roots in root_bundles(rows, state.points, policy="all"):
                if drawable_contains_targets(goal, roots):
                    expected.add(drawable_signature(candidate.drawable))
        self.assertEqual({drawable_signature(item.drawable) for item in helpers}, expected)
        self.assertIn("Line(A, B)", [item.describe() for item in helpers])
        self.assertTrue(report["complete"])
        self.assertTrue(all(item.first_id in {p.point_id for p in points} and
                            item.second_id in {p.point_id for p in points} for item in helpers))

    @unittest.skipIf(numpy is None, "目标见证剪枝需要 numpy")
    def test_coincident_virtual_object_keeps_all_first_helpers(self):
        points = tuple(PointBundle(name, (point,) * 3) for name, point in (
            ("O", (0.0, 0.0)), ("R", (2.0, 0.0)), ("A", (0.0, 1.0)), ("B", (3.0, 1.0))))
        old = DrawableBundle("input_circle", "circle", ((0.0, 0.0, 4.0),) * 3)
        state = State(points, (old,))
        helpers, report = witness_helpers(state, [old], 1, monotonic() + 10)
        self.assertEqual({drawable_signature(item.drawable) for item in helpers},
                         {drawable_signature(item.drawable) for item in generate_candidates(state, 1)})
        self.assertEqual(report["coincident_goal_old_pairs"], 1)


if __name__ == "__main__":
    unittest.main()
