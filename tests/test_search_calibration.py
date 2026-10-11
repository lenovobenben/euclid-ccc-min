"""已证明的几何实例和穷举阳性对照，验证新的搜索诊断。"""

from __future__ import annotations

import sys
import unittest
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
from search_mannheim_root_center_2e import FIXTURES, drawable_contains_targets
from search_parallel_3e import Candidate, DrawableBundle, PointBundle, State, drawable_signature, generate_candidates


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


if __name__ == "__main__":
    unittest.main()
