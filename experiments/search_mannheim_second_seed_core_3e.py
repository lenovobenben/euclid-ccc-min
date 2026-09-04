"""筛查第二条 Mannheim 种子接触弦的 3 E 局部替换。

正规程序先用 5 E 构造第一条种子接触弦；相邻方向类共享一条 ``K'``
定义弦，第二条种子接触弦还需三条定位线和最终接触弦，共 4 E。

本脚本对四种有共享弦的有序种子对，枚举第二块三条标准定位线中的任意
两条，再枚举第三步所有可画直线或圆。若第三步对象经过第二方向类在
第三输入圆上的两个接触点，它便能直接绑定这对点，使两种子核心从
``5 + 4 = 9 E`` 降为 ``5 + 3 = 8 E``。

搜索同步使用三个严格正规 ``D8`` 夹具，只使用已有对象的确定有限实
交点。命中必须另作精确重放和符号证明；零结果只覆盖“两条标准定位线
+ 任意接触点载体”的 3 E 形状，不是 3 E 下界。
"""

from __future__ import annotations

from itertools import combinations

from scan_mannheim_degeneracies import is_d8
from replay_mannheim_center_locus_dependencies import build_seed_tau
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from replay_mannheim_ordered_branches import (
    OTHER_BATCH,
    P0_BATCH,
    SEED_PROFILES,
    TAIL_PROFILES,
    branch_data,
)
from replay_mannheim_three_block_dependencies import collapse_point
from search_mannheim_center_locus_2e import object_value
from search_mannheim_double_kp_global_sequential_2e import float_scalar
from search_mannheim_kp_parallel_center_locus_2e import (
    candidates,
    on_drawable,
    same_drawable,
)
from search_mannheim_root_center_2e import FIXTURES, add_known_drawable
from search_parallel_3e import DrawableBundle, PointBundle, State, apply_candidate


PROFILES = SEED_PROFILES + TAIL_PROFILES
SHARED_SEED_PAIRS = (
    ("P0", "P2"),
    ("P2", "P0"),
    ("P0", "P3"),
    ("P3", "P0"),
    ("P1", "P2"),
    ("P2", "P1"),
    ("P1", "P3"),
    ("P3", "P1"),
)


class TwoSeedReplay(KpCenterLocusReplay):
    """构造两个正规种子块，并记录第二块新增的四个对象。"""

    seed_profiles = ("P0", "P2")

    def run(self):
        if not is_d8(self.centers, self.radii):
            raise AssertionError(f"{self.fixture_name} 不属于严格 D8")
        data = {profile: branch_data(self, profile) for profile in PROFILES}
        if {row["kind"] for row in data.values()} != {"regular"}:
            raise AssertionError("第二种子搜索只接受四方向均正规的夹具")

        self.build_prefix()
        for key in P0_BATCH + OTHER_BATCH:
            self.draw_batch(key)
        if len(self.objects.graph.paid_order) != 13:
            raise AssertionError("正规公共前缀必须恰有 13 E")

        first_profile, second_profile = self.seed_profiles
        _, targets_completed = build_seed_tau(self, first_profile, data[first_profile])
        if targets_completed or len(self.objects.graph.paid_order) != 18:
            raise AssertionError("第一种子接触弦必须恰好增加 5 E")
        self.base_paid_order = tuple(self.objects.graph.paid_order)

        (second_tau, _), targets_completed = build_seed_tau(
            self,
            second_profile,
            data[second_profile],
        )
        if targets_completed:
            raise AssertionError("正规第二种子块不应提前完成目标")
        self.second_paid_ids = tuple(
            self.objects.graph.paid_order[len(self.base_paid_order) :]
        )
        if len(self.second_paid_ids) != 4:
            raise AssertionError("共享定义弦的第二种子块必须恰好增加 4 E")
        if self.second_paid_ids[-1] != self.objects.resolve(
            f"{second_profile}_tau"
        ):
            raise AssertionError("第二种子块的最后对象不是接触弦")
        self.second_locator_ids = self.second_paid_ids[:-1]
        self.second_tau = second_tau
        self.core_data = data
        return {
            "base_paid": len(self.base_paid_order),
            "second_paid": self.second_paid_ids,
        }


def build_state(seed_profiles):
    replays = []
    for index, (centers, radii) in enumerate(FIXTURES):
        replay = TwoSeedReplay(
            f"second_seed_core_{index}",
            centers=centers,
            radii=radii,
            emit=False,
        )
        replay.seed_profiles = seed_profiles
        replay.run()
        replays.append(replay)

    base_sequences = tuple(replay.base_paid_order for replay in replays)
    if any(sequence != base_sequences[0] for sequence in base_sequences[1:]):
        raise AssertionError("三个夹具的第一种子对象顺序不一致")
    locator_sequences = tuple(replay.second_locator_ids for replay in replays)
    if any(sequence != locator_sequences[0] for sequence in locator_sequences[1:]):
        raise AssertionError("三个夹具的第二种子定位线顺序不一致")

    points = tuple(
        PointBundle(
            f"O{center_index + 1}",
            tuple(
                tuple(float(coordinate) for coordinate in centers[center_index])
                for centers, _ in FIXTURES
            ),
        )
        for center_index in range(3)
    )
    input_drawables = tuple(
        DrawableBundle(
            node_id,
            "circle",
            tuple(
                object_value(replay.objects, node_id)[1] for replay in replays
            ),
        )
        for node_id in ("Gamma1", "Gamma2", "Gamma3")
    )
    state = State(points, input_drawables)
    for move, node_ids in enumerate(zip(*base_sequences, strict=True), start=1):
        rows = tuple(
            object_value(replay.objects, node_id)
            for replay, node_id in zip(replays, node_ids, strict=True)
        )
        if {kind for kind, _ in rows} != {rows[0][0]}:
            raise AssertionError(f"第 {move} 个对象种类不一致：{node_ids}")
        state = add_known_drawable(
            state,
            node_ids[0],
            rows[0][0],
            tuple(value for _, value in rows),
            move,
        )

    locator_values = tuple(
        tuple(
            object_value(replay.objects, locator_ids[index])[1]
            for replay, locator_ids in zip(
                replays,
                locator_sequences,
                strict=True,
            )
        )
        for index in range(3)
    )
    second_profile = seed_profiles[1]
    sample_targets = tuple(
        replay.verify_pair(
            second_profile,
            replay.second_tau,
            allow_repeated_physical_signs=True,
        )
        for replay in replays
    )
    keys = tuple(sorted(sample_targets[0]))
    contacts = tuple(
        tuple(
            tuple(
                float_scalar(coordinate)
                for coordinate in collapse_point(targets[key]["contact_3"])
            )
            for targets in sample_targets
        )
        for key in keys
    )
    if len(contacts) != 2:
        raise AssertionError("第二种子方向类没有恰好两个接触点")
    return state, locator_sequences[0], locator_values, contacts


def carries_contacts(drawable, contacts) -> bool:
    first, second = contacts
    return all(
        on_drawable(first_point, drawable.kind, value)
        and on_drawable(second_point, drawable.kind, value)
        for first_point, second_point, value in zip(
            first,
            second,
            drawable.values,
            strict=True,
        )
    )


def carries_contact(drawable, contact) -> bool:
    return all(
        on_drawable(point, drawable.kind, value)
        for point, value in zip(contact, drawable.values, strict=True)
    )


def search_seed_pair(seed_profiles):
    state, locator_ids, locator_values, contacts = build_state(seed_profiles)
    base_candidates = tuple(candidates(state))
    direct_contact_hits = tuple(
        tuple(
            candidate.describe()
            for candidate in base_candidates
            if carries_contact(candidate.drawable, contact)
        )
        for contact in contacts
    )
    locator_candidates = []
    for locator_id, expected_values in zip(
        locator_ids,
        locator_values,
        strict=True,
    ):
        matches = tuple(
            candidate
            for candidate in base_candidates
            if candidate.kind == "line"
            and all(
                same_drawable("line", value, expected)
                for value, expected in zip(
                    candidate.drawable.values,
                    expected_values,
                    strict=True,
                )
            )
        )
        if len(matches) != 1:
            raise AssertionError(
                f"标准定位线 {locator_id} 的候选绑定不是唯一的：{len(matches)}"
            )
        locator_candidates.append(matches[0])

    checked = 0
    hits = []
    subset_reports = []
    for first, second in combinations(locator_candidates, 2):
        state_one = apply_candidate(state, first, len(state.program) + 1)
        state_two = apply_candidate(state_one, second, len(state_one.program) + 1)
        pair_checked = 0
        pair_hits = []
        for third in candidates(state_two):
            checked += 1
            pair_checked += 1
            if carries_contacts(third.drawable, contacts):
                path = (first.describe(), second.describe(), third.describe())
                pair_hits.append(path)
                hits.append(path)
        subset_reports.append(
            {
                "locators": (first.describe(), second.describe()),
                "candidates": pair_checked,
                "hits": len(pair_hits),
            }
        )

    report = {
        "seed_profiles": seed_profiles,
        "base_points": len(state.points),
        "base_drawables": len(state.drawables),
        "base_candidates": len(base_candidates),
        "direct_contact_hits": tuple(len(rows) for rows in direct_contact_hits),
        "locator_subsets": len(subset_reports),
        "third_candidates": checked,
        "hits": len(hits),
    }
    print("second_seed_core_state", report, flush=True)
    for path in hits:
        print("candidate", {"seed_profiles": seed_profiles, "path": path})
    return report


def main() -> None:
    reports = tuple(search_seed_pair(pair) for pair in SHARED_SEED_PAIRS)
    print(
        "mannheim_second_seed_core_3e_search",
        {
            "samples": 3,
            "states": len(reports),
            "base_candidates": sum(
                report["base_candidates"] for report in reports
            ),
            "direct_contact_hits": sum(
                sum(report["direct_contact_hits"]) for report in reports
            ),
            "locator_subsets": sum(
                report["locator_subsets"] for report in reports
            ),
            "third_candidates": sum(
                report["third_candidates"] for report in reports
            ),
            "hits": sum(report["hits"] for report in reports),
        },
    )


if __name__ == "__main__":
    main()
