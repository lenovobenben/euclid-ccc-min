"""搜索先完成两方向类后能否用 2 E 取得余下四个接触点。

正规 49 E 程序先构造四条接触弦，再完成八个目标。本脚本改为先构造
一对共享定义弦的接触弦及根心，立即用当前 6 E 后缀画出这四个目标，
然后检查新增交点能否缩短另外两个方向类的三步接触弦尾部。四种同成本
共享种子对全部检查。

一条新直线或圆若同时经过第三圆上的两个尚缺目标接触点，它与第三圆的
确定交点就能绑定这两点，不要求二者属于同一方向类，也不要求该对象就是
原接触弦。若两步对象能覆盖四个尚缺接触点，它们即可替代当前的
``1 条共享 K' 定位线 + 2 条接触弦``，把正规八解从 49 E 降到 48 E。

脚本同步使用三个严格正规 ``D8`` 夹具，只使用已有对象的有限实交点。
浮点多夹具命中必须另作精确重放和符号证明；零结果不是 2 E 下界。
"""

from __future__ import annotations

from scan_mannheim_degeneracies import is_d8
from check_mannheim_degenerate_fixture import line_intersection, same_line
from replay_mannheim_center_locus_dependencies import build_seed_tau
from replay_mannheim_kp_center_locus_dependencies import KpCenterLocusReplay
from replay_mannheim_ordered_branches import (
    OTHER_BATCH,
    P0_BATCH,
    SEED_PROFILES,
    TAIL_PROFILES,
    branch_data,
)
from replay_mannheim_three_block_dependencies import (
    collapse_point,
)
from search_mannheim_center_locus_2e import object_value
from search_mannheim_double_kp_global_sequential_2e import float_scalar
from search_mannheim_kp_parallel_center_locus_2e import candidates, on_drawable
from search_mannheim_regular_sequential_locus_reuse import same_point
from search_mannheim_root_center_2e import FIXTURES, add_known_drawable
from search_parallel_3e import DrawableBundle, PointBundle, State, apply_candidate


PROFILES = SEED_PROFILES + TAIL_PROFILES
SHARED_SEED_PAIRS = (
    ("P0", "P2"),
    ("P0", "P3"),
    ("P1", "P2"),
    ("P1", "P3"),
)


class SeedPairsFirstReplay(KpCenterLocusReplay):
    """先完成两个根心种子方向类，并在后行接触弦之前停止。"""

    seed_profiles = ("P0", "P2")

    def run(self):
        if not is_d8(self.centers, self.radii):
            raise AssertionError(f"{self.fixture_name} 不属于严格 D8")
        data = {profile: branch_data(self, profile) for profile in PROFILES}
        if {row["kind"] for row in data.values()} != {"regular"}:
            raise AssertionError("种子先行搜索只接受四方向均正规的夹具")

        self.build_prefix()
        for key in P0_BATCH + OTHER_BATCH:
            self.draw_batch(key)
        if len(self.objects.graph.paid_order) != 13:
            raise AssertionError("正规公共前缀必须恰有 13 E")

        seed_results = {}
        for profile in self.seed_profiles:
            seed_results[profile], targets_completed = build_seed_tau(
                self,
                profile,
                data[profile],
            )
            if targets_completed:
                raise AssertionError("正规接触弦核心不应提前完成目标")
        if len(self.objects.graph.paid_order) != 22:
            raise AssertionError("两个正规种子接触弦后必须恰有 22 E")

        first_profile, second_profile = self.seed_profiles
        tau_first, tau_first_id = seed_results[first_profile]
        tau_second, tau_second_id = seed_results[second_profile]
        if same_line(tau_first, tau_second):
            raise AssertionError("两个种子接触弦不能重合")
        mannheim_point = line_intersection(tau_first, tau_second)
        mannheim_point_id = self.objects.point(
            "Mannheim_S_center_locus",
            mannheim_point,
            tau_first_id,
            tau_second_id,
        )

        for profile in self.seed_profiles:
            self.build_pair(
                profile,
                data[profile],
                *seed_results[profile],
                mannheim_point,
                mannheim_point_id,
            )
        if len(self.targets) != 4:
            raise AssertionError("种子先行程序没有恰好完成四个目标")
        if len(self.objects.graph.paid_order) != 34:
            raise AssertionError("两个种子方向类完成后必须恰有 34 E")

        self.core_data = data
        self.seed_results = seed_results
        self.tail_profiles = tuple(
            profile for profile in PROFILES if profile not in self.seed_profiles
        )
        return {
            "paid_objects": len(self.objects.graph.paid_order),
            "targets": len(self.targets),
        }


def build_state(seed_profiles):
    replays = []
    for index, (centers, radii) in enumerate(FIXTURES):
        replay = SeedPairsFirstReplay(
            f"seed_first_tail_contacts_{index}",
            centers=centers,
            radii=radii,
            emit=False,
        )
        replay.seed_profiles = seed_profiles
        replay.run()
        replays.append(replay)

    paid_sequences = tuple(
        tuple(replay.objects.graph.paid_order) for replay in replays
    )
    if any(sequence != paid_sequences[0] for sequence in paid_sequences[1:]):
        raise AssertionError("三个夹具的种子先行对象顺序不一致")

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
    for move, node_ids in enumerate(zip(*paid_sequences, strict=True), start=1):
        rows = tuple(
            object_value(replay.objects, node_id)
            for replay, node_id in zip(replays, node_ids, strict=True)
        )
        kinds = {kind for kind, _ in rows}
        if len(kinds) != 1:
            raise AssertionError(f"第 {move} 个对象种类不一致：{node_ids}")
        state = add_known_drawable(
            state,
            node_ids[0],
            rows[0][0],
            tuple(value for _, value in rows),
            move,
        )

    contacts = {}
    tail_profiles = replays[0].tail_profiles
    if any(replay.tail_profiles != tail_profiles for replay in replays[1:]):
        raise AssertionError("三个夹具的后行方向类不一致")
    for profile in tail_profiles:
        sample_targets = tuple(
            replay.verify_pair(
                profile,
                replay.core_data[profile]["tau"],
                allow_repeated_physical_signs=True,
            )
            for replay in replays
        )
        keys = tuple(sorted(sample_targets[0]))
        if any(tuple(sorted(targets)) != keys for targets in sample_targets[1:]):
            raise AssertionError(f"{profile} 的目标键在夹具间不一致")
        for key in keys:
            contacts[f"{profile}:{key}"] = tuple(
                tuple(
                    float_scalar(coordinate)
                    for coordinate in collapse_point(targets[key]["contact_3"])
                )
                for targets in sample_targets
            )
        if sum(name.startswith(f"{profile}:") for name in contacts) != 2:
            raise AssertionError(f"{profile} 没有恰好两个第三圆接触点")

    return state, contacts, tail_profiles, paid_sequences[0]


def carries_contact(drawable, contact) -> bool:
    return all(
        on_drawable(point, drawable.kind, value)
        for point, value in zip(
            contact,
            drawable.values,
            strict=True,
        )
    )


def contact_hits(drawable, contacts):
    return tuple(
        name
        for name, contact in contacts.items()
        if carries_contact(drawable, contact)
    )


def search_seed_pair(seed_profiles):
    state, contacts, tail_profiles, paid_sequence = build_state(seed_profiles)
    old_contact_points = {
        name: tuple(
            point.point_id
            for point in state.points
            if all(
                same_point(value, expected)
                for value, expected in zip(
                    point.values,
                    contact,
                    strict=True,
                )
            )
        )
        for name, contact in contacts.items()
    }
    old_carriers = tuple(
        (drawable.drawable_id, contact_hits(drawable, contacts))
        for drawable in state.drawables
        if drawable.drawable_id != "Gamma3"
        if len(contact_hits(drawable, contacts)) >= 2
    )

    checked = 0
    hits = []
    for candidate in candidates(state):
        checked += 1
        names = contact_hits(candidate.drawable, contacts)
        if len(names) >= 2:
            hits.append((candidate, names))

    sequential_checked = 0
    sequential_hits = []
    all_contacts = frozenset(contacts)
    for first, first_contacts in hits:
        state_one = apply_candidate(state, first, len(paid_sequence) + 1)
        remaining = all_contacts - frozenset(first_contacts)
        for second in candidates(state_one):
            sequential_checked += 1
            second_contacts = frozenset(contact_hits(second.drawable, contacts))
            if remaining <= second_contacts:
                sequential_hits.append(
                    {
                        "first": first.describe(),
                        "first_contacts": first_contacts,
                        "second": second.describe(),
                        "second_contacts": tuple(sorted(second_contacts)),
                    }
                )

    report = {
        "seed_profiles": seed_profiles,
        "tail_profiles": tail_profiles,
        "paid_seed_first_state": len(paid_sequence),
        "points": len(state.points),
        "drawables": len(state.drawables),
        "old_contact_points": old_contact_points,
        "old_carriers": old_carriers,
        "first_candidates": checked,
        "first_carrier_hits": len(hits),
        "sequential_candidates": sequential_checked,
        "sequential_hits": len(sequential_hits),
    }
    print("seed_first_tail_state", report, flush=True)
    for candidate, names in hits:
        print(
            "first_contact_carrier",
            {"candidate": candidate.describe(), "contacts": names},
        )
    for hit in sequential_hits:
        print("two_step_tail", hit)
    return report


def main() -> None:
    reports = tuple(search_seed_pair(pair) for pair in SHARED_SEED_PAIRS)
    print(
        "mannheim_seed_first_tail_contact_carrier_search",
        {
            "samples": 3,
            "states": len(reports),
            "first_candidates": sum(
                report["first_candidates"] for report in reports
            ),
            "first_carrier_hits": sum(
                report["first_carrier_hits"] for report in reports
            ),
            "sequential_candidates": sum(
                report["sequential_candidates"] for report in reports
            ),
            "sequential_hits": sum(
                report["sequential_hits"] for report in reports
            ),
        },
    )


if __name__ == "__main__":
    main()
