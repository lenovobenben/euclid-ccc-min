"""精确重放由有限 ``K'`` 提前构造的 6 E 双目标后缀。

对 Mannheim 的每个方向类，令 ``K'`` 为圆内接四边形的一组对角弦
交点，并令 ``S`` 为三个输入圆的根心。符号恒等式给出

    J = Line(O3, K') ∩ Line(O1, O2)

位于该方向类的目标圆心线上。因此在 ``K'`` 有限且相关点互异时，先画
``Line(O3,K')`` 和 ``Line(S,J)``，再分别画两个第三圆接触半径与两个
输出圆，一对目标只需 6 E。构造退化时沿用 7 E 圆心线后缀。
"""

from __future__ import annotations

from fractions import Fraction

from check_mannheim_degenerate_fixture import (
    add,
    dot,
    line_intersection,
    line_through,
    same_line,
    subtract,
)
from replay_mannheim_center_locus_dependencies import (
    CenterLocusReplay,
    build_seed_tau,
    build_tail_tau,
    build_targets_7e,
    on_line,
)
from replay_mannheim_ordered_branches import (
    FIXTURES,
    OTHER_BATCH,
    P0_BATCH,
    SEED_PROFILES,
    TAIL_PROFILES,
    branch_data,
    ordered_target_keys,
    paid_count,
    preferred_merged_key,
)
from replay_mannheim_three_block_dependencies import (
    collapse_circle,
    collapse_point,
    point_key,
    role_batch_keys,
)
from scan_mannheim_degeneracies import solve_dot_system


F = Fraction
EXTRA_FIXTURES = {
    "root_on_ell": (
        ((F(0), F(0)), (F(12), F(0)), (F(3), F(9, 2))),
        (F(3), F(2), F(1)),
    ),
    "merge_P0_double_Kp_P2_P3": (
        (
            (F(0), F(0)),
            (F(145, 18), F(0)),
            (F(25, 4), F(7, 4)),
        ),
        (F(5), F(13, 9), F(1)),
    ),
}


def build_targets_6e_from_kp(
    replay,
    profile,
    data,
    tau,
    tau_id,
    mannheim_point,
    mannheim_point_id,
) -> bool:
    """使用有限 ``K'`` 构造一对目标；不可用时返回 ``False``。"""

    kp_data = data.get("Kp")
    ell = line_through(replay.o1, replay.o2)
    if (
        kp_data is None
        or kp_data[0] != "finite"
        or kp_data[2] == replay.o3
        or on_line(mannheim_point, ell)
        or on_line(replay.o3, tau)
    ):
        return False
    kp = kp_data[2]
    assert kp is not None
    kp_id = replay.objects.point_registry.get(point_key(kp))
    if kp_id is None:
        return False

    before = paid_count(replay)
    radial = line_through(replay.o3, kp)
    radial_id = replay.objects.line(
        f"{profile}_O3_Kp",
        radial,
        "O3",
        kp_id,
    )
    j = line_intersection(radial, ell)
    j_id = replay.objects.point(
        f"{profile}_center_locus_J",
        j,
        radial_id,
        replay.objects.resolve("ell"),
    )
    center_locus = line_through(mannheim_point, j)
    center_locus_id = replay.objects.line(
        f"{profile}_center_locus_from_Kp",
        center_locus,
        mannheim_point_id,
        j_id,
    )

    targets = replay.verify_pair(
        profile,
        tau,
        allow_repeated_physical_signs=True,
    )
    contact_ids = {
        key: replay.objects.point(
            f"{profile}_{key}_M3_Kp_locus",
            collapse_point(target["contact_3"]),
            tau_id,
            "Gamma3",
        )
        for key, target in targets.items()
    }
    for key in ordered_target_keys(replay, profile, tau):
        sign = f"{profile}:{key}"
        contact = collapse_point(targets[key]["contact_3"])
        radius = line_through(replay.o3, contact)
        if same_line(center_locus, radius):
            raise AssertionError("非居中有限 K' 分支的定心线不能重合")
        radius_id = replay.objects.line(
            f"{profile}_{sign}_radius_3_from_Kp_locus",
            radius,
            "O3",
            contact_ids[key],
        )
        center = line_intersection(center_locus, radius)
        expected_center = collapse_point(targets[key]["center"])
        if center != expected_center:
            raise AssertionError(f"{profile} 没有由有限 K' 恢复目标圆心")
        center_id = replay.objects.point(
            f"{profile}_{sign}_center_from_Kp_locus",
            center,
            center_locus_id,
            radius_id,
        )
        output_circle = collapse_circle(targets[key]["output_circle"])
        output_id = replay.objects.circle(
            f"target_{sign}",
            output_circle,
            center_id,
            contact_ids[key],
        )
        replay.targets[sign] = {
            "profile": profile,
            "output_id": output_id,
            "circle": output_circle,
            "draw_index": replay.objects.graph.paid_order.index(output_id) + 1,
        }

    if paid_count(replay) - before > 6:
        raise AssertionError(f"{profile} 的有限 K' 双目标后缀超过 6 E")
    return True


class KpCenterLocusReplay(CenterLocusReplay):
    # 本程序针对八解联合祖先；三重外切单解仍沿用独立的 19 E 程序。
    first_ext_limit = 57

    def root_center_on_ell(self):
        e2 = subtract(self.o2, self.o1)
        e3 = subtract(self.o3, self.o1)
        relative = solve_dot_system(
            e2,
            e3,
            (
                dot(e2, e2) - self.r2**2 + self.r1**2
            ) / 2,
            (
                dot(e3, e3) - self.r3**2 + self.r1**2
            ) / 2,
        )
        center = add(self.o1, relative)
        return center if on_line(center, line_through(self.o1, self.o2)) else None

    def run_root_on_ell(self, data, mannheim_point) -> dict:
        before = paid_count(self)
        self.build_prefix()
        for key in P0_BATCH + OTHER_BATCH:
            self.draw_batch(key)
        input_cost = paid_count(self) - before

        before = paid_count(self)
        (tau, tau_id), targets_completed = build_seed_tau(
            self,
            "P0",
            data["P0"],
        )
        if targets_completed or paid_count(self) - before > 5:
            raise AssertionError("根心在线分支的 P0 核心超过 5 E")
        kp = data["P0"].get("Kp", (None, None, None))[2]
        if kp != mannheim_point:
            raise AssertionError("根心在线分支没有满足 P0:K'=S")
        mannheim_point_id = self.objects.point_registry[point_key(kp)]
        seed_costs = {"P0": paid_count(self) - before}

        tau_results = {"P0": (tau, tau_id)}
        tail_core_costs = {}
        for profile in ("P2", "P1", "P3"):
            before = paid_count(self)
            tau_results[profile] = build_tail_tau(
                self,
                profile,
                data[profile],
                mannheim_point,
                mannheim_point_id,
            )
            tail_core_costs[profile] = paid_count(self) - before
            if tail_core_costs[profile] > 3:
                raise AssertionError(f"{profile} 的后行接触弦核心超过 3 E")

        suffix_costs = {}
        for profile in SEED_PROFILES + TAIL_PROFILES:
            before = paid_count(self)
            self.build_pair(
                profile,
                data[profile],
                *tau_results[profile],
                mannheim_point,
                mannheim_point_id,
            )
            suffix_costs[profile] = paid_count(self) - before

        return self.audit_center_locus(
            data,
            input_cost,
            seed_costs,
            tail_core_costs,
            suffix_costs,
        )

    def run(self) -> dict:
        data = {
            profile: branch_data(self, profile)
            for profile in SEED_PROFILES + TAIL_PROFILES
        }
        simple_seeds = tuple(
            profile
            for profile in SEED_PROFILES
            if data[profile]["kind"] == "simple_merge"
        )
        root_on_ell = self.root_center_on_ell()
        if (
            not simple_seeds
            and root_on_ell is not None
            and data["P0"].get("Kp", (None, None, None))[0] == "finite"
        ):
            return self.run_root_on_ell(data, root_on_ell)
        if not simple_seeds:
            return super().run()

        selected = simple_seeds[0]
        preferred_key = preferred_merged_key(
            selected,
            data[selected]["roles"],
        )
        selected_keys = role_batch_keys(selected)
        merged_names = next(
            names
            for names in (("x", "z"), ("y", "w"))
            if data[selected]["roles"][names[0]]
            == data[selected]["roles"][names[1]]
        )
        merged_keys = tuple(selected_keys[name] for name in merged_names)
        redundant_key = next(key for key in merged_keys if key != preferred_key)
        batch_order = (preferred_key,) + tuple(
            key
            for key in P0_BATCH + OTHER_BATCH
            if key not in {preferred_key, redundant_key}
        )

        before = paid_count(self)
        self.build_prefix()
        for key in batch_order:
            self.draw_batch(key)
        self.batch_points[redundant_key] = self.batch_points[preferred_key]
        self.batch_point_ids[redundant_key] = self.batch_point_ids[preferred_key]
        input_cost = paid_count(self) - before
        if input_cost > 12:
            raise AssertionError("简单合并分支的裁剪前缀超过 12 E")

        before = paid_count(self)
        (tau, tau_id), targets_completed = build_seed_tau(
            self,
            selected,
            data[selected],
        )
        if not targets_completed:
            raise AssertionError("选定的简单合并块没有独立完成两个目标")

        selected_targets = tuple(
            target
            for target in self.targets.values()
            if target["profile"] == selected
        )
        if len(selected_targets) != 2:
            raise AssertionError("简单合并块没有恰好两个目标圆")
        centers = tuple(target["circle"][0] for target in selected_targets)
        center_ids = tuple(
            self.objects.point_registry[point_key(center)] for center in centers
        )
        center_locus = line_through(*centers)
        if same_line(center_locus, tau):
            raise AssertionError("有限简单合并的圆心线不能等于接触弦")
        center_locus_id = self.objects.line(
            f"{selected}_completed_pair_center_locus",
            center_locus,
            *center_ids,
        )
        mannheim_point = line_intersection(tau, center_locus)
        mannheim_point_id = self.objects.point(
            "Mannheim_S_from_completed_merge_pair",
            mannheim_point,
            tau_id,
            center_locus_id,
        )
        seed_costs = {selected: paid_count(self) - before}
        if seed_costs[selected] > 14:
            raise AssertionError("简单合并块和根心圆心线超过 14 E")

        tau_results = {selected: (tau, tau_id)}
        tail_core_costs = {}
        for profile in SEED_PROFILES + TAIL_PROFILES:
            if profile == selected:
                continue
            before = paid_count(self)
            tau_results[profile] = build_tail_tau(
                self,
                profile,
                data[profile],
                mannheim_point,
                mannheim_point_id,
            )
            tail_core_costs[profile] = paid_count(self) - before
            if tail_core_costs[profile] > 3:
                raise AssertionError(f"{profile} 的后行接触弦核心超过 3 E")

        suffix_costs = {}
        for profile in SEED_PROFILES + TAIL_PROFILES:
            if profile == selected:
                continue
            before = paid_count(self)
            self.build_pair(
                profile,
                data[profile],
                *tau_results[profile],
                mannheim_point,
                mannheim_point_id,
            )
            suffix_costs[profile] = paid_count(self) - before

        return self.audit_center_locus(
            data,
            input_cost,
            seed_costs,
            tail_core_costs,
            suffix_costs,
        )

    def build_pair(
        self,
        profile,
        data,
        tau,
        tau_id,
        mannheim_point,
        mannheim_point_id,
    ) -> None:
        if build_targets_6e_from_kp(
            self,
            profile,
            data,
            tau,
            tau_id,
            mannheim_point,
            mannheim_point_id,
        ):
            return
        build_targets_7e(
            self,
            profile,
            tau,
            tau_id,
            mannheim_point,
            mannheim_point_id,
        )


def main() -> None:
    reports = {
        fixture_name: KpCenterLocusReplay(fixture_name).run()
        for fixture_name in FIXTURES
    }
    reports.update(
        {
            fixture_name: KpCenterLocusReplay(
                fixture_name,
                centers=centers,
                radii=radii,
            ).run()
            for fixture_name, (centers, radii) in EXTRA_FIXTURES.items()
        }
    )
    expected_costs = {
        "regular": 49,
        "centered_P0": 48,
        "single_merge_P0": 51,
        "double_merge_P2": 46,
        "three_merges": 42,
        "single_K_P2": 47,
        "single_Kp_P0": 52,
        "double_Kp_P2_P3": 55,
        "merge_P1_parallel_P2_K": 48,
        "merge_P1_parallel_P3_K_alias": 46,
        "merge_P1_parallel_P3_K_generic": 50,
        "merge_P2_parallel_P3_Kp": 54,
        "root_on_ell": 55,
        "merge_P0_double_Kp_P2_P3": 55,
    }
    if any(
        reports[name]["all_targets"] != cost
        for name, cost in expected_costs.items()
    ):
        raise AssertionError("有限 K' 圆心线夹具的精确 E 分数发生变化")
    print(
        "kp_center_locus_summary",
        {
            "fixtures": len(reports),
            "maximum_trace": max(report["trace"] for report in reports.values()),
            "maximum_all_targets": max(
                report["all_targets"] for report in reports.values()
            ),
            "costs": {
                name: report["all_targets"] for name, report in reports.items()
            },
        },
    )


if __name__ == "__main__":
    main()
