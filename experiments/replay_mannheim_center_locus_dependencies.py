"""精确重放经过根心的目标圆心线及 7 E 双目标后缀。

固定一个 Mannheim 方向类时，两个有向根的圆心都形如 ``S + rho*u``，
其中 ``S`` 是三个输入圆的根心。因此先用标准后缀完成第一个目标，再画
它的圆心与 ``S`` 的连线；第二个目标圆心只需由这条线与第三圆半径线
相交得到。一对目标至多使用 5 条线和 2 个圆，即 7 E。

程序先只构造 ``P0``、``P2`` 的接触弦以取得 ``S``，再完成尚未输出的
目标。有限简单合并先行块保留原 13 E 完整模块，因为其接触弦修复本身
较贵；其它先行块只画接触弦核心。计入简单合并必然产生的冗余批量线
后，统一 E 上界由 61 降为 57。
"""

from __future__ import annotations

from check_mannheim_degenerate_fixture import (
    add,
    determinant,
    line_intersection,
    line_through,
    multiply,
    same_line,
    subtract,
)
from replay_mannheim_ordered_branches import (
    FIXTURES,
    OTHER_BATCH,
    P0_BATCH,
    SEED_PROFILES,
    TAIL_PROFILES,
    OrderedBranchReplay,
    branch_data,
    merge_pairs,
    ordered_target_keys,
    paid_count,
    preferred_merged_key,
)
from replay_mannheim_three_block_dependencies import (
    circle,
    collapse_circle,
    collapse_point,
    on_circle,
    role_batch_keys,
)
from scan_mannheim_degeneracies import is_d8


def on_line(point, line) -> bool:
    return line[0] * point[0] + line[1] * point[1] + line[2] == 0


def build_targets_7e(
    replay: OrderedBranchReplay,
    profile: str,
    tau,
    tau_id: str,
    mannheim_point,
    mannheim_point_id: str,
) -> None:
    before = paid_count(replay)
    targets = replay.verify_pair(
        profile, tau, allow_repeated_physical_signs=True
    )
    contact_ids = {
        key: replay.objects.point(
            f"{profile}_{key}_M3_center_locus",
            collapse_point(target["contact_3"]),
            tau_id,
            "Gamma3",
        )
        for key, target in targets.items()
    }
    target_order = ordered_target_keys(replay, profile, tau)

    if on_line(replay.o3, tau):
        for key in target_order:
            replay.build_target(
                profile,
                f"{profile}:{key}",
                targets[key],
                contact_ids[key],
            )
        if paid_count(replay) - before > 6:
            raise AssertionError("居中接触弦的两个目标超过 6 E")
        return

    first_key, second_key = target_order
    first_sign = f"{profile}:{first_key}"
    replay.build_target(
        profile,
        first_sign,
        targets[first_key],
        contact_ids[first_key],
    )
    first_center = collapse_point(targets[first_key]["center"])
    second_center = collapse_point(targets[second_key]["center"])
    if determinant(
        subtract(first_center, mannheim_point),
        subtract(second_center, mannheim_point),
    ) != 0:
        raise AssertionError(f"{profile} 的两个目标圆心与根心不共线")
    if first_center == mannheim_point:
        raise AssertionError("非零目标圆的圆心不能等于根心")

    first_center_id = replay.objects.resolve(f"{profile}_{first_sign}_center")
    center_locus = line_through(mannheim_point, first_center)
    center_locus_id = replay.objects.line(
        f"{profile}_center_locus",
        center_locus,
        mannheim_point_id,
        first_center_id,
    )

    second_sign = f"{profile}:{second_key}"
    second_contact = collapse_point(targets[second_key]["contact_3"])
    second_radius = line_through(replay.o3, second_contact)
    if same_line(center_locus, second_radius):
        raise AssertionError("非居中分支的圆心线不能与第二半径线重合")
    second_radius_id = replay.objects.line(
        f"{profile}_{second_sign}_radius_3_from_locus",
        second_radius,
        "O3",
        contact_ids[second_key],
    )
    recovered_center = line_intersection(center_locus, second_radius)
    if recovered_center != second_center:
        raise AssertionError(f"{profile} 没有从圆心线恢复第二目标圆心")
    second_center_id = replay.objects.point(
        f"{profile}_{second_sign}_center_from_locus",
        second_center,
        center_locus_id,
        second_radius_id,
    )
    output_circle = collapse_circle(targets[second_key]["output_circle"])
    output_id = replay.objects.circle(
        f"target_{second_sign}",
        output_circle,
        second_center_id,
        contact_ids[second_key],
    )
    replay.targets[second_sign] = {
        "profile": profile,
        "output_id": output_id,
        "circle": output_circle,
        "draw_index": replay.objects.graph.paid_order.index(output_id) + 1,
    }
    if paid_count(replay) - before > 7:
        raise AssertionError(f"{profile} 的圆心线双目标后缀超过 7 E")


def build_regular_tau_core(replay, profile, data):
    roles = data["roles"]
    keys = role_batch_keys(profile)
    point_ids = {
        name: replay.batch_point_ids[key] for name, key in keys.items()
    }
    diagonal_points = {}
    for diagonal_name, pairs in (
        ("K", (("x", "w"), ("y", "z"))),
        ("Kp", (("x", "y"), ("z", "w"))),
    ):
        line_ids = tuple(
            replay.objects.line(
                f"{profile}_{diagonal_name}_core_line_{index}",
                line_through(roles[first], roles[second]),
                point_ids[first],
                point_ids[second],
            )
            for index, (first, second) in enumerate(pairs, start=1)
        )
        point = data[diagonal_name][2]
        diagonal_points[diagonal_name] = (
            point,
            replay.objects.point(
                f"{profile}_{diagonal_name}_core_point", point, *line_ids
            ),
        )
    tau_id = replay.objects.line(
        f"{profile}_tau",
        data["tau"],
        diagonal_points["K"][1],
        diagonal_points["Kp"][1],
    )
    return data["tau"], tau_id


def build_parallel_tau_core(replay, profile, data):
    roles = data["roles"]
    keys = role_batch_keys(profile)
    point_ids = {
        name: replay.batch_point_ids[key] for name, key in keys.items()
    }
    infinite = "K" if data["kind"] == "K_parallel" else "Kp"
    finite = "Kp" if infinite == "K" else "K"
    _, finite_lines, finite_point = data[finite]
    assert finite_lines is not None and finite_point is not None
    finite_pairs = (
        (("x", "y"), ("z", "w"))
        if finite == "Kp"
        else (("x", "w"), ("y", "z"))
    )
    finite_line_ids = tuple(
        replay.objects.line(
            f"{profile}_{finite}_parallel_core_line_{index}",
            line,
            *(point_ids[name] for name in pair),
        )
        for index, (line, pair) in enumerate(
            zip(finite_lines, finite_pairs, strict=True), start=1
        )
    )
    finite_point_id = replay.objects.point(
        f"{profile}_{finite}_parallel_core_point",
        finite_point,
        *finite_line_ids,
    )
    first_name, second_name = (
        ("x", "y") if infinite == "Kp" else ("x", "w")
    )
    first = roles[first_name]
    second = roles[second_name]
    first_circle = circle(first, finite_point)
    first_circle_id = replay.objects.circle(
        f"{profile}_{infinite}_parallel_core_circle_1",
        first_circle,
        point_ids[first_name],
        finite_point_id,
    )
    reflected = add(first, subtract(first, finite_point))
    reflected_id = replay.objects.point(
        f"{profile}_{infinite}_parallel_core_reflected",
        reflected,
        first_circle_id,
        finite_line_ids[0],
    )
    second_circle = circle(second, reflected)
    second_circle_id = replay.objects.circle(
        f"{profile}_{infinite}_parallel_core_circle_2",
        second_circle,
        point_ids[second_name],
        reflected_id,
    )
    q = add(finite_point, subtract(first, second))
    if not on_circle(q, first_circle) or not on_circle(q, second_circle):
        raise AssertionError("平行核心的第二辅助圆交点错误")
    q_id = replay.objects.point(
        f"{profile}_{infinite}_parallel_core_Q",
        q,
        first_circle_id,
        second_circle_id,
    )
    tau_id = replay.objects.line(
        f"{profile}_tau",
        data["tau"],
        finite_point_id,
        q_id,
    )
    return data["tau"], tau_id


def build_short_tau_core(replay, profile, data):
    roles = data["roles"]
    keys = role_batch_keys(profile)
    kind = data["kind"]
    if kind == "centered_Kp_parallel":
        x, z = roles["x"], roles["z"]
        circle_x = circle(x, replay.o3)
        circle_z = circle(z, replay.o3)
        circle_x_id = replay.objects.circle(
            f"{profile}_centered_core_circle_x",
            circle_x,
            replay.batch_point_ids[keys["x"]],
            "O3",
        )
        circle_z_id = replay.objects.circle(
            f"{profile}_centered_core_circle_z",
            circle_z,
            replay.batch_point_ids[keys["z"]],
            "O3",
        )
        q = add(add(x, z), multiply(-1, replay.o3))
        q_id = replay.objects.point(
            f"{profile}_centered_core_Q", q, circle_x_id, circle_z_id
        )
        dependencies = ("O3", q_id)
    elif kind == "double_merge":
        dependencies = (
            replay.batch_point_ids[keys["x"]],
            replay.batch_point_ids[keys["y"]],
        )
    elif kind == "centered_merge":
        merged_name = merge_pairs(roles)[0][0]
        dependencies = (
            "O3",
            replay.batch_point_ids[keys[merged_name]],
        )
    else:
        raise AssertionError(f"{kind} 不是短接触弦核心")
    tau_id = replay.objects.line(
        f"{profile}_tau", data["tau"], *dependencies
    )
    return data["tau"], tau_id


def build_seed_tau(replay, profile, data):
    kind = data["kind"]
    if kind == "regular":
        return build_regular_tau_core(replay, profile, data), False
    if kind in {"K_parallel", "Kp_parallel"}:
        return build_parallel_tau_core(replay, profile, data), False
    if kind in {"centered_Kp_parallel", "centered_merge", "double_merge"}:
        return build_short_tau_core(replay, profile, data), False
    if kind == "simple_merge":
        before = paid_count(replay)
        tau, tau_id = replay_mannheim_full_merge_seed(replay, profile, data)
        if paid_count(replay) - before > 13:
            raise AssertionError("有限简单合并先行块超过 13 E")
        return (tau, tau_id), True
    raise AssertionError(f"未实现的先行接触弦分支 {kind}")


def replay_mannheim_full_merge_seed(replay, profile, data):
    target_order = ordered_target_keys(replay, profile, data["tau"])
    preferred_key = preferred_merged_key(profile, data["roles"])
    replay.build_merge_pair(
        profile,
        target_order,
        preferred_merged_key=preferred_key,
        allow_repeated_physical_signs=True,
    )
    tau_id = replay.objects.line(
        f"{profile}_tau_for_S",
        data["tau"],
        replay.batch_point_ids[preferred_key],
        replay.objects.resolve(f"{profile}_second_contact"),
    )
    return data["tau"], tau_id


def build_tail_tau(
    replay,
    profile,
    data,
    mannheim_point,
    mannheim_point_id,
):
    roles = data["roles"]
    keys = role_batch_keys(profile)
    if data["kind"] in {"simple_merge", "centered_merge", "double_merge"}:
        merged_name = merge_pairs(roles)[0][0]
        tau_id = replay.objects.line(
            f"{profile}_tau_from_S_merge",
            data["tau"],
            mannheim_point_id,
            replay.batch_point_ids[keys[merged_name]],
        )
        return data["tau"], tau_id

    chosen = None
    for diagonal_name in ("Kp", "K"):
        state, lines, point = data[diagonal_name]
        if state == "finite" and point != mannheim_point:
            chosen = diagonal_name, lines, point
            break
    if chosen is None:
        raise AssertionError(f"{profile} 没有可用有限对角点")
    diagonal_name, lines, point = chosen
    assert lines is not None and point is not None
    pairs = (
        (("x", "y"), ("z", "w"))
        if diagonal_name == "Kp"
        else (("x", "w"), ("y", "z"))
    )
    line_ids = tuple(
        replay.objects.line(
            f"{profile}_{diagonal_name}_locus_line_{index}",
            line,
            *(replay.batch_point_ids[keys[name]] for name in pair),
        )
        for index, (line, pair) in enumerate(
            zip(lines, pairs, strict=True), start=1
        )
    )
    point_id = replay.objects.point(
        f"{profile}_{diagonal_name}_locus_point", point, *line_ids
    )
    tau_id = replay.objects.line(
        f"{profile}_tau_from_S",
        data["tau"],
        mannheim_point_id,
        point_id,
    )
    return data["tau"], tau_id


class CenterLocusReplay(OrderedBranchReplay):
    def run(self) -> dict:
        if not is_d8(self.centers, self.radii):
            raise AssertionError(f"{self.fixture_name} 不属于严格 D8")
        data = {
            profile: branch_data(self, profile)
            for profile in SEED_PROFILES + TAIL_PROFILES
        }
        p0_batch = P0_BATCH
        if data["P0"]["kind"] == "simple_merge":
            preferred_key = preferred_merged_key("P0", data["P0"]["roles"])
            p0_batch = (preferred_key,) + tuple(
                key for key in P0_BATCH if key != preferred_key
            )

        before = paid_count(self)
        self.build_prefix()
        for key in p0_batch:
            self.draw_batch(key)
        for key in OTHER_BATCH:
            self.draw_batch(key)
        input_cost = paid_count(self) - before
        if input_cost > 13:
            raise AssertionError("公共前缀和批量线超过 13 E")

        seed_results = {}
        completed = set()
        seed_costs = {}
        for profile in SEED_PROFILES:
            before = paid_count(self)
            result, targets_completed = build_seed_tau(self, profile, data[profile])
            seed_results[profile] = result
            seed_costs[profile] = paid_count(self) - before
            if targets_completed:
                completed.add(profile)
            elif seed_costs[profile] > 5:
                raise AssertionError(f"{profile} 的先行接触弦核心超过 5 E")

        tau0, tau0_id = seed_results["P0"]
        tau2, tau2_id = seed_results["P2"]
        if same_line(tau0, tau2):
            raise AssertionError("P0、P2 接触弦不能重合")
        mannheim_point = line_intersection(tau0, tau2)
        mannheim_point_id = self.objects.point(
            "Mannheim_S_center_locus",
            mannheim_point,
            tau0_id,
            tau2_id,
        )

        tau_results = dict(seed_results)
        tail_core_costs = {}
        for profile in TAIL_PROFILES:
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
            if profile in completed:
                continue
            tau, tau_id = tau_results[profile]
            before = paid_count(self)
            build_targets_7e(
                self,
                profile,
                tau,
                tau_id,
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

    def audit_center_locus(
        self,
        data,
        input_cost,
        seed_costs,
        tail_core_costs,
        suffix_costs,
    ):
        if len(self.targets) != 8:
            raise AssertionError("圆心线程序没有得到八个有向根")
        if len({target["circle"] for target in self.targets.values()}) != 8:
            raise AssertionError("圆心线程序的八个输出圆不是两两不同")
        graph = self.objects.graph
        ancestor_sets = {
            key: graph.paid_ancestors(target["output_id"])
            for key, target in self.targets.items()
        }
        union = frozenset().union(*ancestor_sets.values())
        if len(union) > 57 or len(graph.paid_order) > 59:
            raise AssertionError("圆心线程序超过 57 E 联合上界")
        line_count = sum(graph.paid_kinds[node] == "line" for node in union)
        circle_count = sum(graph.paid_kinds[node] == "circle" for node in union)
        ext_keys = tuple(
            key for key in self.targets if key.startswith("P0:+++@")
        )
        first_ext = min(len(ancestor_sets[key]) for key in ext_keys)
        if first_ext > 19:
            raise AssertionError("圆心线程序的三重外切目标超过 19 E")
        report = {
            "fixture": self.fixture_name,
            "branches": {
                profile: data[profile]["kind"]
                for profile in SEED_PROFILES + TAIL_PROFILES
            },
            "input": input_cost,
            "seed_cores": seed_costs,
            "tail_cores": tail_core_costs,
            "suffixes": suffix_costs,
            "trace": len(graph.paid_order),
            "union_lines": line_count,
            "union_circles": circle_count,
            "all_targets": len(union),
            "first_ext": first_ext,
            "non_ancestors": sorted(set(graph.paid_order) - set(union)),
            "aliases": len(self.objects.paid_aliases),
        }
        if self.emit:
            print("center_locus_replay", report)
        return report


def main() -> None:
    reports = {
        fixture_name: CenterLocusReplay(fixture_name).run()
        for fixture_name in FIXTURES
    }
    expected_costs = {
        "regular": 53,
        "centered_P0": 52,
        "single_merge_P0": 55,
        "double_merge_P2": 46,
        "three_merges": 47,
        "single_K_P2": 52,
        "single_Kp_P0": 55,
        "double_Kp_P2_P3": 57,
        "merge_P1_parallel_P2_K": 52,
        "merge_P1_parallel_P3_K_alias": 49,
        "merge_P1_parallel_P3_K_generic": 53,
        "merge_P2_parallel_P3_Kp": 57,
    }
    if any(
        reports[name]["all_targets"] != cost
        for name, cost in expected_costs.items()
    ):
        raise AssertionError("圆心线夹具套件的精确 E 分数发生变化")
    print(
        "center_locus_summary",
        {
            "fixtures": len(reports),
            "maximum_trace": max(report["trace"] for report in reports.values()),
            "maximum_all_targets": max(
                report["all_targets"] for report in reports.values()
            ),
        },
    )


if __name__ == "__main__":
    main()
