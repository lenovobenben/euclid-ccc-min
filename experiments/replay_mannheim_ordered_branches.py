"""统一重放 Mannheim 共点程序的有序分支证书。

程序固定先构造 ``P0``、``P2`` 的接触弦，再以它们的交点取得根心
``S``，最后构造 ``P1``、``P3``。先行块覆盖正规、一般平行、居中平行、
简单合并、合并后的居中极线和双对合并；后行块只用一个有限对角点或
一个合并接触点连接 ``S``。所有选择都由已经构造的角色点之间的相等、
重合和平行关系决定。

这份脚本的任务不是以有限夹具代替连续参数证明。连续覆盖来自文档第
8.10、8.12、8.14、8.15 节；这里逐对象检查分支顺序、合法依赖、几何
对象去重、八个输出圆和 ``61 E`` 台账。夹具集合覆盖仓库中已经发现的
正规、单/双平行、单/双/三合并以及合并和平行同现类型。
"""

from __future__ import annotations

from fractions import Fraction

from check_mannheim_degenerate_fixture import (
    add,
    determinant,
    dot,
    line_intersection,
    line_through,
    multiply,
    same_line,
    subtract,
)
from replay_mannheim_centered_parallel_repair import build_roles
from replay_mannheim_three_block_dependencies import (
    ThreeBlockReplay,
    circle,
    collapse_point,
    on_circle,
    role_batch_keys,
)
from scan_mannheim_degeneracies import is_d8


F = Fraction
Point = tuple[Fraction, Fraction]
Line = tuple[Fraction, Fraction, Fraction]

P0_BATCH = ("alphaA", "aB", "a1A", "alpha1B")
OTHER_BATCH = ("alphaB", "aA", "a1B", "alpha1A")
SEED_PROFILES = ("P0", "P2")
TAIL_PROFILES = ("P1", "P3")

FIXTURES = {
    "regular": (
        ((F(0), F(0)), (F(13), F(0)), (F(4), F(15))),
        (F(4), F(2), F(1)),
    ),
    "centered_P0": (
        ((F(0), F(0)), (F(9), F(0)), (F(6), F(4))),
        (F(3), F(2), F(1)),
    ),
    "single_merge_P0": (
        ((F(0), F(0)), (F(7), F(0)), (F(3), F(5))),
        (F(3), F(2), F(1)),
    ),
    "double_merge_P2": (
        ((F(0), F(0)), (F(16), F(0)), (F(12), F(6))),
        (F(5), F(3), F(1)),
    ),
    "three_merges": (
        ((F(0), F(0)), (F(253), F(0)), (F(155), F(120))),
        (F(8), F(5), F(3)),
    ),
    "single_K_P2": (
        ((F(0), F(0)), (F(7), F(0)), (F(8), F(4))),
        (F(3), F(2), F(1)),
    ),
    "single_Kp_P0": (
        ((F(0), F(0)), (F(28), F(0)), (F(13), F(27, 2))),
        (F(5), F(4), F(1)),
    ),
    "double_Kp_P2_P3": (
        ((F(0), F(0)), (F(45), F(0)), (F(30), F(21))),
        (F(4), F(2), F(1)),
    ),
    "merge_P1_parallel_P2_K": (
        ((F(0), F(0)), (F(12), F(0)), (F(13), F(6))),
        (F(4), F(2), F(1)),
    ),
    "merge_P1_parallel_P3_K_alias": (
        ((F(0), F(0)), (F(9), F(0)), (F(-2), F(8))),
        (F(4), F(3), F(2)),
    ),
    "merge_P1_parallel_P3_K_generic": (
        ((F(0), F(0)), (F(16), F(0)), (F(10), F(12))),
        (F(5), F(2), F(1)),
    ),
    "merge_P2_parallel_P3_Kp": (
        ((F(0), F(0)), (F(18), F(0)), (F(13), F(7))),
        (F(5), F(3), F(1)),
    ),
}

EXPECTED_RESULTS = {
    "regular": (("regular", "regular", "regular", "regular"), 57),
    "centered_P0": (
        ("centered_Kp_parallel", "regular", "regular", "regular"),
        55,
    ),
    "single_merge_P0": (
        ("simple_merge", "regular", "regular", "regular"),
        58,
    ),
    "double_merge_P2": (
        ("regular", "double_merge", "regular", "regular"),
        49,
    ),
    "three_merges": (
        ("simple_merge", "simple_merge", "regular", "simple_merge"),
        47,
    ),
    "single_K_P2": (
        ("regular", "K_parallel", "regular", "regular"),
        56,
    ),
    "single_Kp_P0": (
        ("Kp_parallel", "regular", "regular", "regular"),
        59,
    ),
    "double_Kp_P2_P3": (
        ("regular", "Kp_parallel", "regular", "Kp_parallel"),
        61,
    ),
    "merge_P1_parallel_P2_K": (
        ("regular", "K_parallel", "simple_merge", "regular"),
        54,
    ),
    "merge_P1_parallel_P3_K_alias": (
        ("regular", "regular", "simple_merge", "K_parallel"),
        52,
    ),
    "merge_P1_parallel_P3_K_generic": (
        ("regular", "regular", "simple_merge", "K_parallel"),
        55,
    ),
    "merge_P2_parallel_P3_Kp": (
        ("regular", "simple_merge", "regular", "Kp_parallel"),
        60,
    ),
}


def paid_count(replay: ThreeBlockReplay) -> int:
    return len(replay.objects.graph.paid_order)


def merge_pairs(roles: dict[str, Point]) -> tuple[tuple[str, str], ...]:
    return tuple(
        (first, second)
        for first, second in (("x", "z"), ("y", "w"))
        if roles[first] == roles[second]
    )


def preferred_merged_key(profile: str, roles: dict[str, Point]) -> str:
    keys = role_batch_keys(profile)
    names = merge_pairs(roles)[0]
    candidates = tuple(keys[name] for name in names)
    return next(
        key for key in candidates if key.startswith(("a1", "alpha1"))
    )


def diagonal(
    roles: dict[str, Point], name: str
) -> tuple[str, tuple[Line, Line] | None, Point | None]:
    pairs = (
        (("x", "w"), ("y", "z"))
        if name == "K"
        else (("x", "y"), ("z", "w"))
    )
    if any(roles[first] == roles[second] for first, second in pairs):
        return "undefined", None, None
    lines = tuple(
        line_through(roles[first], roles[second]) for first, second in pairs
    )
    if same_line(*lines):
        return "coincident", lines, None
    if determinant(lines[0][:2], lines[1][:2]) == 0:
        return "infinite", lines, None
    return "finite", lines, line_intersection(*lines)


def simple_merge_tau(
    roles: dict[str, Point], o3: Point, r3: Fraction
) -> tuple[Line, bool]:
    pairs = merge_pairs(roles)
    if len(pairs) != 1:
        raise AssertionError("简单合并极线需要恰好一个对向合并")
    first, _ = pairs[0]
    merged = roles[first]
    other_names = ("y", "w") if first == "x" else ("x", "z")
    other = line_through(*(roles[name] for name in other_names))
    tangent = (
        merged[0] - o3[0],
        merged[1] - o3[1],
        -(
            (merged[0] - o3[0]) * merged[0]
            + (merged[1] - o3[1]) * merged[1]
        ),
    )
    if determinant(tangent[:2], other[:2]) == 0:
        return line_through(o3, merged), True
    center = line_intersection(tangent, other)
    radius_squared = dot(subtract(center, merged), subtract(center, merged))
    tau = (
        2 * (center[0] - o3[0]),
        2 * (center[1] - o3[1]),
        o3[0] ** 2
        + o3[1] ** 2
        - r3**2
        - center[0] ** 2
        - center[1] ** 2
        + radius_squared,
    )
    return tau, False


def branch_data(replay: ThreeBlockReplay, profile: str) -> dict:
    roles = build_roles(replay.centers, replay.radii, profile)
    merges = merge_pairs(roles)
    if len(merges) == 2:
        return {
            "kind": "double_merge",
            "roles": roles,
            "tau": line_through(roles["x"], roles["y"]),
        }
    if len(merges) == 1:
        tau, centered = simple_merge_tau(roles, replay.o3, replay.r3)
        return {
            "kind": "centered_merge" if centered else "simple_merge",
            "roles": roles,
            "tau": tau,
        }

    k_state, k_lines, k = diagonal(roles, "K")
    kp_state, kp_lines, kp = diagonal(roles, "Kp")
    if {k_state, kp_state} - {"finite", "infinite"}:
        raise AssertionError(f"{profile} 的四点分支含未分类对角退化")
    if k_state == kp_state == "infinite":
        raise AssertionError(f"{profile} 的两个对角点不能同时位于无穷远")
    if kp_state == "infinite":
        if k == replay.o3:
            kind = "centered_Kp_parallel"
        else:
            kind = "Kp_parallel"
        assert k is not None and kp_lines is not None
        tau = (kp_lines[0][0], kp_lines[0][1], 0)
        tau = (tau[0], tau[1], -(tau[0] * k[0] + tau[1] * k[1]))
    elif k_state == "infinite":
        kind = "K_parallel"
        assert kp is not None and k_lines is not None
        tau = (k_lines[0][0], k_lines[0][1], 0)
        tau = (tau[0], tau[1], -(tau[0] * kp[0] + tau[1] * kp[1]))
    else:
        kind = "regular"
        assert k is not None and kp is not None
        tau = line_through(k, kp)
    return {
        "kind": kind,
        "roles": roles,
        "tau": tau,
        "K": (k_state, k_lines, k),
        "Kp": (kp_state, kp_lines, kp),
    }


def ordered_target_keys(
    replay: ThreeBlockReplay, profile: str, tau: Line
) -> tuple[str, str]:
    targets = replay.verify_pair(
        profile, tau, allow_repeated_physical_signs=True
    )
    keys = tuple(targets)
    if profile != "P0":
        return keys
    return tuple(
        sorted(
            keys,
            key=lambda key: targets[key]["physical_sign"] != "+++",
        )
    )  # type: ignore[return-value]


def build_targets_from_tau(
    replay: ThreeBlockReplay,
    profile: str,
    tau: Line,
    tau_id: str,
) -> None:
    targets = replay.verify_pair(
        profile, tau, allow_repeated_physical_signs=True
    )
    contact_ids = {
        key: replay.objects.point(
            f"{profile}_{key}_M3_ordered",
            collapse_point(target["contact_3"]),
            tau_id,
            "Gamma3",
        )
        for key, target in targets.items()
    }
    keys = ordered_target_keys(replay, profile, tau)
    for key in keys:
        replay.build_target(
            profile,
            f"{profile}:{key}",
            targets[key],
            contact_ids[key],
        )


def build_centered_seed(
    replay: ThreeBlockReplay, profile: str, data: dict
) -> str:
    roles = data["roles"]
    keys = role_batch_keys(profile)
    x, z = roles["x"], roles["z"]
    circle_x = circle(x, replay.o3)
    circle_z = circle(z, replay.o3)
    circle_x_id = replay.objects.circle(
        f"{profile}_centered_circle_x",
        circle_x,
        replay.batch_point_ids[keys["x"]],
        "O3",
    )
    circle_z_id = replay.objects.circle(
        f"{profile}_centered_circle_z",
        circle_z,
        replay.batch_point_ids[keys["z"]],
        "O3",
    )
    q = add(add(x, z), multiply(-1, replay.o3))
    if q == replay.o3 or not on_circle(q, circle_x) or not on_circle(q, circle_z):
        raise AssertionError("居中平行分支的第二圆交点错误")
    q_id = replay.objects.point(
        f"{profile}_centered_Q", q, circle_x_id, circle_z_id
    )
    tau_id = replay.objects.line(
        f"{profile}_tau",
        data["tau"],
        "O3",
        q_id,
    )
    build_targets_from_tau(replay, profile, data["tau"], tau_id)
    return tau_id


def build_direct_seed(
    replay: ThreeBlockReplay, profile: str, data: dict
) -> str:
    roles = data["roles"]
    keys = role_batch_keys(profile)
    if data["kind"] == "double_merge":
        point_names = ("x", "y")
    else:
        point_names = merge_pairs(roles)[0]
        point_names = (point_names[0],)
    if len(point_names) == 2:
        dependencies = tuple(
            replay.batch_point_ids[keys[name]] for name in point_names
        )
    else:
        dependencies = ("O3", replay.batch_point_ids[keys[point_names[0]]])
    tau_id = replay.objects.line(
        f"{profile}_tau",
        data["tau"],
        *dependencies,
    )
    build_targets_from_tau(replay, profile, data["tau"], tau_id)
    return tau_id


def build_seed(
    replay: ThreeBlockReplay, profile: str, data: dict
) -> tuple[Line, str]:
    kind = data["kind"]
    target_order = ordered_target_keys(replay, profile, data["tau"])
    options = {"allow_repeated_physical_signs": True}
    if kind == "regular":
        replay.build_regular_pair(profile, target_order, **options)
        tau_id = replay.objects.resolve(f"{profile}_tau")
    elif kind in {"K_parallel", "Kp_parallel"}:
        replay.build_parallel_pair(
            profile,
            "K" if kind == "K_parallel" else "Kp",
            target_order,
            **options,
        )
        tau_id = replay.objects.resolve(f"{profile}_tau")
    elif kind == "centered_Kp_parallel":
        tau_id = build_centered_seed(replay, profile, data)
    elif kind == "simple_merge":
        preferred_key = preferred_merged_key(profile, data["roles"])
        replay.build_merge_pair(
            profile,
            target_order,
            preferred_merged_key=preferred_key,
            **options,
        )
        roles = data["roles"]
        tau_id = replay.objects.line(
            f"{profile}_tau_for_S",
            data["tau"],
            replay.batch_point_ids[preferred_key],
            replay.objects.resolve(f"{profile}_second_contact"),
        )
    elif kind in {"centered_merge", "double_merge"}:
        tau_id = build_direct_seed(replay, profile, data)
    else:
        raise AssertionError(f"未实现的先行分支 {kind}")
    return data["tau"], tau_id


def build_tail(
    replay: ThreeBlockReplay,
    profile: str,
    data: dict,
    mannheim_point: Point,
    mannheim_point_id: str,
) -> str:
    roles = data["roles"]
    keys = role_batch_keys(profile)
    if data["kind"] in {"simple_merge", "centered_merge", "double_merge"}:
        merged_name = merge_pairs(roles)[0][0]
        tau = line_through(mannheim_point, roles[merged_name])
        if not same_line(tau, data["tau"]):
            raise AssertionError(f"{profile} 的合并点不在共点接触弦上")
        tau_id = replay.objects.line(
            f"{profile}_tau_from_S_merge",
            tau,
            mannheim_point_id,
            replay.batch_point_ids[keys[merged_name]],
        )
        build_targets_from_tau(replay, profile, tau, tau_id)
        return "merge_point"

    chosen = None
    for diagonal_name in ("Kp", "K"):
        state, lines, point = data[diagonal_name]
        if state == "finite" and point != mannheim_point:
            chosen = diagonal_name, lines, point
            break
    if chosen is None:
        raise AssertionError(f"{profile} 没有可与 S 相连的有限对角点")
    diagonal_name, lines, point = chosen
    assert lines is not None and point is not None
    pairs = (
        (("x", "y"), ("z", "w"))
        if diagonal_name == "Kp"
        else (("x", "w"), ("y", "z"))
    )
    line_ids = tuple(
        replay.objects.line(
            f"{profile}_{diagonal_name}_ordered_line_{index}",
            line,
            *(replay.batch_point_ids[keys[name]] for name in pair),
        )
        for index, (line, pair) in enumerate(
            zip(lines, pairs, strict=True), start=1
        )
    )
    point_id = replay.objects.point(
        f"{profile}_{diagonal_name}_ordered_point", point, *line_ids
    )
    tau = line_through(mannheim_point, point)
    if not same_line(tau, data["tau"]):
        raise AssertionError(f"{profile} 的有限对角点没有恢复接触弦")
    tau_id = replay.objects.line(
        f"{profile}_tau_from_S",
        tau,
        mannheim_point_id,
        point_id,
    )
    build_targets_from_tau(replay, profile, tau, tau_id)
    return diagonal_name


class OrderedBranchReplay(ThreeBlockReplay):
    def __init__(
        self,
        fixture_name: str,
        centers=None,
        radii=None,
        *,
        emit: bool = True,
    ) -> None:
        if centers is None or radii is None:
            centers, radii = FIXTURES[fixture_name]
        super().__init__(centers, radii)
        self.fixture_name = fixture_name
        self.emit = emit
        self.accounting: dict[str, int] = {}

    def measure(self, name: str, operation) -> None:
        before = paid_count(self)
        operation()
        self.accounting[name] = paid_count(self) - before

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

        self.measure("prefix", self.build_prefix)
        self.measure("P0_batch", lambda: [self.draw_batch(key) for key in p0_batch])

        seed_results: dict[str, tuple[Line, str]] = {}
        self.measure(
            "P0",
            lambda: seed_results.setdefault("P0", build_seed(self, "P0", data["P0"])),
        )
        self.measure(
            "other_batch", lambda: [self.draw_batch(key) for key in OTHER_BATCH]
        )
        self.measure(
            "P2",
            lambda: seed_results.setdefault("P2", build_seed(self, "P2", data["P2"])),
        )

        tau0, tau0_id = seed_results["P0"]
        tau2, tau2_id = seed_results["P2"]
        if same_line(tau0, tau2):
            raise AssertionError("P0、P2 的接触弦在 D8 中不能重合")
        mannheim_point = line_intersection(tau0, tau2)
        mannheim_point_id = self.objects.point(
            "Mannheim_S_ordered",
            mannheim_point,
            tau0_id,
            tau2_id,
        )

        tail_choices: dict[str, str] = {}
        for profile in TAIL_PROFILES:
            self.measure(
                profile,
                lambda profile=profile: tail_choices.setdefault(
                    profile,
                    build_tail(
                        self,
                        profile,
                        data[profile],
                        mannheim_point,
                        mannheim_point_id,
                    ),
                ),
            )
        return self.audit(data, tail_choices)

    def audit(self, data: dict, tail_choices: dict[str, str]) -> dict:
        if len(self.targets) != 8:
            raise AssertionError("统一分支程序没有得到八个有向根")
        if len({target["circle"] for target in self.targets.values()}) != 8:
            raise AssertionError("统一分支程序的八个输出圆不是两两不同")

        input_cost = sum(
            self.accounting[name] for name in ("prefix", "P0_batch", "other_batch")
        )
        if input_cost > 13:
            raise AssertionError("公共前缀与八条批量线超过 13 E")
        for profile in SEED_PROFILES:
            if self.accounting[profile] > 13:
                raise AssertionError(f"先行 {profile} 超过 13 E")
        for profile in TAIL_PROFILES:
            if self.accounting[profile] > 11:
                raise AssertionError(f"后行 {profile} 超过 11 E")

        graph = self.objects.graph
        ancestor_sets = {
            key: graph.paid_ancestors(target["output_id"])
            for key, target in self.targets.items()
        }
        union = frozenset().union(*ancestor_sets.values())
        trace = len(graph.paid_order)
        if trace > 61 or len(union) > 61:
            raise AssertionError("统一分支程序超过完整 61 E 上界")
        line_count = sum(graph.paid_kinds[node] == "line" for node in union)
        circle_count = sum(graph.paid_kinds[node] == "circle" for node in union)
        ext_keys = tuple(
            key for key in self.targets if key.startswith("P0:+++@")
        )
        first_ext = min(len(ancestor_sets[key]) for key in ext_keys)
        first_ext_trace = min(self.targets[key]["draw_index"] for key in ext_keys)
        if first_ext > 19:
            raise AssertionError(
                f"三重外切圆没有在 19 E 内完成：{first_ext} E"
            )

        report = {
            "fixture": self.fixture_name,
            "branches": {
                profile: data[profile]["kind"]
                for profile in SEED_PROFILES + TAIL_PROFILES
            },
            "tail_choices": tail_choices,
            "accounting": self.accounting,
            "trace": trace,
            "union_lines": line_count,
            "union_circles": circle_count,
            "all_targets": len(union),
            "first_ext": first_ext,
            "first_ext_trace": first_ext_trace,
            "non_ancestors": sorted(set(graph.paid_order) - set(union)),
            "aliases": len(self.objects.paid_aliases),
        }
        if self.emit:
            print("ordered_branch_replay", report)
        return report


def main() -> None:
    reports = {
        fixture_name: OrderedBranchReplay(fixture_name).run()
        for fixture_name in FIXTURES
    }
    for fixture_name, (expected_branches, expected_cost) in EXPECTED_RESULTS.items():
        report = reports[fixture_name]
        branches = tuple(
            report["branches"][profile]
            for profile in ("P0", "P2", "P1", "P3")
        )
        if branches != expected_branches or report["all_targets"] != expected_cost:
            raise AssertionError(f"{fixture_name} 的有序分支证书发生变化")
    saturated = reports["double_Kp_P2_P3"]
    if saturated["trace"] != 61 or saturated["all_targets"] != 61:
        raise AssertionError("双 K' 夹具没有同时饱和有序 61 E 台账")
    print(
        "ordered_branch_summary",
        {
            "fixtures": len(reports),
            "maximum_trace": max(report["trace"] for report in reports.values()),
            "maximum_all_targets": max(
                report["all_targets"] for report in reports.values()
            ),
            "saturated_fixture": "double_Kp_P2_P3",
        },
    )


if __name__ == "__main__":
    main()
