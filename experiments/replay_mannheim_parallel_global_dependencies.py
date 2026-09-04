"""审计 Mannheim 平行退化与其余方向类的全局 E 复用。

脚本覆盖一个 ``K`` 无穷远块、一个 ``K'`` 无穷远块，以及两个 ``K'``
无穷远块同时出现的严格 ``D8`` 夹具。每个退化块使用第 8.8 节的两圆
一线修复，其余方向类使用正规五线块。
"""

from __future__ import annotations

from fractions import Fraction

from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
FIXTURES = {
    "single_K": (
        ((F(0), F(0)), (F(7), F(0)), (F(8), F(4))),
        (F(3), F(2), F(1)),
        {"P2:parallel:K"},
    ),
    "single_Kp": (
        ((F(0), F(0)), (F(28), F(0)), (F(13), F(27, 2))),
        (F(5), F(4), F(1)),
        {"P0:parallel:Kp"},
    ),
    "double_Kp_worst": (
        ((F(0), F(0)), (F(45), F(0)), (F(30), F(21))),
        (F(4), F(2), F(1)),
        {"P2:parallel:Kp", "P3:parallel:Kp"},
    ),
}

EXPECTED_COUNTS = {
    "single_K": (48, 12, 60),
    "single_Kp": (51, 12, 63),
    "double_Kp_worst": (51, 14, 65),
}


class ParallelGlobalReplay(ThreeBlockReplay):
    def __init__(self, fixture_name: str) -> None:
        centers, radii, expected_events = FIXTURES[fixture_name]
        super().__init__(centers, radii)
        self.fixture_name = fixture_name
        self.expected_events = expected_events

    def run(self) -> dict:
        if not is_d8(self.centers, self.radii):
            raise AssertionError(f"{self.fixture_name} 不属于 D8")
        events = analyze_fixture(self.centers, self.radii)
        parallel_events = {
            event for event in events if ":parallel:K" in event
        }
        if parallel_events != self.expected_events:
            raise AssertionError(
                f"{self.fixture_name} 的退化型错误：{sorted(events)}"
            )

        self.build_prefix()
        for key in (
            "alphaA",
            "aB",
            "a1A",
            "alpha1B",
            "alphaB",
            "aA",
            "a1B",
            "alpha1A",
        ):
            self.draw_batch(key)
        for profile in ("P0", "P1", "P2", "P3"):
            if f"{profile}:parallel:K" in events:
                self.build_parallel_pair(
                    profile,
                    "K",
                    None,
                    allow_repeated_physical_signs=True,
                )
            elif f"{profile}:parallel:Kp" in events:
                self.build_parallel_pair(
                    profile,
                    "Kp",
                    None,
                    allow_repeated_physical_signs=True,
                )
            else:
                self.build_regular_pair(
                    profile,
                    None,
                    allow_repeated_physical_signs=True,
                )
        return self.audit_parallel()

    def audit_parallel(self) -> dict:
        graph = self.objects.graph
        if len(self.targets) != 8:
            raise AssertionError("没有恢复八个目标圆")
        if len({target["circle"] for target in self.targets.values()}) != 8:
            raise AssertionError("八个目标圆不是两两不同")
        ancestor_sets = {
            sign: graph.paid_ancestors(target["output_id"])
            for sign, target in self.targets.items()
        }
        union = frozenset().union(*ancestor_sets.values())
        line_count = sum(kind == "line" for kind in graph.paid_kinds.values())
        circle_count = sum(
            kind == "circle" for kind in graph.paid_kinds.values()
        )
        union_lines = sum(graph.paid_kinds[node] == "line" for node in union)
        union_circles = sum(graph.paid_kinds[node] == "circle" for node in union)
        report = {
            "fixture": self.fixture_name,
            "trace": len(graph.paid_order),
            "trace_lines": line_count,
            "trace_circles": circle_count,
            "union": len(union),
            "union_lines": union_lines,
            "union_circles": union_circles,
            "non_ancestors": sorted(set(graph.paid_order) - set(union)),
            "aliases": dict(sorted(self.objects.paid_aliases.items())),
            "per_target": {
                sign: len(ancestor_sets[sign]) for sign in sorted(ancestor_sets)
            },
        }
        expected_lines, expected_circles, expected_total = EXPECTED_COUNTS[
            self.fixture_name
        ]
        if (union_lines, union_circles, len(union)) != (
            expected_lines,
            expected_circles,
            expected_total,
        ):
            raise AssertionError(
                f"{self.fixture_name} 的八目标联合成本发生变化"
            )
        if self.fixture_name == "double_Kp_worst" and self.objects.paid_aliases:
            raise AssertionError("双 K' 最坏夹具意外出现计费对象复用")
        print("parallel_global", report)
        return report


def main() -> None:
    reports = {
        fixture_name: ParallelGlobalReplay(fixture_name).run()
        for fixture_name in FIXTURES
    }
    print(
        "maximum_pruned_all_targets",
        max(report["union"] for report in reports.values()),
    )


if __name__ == "__main__":
    main()
