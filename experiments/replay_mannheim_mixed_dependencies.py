"""完整重放 Mannheim 简单合并与平行块同现的整数夹具。

第 8.13 节的无合并分析不能单独排除两类退化同时发生。本脚本覆盖
有界整数扫描实际找到的四个同现夹具，包括 ``K`` 与 ``K'`` 平行，
并按几何对象去重、按八目标联合祖先裁剪。
"""

from __future__ import annotations

from fractions import Fraction

from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
FIXTURES = {
    "merge_P1_parallel_P2_K": {
        "centers": ((F(0), F(0)), (F(12), F(0)), (F(13), F(6))),
        "radii": (F(4), F(2), F(1)),
        "merge": {"P1:merge:alpha2=alpha2_prime"},
        "parallel": {"P2:parallel:K"},
        "expected": (45, 16, 61),
    },
    "merge_P1_parallel_P3_K_alias": {
        "centers": ((F(0), F(0)), (F(9), F(0)), (F(-2), F(8))),
        "radii": (F(4), F(3), F(2)),
        "merge": {"P1:merge:alpha2=alpha2_prime"},
        "parallel": {"P3:parallel:K"},
        "expected": (44, 16, 60),
    },
    "merge_P1_parallel_P3_K_generic": {
        "centers": ((F(0), F(0)), (F(16), F(0)), (F(10), F(12))),
        "radii": (F(5), F(2), F(1)),
        "merge": {"P1:merge:alpha2=alpha2_prime"},
        "parallel": {"P3:parallel:K"},
        "expected": (46, 16, 62),
    },
    "merge_P2_parallel_P3_Kp": {
        "centers": ((F(0), F(0)), (F(18), F(0)), (F(13), F(7))),
        "radii": (F(5), F(3), F(1)),
        "merge": {"P2:merge:a2=a2_prime"},
        "parallel": {"P3:parallel:Kp"},
        "expected": (47, 16, 63),
    },
}

BATCH_KEYS = (
    "alphaA",
    "aB",
    "a1A",
    "alpha1B",
    "alphaB",
    "aA",
    "a1B",
    "alpha1A",
)


class MixedReplay(ThreeBlockReplay):
    def __init__(self, fixture_name: str) -> None:
        fixture = FIXTURES[fixture_name]
        super().__init__(fixture["centers"], fixture["radii"])
        self.fixture_name = fixture_name
        self.expected_merge = fixture["merge"]
        self.expected_parallel = fixture["parallel"]
        self.expected_count = fixture["expected"]

    def run(self) -> dict:
        if not is_d8(self.centers, self.radii):
            raise AssertionError(f"{self.fixture_name} 不属于 D8")
        events = analyze_fixture(self.centers, self.radii)
        merge_events = {event for event in events if ":merge:" in event}
        parallel_events = {
            event for event in events if ":parallel:K" in event
        }
        if merge_events != self.expected_merge:
            raise AssertionError(f"{self.fixture_name} 的合并型错误")
        if parallel_events != self.expected_parallel:
            raise AssertionError(f"{self.fixture_name} 的平行型错误")

        self.build_prefix()
        for key in BATCH_KEYS:
            self.draw_batch(key)
        for profile in ("P0", "P1", "P2", "P3"):
            options = {"allow_repeated_physical_signs": True}
            if any(event.startswith(f"{profile}:merge:") for event in events):
                self.build_merge_pair(profile, None, **options)
            elif f"{profile}:parallel:K" in events:
                self.build_parallel_pair(profile, "K", None, **options)
            elif f"{profile}:parallel:Kp" in events:
                self.build_parallel_pair(profile, "Kp", None, **options)
            else:
                self.build_regular_pair(profile, None, **options)
        return self.audit_mixed()

    def audit_mixed(self) -> dict:
        if len(self.targets) != 8:
            raise AssertionError("没有恢复八个目标根")
        if len({target["circle"] for target in self.targets.values()}) != 8:
            raise AssertionError("八个目标圆不是两两不同")

        graph = self.objects.graph
        ancestor_sets = {
            key: graph.paid_ancestors(target["output_id"])
            for key, target in self.targets.items()
        }
        union = frozenset().union(*ancestor_sets.values())
        union_lines = sum(graph.paid_kinds[node] == "line" for node in union)
        union_circles = sum(
            graph.paid_kinds[node] == "circle" for node in union
        )
        if (union_lines, union_circles, len(union)) != self.expected_count:
            raise AssertionError(f"{self.fixture_name} 的八目标成本发生变化")

        report = {
            "fixture": self.fixture_name,
            "trace": len(graph.paid_order),
            "union_lines": union_lines,
            "union_circles": union_circles,
            "union": len(union),
            "non_ancestors": sorted(set(graph.paid_order) - set(union)),
            "aliases": dict(sorted(self.objects.paid_aliases.items())),
        }
        print("mixed_global", report)
        return report


def main() -> None:
    reports = {
        fixture_name: MixedReplay(fixture_name).run()
        for fixture_name in FIXTURES
    }
    print(
        "maximum_pruned_all_targets",
        max(report["union"] for report in reports.values()),
    )


if __name__ == "__main__":
    main()
