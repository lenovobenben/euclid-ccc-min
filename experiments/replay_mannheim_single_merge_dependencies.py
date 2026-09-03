"""完整审计一个仅含单个简单合并块的 Mannheim 八解程序。

这个夹具用于检查当前 66 E 最坏分支是否必然能与另外三个正规方向类
复用对象。几何和依赖图实现复用
``replay_mannheim_three_block_dependencies`` 中的精确对象注册器。
"""

from __future__ import annotations

from fractions import Fraction

from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
CENTERS = ((F(0), F(0)), (F(7), F(0)), (F(3), F(5)))
RADII = (F(3), F(2), F(1))


class SingleMergeReplay(ThreeBlockReplay):
    def __init__(self) -> None:
        super().__init__(CENTERS, RADII)

    def run(self) -> None:
        if not is_d8(self.centers, self.radii):
            raise AssertionError("单合并夹具不属于 D8")
        events = analyze_fixture(self.centers, self.radii)
        expected_events = {
            "P0:collapse:K=Kp",
            "P0:merge:alpha2=alpha2_prime",
        }
        if events != expected_events:
            raise AssertionError(f"单合并夹具的退化型错误：{sorted(events)}")

        self.build_prefix()
        # P0 的 y=w。优先用 alpha1B 表示合并点，使 H23_ext 已经可用；
        # 与它重合的 alphaA 稍后仍服务 P3。
        for key in ("aB", "a1A", "alpha1B"):
            self.draw_batch(key)
        self.build_merge_pair(
            "P0",
            ("+++", "---"),
            preferred_merged_key="alpha1B",
        )

        for key in ("alphaA", "alphaB", "aA", "a1B", "alpha1A"):
            self.draw_batch(key)
        self.build_regular_pair("P1", ("++-", "--+"))
        self.build_regular_pair("P2", ("+--", "-++"))
        self.build_regular_pair("P3", ("+-+", "-+-"))
        self.audit_single_merge()

    def audit_single_merge(self) -> None:
        expected_signs = {
            "+++",
            "---",
            "++-",
            "--+",
            "+--",
            "-++",
            "+-+",
            "-+-",
        }
        if set(self.targets) != expected_signs:
            raise AssertionError("八个物理切向符号没有全部出现")
        if len({target["circle"] for target in self.targets.values()}) != 8:
            raise AssertionError("八个目标圆不是两两不同")

        graph = self.objects.graph
        ancestor_sets = {
            sign: graph.paid_ancestors(target["output_id"])
            for sign, target in self.targets.items()
        }
        union = frozenset().union(*ancestor_sets.values())
        missing = set(graph.paid_order) - set(union)
        line_count = sum(kind == "line" for kind in graph.paid_kinds.values())
        circle_count = sum(
            kind == "circle" for kind in graph.paid_kinds.values()
        )
        union_lines = sum(graph.paid_kinds[node] == "line" for node in union)
        union_circles = sum(graph.paid_kinds[node] == "circle" for node in union)
        reuse = sum(len(items) for items in ancestor_sets.values()) - len(union)
        expected_aliases = {
            "P0_---_contact_line": "batch_alpha1B",
            "P0_---_radius_2": "ell",
            "P0_---_radius_3": "P0_merged_radius",
            "P2_Kp_line_1": "P1_Kp_line_1",
            "P3_Kp_line_2": "P1_Kp_line_2",
        }
        if self.objects.paid_aliases != expected_aliases:
            raise AssertionError("单合并程序的计费对象复用关系错误")
        if (line_count, circle_count, len(graph.paid_order)) != (48, 14, 62):
            raise AssertionError("单合并去重轨迹不是 48 线、14 圆、62 E")
        if self.targets["+++"]["draw_index"] != 19:
            raise AssertionError("三重外切圆没有在第 19 E 完成")
        if (union_lines, union_circles, len(union)) != (47, 14, 61):
            raise AssertionError("单合并八目标祖先不是 47 线、14 圆、61 E")
        if missing != {"batch_alphaA"}:
            raise AssertionError("单合并程序的非祖先对象错误")
        if reuse != 73:
            raise AssertionError("单合并八目标复用量错误")

        print(
            "score",
            {
                "trace_lines": line_count,
                "trace_circles": circle_count,
                "deduplicated_trace": len(graph.paid_order),
                "first_ext": self.targets["+++"]["draw_index"],
                "pruned_ext": len(ancestor_sets["+++"]),
                "pruned_all_targets": len(union),
                "pruned_lines": union_lines,
                "pruned_circles": union_circles,
            },
        )
        print(
            "dependencies",
            {
                "per_target": {
                    sign: len(ancestor_sets[sign]) for sign in sorted(ancestor_sets)
                },
                "union": len(union),
                "reuse": reuse,
                "non_ancestors": sorted(missing),
            },
        )
        print("paid_aliases", dict(sorted(expected_aliases.items())))


def main() -> None:
    replay = SingleMergeReplay()
    replay.run()


if __name__ == "__main__":
    main()
