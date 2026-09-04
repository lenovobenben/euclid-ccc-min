"""精确审计 Mannheim 四个正规方向类的全局对象复用。

固定使用 ``replay_mannheim_fixed.py`` 的严格 ``D8`` 校准实例，但通过
几何对象注册器执行真正的 E 去重。四个五线块的八条 ``K'`` 定义弦只
有四条不同直线，因此合法正规程序为 61 E，而不是逐块相加所得的 65 E。
"""

from __future__ import annotations

from fractions import Fraction

from replay_mannheim_three_block_dependencies import ThreeBlockReplay
from scan_mannheim_degeneracies import analyze_fixture, is_d8


F = Fraction
CENTERS = ((F(0), F(0)), (F(13), F(0)), (F(4), F(15)))
RADII = (F(4), F(2), F(1))


class RegularReplay(ThreeBlockReplay):
    def __init__(self) -> None:
        super().__init__(CENTERS, RADII)

    def run(self) -> None:
        if not is_d8(self.centers, self.radii):
            raise AssertionError("正规校准实例不属于 D8")
        if analyze_fixture(self.centers, self.radii):
            raise AssertionError("正规校准实例含有五线块退化")

        self.build_prefix()
        for key in ("alphaA", "aB", "a1A", "alpha1B"):
            self.draw_batch(key)
        self.build_regular_pair("P0", ("+++", "---"))

        for key in ("alphaB", "aA", "a1B", "alpha1A"):
            self.draw_batch(key)
        self.build_regular_pair("P2", ("+--", "-++"))
        self.build_regular_pair("P1", ("++-", "--+"))
        self.build_regular_pair("P3", ("+-+", "-+-"))
        self.audit_regular()

    def audit_regular(self) -> None:
        graph = self.objects.graph
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

        expected_aliases = {
            "P2_Kp_line_2": "P0_Kp_line_2",
            "P1_Kp_line_1": "P2_Kp_line_1",
            "P3_Kp_line_1": "P0_Kp_line_1",
            "P3_Kp_line_2": "P1_Kp_line_2",
        }
        if self.objects.paid_aliases != expected_aliases:
            raise AssertionError("四个正规五线块的公共弦关系错误")

        line_count = sum(kind == "line" for kind in graph.paid_kinds.values())
        circle_count = sum(
            kind == "circle" for kind in graph.paid_kinds.values()
        )
        if (line_count, circle_count, len(graph.paid_order)) != (51, 10, 61):
            raise AssertionError("正规八解程序不是 51 线、10 圆、61 E")
        if self.targets["+++"]["draw_index"] != 18:
            raise AssertionError("三重外切圆没有在第 18 E 完成")

        ancestor_sets = {
            sign: graph.paid_ancestors(target["output_id"])
            for sign, target in self.targets.items()
        }
        if {len(items) for items in ancestor_sets.values()} != {18}:
            raise AssertionError("正规实例每个目标应有 18 个计费祖先")
        union = frozenset().union(*ancestor_sets.values())
        if len(union) != 61 or set(union) != set(graph.paid_order):
            raise AssertionError("正规八目标联合祖先错误")
        reuse = sum(len(items) for items in ancestor_sets.values()) - len(union)
        if reuse != 83:
            raise AssertionError("正规八目标复用量错误")

        print(
            "score",
            {
                "lines": line_count,
                "circles": circle_count,
                "first_ext": self.targets["+++"]["draw_index"],
                "all_targets": len(union),
                "saved_shared_chords": 4,
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
            },
        )
        print("paid_aliases", expected_aliases)


def main() -> None:
    replay = RegularReplay()
    replay.run()


if __name__ == "__main__":
    main()
