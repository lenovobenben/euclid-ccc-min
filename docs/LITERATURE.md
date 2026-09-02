# 一般 CCC 最小步数文献台账（第一轮）

## 1. 本轮结论

检索日期：2026-09-03。

针对本项目的精确问题——三个两两外离、半径不同、圆心不共线的给定圆，只画出唯一一个与三圆均外切的圆，使用无任意点的直尺与塌圆规，按 `Line + Circle` 计 E——本轮**没有找到公开发表的最小已知 E 数，也没有找到最优性证明**。

当前能可靠支持的结论分为三层：

1. Lemoine 1892 对一般 CCC 的四种经典方案作了逐项计数，但目标是画出全部八个阿波罗尼乌斯圆，工具模型还允许任意点和异地搬运圆规开度。
2. 在 Lemoine 比较的四种全部八解方案中，Mannheim 方案最少：54 条直线和 10 个圆，共 **64 个已画对象**。这是历史比较中的最小值，不是本项目的 64 E 结果，也没有最优性证明。
3. 从 Lemoine 给出的 Mannheim 明细裁剪到同号分支中的一个目标圆，可得到 **21 个已画对象**的候选前缀。这个 21 是本项目下一阶段最值得严格化的线索，但它不是文献原报数字，也不是合法 E 分数；其中的平行线宏隐藏了任意点和距离搬运。

另有与本项目计数规则高度相近的 7 E 和 11 E 结果，但它们只适用于三个输入圆已经两两相切、且免费给出切点和若干圆心连线的特殊题，不能用作一般 CCC 的记录。

## 2. 比较口径

### 2.1 本项目的目标口径

本台账以 [CCC 研究设计](DESIGN.md) 为准：

- 三个输入圆的圆心和圆本身免费；
- 圆周上不免费给点；
- 不允许任意点；
- 每一步必须是 `Line(P,Q)` 或 `Circle(center=P, through=Q)`；
- 圆规为塌圆规，不能把 `|AB|` 直接搬到第三点 `C`；
- 只要求实际画出唯一的三重外切圆；
- 另外七个含内切关系的解不计入目标。

### 2.2 Lemoine 计数与 E 的数值关系

Lemoine 把一份构造写成

\[
\operatorname{Op}=
M_1R_1+M_2R_2+N_1C_1+N_2C_2+N_3C_3.
\]

其中：

- `R1`、`C1`、`C2` 统计放置直尺或圆规端点等动作；
- `R2` 是实际画一条直线；
- `C3` 是实际画一个圆；
- 简洁度 `S = M1 + M2 + N1 + N2 + N3`；
- 精确度 `X = M1 + N1 + N2`。

所以

\[
S-X=M_2+N_3
=\#\text{直线}+\#\text{圆}.
\]

这与 E 的**数值目标函数**一致，但两套工具的合法程序集合不一致。Lemoine 可以：

- 任选平面点或圆上点；
- 任取“足够大”的辅助半径；
- 夹取 `|AB|` 后，以另一个点为圆心画该半径的圆。

因此本文把 Lemoine 的换算值称为“已画对象数”或 `E_trace`。只有全部宏被改写成项目合法程序后，才能称为 E。

## 3. Lemoine 1892：一般 CCC 的原始计数

原始来源：Émile Lemoine，*Application d'une méthode d'évaluation de la simplicité des constructions à la comparaison de quelques solutions du problème d'Apollonius*，*Nouvelles annales de mathématiques*，第 3 辑第 11 卷，1892，453–474 页。

- [Numdam 书目页](https://www.numdam.org/item/NAM_1892_3_11__453_1/)
- [Numdam PDF](https://www.numdam.org/item/NAM_1892_3_11__453_1.pdf)

论文明确按一般半径不等的三个圆讨论，并在总表中比较画出全部八解的四种路线。

| 方案 | Lemoine 符号 | 简洁度 | 精确度 | 直线 | 圆 | `E_trace` | 输出 |
|---|---:|---:|---:|---:|---:|---:|---|
| Mannheim | `108R1 + 54R2 + 20C1 + 10C3` | 192 | 128 | 54 | 10 | **64** | 全部八解 |
| Fouché | `112R1 + 56R2 + 53C1 + 26C3` | 247 | 165 | 56 | 26 | **82** | 全部八解 |
| Viète | `52R1 + 26R2 + 98C1 + C2 + 58C3` | 235 | 151 | 26 | 58 | **84** | 全部八解 |
| Bobillier–Gergonne | `120R1 + 60R2 + 104C1 + 72C3` | 356 | 224 | 60 | 72 | **132** | 全部八解 |

Lemoine 的结论是 Mannheim 在这四种方案中明显最好。他没有声称 64 是所有方案中的下界。

现代史学复核也抄录了同一张表：Sylviane R. Schwer，*La géométrographie*，表 12.1。

- [Publimath PDF](https://bibnum.publimath.fr/IWH/IWH22013.pdf)

### 3.1 Fouché 首对解：37 个已画对象

Lemoine 在第 469–470 页说明：Fouché 构造先得到一对圆，其中一个与三圆均外切，另一个与三圆均内切；再重复三个后缀，补齐其余六解。

全部八解的符号为

\[
T=112R_1+56R_2+53C_1+26C_3.
\]

论文给出每个后续解对的重复后缀

\[
B=20R_1+10R_2+11C_1+5C_3.
\]

因此第一对解的前缀为

\[
T-3B=52R_1+26R_2+20C_1+11C_3,
\]

即 26 条直线加 11 个圆，共 **37 个已画对象**。

这是从论文总式直接相减得到的可复核数字。它仍然：

- 同时画出外切和内切两个圆；
- 预先构造了服务于全部八解的对象；
- 使用 Lemoine 的任意点和圆规搬运能力。

所以 37 既不是单目标最短数，也不是项目 E 上界。

### 3.2 Mannheim 单目标裁剪候选：21 个已画对象

Mannheim 原始构造见 1885 年的 *Note de géométrie*，105–109 页；一般三圆构造在 107–109 页。

- [Numdam 书目页](https://www.numdam.org/item/NAM_1885_3_4__105_0/)
- [Numdam PDF](https://www.numdam.org/item/NAM_1885_3_4__105_0.pdf)

Lemoine 1892 第 470–472 页把 Mannheim 方案完全计数。若只保留同号分支的一对接触点，并只画其中一个目标圆，原文明细给出：

| 阶段 | 已画对象 |
|---|---:|
| 画圆心连线 `oo1` | 1 条直线 |
| 过第三圆心作平行线 | 1 条直线 + 2 个圆 = 3 |
| 将四个已得圆周点分别连到 `A,B` | 8 条直线 |
| 为所选解对得到第三圆上的两个接触点 | 5 条直线 |
| 由一个接触点恢复另外接触点和目标圆心 | 3 条直线 |
| 实际画目标圆 | 1 个圆 |
| **合计** | **21** |

最后两行也可由论文的八个输出圆总式核对：

\[
48R_1+24R_2+16C_1+8C_3,
\]

平均每个输出圆需要 3 条定心直线和 1 个目标圆。

这个 21 是本项目的**首选严格化候选**，不是当前记录，原因包括：

1. Lemoine 的平行线构造从任意点 `λ` 开始；
2. 它先夹取 `|OO''|`，再以另一个点为圆心画同半径圆，属于本项目禁止的距离直接搬运；
3. 同号分支中哪一个交点对应三重外切，需要写成一般输入上稳定的精确分支条件；
4. 还没有证明这 21 个对象全部属于单目标依赖祖先，也没有证明其它批量步骤不能继续删减；
5. 替换非法宏后，项目 E 可能高于 21；若上下文中已有对象能承担替换，也可能出现新的共享优化。

因此当前只能记为：

\[
E_{\text{trace,candidate}}=21,
\qquad
E_{\text{project}}\ \text{未知}.
\]

## 4. 现代来源

### 4.1 Gergonne 一般方案

下列来源给出或解释一般三圆的 Gergonne 路线：相似中心、四条相似轴、根心、极点和切点。

- David Gisch、Jason M. Ribando，*Apollonius' Problem: A Study of Solutions and Their Connections*，2004，[DOI](https://doi.org/10.33697/ajur.2004.010)，[PDF](https://www.ajuronline.org/uploads/Volume%203/Issue%201/31D-GischArt.pdf)
- [Cut-the-Knot：Gergonne CCC](https://www.cut-the-knot.org/Curriculum/Geometry/GeoGebra/CCC_Gergonne.shtml)
- Tom Davis，*Geometry*，[PDF](http://www.geometer.org/geometer/geometry.pdf)

这些来源足以恢复宏观几何路线，但没有按 `Line + Circle` 展开，也没有单独统计三重外切目标。Cut-the-Knot 页面上的“五步”全部是高级宏，不能作为 5 E。

Gisch–Ribando 还明确区分：

- Gergonne 处理一般位置三圆并产生八解；
- Eppstein、Soddy 等短方案处理三个圆已经两两相切的特殊题。

### 4.2 Viète 半径平移路线

Tom Davis 第 9.15.5 节描述：把最小圆缩成圆心点，同时对另两个圆的半径作相应加减，再解 PCC，最后恢复半径。

该路线适用于一般 CCC，但没有 E 计数，并隐藏：

- 半径加减；
- 塌圆规下的距离搬运；
- PCC 的反演或公切线构造；
- 逆变换后的目标圆恢复。

它仍是后续形式化所需的第二条独立基线。

### 4.3 现代反演方案

Sabihi 2019：

- [arXiv 摘要](https://arxiv.org/abs/1906.00068)
- [HTML](https://arxiv.org/html/1906.00068v1)

其 4.1–4.10 是十个宏步骤，包含任意辅助圆、根轴、反演圆像和长度乘积；没有基础作图数。

Azizov–Litvinov 2025：

- [arXiv HTML](https://arxiv.org/html/2409.17153)

该文解决更广的三圆定角相交问题，切圆是零角特例。它使用把两个不交圆反演成同心圆的路线，但仍含任意点、平行线、公切线、根轴、任意反演半径以及乘除和平方根宏，没有 E 统计。

两者均可作为反演 gadget 的来源，不能作为当前最短记录。

### 4.4 三圆已经两两相切的特殊题

Baragar–Kontorovich 2020：

- Arthur Baragar、Alex Kontorovich，*Efficiently Constructing Tangent Circles*，[PDF](https://par.nsf.gov/servlets/purl/10143905)

该文的计数模型与项目很接近：画一条直线或一个以已有点为圆心、过已有点的圆，各算一步，交点免费，使用塌圆规。

但输入免费给出：

- 三个已经两两相切的圆；
- 三个圆心；
- 三个已知切点；
- 若干圆心连线。

在这个特殊输入下：

- Baragar–Kontorovich 构造一个目标圆用 7 E；
- 论文引用的 Eppstein 旧方案用 11 E。

论文没有证明 7 E 最优。更关键的是，本项目的三个输入圆两两外离，既无切点也无免费圆心连线，所以 7 和 11 都不能进入一般 CCC 排名。

## 5. 未提供一般 CCC 记录的数据库

### 5.1 Labelle 构造数据库

Paul-Émile Labelle 的旧 Java 构造数据库采用与 E 相近的“每条线或圆一步”模型。通过 Internet Archive 恢复页面和 class 文件后，问题列表只有 37 个传统点输入问题和 4 个不可作问题，没有 CCC。

- [Internet Archive CDX 索引](https://web.archive.org/cdx/search/cdx?url=www.cs.mcgill.ca/~sqrt/cons/*)

因此该数据库不能给出 CCC 历史分数。

## 6. 需要排除的网络数字

### 6.1 “154 operations”

英文 Wikipedia 的 Émile Lemoine 条目声称 Lemoine 把阿波罗尼乌斯问题降到 154 次操作：

- [Émile Lemoine — Wikipedia](https://en.wikipedia.org/wiki/%C3%89mile_Lemoine)

该条目没有链接到一份含 154 明细的原始构造。154 即便属实，也指 Lemoine 的五类手工动作总数，不是直线加圆；输出一个、两个还是八个圆也不明确。本轮不采用。

### 6.2 “199 operations”

MacTutor 声称通常方案超过 400 次 Lemoine 操作，而 Lemoine 降到 199：

- [MacTutor：Émile Lemoine](https://mathshistory.st-andrews.ac.uk/Biographies/Lemoine/)

该页面同样没有给出可重放的 199 操作程序。Lemoine 1892 对自己 1888 年比较的回顾给出的却是 Viète 335、Bobillier–Gergonne 500，随后 1892 年表中最好的 Mannheim 是 192。不同数字很可能混合了不同方案、输出范围或修订年代。

在找到原始逐项符号和明确输出目标之前，154、199、400 都不能换算成项目 E。

## 7. 第一轮检索式与筛选

本轮使用的主要检索式：

```text
"Apollonius problem" straightedge compass minimum steps
"problème d'Apollonius" "nombre" "construction" "simplicité"
Mannheim 1885 Apollonius construction Nouvelles annales 108
Lemoine Apollonius "154" operations
Emile Lemoine Géométrographie 1902 pdf Apollonius
```

纳入标准：

- 能定位原始论文、稳定 PDF 或足以复原路线的技术资料；
- 明确输入是一般 CCC 还是三圆两两相切的特殊题；
- 明确输出一个解、一对解还是全部八解；
- 对数字说明原始指标及其与 E 的关系。

排除标准：

- 把相似中心、根轴、反演或公切线各算“一步”的网页宏计数；
- 没有构造明细的二手数字；
- 三圆已两两相切而未注明额外免费对象的短方案；
- 只给代数解或数值算法，没有尺规程序的来源。

## 8. 当前候选排序

| 地位 | 数字 | 适用范围 | 当前判断 |
|---|---:|---|---|
| 一般 CCC、单个三重外切圆、项目规则 | 未知 | 本项目精确问题 | 尚无公开记录 |
| Mannheim 单目标裁剪 | 21 个已画对象 | 一般 CCC 的同号分支 | 最优先严格化；含非法宏，非 E |
| Fouché 首对解 | 37 个已画对象 | 一般 CCC，外切圆 + 内切圆 | 原文可换算；非单目标，非 E |
| Mannheim 全部八解 | 64 个已画对象 | 一般 CCC | Lemoine 四方案中最好；非 E，无最优性证明 |
| Baragar–Kontorovich | 7 E | 三圆已两两相切并免费给切点等 | 同规则特殊题，不可比较 |

下一阶段应先严格展开 Mannheim 的 21 对象候选。第一处硬缺口是“过第三圆心作平行线”的 3 对象宏：需要在无任意点、塌圆规的当前状态中找到合法替代并计算真实边际 E。完成后再与 Viète/反演路线按同一 profile 比较。
