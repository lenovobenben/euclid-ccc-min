# 18 E 单目标动画

两段动画分别从输入三圆及其圆心开始，各画一个目标圆：

| 动画 | 目标 | 4K60 成片 |
|---|---|---|
| `E18External` | 与三个输入圆都外切 | [下载](media/E18External_4k60.mp4?raw=true) |
| `E18Internal` | 包住三个输入圆，并与它们分别内切 | [下载](media/E18Internal_4k60.mp4?raw=true) |

两条程序各为 **15 条直线 + 3 个圆 = 18 E**，共用前 14 E 的构造，
第 15–18 E 分别使用不同的第三圆接触点恢复圆心并画出目标圆。
每段视频都是一条独立程序；全内切动画不使用另一段已经画出的全外切圆。

本目录演示 Mannheim 的严格正规一般位置。两段动画的 18 E 都不是
最优性证明；这里没有把全外切单解的全域 18 E 覆盖结论移用到全内切单解。
一般参数推导与项目正式覆盖口径见 [Mannheim 文档](../../docs/MANNHEIM.md)。

## 几何与计数

圆心为 `(0,0)、(11,0)、(6,12)`，半径为 `9/2、4、2`。
沿用 49 E 动画的圆心布局，将第二圆半径调整为 4，使恢复圆心时用到的
外相似中心 `H = (1,24)` 更靠近主体；此夹具仍通过严格 `D8` 和正规性检查。

`export_e18_geometry.py` 调用正式 `ThreeBlockReplay` 的平行前缀、四条
批量线和 `P0` 正规分支，以选定目标圆的全部计费祖先裁剪并重新编号。
目标相切方向、包围关系和最终三个切点均在精确二次域中检查。
只有导给 Manim 的显示坐标转换为浮点数。

步骤安排：

1. 第 1–5 E：圆心连线和合法 4 E 平行线构造。
2. 第 6–9 E：取得第三圆上的四个批量点。
3. 第 10–14 E：得到 `K`、`K′`，画出接触弦；两个候选接触点免费。
4. 第 15–17 E：从本片选中的接触点恢复另一个切点及目标圆心。
5. 第 18 E：实际画出目标圆；结尾标出的三个切点是免费交点。

继承 49 E 动画的细线、小定位点、颜色高亮、镜头尺寸补偿和视野内匀速
画线方式。第 14 E 展示两个接触点并突出所选分支；最终镜头完整容纳
输入三圆和目标圆。

## 导出、检查与渲染

从仓库根目录执行：

```sh
python3 animations/e18/export_e18_geometry.py
python3 -m unittest discover -v
```

镜头测试需要 Manim 环境。沿用现有镜像：

```sh
docker build -t euclid-min-manim:0.21.0 animations/e49
docker run --rm -v "$PWD:/manim" -w /manim euclid-min-manim:0.21.0 python -m unittest discover -v
```

720p60 预览：

```sh
docker run --rm -v "$PWD:/manim" -w /manim euclid-min-manim:0.21.0 \
  manim --config_file animations/e18/manim.cfg -r 1280,720 --fps 60 \
  animations/e18/e18_progress.py E18External E18Internal
```

4K60 母片：

```sh
docker run --rm -v "$PWD:/manim" -w /manim euclid-min-manim:0.21.0 \
  manim --config_file animations/e18/manim.cfg --fps 60 \
  animations/e18/e18_progress.py E18External E18Internal
```

母片位于 `animations/e18/media/videos/e18_progress/2160p60/`。
用本机 FFmpeg 对整片重编码为恒定 60 fps、两秒关键帧间隔，并将 MP4
索引移到开头，再将这两份正式播放版提交到仓库：

```sh
for scene in E18External E18Internal; do
  ffmpeg -n -i "animations/e18/media/videos/e18_progress/2160p60/${scene}.mp4" \
    -vf fps=60 -c:v libx264 -preset fast -crf 18 -maxrate 24M -bufsize 48M \
    -g 120 -keyint_min 120 -sc_threshold 0 -bf 2 -refs 2 \
    -pix_fmt yuv420p -fps_mode cfr -movflags +faststart -an \
    "animations/e18/media/${scene}_4k60.mp4"
done
```

## 许可

源代码和配置采用 [Apache License 2.0](../../LICENSE)；文档、几何数据
和视频成品采用 [CC BY 4.0](../../LICENSE-CONTENT)。详见
[许可范围](../../LICENSE-SCOPE.md)。
