# E49 Manim 演示动画

[下载 4K60 成片](https://raw.githubusercontent.com/lenovobenben/euclid-ccc-min/main/animations/e49/media/E49Progress_4k60.mp4)
（3840 × 2160，60 fps，约 1 分 39 秒）。仓库只跟踪这份正式播放版，
其余渲染中间文件和预览仍由 `.gitignore` 排除。

本目录把已验证的 Mannheim 正规 49 E 八解程序制作成中文 Manim 动画。
视觉风格与 Euclid-Min 的正十七边形动画一致：4K、60 fps、暗色全屏几何、
右上角单一 E 计数器、落笔前定位点高亮，以及输入圆、辅助对象和目标圆的
分层配色。

动画中的数值坐标不参与正确性判定。`geometry.json` 由正式程序对圆心
`(0,0)、(11,0)、(6,12)`、半径 `9/2、5/2、2` 的严格正规夹具做精确重放
后导出；Manim 只读取快照并负责显示。这个夹具增大了第三个输入圆，仍走
四个正规分支。正式轨迹为 39 条直线和 10 个圆，共 49 E，八个目标圆分别
在第 29、31、35、37、41、43、47、49 E 画出。

第三圆附近的接触弦阶段使用分段局部镜头。渲染时逐步检查两个作图定位点
都位于镜头安全区内。线条、定位点和标签保持稳定的屏幕尺寸；定位点使用
小实心点加细空心环，当前对象主要靠颜色高亮区分。镜头按缩放比例平滑
过渡，并给定位、画线和画圆留出停留时间。直线先在当前视野内匀速绘制，
完成后无缝恢复完整直线，避免局部放大时可见部分瞬间扫过。
镜头静止时不运行逐对象的样式更新，计数器直接更新整数，不进行数字形变。

这里的 49 E 是 `CCC-ALL-8 / gen` 的当前已验证上界，不是最优性证明或
世界纪录。项目独立的 `CCC-EXT-1 / gen = 18 E` 使用另一条程序，不能解释为
这段 49 E 动画的前 18 步。

## 文件

- `STORYBOARD.md`：表述口径、画面规则和 49 E 分镜；
- `export_geometry.py`：从正式精确重放导出动画快照；
- `geometry.json`：提交到仓库的确定性动画输入；
- `e49_progress.py`：Manim Community v0.21.0 场景；
- `Dockerfile`：固定 Manim 镜像并安装中文字体；
- `manim.cfg`：4K、60 fps 和输出目录。

## 导出几何数据

```powershell
docker run --rm `
  -v "${PWD}:/workspace" `
  -w /workspace `
  -e PYTHONPATH=/workspace/experiments `
  euclid-min-manim:0.21.0 `
  python animations/e49/export_geometry.py
```

## 构建渲染环境

若本机尚无参考项目使用的镜像：

```powershell
docker build -t euclid-min-manim:0.21.0 animations/e49
```

## 渲染

720p60 预览（保持与成片相同的帧率）：

```powershell
docker run --rm `
  -v "${PWD}:/manim" `
  -w /manim `
  euclid-min-manim:0.21.0 `
  manim --config_file animations/e49/manim.cfg -r 1280,720 --fps 60 `
  animations/e49/e49_progress.py E49Progress
```

4K 成片：

```powershell
docker run --rm `
  -v "${PWD}:/manim" `
  -w /manim `
  euclid-min-manim:0.21.0 `
  manim --config_file animations/e49/manim.cfg --fps 60 `
  animations/e49/e49_progress.py E49Progress
```

成片位于
`animations/e49/media/videos/e49_progress/2160p60/E49Progress.mp4`。

## 播放版编码

Manim 会拼接许多独立动画片段。最终播放版应对整片重新编码，统一为恒定
60 fps 和两秒关键帧间隔，并将 MP4 索引移到文件开头。下面的命令需要
本机安装 FFmpeg，从仓库根目录执行；它另存为较轻的 1080p60 文件：

```sh
ffmpeg -i animations/e49/media/videos/e49_progress/2160p60/E49Progress.mp4 -vf "fps=60,scale=1920:1080:flags=lanczos" -c:v libx264 -preset fast -crf 20 -maxrate 8M -bufsize 16M -g 120 -keyint_min 120 -sc_threshold 0 -bf 2 -refs 2 -pix_fmt yuv420p -fps_mode cfr -movflags +faststart -an animations/e49/media/E49Progress_1080p60.mp4
```

需要保留 4K 分辨率时，可另存 4K60 播放版：

```sh
ffmpeg -i animations/e49/media/videos/e49_progress/2160p60/E49Progress.mp4 -vf fps=60 -c:v libx264 -preset fast -crf 18 -maxrate 24M -bufsize 48M -g 120 -keyint_min 120 -sc_threshold 0 -bf 2 -refs 2 -pix_fmt yuv420p -fps_mode cfr -movflags +faststart -an animations/e49/media/E49Progress_4k60.mp4
```

播放卡顿时优先比较 1080p60 与 4K60 播放版，以区分动画节奏和播放器解码负担。

## 许可

Manim 源代码、导出程序和配置采用 [Apache License 2.0](../../LICENSE)；
几何导出数据、说明、分镜和视频成品采用 [CC BY 4.0](../../LICENSE-CONTENT)。
转载或改编视频时请按 [许可范围](../../LICENSE-SCOPE.md)保留署名并注明修改。
