# E49 Manim 演示动画

本目录把已验证的 Mannheim 正规 49 E 八解程序制作成中文 Manim 动画。
视觉风格与 Euclid-Min 的正十七边形动画一致：4K、30 fps、暗色全屏几何、
右上角单一 E 计数器、落笔前定位点高亮，以及输入圆、辅助对象和目标圆的
分层配色。

动画中的数值坐标不参与正确性判定。`geometry.json` 由
`KpCenterLocusReplay("regular")` 对严格正规夹具做精确重放后导出；Manim
只读取快照并负责显示。正式轨迹为 39 条直线和 10 个圆，共 49 E，八个
目标圆分别在第 29、31、35、37、41、43、47、49 E 画出。

这里的 49 E 是 `CCC-ALL-8 / gen` 的当前已验证上界，不是最优性证明或
世界纪录。项目独立的 `CCC-EXT-1 / gen = 18 E` 使用另一条程序，不能解释为
这段 49 E 动画的前 18 步。

## 文件

- `STORYBOARD.md`：表述口径、画面规则和 49 E 分镜；
- `export_geometry.py`：从正式精确重放导出动画快照；
- `geometry.json`：提交到仓库的确定性动画输入；
- `e49_progress.py`：Manim Community v0.21.0 场景；
- `Dockerfile`：固定 Manim 镜像并安装中文字体；
- `manim.cfg`：4K、30 fps 和输出目录。

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

快速预览：

```powershell
docker run --rm `
  -v "${PWD}:/manim" `
  -w /manim `
  euclid-min-manim:0.21.0 `
  manim --config_file animations/e49/manim.cfg -ql `
  animations/e49/e49_progress.py E49Progress
```

4K 成片：

```powershell
docker run --rm `
  -v "${PWD}:/manim" `
  -w /manim `
  euclid-min-manim:0.21.0 `
  manim --config_file animations/e49/manim.cfg --fps 30 `
  animations/e49/e49_progress.py E49Progress
```

成片位于
`animations/e49/media/videos/e49_progress/2160p30/E49Progress.mp4`。
