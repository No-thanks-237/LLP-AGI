# 结构整理后的重跑与一致性

2026-10-25 整理 work 目录（work/common -> work/analysis、work/model/audit -> work/review_check，
见 logs/restructure.md）后，清空响应表缓存与逐质量扫描，按新结构**完整重跑**：
foresee（66 质量，B=6000；端点 4.55151/4.784 用 B=20000）、paper_mb（64+2）、exclusive（64+2）。
全部 0 失败。

## 与整理前对比（同一代码物理、同种子、同 nsamp）

用 work/analysis/compare_scans.py 比较整理前快照 temp/2hdm_pre_restructure/
与重跑后的 outputs/：

```
points: A=2640 B=2640 common=2640 onlyA=0 onlyB=0
 L=25  n=2617  max|ds/s|=0.000e+00  median=0.000e+00
 L=300 n=2617  max|ds/s|=0.000e+00  median=0.000e+00
boundary entries compared: 264, mismatches(>1e-6): 0
```

即 **s(mH,tanb) 与全部交点逐位相同**（最大相对差 0），说明整理只动了路径/目录，
没有改变任何物理结果。边界表、主图、显著度图、变体图均已按新数据重新生成。

## 主图风格

按用户要求改为 BC5 风格：只画曲线（实线=短寿交点，虚线=长寿边；蓝=25 fb^-1，红=300 fb^-1），
不加阈值竖线/区域阴影/MC 带；MC 误差与 EFT 有效性标记保留在
outputs/data/exclusion_boundary.csv、exclusion_crossings_mc.png 和 run_report.md 中。

## 其它清理

* 删除 work/magnus/jobs 下可重建的大中间产物（pik_out.tar.gz/pik_extract 等，约 1.1 GB）；
  样本 gz 的 sha256 已核对（minbias_pik.txt.gz = cc157eee…9de41f）。
* 删除 work/common、work/model 空目录；删除 __pycache__；删除被 final_analysis.py 取代的
  make_curves.py 与测试脚本。
