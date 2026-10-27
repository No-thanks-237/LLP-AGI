# 工作区整理记录（2026-10-25 23:05）

用户要求：把 run 目录整理成接近 BC5 v3.6 的结构，衰变与产生代码各占 work/ 下独立文件夹；
整理后重新扫描；主排除图去掉竖线、只画曲线。

## 结构调整（物理内容未改）

| 原路径 | 新路径 |
|---|---|
| work/common/{launch_scan,aggregate_scan,final_analysis,final_figure,signal_map}.py | work/analysis/ |
| work/common/{mass_grid,scan_mass_grid,mass_list_*,scan_points_*,widths_img_*} | work/analysis/grids/ |
| work/model/audit/ | work/review_check/ |
| work/common（空）; work/model（空） | 删除 |
| work/common/make_curves.py、_test_*、_synth_scan | 删除（被 final_analysis.py 取代） |
| work/magnus/jobs/{pik_out.tar.gz,pik_extract,smoke_*,hadronlib_*} | 删除（样本 gz 已核 sha256=cc157eee…9de41f；结果 json 保留） |

最终结构：
    inputs/  logs/  outputs/{data,figures}
    work/decay/      作者衰变程序（独立）
    work/production/ 解析产生模块（独立）
    work/response/   探测器内核 + 逐通道响应库
    work/analysis/   引擎、扫描、聚合、曲线、作图、验证（含 grids/ 与 cache/）
    work/magnus/     Pythia 样本、强子化库、作业记录
    work/review_check/ 独立解析审计脚本
与 BC5 v3.6 的 work/{widths->decay, production, response, analysis, magnus, review_check} 一一对应。

## 路径更新

* work/production/twohdm_production.py：mass_grid() 读 ../analysis/grids/mass_grid.txt；
  MODEL_TAG 1.0.0 -> 1.0.1（仅路径重构，物理未变）。
* work/decay/make_grids.py、work/analysis/launch_scan.py、work/magnus/{validate_samples,
  check_lib_interface,_finalise}.py、templates/script_hadronlib.sh 的网格路径同步更新；
  work/analysis/grids/ 已用 make_grids.py 重新生成（内容与移动前一致，sha256 见
  work/decay/grids_provenance.json）。
* work/production/validate_production.py 重跑，validation.json 记录了新路径与新的源码 sha256。
* README.md / DESIGN.md / logs/run_report.md 的结构与命令已更新；历史验证日志
  （decay_validation.md、production_validation.md、samples.md、model_audit.md、final_review.md）
  保留原样作为当时证据，其中出现的 work/common、work/model/audit 为整理前路径。

## 重跑与一致性

清空 work/analysis/{cache,scan_2hdm*} 后按新结构完整重跑（foresee/paper_mb/exclusive）并与
整理前的 outputs/data/signal_scan.json 逐点比较，结果见 logs/restructure_reproducibility.md。
