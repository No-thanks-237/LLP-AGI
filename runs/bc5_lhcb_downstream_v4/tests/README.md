# tests — 单元与回归测试

```powershell
cd runs\bc5_lhcb_downstream_v4
python -m pytest tests -q        # 22 项，约 8 s
```

| 文件 | 覆盖 |
|---|---|
| `test_kinematics.py` | 二体阈值、四动量守恒、无质量四体 ⟨E_i/M⟩ = 1/4、有质量三体 ⟨s23⟩、boost 往返 |
| `test_lifetime.py` | A(z)=1 的解析衰变窗积分、标量/向量一致、出生点在窗外为 0 |
| `test_widths.py` | BR 和 = 1、阈值、ctau ∝ 1/θ²、轻子宽度闭式、数字化表已装载 |
| `test_production.py` | Boiarska tab:BR 逐道、B⁰/B± = 0.9275、K 道、四点道关闭、α(h→SS) |
| `test_curves.py` | 从已提交的 `signal_scan.json` 重算边界，与已提交 CSV **逐字节**比较 |
| `test_kernel.py` | 调用 `kernel_suite.py`，断言只有已知的 a2 冲突失败 |
| `kernel_suite.py` | 33 项探测器内核自检（独立参考实现 + 边界/符号/确定性） |

**已知 a2**：`x = y = 0` 的径迹正好在束流轴上，内核按 `holeRadius = 0.09 m` 正确判为
`BeamHole`；测试里的旧期望是"接受"，属于规格冲突，不是内核缺陷。`kernel_suite.py`
把它记为已知失败，退出码仍为 0。

测试不读取母样本（除 `test_curves` 读 `outputs/data/signal_scan.json`），可在无 Magnus
环境的机器上跑。
