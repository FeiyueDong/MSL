# MSL 更新记录

本文件记录 MSL 每轮修改的时间、范围、主要内容与验证结果，按时间倒序排列。
新增修改时，在最上方追加一节即可。

---

## 2026-09-11 23:24 · 收尾修复（P1–P3）

**范围**：`msl::difference`、`msl::interp`、`msl::signal`（Butterworth）

**修改内容**

1. **Difference**
   - 矩阵非均匀 `central_gradient(mat, x, axis)` 内部点改用与向量版一致的
     三点加权中心差分公式（`axis=0`、`axis=1` 均修复）。
   - 分配型重载（`diff(y)`、`forward_gradient(y, ...)`）先校验长度，
     空输入不再发生 `size()-1` 无符号下溢。
   - `dx <= 0`、坐标非严格递增（重复/递减）统一抛
     `std::invalid_argument`；新增内部 helper 复用校验。
2. **Interpolation**
   - `find_interval()` 明确返回“区间编号”，钳制到 `[0, n-2]`；
     修复 `CubicSpline::derivative()` / `second_derivative()` 在
     `x == x_.back()` 处访问区间系数越界。
3. **Signal / Butterworth**
   - `poles_to_polynomial()` 检查共轭极点产生的虚部残差，超容差抛异常。
   - `normalize_gain()` 校验增益有限且大于容差，不再静默跳过归一化。

**测试**：新增非均匀矩阵梯度（解析场）、CubicSpline 端点一/二阶导、
Difference 非法输入（空/单点/`dx=0`/重复与递减坐标）等回归用例；
已确认这些用例在旧实现下失败。

**提交**：`1ba1724`、`0cc981d`、`755e6df`

**验证**：C++ 全量 8/8 通过（difference 151、interp 89 项检查）；
MATLAB `run_all` 8/8 PASS。本轮只新增校验与修正越界，未改变既有算法数值。

---

## 2026-09-11 19:55–21:57 · Correctness & API Cleanup（Stage A–E）

**范围**：Matrix、Signal、Integral、Interp、Polynomial、Difference，
以及测试与构建基础设施。

**修改内容**

1. **正确性修复（Stage A/B）**
   - LU 保留 partial pivoting：返回 `lu_result`，语义改为 `PA = LU`；
   - row-wise `ifft_rows` / `ifft_rows_real` 写回 `output`；
   - cumulative Simpson 奇数点改为三点抛物线半区间积分（均匀/非均匀一致）；
   - wide 矩阵 economy QR 尺寸改为 `k = min(m, n)`；
   - SVD 测试重建改用 `U·S·Vᵀ`（复数 `Vᴴ`）；
   - Butterworth bandpass 增益中心反 warp，bandstop 极点变换与零点修正。
2. **API / 语义清理（Stage C 及 P1）**
   - `get_row/get_column/submatrix` → `row_copy/column_copy/submatrix_copy`；
   - 移除 `adjoint()`，统一 `transpose()` / `conjugate_transpose()`；
   - 公开矩阵操作以 `std::invalid_argument` / `std::out_of_range`
     做运行时校验（替代 `assert`）；
   - 单帧 `psd/cpsd` → `power_spectrum/cross_power_spectrum`；
   - `signal.hpp` 聚合补齐 `filter/filtfilt`；
   - Interp 文档与空状态处理、Polynomial QR 三角回代与秩亏检测；
   - Difference `savgol_gradient` 修正为真正的 Savitzky–Golay 导数。
3. **MATLAB 交叉验证（Stage D）**
   - Matrix：LU / QR / SVD / eig / 广义 eig / pinv / inverse / det /
     truncated SVD；
   - Signal：行/列 FFT-IFFT、Welch PSD/CPSD、xcov、detrend、
     Butterworth / filter / filtfilt；
   - Integral：cumulative Simpson（解析原函数）；
   - Interp：非线性非均匀数据下的 linear/spline/pchip/makima/polynomial；
   - Difference：解析场梯度/二阶导/Laplacian/divergence/curl 与 S-G。
4. **基础设施（Stage E）**
   - 共享测试工具 `src/tests/test_utils.hpp` 与统一 `project_root()`；
   - MATLAB 统一 runner `matlab/validation/run_all.m`；
   - xmake header-only `msl` target、依赖路径清理；
   - GitHub Actions CI（Linux GCC / Linux Clang / Windows MSVC）。
5. **配套维护**
   - 新增接口变更说明 `docs/dev/qrest-msl-api-migration.md`；
   - `MSL.sln` 与测试 include 修正。

**破坏性接口变更**：详见 `docs/dev/qrest-msl-api-migration.md`。

**提交**：`dbb5690` … `65b6e22`（共 24 个），配套维护 `3d700bc`、`edb3458`

**验证**：C++ 全量 8/8 通过；MATLAB `run_all` 8/8 PASS，
其中 Butterworth 与 MATLAB `freqz` 误差 3.2e-14、`filtfilt` 5.2e-13、
Welch PSD 2.8e-17。
