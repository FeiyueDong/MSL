# MSL API 变更说明（qREST 对接用）

本文档汇总本轮 Correctness & API Cleanup 中 **会影响 qREST 调用的接口变更**，
用于在 qREST 中检索并修改被破坏的调用点。

- 变更范围：Stage A–E（提交 `dbb5690` … `65b6e22`，分支 `opencode-fix`）
- MSL 无兼容层，旧名称/旧签名直接移除或替换
- 文档中的“行为/数值变化”不破坏编译，但会改变结果，需重设期望值

---

## 1. 编译期破坏性变更（必须改代码）

### 1.1 `matrix::lu()` 返回类型改变

旧接口丢失 partial pivoting 的 permutation，数学上错误（曾声明 `A = LU`）。

| | 旧 | 新 |
|---|---|---|
| 返回类型 | `std::array<real_matrix_owned, 2>` | `msl::matrix::lu_result` |
| 复数返回 | `std::array<complex_matrix_owned, 2>` | `msl::matrix::complex_lu_result` |
| 语义 | `A = L * U`（错误） | `P * A = L * U` |

```cpp
// 旧
auto [L, U] = msl::matrix::lu(A);

// 新
auto result = msl::matrix::lu(A);
result.L;            // unit lower triangular
result.U;            // upper triangular
result.permutation;  // std::vector<size_t>，行置换
```

`permutation` 约定：`(P*A)(i,j) = A(permutation[i], j)`。
求解 `A*x = b` 时需先按该向量重排 `b`（`LU*x = P*b`）。

### 1.2 Matrix 拷贝访问器改名

旧名称返回 owning copy，容易误当作 view 使用（曾直接导致 row-IFFT bug）。

| 旧 | 新 |
|---|---|
| `A.get_row(i)` | `A.row_copy(i)` |
| `A.get_column(j)` | `A.column_copy(j)` |
| `A.submatrix(r0, r1, c0, c1)` | `A.submatrix_copy(r0, r1, c0, c1)` |

`A.column(j)` 保持不变，返回 `std::span`（零拷贝 view）。

### 1.3 `matrix::adjoint()` 移除

旧实现是 Eigen 的 conjugate transpose，却曾被文档描述为 classical adjoint
（adjugate），且强制要求方阵。

| 场景 | 旧 | 新 |
|---|---|---|
| 实矩阵 | `matrix::adjoint(A)` | `matrix::transpose(A)` |
| 复矩阵 | `matrix::adjoint(A)` | `matrix::conjugate_transpose(A)` |

`conjugate_transpose` 对矩形矩阵同样适用；如需真正的 adjugate，请另行实现。

### 1.4 单帧 `psd` / `cpsd` 改名

这两个函数没有 PSD density 归一化，因此按语义改名；Welch 版本名称不变。

| 旧 | 新 |
|---|---|
| `signal::psd(x, nfft)` | `signal::power_spectrum(x, nfft)` |
| `signal::psd(x, out, nfft)` | `signal::power_spectrum(x, out, nfft)` |
| `signal::cpsd(x, y, nfft)` | `signal::cross_power_spectrum(x, y, nfft)` |
| `signal::cpsd(x, y, out, nfft)` | `signal::cross_power_spectrum(x, y, out, nfft)` |
| `signal::psd_welch` / `signal::cpsd_welch` | 不变 |

注意：`msl::signal` 中另有来自 `fft.hpp` 的
`power_spectrum(std::span<const std::complex<double>>)`（对 FFT bin 求 `|X|²`），
按参数类型重载，不受影响。

### 1.5 检索命令

```powershell
# 在 qREST 仓库根目录执行
rg -n -g "*.hpp" -g "*.cpp" "get_row|get_column|\.submatrix\(" .
rg -n -g "*.hpp" -g "*.cpp" "matrix::lu\(|matrix::adjoint\(" .
rg -n -g "*.hpp" -g "*.cpp" "signal::psd\(|signal::cpsd\(" .
```

注意 `\blu\(` 会误命中 Eigen 的 `.lu()`，请使用 `matrix::lu(`。
`submatrix(` 也可能命中非矩阵代码，按上下文判断。

---

## 2. 行为/数值变化（签名不变，需重新基线）

以下接口调用可以原样编译，但数值或异常行为会变。

### 2.1 `integral::cumsimpson()` 奇数点公式修正

| 位置 | 旧 | 新 |
|---|---|---|
| 均匀奇数点 | `I[i-1] = (I[i-2] + I[i]) / 2` | `I[i-1] = I[i-2] + dx*(5y[i-2] + 8y[i-1] - y[i])/12` |
| 非均匀奇数点 | 前一段梯形 | 同一三点抛物线在左半区间积分 |
| 偶数长度末点 | 梯形 | 不变 |

影响：所有偶数噪声/多项式数据的奇数累计点数值都会改变（这是修复项）。
qREST 若有基于旧值的 expected 数据需更新。

### 2.2 Butterworth bandpass / bandstop 系数修正

函数名与签名（`ButterworthFilter::bandpass/bandstop`、
`butterworth_bandpass_design/bandstop_design` 等）均未变，但系数改变：

- **Bandpass**：增益归一化中心由 `(f1+f2)/2` 改为反 warp 中心
  `f0 = (2/π)·atan(Ω0/2)`，其中 `Ω0 = sqrt(tan(π f1/2)·tan(π f2/2))`。
- **Bandstop**：极点变换由错误的 BP 变换改为 `s_L = B·s/(s²+Ω0²)`
  （即 `α = B/(2p)`），零点也落在反 warp 后的真实 notch 频率。

这两个修复会明显改变 `[0.004, 0.8]` 等频带的系数与滤波结果
（RR 位移的 `RelativeL2 ≈ 5.5e-3` 应随之消失），需重新基线。

### 2.3 `difference::savgol_gradient()` 现在是真正的 Savitzky–Golay

旧实现只是“局部差分平均”，不是 S-G。新实现为局部多项式最小二乘导数
（内部闭式权重 + 与 MATLAB 一致的边界拟合），数值结果会改变。
函数签名不变：

```cpp
auto g = msl::difference::savgol_gradient(y, window_size, poly_order, dx);
```

### 2.4 `polynomial::polyfit` / `Polynomial::fit`

- 秩亏（重复 x / 线性相关采样）现在抛 `std::runtime_error`；
  旧实现会通过显式 `inverse(R)` 给出无意义结果。
- 采样数不足仍为 `std::invalid_argument`。
- 良态问题使用 QR + 三角回代，系数与旧值差异通常在 1e-12 量级。

### 2.5 `matrix::qr()` wide 矩阵尺寸修正

- 旧 thin QR：Q `m×n`、R `n×n`（wide 时越界/错误）。
- 新 thin QR：`k = min(m,n)`，Q `m×k`、R `k×n`。
- tall / square 维度不变；**wide 矩阵返回值尺寸会变**，并按 `k` 截断。
- full QR（`qr(A, true)`）不变：Q `m×m`、R `m×n`。

若 qREST 只在多项式拟合（tall Vandermonde）中使用 `qr`，则无影响。

### 2.6 运行时校验：`assert` 改为异常

Release 下 `assert` 会消失并变成 UB；现在公开操作主动抛异常。

- `std::invalid_argument`：shape/size 不匹配
  - 矩阵构造（span / initializer_list 尺寸不符）
  - `operator+`、`operator-`、`operator*`、`+=`、`-=`、
    `msl::matrix::eigen_interface::matmul/matvec/solve`
  - `matrix::trace` / 成员 `trace()`（非方阵）
  - view 的 `copy_from`、`+=`、`-=`、`subview`（非整列子视图）
- `std::out_of_range`：索引/范围越界
  - `row_copy`、`column_copy`、`submatrix_copy`、`column`、view `subview`

`operator()(i,j)` 仍然不做检查（与容器一致）。
另外 `column(j)` 不再是 `noexcept`（新增边界检查）。

如果 qREST 在调用前已自行校验维度，通常无需改动；
若之前依赖 release 下静默越界，现在会收到异常，请补捕获或修正输入。

### 2.7 插值器空状态

对 default-constructed / 未 `set_data()` 的插值器调用
`operator()`、`in_range()`、`range()` 现在抛 `std::runtime_error`，
不再对空 vector 调用 `front()/back()`（旧行为是 UB）。

CubicSpline 默认边界仍为 Not-a-knot，Akima 默认仍为 modified Akima
（MAKIMA），与 MATLAB `interp1(...,'spline')` / `makima` 对齐；
仅文档被修正，调用语义未变。

---

## 3. 未受影响的接口（无需检索）

以下接口本轮未改签名、未改语义：

- `matrix::svd`、`eig`、`pinv`、`truncated_svd`、`inverse`、`determinant`
- `matrix::transpose`、`conjugate_transpose`
- `signal::fft/ifft/ifft_real`、`fft_rows/ifft_rows/ifft_rows_real`、`fft_columns/ifft_columns/ifft_columns_real`
- `signal::filter`、`filtfilt`、`filter_columns`、`filtfilt_columns`
- `signal::psd_welch`、`cpsd_welch`、`xcov_unbiased`、`detrend`、窗函数
- `integral::trapz/cumtrapz/simpson/romberg/adaptive_simpson/quad`
- `polynomial::polyfit/polyval`（签名不变）
- `interp::interp1_*` 全部自由函数（签名不变）
- `equation::*`、`ode::*`

补充：`#include <msl/signal.hpp>` 现在显式包含 `filter.hpp` / `filtfilt.hpp`，
可以移除 qREST 中的重复 include（不删也不影响）。

---

## 4. 迁移示例

```cpp
// 1) LU
auto lu = msl::matrix::lu(A);
msl::matrix::matrixd P(A.rows(), A.cols(), 0.0);
for (size_t i = 0; i < A.rows(); ++i) {
    P(i, lu.permutation[i]) = 1.0;   // P * A == L * U
}

// 2) copy 访问器
auto row = A.row_copy(i);
auto col = A.column_copy(j);
auto sub = A.submatrix_copy(r0, r1, c0, c1);

// 3) adjoint
auto Ah = msl::matrix::conjugate_transpose(A);  // 复数；实数用 transpose

// 4) 单帧谱
auto p  = msl::signal::power_spectrum(x, nfft);          // |X|^2（未归一化）
auto pxy = msl::signal::cross_power_spectrum(x, y, nfft); // X * conj(Y)
```

---

## 5. 验证建议

```sh
# C++ 全量
xmake run test_matrix
xmake run test_signal
xmake run test_integral
xmake run test_difference
xmake run test_interp
xmake run test_polynomial
xmake run test_equation
xmake run test_ode

# MATLAB 独立交叉验证（8/8 PASS）
matlab -batch "run('matlab/validation/run_all.m')"
```

qREST 完成替换后，建议先跑自身与 ODE/Butterworth 相关的回归，
确认 `[0.004, 0.8]` bandpass 的固定比例误差消失，再更新数值基线。
