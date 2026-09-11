# MSL Correctness & API Cleanup 完整修复计划指南

## 1. 背景与目标

MSL 当前已经形成较完整的小型 C++ 数值计算库，主要包括：

- Matrix
- Signal
- Integral
- Difference
- Interpolation
- Polynomial
- Equation
- ODE

等基础模块。

整体架构是合理的，header-only + C++20 的形式也适合作为 qREST 等上层工程的基础数值层。本轮改造**不以重写 MSL 为目标**，而是针对现有实现进行一次系统性的：

> **Correctness Cleanup + API Cleanup + Validation Hardening**

核心目标包括：

1. 修复已经确认的数学/逻辑错误；
2. 清理容易导致误用的 API 语义；
3. 加强关键数值算法的数学不变量测试；
4. 扩展 MATLAB 独立交叉验证；
5. 整理 signal/Butterworth 的正确性；
6. 改善 build/test 基础设施；
7. 在不显著增加库复杂度的前提下，使 MSL 达到可长期作为 qREST 基础数值库使用的稳定程度。

本轮原则上：

> **优先修正确性，再清 API，再补测试，最后处理工程结构。**

不应为了“顺便优化”而引入大规模架构重构。

---

# 2. 本轮修改原则

## 2.1 正确性优先于兼容错误行为

对于已经确认存在数学错误的接口，例如：

- LU decomposition；
- row-wise IFFT；
- cumulative Simpson；
- wide QR；
- Butterworth bandpass/bandstop；

应直接修正正确行为。

如果修改导致历史错误结果发生变化，这是合理且必要的。

---

## 2.2 每个 bug 必须同时加入 regression test

禁止采用：

```text
先修代码
→ 后续统一补测试
```

应采用：

```text
确认 bug
→ 添加能够复现 bug 的测试
→ 修复实现
→ 测试通过
```

因此即使“测试体系扩展”整体属于后续阶段，每个 P0 correctness bug 的测试都必须与修复同时完成。

---

## 2.3 C++ Test 与 MATLAB Validation 职责分离

建立两层验证体系：

### Level 1 — C++ Unit / Invariant Test

负责：

- 参数检查；
- shape/dimension；
- exception/status；
- copy/view；
- 边界情况；
- 数学不变量；
- round-trip；
- regression。

例如：

\[
PA\approx LU
\]

\[
A\approx QR
\]

\[
A\approx U\Sigma V^T
\]

\[
AV\approx VD
\]

\[
IFFT(FFT(x))\approx x
\]

---

### Level 2 — MATLAB Cross Validation

负责：

> 将 MSL 的数值结果与独立成熟实现或解析解比较。

优先使用：

- MATLAB 官方函数；
- MATLAB Toolbox；
- 解析解。

避免自己在 MATLAB 中重新写一遍与 C++ 相同的算法。

现有 MATLAB validation 已覆盖全部 8 个主要模块，但目前更多属于模块级 smoke validation，算法/API 级覆盖仍不足。

---

# 3. 修复优先级

建议分为：

## P0 — 必须修复的正确性错误

1. LU permutation 缺失；
2. `ifft_rows()` / `ifft_rows_real()` 写入 row copy；
3. uniform `cumsimpson()` 奇数位置公式错误；
4. wide matrix thin QR 尺寸错误；
5. SVD test reconstruction 使用 `V` 而不是 `V^T`；
6. Butterworth bandpass center gain normalization 错误；
7. Butterworth bandstop pole transformation 错误；
8. Butterworth bandstop zero center / gain consistency。

---

## P1 — API 与数值语义清理

1. matrix copy/view API 命名；
2. `adjoint()` 语义及文档；
3. matrix public API runtime validation；
4. PSD / power spectrum 命名和归一化语义；
5. Butterworth 全部四类滤波器一致性审查；
6. interpolation 文档与实际默认行为统一；
7. polynomial fit 数值稳定性；
8. signal aggregator 完整性。

---

## P2 — 测试、构建与工程基础设施

1. MATLAB validation 扩充；
2. `run_all.m`；
3. test helper/project root；
4. xmake portability；
5. header-only library target；
6. CI；
7. 公共测试工具整理。

---

# 4. Phase 0：建立修改基线

正式修改前建议先：

1. 保留当前 `main` 作为 baseline；
2. 新建 correctness cleanup 分支；
3. 确认所有现有 C++ tests 当前状态；
4. 保存当前 MATLAB validation 输出；
5. 保存几个实际 qREST regression case。

特别建议保留 Kunming/RR 中实际使用的 Butterworth：

\[
f=[0.004,\;0.8]
\]

作为长期 regression case。该频带对 bilinear frequency warping 很敏感，是非常有效的错误检测样本。

---

# 5. Phase 1：Matrix Correctness

## 5.1 LU decomposition

### 当前问题

实现使用：

```cpp
Eigen::PartialPivLU
```

但 API 只返回：

```text
L
U
```

并声明：

\[
A=LU
\]

实际 partial-pivot LU 应包含 permutation。

正确关系应按接口约定明确为类似：

\[
PA=LU
\]

当前忽略 permutation 会导致需要 pivot 的矩阵得到错误分解。

---

### 修改方案

不建议取消 pivoting。

建议引入明确结果类型：

```cpp
struct LUResult {
    matrixd P;
    matrixd L;
    matrixd U;
};
```

或者保存 permutation vector：

```cpp
struct LUResult {
    permutation_type permutation;
    matrixd L;
    matrixd U;
};
```

优先推荐 permutation representation，避免无必要构造 dense P。

API 文档必须明确关系。

---

### Regression case

至少：

```text
A =
[0 1
 1 1]
```

验证：

\[
\|PA-LU\| < tolerance
\]

同时覆盖：

- 不需要 pivot；
- 必须 pivot；
- rectangular matrix（若 API 支持）。

---

## 5.2 QR decomposition

### 当前问题

thin QR 当前实际假设：

\[
m\ge n
\]

对 wide matrix：

\[
m<n
\]

仍试图产生：

```text
Q: m × n
R: n × n
```

这是错误的。

---

### 正确规则

定义：

\[
k=\min(m,n)
\]

economy QR：

\[
Q\in\mathbb R^{m\times k}
\]

\[
R\in\mathbb R^{k\times n}
\]

并满足：

\[
A\approx QR
\]

\[
Q^TQ\approx I_k
\]

---

### Tests

必须覆盖：

- tall；
- square；
- wide。

不能只测试 polynomial fit 中天然满足 `m>=n` 的情况。

---

## 5.3 SVD

### 当前问题

算法返回的是：

```text
U, S, V
```

但 real SVD test 中使用：

\[
USV
\]

进行 reconstruction。

正确为：

\[
USV^T
\]

complex：

\[
USV^H
\]

---

### 修改

统一变量命名：

```cpp
auto [U, S, V] = svd(A);
```

real：

```cpp
U * S * transpose(V)
```

complex：

```cpp
U * S * conjugate_transpose(V)
```

不要把 `V` 命名成 `Vt`。

---

### Tests

检查：

\[
\frac{\|A-U\Sigma V^T\|}{\|A\|}
\]

以及 singular values 非负、降序。

MATLAB validation 应调用真正的：

```matlab
[U,S,V] = svd(A)
```

而不仅仅检查 C++ 自己 reconstruction 的结果。

---

# 6. Phase 2：Signal Correctness

## 6.1 Row-wise IFFT

### 当前问题

当前逻辑类似：

```cpp
auto row_ifft = output.get_row(i);
fft_engine.inv(row_ifft.data(), ...);
```

而：

```cpp
get_row()
```

返回的是 owning copy。

因此 IFFT 写入 temporary matrix，而不是 `output`。

---

### 修改方案

短期最安全的方法：

```text
计算 row IFFT 到临时 buffer
→ 显式写 output(i,j)
```

暂时不要为了这个 bug 引入复杂 strided row-view。

---

### Tests

新增：

\[
IFFT_{rows}(FFT_{rows}(A))\approx A
\]

包括：

- real；
- complex；
- zero-padding；
- multiple rows。

同时检查 columns round-trip。

---

# 7. Phase 3：Integral Correctness

## 7.1 `cumsimpson()`

### 当前问题

uniform cumulative Simpson 当前对 odd index 使用类似：

\[
I_{i-1}=\frac{I_{i-2}+I_i}{2}
\]

该公式没有正确的 Simpson 数学依据。

例如：

\[
y=x^2,\quad x=0,1,2
\]

精确结果：

\[
\int_0^1x^2dx=\frac13
\]

当前算法会得到明显错误结果。

---

### 修改原则

不要简单把 odd point 当成两个 Simpson 点积分的算术平均。

应实现明确的 cumulative Simpson formulation。

若希望第一版保持简单，可以：

- even point 使用 Simpson；
- odd point 使用局部合理 quadrature；

但必须明确数学定义并测试。

---

### Regression test

使用：

\[
f(x)=x^2
\]

解析 primitive：

\[
F(x)=\frac{x^3}{3}
\]

检查**每一个 cumulative point**。

不能再像现有测试一样只检查首尾/偶数位置。

---

# 8. Phase 4：Butterworth Correctness Audit

Butterworth 本轮应作为正式 P0/P1 工作，而不是未来 enhancement。

实际工程中已经观察到 RR 位移存在稳定约：

```text
RelativeL2 ≈ 5.5e-3
```

且 MATLAB/C++ 位移基本表现为固定比例关系。进一步定位表明误差主要来自 Butterworth bandpass gain design。

---

## 8.1 本轮 Butterworth 修改边界

本轮完成：

```text
prototype
↓
prewarp
↓
frequency transform
↓
bilinear transform
↓
zeros/poles
↓
numerator/denominator
↓
gain normalization
↓
frequency-response validation
```

但**暂不全面重构为通用 ZPK pipeline**。

长期可以考虑：

```text
ZPK
↓
analog frequency transform
↓
bilinear
↓
ZPK → BA
```

这种形式理论上更稳定，也更不容易遗漏 zeros/gain，但不应扩大本轮 scope。

---

# 9. Butterworth Lowpass / Highpass

当前基本路线：

```text
digital cutoff
→ prewarp
→ analog prototype
→ frequency scaling
→ bilinear
→ digital poles
→ numerator
→ gain normalization
```

整体设计合理。

当前：

- lowpass 在 DC 归一化；
- highpass 在 Nyquist 归一化；

也是合理定义。

本轮主要进行 consistency audit，不需要结构重写。

---

### Tests

Lowpass：

\[
|H(0)|\approx1
\]

\[
|H(\omega_c)|\approx\frac1{\sqrt2}
\]

Highpass：

\[
|H(\pi)|\approx1
\]

\[
|H(\omega_c)|\approx\frac1{\sqrt2}
\]

同时：

\[
|p_i|<1
\]

---

# 10. Butterworth Bandpass

## 10.1 已确认 bug：中心频率错误

当前正确执行：

\[
\Omega_1=2\tan\frac{\pi f_1}{2}
\]

\[
\Omega_2=2\tan\frac{\pi f_2}{2}
\]

以及：

\[
\Omega_0=\sqrt{\Omega_1\Omega_2}
\]

但 gain normalization 又使用：

\[
f_c=\frac{f_1+f_2}{2}
\]

这是不一致的。

---

## 10.2 正确中心频率

应该从 analog center inverse-warp：

\[
f_0=
\frac{2}{\pi}
\arctan\frac{\Omega_0}{2}
\]

然后：

```cpp
normalize_gain(f0);
```

---

## 10.3 需要同时验证

- LP prototype pole；
- prewarp；
- LP→BP pole transform；
- zeros at \(z=\pm1\)；
- bilinear；
- gain；
- cutoff；
- center response。

---

### Tests

至少包括：

- narrow bandpass；
- medium bandpass；
- wide bandpass；
- `[0.004, 0.8]` regression。

验证：

\[
|H(f_0)|\approx1
\]

两个 cutoff：

\[
|H(f_1)|\approx|H(f_2)|\approx1/\sqrt2
\]

---

# 11. Butterworth Bandstop

这是本轮需要重点修正的另一项 P0。

## 11.1 当前 pole transform 已确认错误

当前 bandpass 与 bandstop 共用相同的 pole generation。

Bandpass 使用：

\[
s_L=\frac{s^2+\Omega_0^2}{Bs}
\]

Bandstop 应使用：

\[
s_L=\frac{Bs}{s^2+\Omega_0^2}
\]

两者不能使用同一 quadratic transformation。

因此必须将：

```text
design_bandpass_bandstop()
```

内部真正拆分 transform。

可以保留共同流程，但：

```text
BP analog pole transform
BS analog pole transform
```

必须是独立函数。

例如：

```cpp
transform_lowpass_to_bandpass(...)
transform_lowpass_to_bandstop(...)
```

---

## 11.2 Bandstop zero center

bandstop notch frequency 同样不能使用：

\[
\frac{f_1+f_2}{2}
\]

应由：

\[
\Omega_0=\sqrt{\Omega_1\Omega_2}
\]

inverse warp 得到真正的 digital center。



---

## 11.3 Gain normalization

Bandstop 推荐明确采用 passband reference，例如：

- DC；
- Nyquist；

并验证两端 passband gain。

不要通过错误的 notch center 做 unity normalization。

---

### Tests

验证：

\[
|H(f_0)|\approx0
\]

以及：

\[
|H(0)|\approx1
\]

和/或：

\[
|H(1)|\approx1
\]

并检查两个 cutoff。

---

# 12. Butterworth MATLAB Cross Validation

C++ 测试完成数学性质。

MATLAB 负责独立 reference。

建议新增：

```matlab
[b,a] = butter(...)
[h,w] = freqz(...)
```

比较：

- coefficients；
- frequency response；
- cutoff；
- center；
- notch；
- filtered waveform。

对于实际 filtering：

```matlab
filter(...)
filtfilt(...)
```

分别比较。

---

## 12.1 MATLAB Toolbox

Butterworth、filtfilt、freqz 等通常依赖 Signal Processing Toolbox。

因此 MATLAB validation 应支持：

```text
Toolbox available → run
Toolbox missing   → SKIPPED
```

不能因为开发机没有 toolbox 而导致整个 MSL 测试失败。

---

# 13. Matrix API Cleanup

完成 P0 后处理。

## 13.1 Copy / View 命名

当前：

```cpp
get_row()
get_column()
submatrix()
column()
```

copy/view 语义不够明确。

这已经直接导致 row-IFFT bug。

建议长期命名明确化：

```cpp
row_copy()
column_copy()
submatrix_copy()

column_view()
submatrix_view()
```

如果 row 无法提供 contiguous span，则不要伪装成普通 view。

---

## 13.2 `adjoint()`

当前实现使用 Eigen：

```cpp
.adjoint()
```

数学上表示：

- real：transpose；
- complex：conjugate transpose。

但文档曾将其描述为：

> classical adjoint / adjugate

这是完全不同的数学概念。

此外 conjugate transpose 对 rectangular matrix 也成立，不应要求 square。

建议：

### Preferred

弃用模糊的：

```cpp
adjoint()
```

统一：

```cpp
transpose()
conjugate_transpose()
```

如果未来真正需要：

```cpp
adjugate()
```

则单独实现。

---

# 14. Public Runtime Validation

matrix 底层目前大量使用：

```cpp
assert(...)
```

Release build 中 assert 会消失。

因此结构性错误可能转变为：

- out-of-bounds；
- UB；
- silent wrong result。

建议：

### 保留 unchecked 的部分

性能敏感的：

```cpp
operator()(i,j)
```

可以保持轻量。

### Public high-level operation

如：

```text
matmul
solve
submatrix
decomposition
resize-related
```

应使用：

```cpp
std::invalid_argument
std::out_of_range
```

或者明确 Result/Status。

Equation/ODE 当前的 result/status 模式已经比 matrix 更健壮，可以作为参考。

---

# 15. Signal PSD / Spectrum 语义

当前：

```cpp
cpsd(x,y,nfft)
psd(x,nfft)
```

本质更接近：

\[
X\overline Y
\]

而不一定具有真正 PSD density：

\[
power/Hz
\]

的归一化。

相反 Welch 版本已经明确采用：

\[
\frac{1}{f_s\sum w^2}
\]

归一化。

建议：

```text
power_spectrum()
cross_power_spectrum()

psd_welch()
cpsd_welch()
```

避免把未经 density normalization 的函数称为 PSD。

同时统一：

```text
nfft > input length
```

时的 zero-padding 规则。

---

# 16. Interpolation Cleanup

## 16.1 测试数据

当前 MATLAB interpolation validation 使用严格线性数据：

\[
y=2x+1
\]

这无法真正区分：

- linear；
- spline；
- PCHIP；
- MAKIMA。

应换成 nonlinear、nonuniform dataset。

例如：

```text
x = [0.0, 0.7, 1.8, 3.0, 5.0]
y = [0.0, 1.2, 0.8, 2.5, 1.0]
```

再建立 dense `xq`。

---

## 16.2 Cubic 文档

当前实现默认实际为：

```text
NotAKnot
```

并与 MATLAB `interp1(...,'spline')` 的 reference 方向一致。

但源码注释/文档部分仍称其为：

```text
Natural cubic spline
```

需要统一。

---

## 16.3 Akima

当前实现默认采用 modified Akima，即目标 MATLAB：

```matlab
makima
```

因此 MATLAB reference 合理。

文档应明确：

```text
default = modified Akima / MAKIMA
```

如果支持 original Akima，则单独测试。

---

## 16.4 Default / empty state

还需检查：

- default-constructed interpolator；
- empty `x_ / y_`；
- `range()`；
- `in_range()`；
- evaluation before `set_data()`。

应避免使用空 vector 的 `front()/back()`。

---

# 17. Polynomial Cleanup

当前：

```text
Vandermonde
→ QR
→ inverse(R)
→ ...
```

不建议显式：

\[
R^{-1}
\]

应使用 triangular solve：

\[
Rx=Q^Ty
\]

理由：

- 更稳定；
- 更快；
- 避免显式 inverse。

还应检查：

- rank deficiency；
- duplicate x；
- degree ≥ number of points；
- noisy least-squares。

---

# 18. Difference Tests

Difference 模块本身暂未确认严重 correctness bug，但 MATLAB coverage 偏低。

应该补：

- nonuniform forward gradient；
- nonuniform central gradient；
- second derivative；
- `gradient2d`；
- Laplacian；
- divergence；
- curl；
- Savitzky-Golay。

优先使用解析场。

例如：

\[
f=x^2+y^2
\]

则：

\[
\nabla^2f=4
\]

取：

\[
F=(x,y)
\]

则：

\[
\nabla\cdot F=2
\]

取：

\[
F=(-y,x)
\]

则：

\[
\nabla\times F=2
\]

这样比依赖 MATLAB 不同边界定义更稳定。

---

# 19. Equation / ODE

当前整体状态较好。

## Equation

MATLAB 已使用：

```matlab
fzero
```

验证多种 root solver。

继续由 C++ 测试：

- invalid bracket；
- derivative failure；
- max iteration；
- residual；
- status。

MATLAB 不需要重复异常测试。

---

## ODE

现有体系较成熟：

- Euler vs analytic；
- Heun vs analytic；
- RK4 vs analytic；
- oscillator；
- MSL ode45 vs MATLAB ode45。

后续可增加：

- non-autonomous ODE；
- coupled state；
- convergence order；
- tolerance sensitivity。

但不是本轮 P0。

---

# 20. MATLAB Validation 完整规划

当前 validation 目录已经实现：

```text
difference
equation
integral
interp
matrix
ode
polynomial
signal
```

即模块层面覆盖良好。

本轮重点提升“算法覆盖”。

---

## 20.1 Matrix MATLAB

增加：

```text
LU
QR
SVD
eig
generalized eig
pinv
inverse
determinant
truncated SVD
```

注意不要直接逐元素比较 eigenvector/SVD vector。

使用 residual / invariant。

---

## 20.2 Signal MATLAB

增加：

```text
fft_rows
ifft_rows
fft_columns
ifft_columns

Butterworth LP
Butterworth HP
Butterworth BP
Butterworth BS

filter
filtfilt

pwelch
cpsd
xcov
detrend
```

---

## 20.3 Integral MATLAB

增加：

```text
cumsimpson
matrix integration
nonuniform Simpson
edge cases
```

cumsimpson 优先使用 analytic primitive。

---

## 20.4 Interpolation MATLAB

增加：

```text
nonlinear dataset
nonuniform x
dense xq

linear
spline
pchip
makima
nearest
polynomial
```

extrapolation mode 若 MATLAB 没有完全等价语义，可由 C++ invariant tests 负责。

---

# 21. MATLAB Comparison Rules

不要统一使用：

```matlab
abs(a-b) < 1e-12
```

建议建立：

\[
|a-b|
\le
atol+rtol\cdot scale
\]

例如 helper：

```matlab
assert_close(actual, expected, rtol, atol)
```

matrix：

\[
\frac{\|A-B\|}{\max(\|B\|,\epsilon)}
\]

---

## 特殊算法

### Eigen

允许：

- ordering 不同；
- sign 不同；
- complex phase 不同。

检查：

\[
AV-VD
\]

### SVD

检查：

- singular values；
- reconstruction；
- orthogonality。

### Repeated eigenvalue/singular value

应比较 subspace，而不是单个 vector。

---

# 22. MATLAB Runner

增加：

```text
matlab/validation/run_all.m
```

输出类似：

```text
Matrix       PASS
Signal       PASS
Integral     PASS
Difference   PASS
Interp       PASS
Polynomial   PASS
Equation     PASS
ODE          PASS

Signal Toolbox:
Butterworth  PASS
filtfilt     PASS
Welch        PASS

Total: XX passed, XX skipped, XX failed
```

---

# 23. Signal Aggregator

当前：

```cpp
msl/signal.hpp
```

没有完整聚合：

```text
filter.hpp
filtfilt.hpp
```

因此 tests 仍需要显式 include `signal/filter.hpp`。

应修复 aggregator。

原则：

```cpp
#include <msl/signal.hpp>
```

应该能够访问 signal public API。

---

# 24. Test Infrastructure Cleanup

## 24.1 project_root()

多个测试重复实现：

```cpp
if (exists(path / "msl") &&
    exists(path / "tests"))
```

但当前结构已经变为：

```text
src/msl
src/tests
```

因此 helper 已过时。

建议建立：

```text
src/tests/test_utils.hpp
```

统一：

- project root；
- result directory；
- tolerance helper；
- matrix compare；
- complex compare；
- file export。

避免每个 test 复制一遍。

---

# 25. xmake Cleanup

当前 root `xmake.lua` 存在环境相关路径，例如 Windows MinGW SDK / vcpkg include 的 hard-coded path。

本轮建议建立真正的 header-only target：

```text
target("msl")
    set_kind("headeronly")
```

tests：

```text
add_deps("msl")
```

统一：

```text
src/msl
```

include/export。

避免测试 target 自己隐式拼 include path。

---

## Dependency

Eigen 依赖尽量统一为 xmake package/system package 机制。

不要把个人开发环境路径写入仓库。

---

# 26. CI

本轮 correctness cleanup 完成后再加 CI。

建议至少：

```text
Linux GCC
Linux Clang
Windows MSVC
```

普通 CI 只跑：

```text
C++ tests
```

MATLAB validation 不应成为普通开源 CI 的硬依赖。

可作为：

```text
local validation
manual CI
release validation
```

---

# 27. 推荐实施顺序

## Stage A — Critical Correctness

依次修：

```text
A1 LU permutation
A2 row IFFT
A3 cumsimpson
A4 wide QR
A5 SVD tests
A6 Butterworth bandpass
A7 Butterworth bandstop
```

每项：

```text
regression test
→ implementation fix
→ C++ pass
```

---

## Stage B — Butterworth Validation

完成：

```text
LP/HP consistency
BP regression
BS regression
frequency response helper
stability
cutoff
center/notch
MATLAB butter/freqz/filter/filtfilt comparison
```

---

## Stage C — Matrix / Signal API Cleanup

完成：

```text
copy/view naming
adjoint
runtime validation
PSD semantics
signal aggregator
```

注意这一阶段尽量避免与 P0 修改混在同一个 commit。

---

## Stage D — MATLAB Validation Expansion

优先：

```text
Matrix
Signal
Integral
Interp
```

然后：

```text
Difference
Polynomial
```

Equation/ODE 主要保留已有体系。

---

## Stage E — Test / Build Infrastructure

完成：

```text
shared test helper
project_root
run_all.m
xmake headeronly target
dependency cleanup
CI
```

---

# 28. 建议的 Commit 拆分

不建议一次 giant commit。

推荐：

```text
fix(matrix): preserve LU permutation
fix(signal): correct row-wise inverse FFT
fix(integral): correct cumulative Simpson integration
fix(matrix): support economy QR for wide matrices
test(matrix): correct SVD reconstruction semantics

fix(signal): correct Butterworth bandpass center normalization
fix(signal): implement proper Butterworth bandstop transform
test(signal): add Butterworth response regression tests

refactor(matrix): clarify copy and view API
refactor(matrix): clarify transpose and adjoint semantics
refactor(signal): clarify spectrum and PSD naming

test(matlab): expand matrix reference validation
test(matlab): expand signal reference validation
test(matlab): improve interpolation reference cases
test(matlab): add cumulative integration validation

build: clean xmake header-only target
test: centralize test utilities
ci: add C++ test workflow
```

这样出现 regression 时更容易定位。

---

# 29. 暂不在本阶段处理的内容

以下内容建议明确排除：

## Butterworth ZPK 全量重构

虽然值得做，但暂不执行。

## SOS / biquad

高阶 IIR 使用 SOS 会更稳定，但属于后续 filter architecture enhancement。

## 通用 DSP framework

不应把 MSL 扩展成完整 DSP 库。

## Matrix expression templates

目前没有必要。

## SIMD / GPU

先保证 correctness。

## 性能专项优化

除非当前修复明显引入性能问题，否则推迟。

---

# 30. 本轮完成后的验收标准

## Correctness

所有已知 P0 问题有：

```text
明确 regression test
+
正确 implementation
```

---

## Matrix

至少满足：

\[
PA\approx LU
\]

\[
A\approx QR
\]

\[
Q^TQ\approx I
\]

\[
A\approx U\Sigma V^T
\]

\[
AV\approx VD
\]

pinv 满足 Moore-Penrose identities。

---

## Signal

至少满足：

\[
IFFT(FFT(x))\approx x
\]

row/column round-trip 正确。

Welch：

- auto PSD 非负；
- CPSD conjugate symmetry；
- frequency axis 正确。

Butterworth：

- stable poles；
- cutoff response；
- LP DC；
- HP Nyquist；
- BP warped center；
- BS warped notch；
- MATLAB frequency response 一致；
- qREST regression 不再出现原有固定比例误差。

---

## Integral

`cumsimpson()` 在解析多项式上所有 cumulative points 正确。

---

## Interpolation

非线性数据下：

```text
spline
PCHIP
MAKIMA
linear
```

分别与预期 reference 对应。

---

## MATLAB

所有主要模块能够：

```text
PASS / FAIL / SKIPPED
```

统一报告。

---

## Build

至少：

```text
Linux GCC
Linux Clang
Windows
```

不存在个人绝对路径依赖。

---

# 31. 最终目标状态

本轮完成后，MSL 应达到如下定位：

> **一个轻量、可独立测试、数学定义明确、能够被 qREST 稳定依赖的 C++ 基础数值库。**

其可信度来源不再只是“当前例子能够运行”，而应同时来自：

```text
数学不变量
+
边界/异常测试
+
MATLAB 独立 reference
+
真实工程 regression
```

形成：

\[
\boxed{
\text{Implementation}
+
\text{Invariant Test}
+
\text{Independent Validation}
+
\text{Application Regression}
}
\]

四层验证结构。

---

# 32. Codex 执行要求

在使用 Codex 实施本计划时，应遵守以下规则：

1. **先阅读现有实现和测试，再修改。**
2. 不因为某算法“看起来不标准”就直接重写。
3. 每个 correctness change 必须先确认数学定义。
4. 优先做最小正确修改。
5. 不扩大 public API，除非现有 API 无法表达正确结果。
6. 避免无关格式化和目录重构。
7. 每个 P0 修复必须附 regression test。
8. MATLAB reference 只用于独立验证，不应反向决定算法公式。
9. 若 MATLAB 与 MSL 不一致，首先检查：
   - frequency normalization；
   - boundary condition；
   - one-sided/two-sided；
   - scaling；
   - ordering；
   - sign/phase ambiguity；
   - algorithm variant。
10. 只有确认数学语义相同后，才能将差异视为 implementation bug。
11. Butterworth 本阶段只完成正确性修复，不进行完整 ZPK/SOS 架构升级。
12. 完成每一阶段后运行当前全部 C++ tests，避免跨模块 regression。

---

# 33. 推荐的第一批实际任务

如果从现在开始实施，建议 Codex 第一轮只处理：

```text
Task 1
修复 LU permutation + tests

Task 2
修复 ifft_rows / ifft_rows_real + round-trip tests

Task 3
修复 cumsimpson + polynomial cumulative tests

Task 4
修复 economy QR wide matrix + orthogonality/reconstruction tests

Task 5
修正 real SVD test semantics

Task 6
完整审查并修复 Butterworth BP/BS
+ frequency-response tests
+ MATLAB comparison data
```

第一轮完成并稳定后，再进入 API cleanup 和 validation expansion。

这能有效避免一次修改过多模块，使错误来源难以追踪。

---

## 结论

当前 MSL 不需要推倒重写。

本轮最重要的任务是：

> **把已经具备良好雏形的数值库，从“多数功能可用”提升为“核心数学行为经过系统验证”。**

其中 Matrix、Signal 和 Integral 是最高优先级；Butterworth 已经通过实际 qREST 结果证明存在 correctness 问题，因此应与 LU、IFFT、Simpson、QR 等问题一并纳入本轮修复。

Butterworth 修复后的目标也不应定义为“与 MATLAB 完全一致”，而应是：

> **自身数字滤波理论自洽，同时能够通过 MATLAB 等成熟实现的独立交叉验证。**

这也是本轮整个 MSL 改造最适合采用的原则。