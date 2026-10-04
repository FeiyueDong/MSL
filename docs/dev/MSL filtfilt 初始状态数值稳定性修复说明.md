# MSL `filtfilt` 初始状态数值稳定性修复说明

## 1. 问题概述

在 qREST 的 RR 双积分流程中，对 MSL `filtfilt` 与 Python/SciPy reference 进行数值对比时发现：

- 加速度结果基本一致；
- 使用默认四阶 Butterworth 滤波后，位移结果出现约 `3.8e-6` 的最大差异；
- 问题进一步定位到 `filtfilt` 的 IIR 初始状态计算。

当前实现中，初始状态通过显式矩阵求逆计算：

```cpp
zi = inverse(sp) * rhs;
```

数学上等价于求解：

\[
sp\,z_i=rhs,
\]

但对于条件数较高的矩阵，显式求逆会放大数值舍入误差。

本次测试中该矩阵条件数约为：

\[
\kappa(sp)\approx1.5\times10^8.
\]

当前显式求逆得到的线性方程残差约为：

\[
3.4\times10^{-9},
\]

而采用直接线性方程求解后可降低到约：

\[
3.6\times10^{-16}.
\]

修改后 qREST Python/C++ RR 位移差异由约：

\[
3.8\times10^{-6}
\]

下降至约：

\[
5.8\times10^{-11}.
\]

因此问题可以确定为 `filtfilt` 初始状态求解中的数值稳定性问题，而不是滤波公式、padding 或滤波器系数定义不一致。

## 2. 建议修改

涉及文件：

```text
include/msl/signal/filtfilt.hpp
```

当前逻辑：

```cpp
auto sp_eig = matrix::eigen_interface::as_eigen(sp);
Eigen::MatrixXd sp_inv = sp_eig.inverse();

auto rhs_eig =
    Eigen::Map<const Eigen::VectorXd>(rhs.data(), rhs.size());

Eigen::VectorXd zi_eig = sp_inv * rhs_eig;

std::vector<double> zi(
    zi_eig.data(),
    zi_eig.data() + zi_eig.size());

return zi;
```

建议直接使用 MSL 已有线性方程求解接口：

```cpp
return matrix::eigen_interface::solve(sp, rhs);
```

即将：

\[
z_i=sp^{-1}rhs
\]

的显式逆实现改为直接求解：

\[
sp\,z_i=rhs.
\]

## 3. 修改范围

本轮只进行局部数值修复，不修改：

- `filtfilt` public API；
- Butterworth/filter coefficient 定义；
- padding 长度；
- odd reflection 规则；
- forward/backward filtering 流程；
- qREST 上层代码；
- MSL 其他矩阵接口。

不建议在本轮顺带进行全库 `inverse()` 重构。

## 4. 建议测试

建议至少增加一组 regression test。

### 4.1 初始状态线性方程残差

验证：

\[
\|sp\,z_i-rhs\|
\]

达到接近双精度机器精度的水平，例如：

```text
relative residual < 1e-12
```

避免以后重新使用显式求逆。

### 4.2 `filtfilt` 数值回归

使用固定四阶 Butterworth 和固定输入序列，对 `filtfilt` 输出进行回归验证。

建议覆盖当前暴露问题的典型参数：

```text
filter order: 4
sampling frequency: 50 Hz
bandpass: 0.1–20 Hz
zero-phase filtering
```

可使用已经验证的 SciPy / qREST Python reference 作为数值对照。

## 5. 验收条件

修复完成后应满足：

1. MSL 原有 `filtfilt` API 与行为语义保持不变；
2. 初始状态方程残差显著降低；
3. 原有 MSL 测试全部通过；
4. 新增 `filtfilt` regression 通过；
5. qREST Python/C++ ABR direct parity 从当前 `34/35` 恢复为 `35/35`。

本问题属于局部数值稳定性修复，预计不需要架构调整。