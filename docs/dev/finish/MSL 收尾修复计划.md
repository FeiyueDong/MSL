# MSL 收尾修复计划

## 1. 目标

当前 MSL 主体修复已经完成，核心 Matrix / Signal 数值问题基本解决，qREST RR 回归测试也已恢复到数值误差级别。

本轮仅进行一次小范围收尾修复，目标是：

- 修复剩余明确的数值逻辑问题；
- 补齐少量边界输入检查；
- 不进行新的架构调整或大规模 API 修改；
- 完成后结束本轮 Correctness & API Cleanup。

---

## 2. P1：修复非均匀网格 Matrix Central Gradient

文件：

`src/msl/difference/difference.hpp`

### 当前问题

向量版本：

```cpp
central_gradient(x, y)
```

已经使用适用于非均匀网格的三点加权中心差分。

但矩阵版本：

```cpp
central_gradient(mat, x, axis)
```

内部点仍采用：

```cpp
(f[i+1] - f[i-1]) / (x[i+1] - x[i-1])
```

该公式只在均匀网格下等价于标准中心差分。

例如：

```text
x = [0, 1, 3]
f = x²
```

在 `x=1` 处真实导数为 `2`，当前矩阵实现得到 `3`。

### 修改要求

使 Matrix non-uniform overload 与 vector overload 使用相同的三点非均匀中心差分公式。

`axis=0` 和 `axis=1` 均需要处理。

尽量提取一个内部 helper，避免 vector / matrix 两套公式再次发生偏差。

### 测试

增加非均匀坐标下的矩阵测试，例如：

```text
f(x) = x²
x = [0, 1, 3, 4]
```

验证内部节点导数精确等于：

```text
2x
```

---

## 3. P1：修复 CubicSpline 右端点导数越界

文件：

`src/msl/interp/interp_1d_base.hpp`

以及必要时：

`src/msl/interp/interp_1d_cubic.hpp`

### 当前问题

`find_interval(x)` 在：

```cpp
x == x_.back()
```

时可能返回：

```cpp
x_.size() - 1
```

但 CubicSpline 的区间系数：

```cpp
b_
d_
```

只有：

```cpp
x_.size() - 1
```

个元素，因此：

```cpp
derivative(x_.back())
second_derivative(x_.back())
```

可能访问越界。

普通 `operator()` 因为 exact-match 提前返回，所以暂时不会暴露该问题。

### 修改要求

统一规定 `find_interval()` 返回的是“区间编号”，因此合法范围必须始终满足：

```text
0 <= index <= n - 2
```

对于：

```cpp
x == x_.back()
```

应返回最后一个有效区间：

```cpp
n - 2
```

不要只在 CubicSpline derivative 中做特殊补丁，应从 interval helper 的语义上解决。

### 测试

增加：

```cpp
spline.derivative(x.front())
spline.derivative(x.back())

spline.second_derivative(x.front())
spline.second_derivative(x.back())
```

确保：

- 不越界；
- 返回有限值；
- 对已知多项式数据能够得到正确结果。

---

## 4. P2：补齐 Difference 模块边界检查

文件：

`src/msl/difference/difference.hpp`

### 4.1 空输入下的 size_t 下溢

部分 allocating overload 当前先执行：

```cpp
std::vector<double> result(y.size() - 1);
```

再调用实际实现。

当：

```cpp
y.size() == 0
```

时会发生无符号整数下溢，可能尝试分配极大内存，而不是正常抛出参数异常。

涉及至少：

```cpp
diff(y)
forward_gradient(y, ...)
```

### 修改要求

在进行 `size()-1` 分配之前先验证输入长度。

原则上：

```text
diff / forward_gradient:
至少需要 2 个数据点
```

---

### 4.2 检查无效 spacing

对于 uniform gradient：

```cpp
dx
```

至少应拒绝：

```text
dx == 0
```

推荐要求：

```text
dx > 0
```

对于 non-uniform coordinate：

```cpp
x
```

建议统一要求：

```text
严格递增
x[i+1] > x[i]
```

避免：

- 重复坐标；
- 除零；
- 非预期的反向坐标。

建议将坐标验证提取成内部 helper，在：

- `forward_gradient`
- `central_gradient`
- matrix overload

中复用。

---

## 5. P3：少量数值健壮性完善

以下项目不是当前 correctness blocker，可在本轮顺手完成，也可以保留为 TODO。

### Butterworth polynomial realness

`poles_to_polynomial()` 当前直接取：

```cpp
poly[i].real()
```

理论上共轭极点应得到实系数，但建议检查虚部 residual。

若：

```text
abs(imag) > tolerance
```

应报告数值异常，而不是静默丢弃。

### Butterworth gain normalization

当前类似：

```cpp
if (gain == 0.0)
    return;
```

建议改为验证：

```text
isfinite(gain)
gain > tolerance
```

否则抛出异常。

合法 Butterworth 设计的 normalization point 不应出现零增益，因此异常情况不应静默生成未归一化滤波器。

---

## 6. 测试要求

本轮新增测试重点针对本次修改，不需要再次大规模扩展测试体系。

至少增加以下 regression tests：

```text
1. non-uniform matrix central gradient
   - quadratic function
   - axis=0
   - axis=1

2. CubicSpline endpoint derivative
   - left endpoint
   - right endpoint
   - first derivative
   - second derivative

3. Difference invalid input
   - empty input
   - single-point input
   - dx == 0
   - repeated x coordinate
   - non-increasing x coordinate
```

完成后保证现有全部 C++ 测试继续通过。

---

## 7. 本轮不再处理的内容

以下项目已经完成验证或不属于本轮收尾范围：

- LU permutation 修复；
- row-wise IFFT；
- cumulative Simpson；
- wide QR；
- SVD reconstruction；
- Matrix copy/view API；
- PSD / power spectrum 命名；
- Butterworth BP / BS 核心算法；
- qREST RR 对接误差；
- MSVC 本地编译运行；
- ZPK / SOS 高阶 IIR 重构；
- 新增算法或性能优化。

其中高阶 IIR 的 SOS 化可以作为未来独立优化任务，不应继续扩大本次修改范围。

---

## 8. 完成标准

本次收尾完成后应满足：

```text
[ ] Matrix non-uniform central gradient 与 vector 版本一致
[ ] CubicSpline endpoint derivative 无越界
[ ] Difference 空输入不会发生 size_t 下溢
[ ] spacing / coordinate 非法输入得到明确异常
[ ] 新增 regression tests 通过
[ ] 原有全部 C++ tests 通过
[ ] qREST RR 回归结果保持数值误差级别
```

完成上述内容后，可以结束本轮 MSL Correctness & API Cleanup，并将当前开发版本作为后续数值算法开发的稳定基础。