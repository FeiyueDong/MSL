# OMA 算法对 MSL 的数值能力需求

## 1. 目的与边界

本文记录 qREST 将 MATLAB FDD 与 SSI-COV 移植到 C++ 时需要 MSL 提供的基础数值能力。
MSL 只负责通用矩阵分解和信号处理；候选模态筛选、MAC、EFDD 阻尼、稳定极点跟踪、
聚类和 OMA 结果解释仍属于 `qrest_algorithm`。

参考实现：

- `matlab/algorithm/OperationalModalAnalysis/@FrequencyDomainDecomposition/FrequencyDomainDecomposition.m`
- `matlab/algorithm/OperationalModalAnalysis/@SSICOV/SSICOV.m`

## 2. 优先级概览

| 优先级 | 能力 | 用途 |
|---|---|---|
| P0 | 保留复特征值的非对称实矩阵特征分解 | SSI-COV 由离散极点计算频率和阻尼 |
| P0 | 基于 SVD 和容差的矩阵伪逆 | SSI-COV 可观测矩阵移位求解 |
| P1 | Welch CPSD 的窗口长度与 FFT 长度分离 | 复现 FDD 的分段、重叠和零填充 |
| P1 | 薄型或前 `r` 阶 SVD | 降低 SSI-COV block-Toeplitz 分解成本 |
| P1 | FFT 自动零填充 | 用线性相关实现 SSI-COV NExT |
| P2 | 无偏互协方差 | 支持 SSI-COV `methodCOV = 2` |

## 3. P0：保留复特征值

### 需求

对实数非对称方阵 `A` 计算

```text
A * V = V * D
```

时，`V` 和 `D` 都必须保留复数结果。当前 MSL 的普通和广义特征分解把特征值写成
`eig_D(i).real()`，会丢失共轭极点的虚部。SSI-COV 需要计算：

```text
mu   = log(lambda) / dt
fn   = abs(mu) / (2*pi)
zeta = -real(mu) / abs(mu)
```

因此至少应修正普通实矩阵特征分解；广义特征分解建议保持相同返回语义。

### 验收建议

- 对 `[[a, -b], [b, a]]` 返回 `a + ib` 和 `a - ib`；
- 验证每个特征对的相对残差 `||A*v-lambda*v|| / (||A||*||v||)`；
- 结果中不得把非零虚部截断为零。

## 4. P0：矩阵伪逆

### 需求

提供实矩阵 Moore-Penrose 伪逆，输入为 `m x n`，输出为 `n x m`。实现应基于 SVD，
并允许调用方指定容差；未指定时建议采用与 MATLAB `pinv` 接近的规则：

```text
tol = max(m, n) * eps * max(singular_values)
```

小于等于容差的奇异值按零处理，不能通过普通矩阵求逆代替。第一阶段只要求实矩阵版本；
若 MSL 希望保持接口完整，可同时提供复矩阵版本。

### 验收建议

- 覆盖高矩阵、宽矩阵、秩亏矩阵和零矩阵；
- 验证 `A*A+*A ~= A` 与 `A+*A*A+ ~= A+`；
- 使用病态矩阵检查容差改变时的有效秩；
- 与 MATLAB `pinv` 在相同容差下比较。

## 5. P1：Welch 互功率谱

### 需求

扩展 Welch CPSD，使以下参数相互独立：

- `window` / `segment_length`；
- `noverlap`；
- `nfft`，允许 `nfft >= segment_length`；
- `sampling_rate`，用于功率谱密度尺度和频率轴。

当 `nfft > segment_length` 时，每段加窗后应补零到 `nfft` 再执行 FFT。当前 MSL 将
段长和 FFT 长度合并为 `nperseg`，不能表达 MATLAB FDD 使用的“较短窗口 + 较长 FFT”。

若提供单边 CPSD，实信号内部频点应采用 MATLAB `cpsd` 一致的 PSD density 缩放；
DC 和 Nyquist 不加倍，其余正频率点加倍。若继续返回完整双边谱，也应明确归一化定义，
由 qREST 统一截取单边谱。

### 验收建议

- 覆盖 `segment_length == nfft` 和 `segment_length < nfft`；
- 与 MATLAB `cpsd(x,y,window,noverlap,nfft,fs)` 比较频率轴和复 CPSD；
- 自谱应为实数非负，互谱矩阵应满足 Hermitian 共轭关系；
- 检查最后一个完整分段的计数、窗能量归一化和零填充后的频率分辨率。

## 6. P1：薄型或截断 SVD

### 需求

SSI-COV 对 block-Toeplitz 矩阵只使用前 `Nmax` 个左奇异向量和奇异值。当前基准中：

- Kunming 矩阵约为 `2232 x 2232`；
- Wuhan 矩阵约为 `3348 x 3348`；
- `Nmax = 30`。

建议在保留现有完整 SVD 的同时，增加 thin SVD 或前 `r` 阶 SVD。结果至少应包含：

- 按降序排列的前 `r` 个奇异值；
- 对应的 `m x r` 左奇异向量；
- 右奇异向量可按接口选择是否计算。

该项不是数学正确性的硬阻塞，但完整 Jacobi SVD 可能使 SSI-COV 的内存和耗时不可接受。

### 验收建议

- 前 `r` 个奇异值与完整 SVD 一致；
- 验证 `U_r` 的正交性和截断重构误差；
- 使用上述两种矩阵规模做内存和运行时间检查。

## 7. P1：FFT 零填充

### 需求

FFT 接口应允许 `nfft > input.size()`，并按标准语义在输入末尾补零。SSI-COV 的 FFT
相关法需要至少 `2*N-1` 点的线性相关长度，通常取不小于该长度的二次幂。

如果 MSL 不调整现有 FFT 接口，qREST 也可以在调用前显式构造补零数组；但由 MSL
统一实现更不易产生循环相关与线性相关混淆。

### 验收建议

- 补零 FFT 与显式补零数组的 FFT 一致；
- FFT/IFFT 往返误差满足双精度预期；
- 用脉冲和短序列验证互相关不存在循环回绕。

## 8. P2：无偏互协方差

### 需求

为完整支持 SSI-COV `methodCOV = 2`，需要等价于 MATLAB：

```matlab
xcov(x, y, max_lag, 'unbiased')
```

的实信号互协方差。应明确：

- 输出 lag 范围和顺序；
- 是否移除均值；
- 无偏归一化分母为 `N - abs(lag)`。

当前数据集配置使用 `methodCOV = 1`，因此该项可以晚于第一版 SSI-COV。

## 9. 不要求放入 MSL 的内容

以下逻辑保留在 qREST OMA 算法层：

- FDD 候选频率、MAC 去重和 EFDD 钟形谱；
- 峰值选择及阻尼拟合的业务规则；
- SSI-COV 稳定性状态、极点匹配和轨迹；
- MAC 距离、层次聚类、IQR 异常值剔除；
- `ObservationLayout`、X/Y/Rz 振型和 `Other` 成分解释。

## 10. qREST 接入条件

MSL 完成修改后，请提供对应头文件入口和最小使用示例。qREST 接入时将补充：

- MSL 基础函数的针对性数值测试；
- Kunming/Wuhan 的 FDD 与 SSI-COV MATLAB 对比；
- 频率与阻尼容差、相位对齐后的振型 MAC；
- 大型 Toeplitz 矩阵的运行时间和内存记录。
