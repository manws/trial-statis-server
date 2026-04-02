# 临床泊松分布统计分析模块 (Clinical Poisson Distribution Statistics Module)

## 模块概述

临床泊松分布统计分析模块提供全面的泊松分布统计分析功能，用于临床数据中稀有事件发生次数的分析，是临床试验数据分析的重要组成部分。泊松分布在医学研究中广泛应用，如计算某段时间内的疾病发病次数、意外事件发生次数、突变计数等，特别适用于稀有事件的统计推断。

## 模块功能

### 1. 概率计算
- 计算特定事件数的概率
- 计算累积概率

### 2. 临界值计算
- 给定累积概率下的临界事件数
- 计算分位数值

### 3. 区间估计
- 总体均数的区间估计
- 提供多种估计方法

### 4. 假设检验
- 样本与总体均数的比较
- 适合性检验

## 临床应用价值

- **疾病监测**：监测罕见疾病或不良事件的发生频率
- **质量控制**：监控医疗差错或感染事件的发生
- **风险评估**：评估医疗操作中并发症的预期发生数
- **资源规划**：预测急诊科就诊人数或住院患者数

## 统计方法选择指南

### 1. 泊松分布适用条件
- 事件发生的概率很小
- 试验次数很多
- 事件之间相互独立
- 在相同条件下，事件发生的概率恒定

### 2. 临床应用场景
- 某时间段内医院感染病例数
- 某区域内的罕见疾病发生数
- 某项操作的并发症发生数

## 主要函数说明

### 泊松分布概率函数 - probability_poisson

计算泊松分布的概率 P(X=k)。在临床研究中，这用于计算在给定时间内恰好发生k次稀有事件的概率，以评估事件发生的可能性。

**参数**：
- `k`: 事件发生次数
- `lam`: 泊松分布的均数λ（lambda）

**返回值**：
- `float`: 概率值 P(X=k)

### 泊松分布累积概率函数 - cumulative_probability_poisson

计算泊松分布的累积概率 P(X≤k)。在临床研究中，这用于计算在给定时间内最多发生k次稀有事件的概率，以评估事件发生频率的上限。

**参数**：
- `k`: 最大事件发生次数
- `lam`: 泊松分布的均数λ（lambda）

**返回值**：
- `float`: 累积概率值 P(X≤k)

### 临界值计算函数 - calculate_critical_value_poisson

计算泊松分布的临界值。在临床研究中，这用于确定在给定显著性水平alpha下，事件发生的临界次数，以设定临床决策的阈值。

**参数**：
- `alpha`: 显著性水平
- `lam`: 泊松分布的均数λ（lambda）

**返回值**：
- `int`: 临界值

### 概率区间函数 - probability_interval_poisson

计算泊松分布在指定区间内的概率 P(k_start ≤ X ≤ k_end)。在临床研究中，这用于计算事件发生次数在特定区间内的概率，以评估事件发生的可能范围。

**参数**：
- `k_start`: 区间起始值
- `k_end`: 区间结束值
- `lam`: 泊松分布的均数λ（lambda）

**返回值**：
- `float`: 区间概率值

### 结果计算函数 - cal_result_pois1

生成泊松分布统计分析的完整报告字典。此函数整合了泊松分布的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的泊松分布分析结果。

**参数**：
- `param`: BaseParamPois1对象，包含total_avg, total_p, sample_num

**返回值**：
- 包含泊松分布统计分析指标的字典，包括：
  - `parameters`: 参数信息
  - `distribution_properties`: 分布性质
  - `percentiles`: 百分位数
  - `probabilities`: 各事件数的概率列表
  - `interpretation`: 结果解释
  - `remark`: 备注信息

## 结果解读注意事项

1. **概率解释**：P(X=k)表示恰好发生k次事件的概率
2. **累积概率解释**：P(X≤k)表示最多发生k次事件的概率
3. **均数解释**：λ表示单位时间（空间）内事件的平均发生次数
4. **适用性检查**：注意稀有事件和独立性假设
5. **临床意义**：结合实际情况解释统计结果

## 使用示例

```python
from app.stats.base.base_stats_pois1 import cal_result_pois1
from app.schemas.request_data.base_param import BaseParamPois1

# 创建参数对象
param = BaseParamPois1(
    total_avg=2.5,      # 总体均数λ为2.5
    total_p=0.1,        # 总体概率为0.1
    sample_num=100      # 样本数为100
)

# 计算泊松分布统计分析
result = cal_result_pois1(param)

# 查看结果
print(f"参数信息: {result['parameters']}")
print(f"分布性质: {result['distribution_properties']}")
print(f"百分位数: {result['percentiles']}")
print(f"概率列表: {result['probabilities'][:5]}")  # 显示前5个
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范