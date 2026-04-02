# 临床二项分布统计分析模块 (Clinical Binomial Distribution Statistics Module)

## 模块概述

临床二项分布统计分析模块提供全面的二项分布统计分析功能，用于临床数据中阳性率、成功率等二项分布参数的估计和检验，是临床试验数据分析的重要组成部分。二项分布在医学研究中广泛应用，如计算某种疾病的发生率、某种治疗的成功率、某种基因型的比例等。

## 模块功能

### 1. 概率计算
- 计算特定阳性数的概率 P(X=k)
- 计算累积概率 P(X≤k)
- 计算区间概率

### 2. 区间估计
- 计算总体率的置信区间
- 提供多种区间估计方法

### 3. 假设检验
- 样本率与总体率的比较
- 两样本率的比较

### 4. 临界值计算
- 计算给定累积概率对应的阳性数
- 计算检验的临界值

## 临床应用价值

- **疗效评估**：评估某种治疗的成功率是否达到预期
- **疾病发生率**：估计特定人群中的疾病发生率
- **诊断准确性**：评估诊断测试的敏感性或特异性
- **临床决策**：为临床决策提供概率依据

## 统计方法选择指南

### 1. 二项分布适用条件
- 试验只有两种互斥结果（成功/失败）
- 试验次数固定
- 每次试验的成功概率相同
- 试验之间相互独立

### 2. 临床应用场景
- 治疗成功率的区间估计
- 新药有效性的假设检验
- 不同疗法效果的比较

## 主要函数说明

### 二项分布概率计算函数 - probability_binomial

计算二项分布的概率 P(X=k)。在临床研究中，这用于计算在n次独立试验中恰好出现k次阳性（成功）的概率，以评估某种现象的出现可能性。

**参数**：
- `k`: 阳性数（成功次数）
- `n`: 试验总数（样本量）
- `p`: 总体阳性概率（每次试验成功的概率）

**返回值**：
- `float`: 概率值 P(X=k)

### 二项分布累积概率函数 - cumulative_probability_binomial

计算二项分布的累积概率 P(X≤k)。在临床研究中，这用于计算在n次独立试验中最多出现k次阳性（成功）的概率，以评估某种现象出现频率的上限。

**参数**：
- `k`: 最大阳性数
- `n`: 试验总数
- `p`: 总体阳性概率

**返回值**：
- `float`: 累积概率值 P(X≤k)

### 临界值计算函数 - calculate_critical_value_binomial

计算二项分布的临界值。在临床研究中，这用于确定在给定显著性水平alpha下，n次试验中最可能出现的最大阳性数，以设定临床决策的阈值。

**参数**：
- `alpha`: 显著性水平
- `n`: 试验总数
- `p`: 总体阳性概率

**返回值**：
- `int`: 临界值

### 结果计算函数 - cal_result_bino1

生成二项分布统计分析的完整报告字典。此函数整合了二项分布的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的二项分布分析结果。

**参数**：
- `param`: BaseParamBino1对象，包含total_posi_p, sample_size

**返回值**：
- 包含二项分布统计分析指标的字典，包括：
  - `total_positive_probability`: 总体阳性概率
  - `sample_size`: 样本量
  - `expected_value`: 期望值
  - `variance`: 方差
  - `standard_deviation`: 标准差
  - `probabilities`: 各阳性数的概率列表

## 结果解读注意事项

1. **概率解释**：P(X=k)表示恰好k次成功的概率
2. **累积概率解释**：P(X≤k)表示最多k次成功的概率
3. **置信区间解释**：在给定置信水平下，总体率的可能范围
4. **P值解释**：在零假设成立的前提下，观察到当前或更极端结果的概率
5. **临床意义**：统计显著性不等同于临床重要性，需结合实际意义解读

## 使用示例

```python
from app.stats.base.base_stats_bino1 import cal_result_bino1
from app.schemas.request_data.base_param import BaseParamBino1

# 创建参数对象
param = BaseParamBino1(
    total_posi_p=0.3,   # 总体阳性概率为30%
    sample_size=20      # 样本量为20
)

# 计算二项分布统计分析
result = cal_result_bino1(param)

# 查看结果
print(f"总体阳性概率: {result['total_positive_probability']}")
print(f"样本量: {result['sample_size']}")
print(f"期望值: {result['expected_value']}")
print(f"方差: {result['variance']}")
print(f"标准差: {result['standard_deviation']}")
print(f"各阳性数概率: {result['probabilities'][:5]}")  # 显示前5个
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范