# 临床总体率的区间估计统计分析模块 (Clinical Overall Rate Interval Estimation Statistics Module)

## 模块概述

临床总体率的区间估计统计分析模块提供全面的总体率区间估计统计分析功能，用于临床数据中总体率的置信区间计算，是临床试验数据分析的重要组成部分。总体率的区间估计在医学研究中广泛应用，如计算某种疾病的发病率、某种治疗的成功率、某种不良反应的发生率等的置信区间，为临床决策提供精确的概率范围。

## 模块功能

### 1. 区间估计
- 计算总体率的置信区间（Wilson得分法、正态近似法等）
- 提供多种区间估计方法的结果对比

### 2. 精度评估
- 评估估计的精度和可靠性
- 比较不同方法的区间宽度

### 3. 结果解释
- 提供统计和临床意义的解释
- 判断方法的适用性

## 临床应用价值

- **疗效评估**：评估治疗成功率的置信区间
- **疾病发生率**：估计特定人群中疾病发生率的置信区间
- **诊断准确性**：估计诊断测试敏感性/特异性的置信区间
- **临床决策**：为临床决策提供概率范围而非单一估计值

## 统计方法选择指南

### 1. Wilson得分法适用条件
- 适用于小样本或比例接近0或1的情况
- 不受正态近似的限制
- 在边界情况下表现更好

### 2. 正态近似法适用条件
- 样本量足够大
- np ≥ 5 且 n(1-p) ≥ 5
- 比例远离0或1

### 3. 临床应用场景
- 治疗成功率的置信区间估计
- 疾病发生率的区间估计
- 不良反应发生率的区间估计

## 主要函数说明

### Wilson得分区间函数 - wilson_score_interval

计算Wilson得分置信区间。在临床研究中，这用于计算总体率的置信区间，特别适用于小样本或比例接近0或1的情况，以提供更可靠的区间估计。

**参数**：
- `n`: 样本量
- `x`: 样本阳性数
- `confidence_level`: 置信水平

**返回值**：
- `Dict[str, float]`: 包含置信区间上下限的字典

### 正态近似区间函数 - normal_approximation_interval

计算正态近似置信区间。在临床研究中，这用于计算总体率的置信区间，适用于大样本且比例远离0或1的情况，是最常用的区间估计方法之一。

**参数**：
- `n`: 样本量
- `x`: 样本阳性数
- `confidence_level`: 置信水平

**返回值**：
- `Dict[str, float]`: 包含置信区间上下限的字典

### Agresti-Coull区间函数 - agresti_coull_interval

计算Agresti-Coull置信区间。在临床研究中，这用于计算总体率的置信区间，是对正态近似法的一种改进，通过增加2个成功和2个失败来提高估计精度。

**参数**：
- `n`: 样本量
- `x`: 样本阳性数
- `confidence_level`: 置信水平

**返回值**：
- `Dict[str, float]`: 包含置信区间上下限的字典

### 样本量计算函数 - calculate_sample_size_for_ci

计算达到特定精度所需的样本量。在临床研究设计中，这用于确定获得指定精度的置信区间所需的样本量，有助于研究设计阶段的样本量估算。

**参数**：
- `confidence_level`: 置信水平
- `precision`: 所需精度（置信区间半宽）
- `p_estimate`: 比例预估值（如果未知，默认为0.5以获得最大样本量）

**返回值**：
- `int`: 所需样本量

### 结果计算函数 - cal_result_bino2

生成总体率区间估计统计分析的完整报告字典。此函数整合了总体率区间估计的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的区间估计分析结果。

**参数**：
- `param`: BaseParamBino2对象，包含confidence_level, sample_size, sample_posi_num

**返回值**：
- 包含总体率区间估计统计分析指标的字典，包括：
  - `sample_size`: 样本量
  - `sample_positive_number`: 样本阳性数
  - `sample_rate`: 样本率
  - `confidence_level`: 置信水平
  - `wilson_ci`: Wilson得分法置信区间
  - `normal_ci`: 正态近似法置信区间
  - `agresti_coull_ci`: Agresti-Coull法置信区间
  - `interval_widths`: 各方法的区间宽度
  - `normal_approximation_valid`: 正态近似法是否适用

## 结果解读注意事项

1. **置信区间解释**：95%置信区间表示真实总体率有95%的可能性落在该区间内
2. **区间宽度解释**：区间越窄，估计越精确
3. **临床意义**：置信区间应结合临床实际意义进行解读
4. **样本量影响**：样本量越大，置信区间越窄
5. **临床决策**：置信区间可用于临床决策的参考

## 使用示例

```python
from app.stats.base.base_stats_bino2 import cal_result_bino2
from app.schemas.request_data.base_param import BaseParamBino2

# 创建参数对象
param = BaseParamBino2(
    confidence_level=0.95,  # 95%置信水平
    sample_size=100,        # 样本量为100
    sample_posi_num=25      # 样本阳性数为25
)

# 计算总体率区间估计
result = cal_result_bino2(param)

# 查看结果
print(f"样本量: {result['sample_size']}")
print(f"样本阳性数: {result['sample_positive_number']}")
print(f"样本率: {result['sample_rate']}")
print(f"Wilson置信区间: [{result['wilson_ci']['lower_bound']:.4f}, {result['wilson_ci']['upper_bound']:.4f}]")
print(f"正态近似置信区间: [{result['normal_ci']['lower_bound']:.4f}, {result['normal_ci']['upper_bound']:.4f}]")
print(f"Agresti-Coull置信区间: [{result['agresti_coull_ci']['lower_bound']:.4f}, {result['agresti_coull_ci']['upper_bound']:.4f}]")
print(f"各方法区间宽度: {result['interval_widths']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范