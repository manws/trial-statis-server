# 临床泊松分布两样本比较统计分析模块 (Clinical Two-Sample Poisson Distribution Comparison Statistics Module)

## 模块概述

临床泊松分布两样本比较统计分析模块提供全面的泊松分布两样本比较统计分析功能，用于临床数据中两个样本均数的差异检验，是临床试验数据分析的重要组成部分。两样本均数比较在医学研究中广泛应用，如比较两组患者的事件发生率、两种治疗方法的不良事件发生率、不同时期的疾病发病率等。

## 模块功能

### 1. 两样本均数比较
- 计算两个样本均数的差异及其统计显著性
- 使用适当的检验方法

### 2. 比率差计算
- 计算两个样本率的差值及置信区间
- 提供精确的区间估计

### 3. 比值比计算
- 计算相对风险比及其置信区间
- 评估相对风险倍数

### 4. 精确检验
- 提供泊松分布的精确检验方法
- 适用于小样本数据

### 5. 效应量计算
- 计算差异的效应量指标
- 评估差异的实际意义

## 临床应用价值

- **疾病监测**：比较不同时期的疾病发生率
- **安全评估**：比较不同治疗方法的不良事件发生率
- **政策评估**：比较不同政策实施前后的事件发生率
- **质量改进**：比较不同干预措施的效果

## 统计方法选择指南

### 1. 精确检验适用条件
- 适用于所有样本量
- 特别适用于小样本或事件数较少的情况

### 2. 正态近似法适用条件
- 两个样本的事件数都足够大
- 通常要求事件数 ≥ 5

### 3. 临床应用场景
- 比较新旧两种治疗方法的安全性
- 比较不同手术方式的并发症发生率
- 比较不同时期的感染率

## 主要函数说明

### 两样本泊松检验函数 - two_sample_poisson_test

两样本泊松分布比较检验。在临床研究中，这用于比较两组数据中事件发生率的差异，以评估其统计学显著性。

**参数**：
- `sample1_num`: 样本1的观察单位数
- `sample1_tick`: 样本1的事件数
- `sample2_num`: 样本2的观察单位数
- `sample2_tick`: 样本2的事件数
- `alpha`: 显著性水平 (默认0.05)

**返回值**：
- `Dict[str, Any]`: 包含检验统计量和结论的字典

### 率比置信区间函数 - calculate_confidence_interval_for_rate_ratio

计算两个泊松率比的置信区间。在临床研究中，这用于估计两个样本率比值的置信区间，以提供相对风险的精确范围估计。

**参数**：
- `events1`: 样本1的事件数
- `units1`: 样本1的观察单位数
- `events2`: 样本2的事件数
- `units2`: 样本2的观察单位数
- `confidence_level`: 置信水平 (默认0.95)

**返回值**：
- `Dict[str, float]`: 包含率比及其置信区间的字典

### 率差置信区间函数 - calculate_confidence_interval_for_rate_diff

计算两个泊松率差的置信区间。在临床研究中，这用于估计两个样本率差值的置信区间，以提供绝对风险差的精确范围估计。

**参数**：
- `events1`: 样本1的事件数
- `units1`: 样本1的观察单位数
- `events2`: 样本2的事件数
- `units2`: 样本2的观察单位数
- `confidence_level`: 置信水平 (默认0.95)

**返回值**：
- `Dict[str, float]`: 包含率差及其置信区间的字典

### 结果计算函数 - cal_result_pois4

生成泊松分布两样本比较统计分析的完整报告字典。此函数整合了泊松分布两样本比较的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的泊松分布两样本比较分析结果。

**参数**：
- `param`: BaseParamPois4对象，包含sample1_num, sample1_tick, sample2_num, sample2_tick

**返回值**：
- 包含泊松分布两样本比较统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `test_results`: 检验结果
  - `rate_ratio_ci`: 率比的置信区间
  - `rate_diff_ci`: 率差的置信区间
  - `interpretation`: 结果的专业解释

## 结果解读注意事项

1. **P值解释**：在零假设成立的前提下，观察到当前或更极端结果的概率
2. **置信区间解释**：均数差或比值的置信区间
3. **显著性水平**：通常使用α=0.05作为判断标准
4. **临床意义**：统计显著性不等同于临床重要性
5. **混杂因素**：注意潜在混杂因素的影响

## 使用示例

```python
from app.stats.base.base_stats_pois4 import cal_result_pois4
from app.schemas.request_data.base_param import BaseParamPois4

# 创建参数对象
param = BaseParamPois4(
    sample1_num=100,       # 样本1观察单位数
    sample1_tick=25,       # 样本1事件数
    sample2_num=120,       # 样本2观察单位数
    sample2_tick=35        # 样本2事件数
)

# 计算泊松分布两样本比较
result = cal_result_pois4(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"检验结果: {result['test_results']}")
print(f"率比置信区间: {result['rate_ratio_ci']}")
print(f"率差置信区间: {result['rate_diff_ci']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范