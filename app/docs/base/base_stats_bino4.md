# 临床二项分布两样本比较统计分析模块 (Clinical Two-Sample Binomial Distribution Comparison Statistics Module)

## 模块概述

临床二项分布两样本比较统计分析模块提供全面的二项分布两样本比较统计分析功能，用于临床数据中两个样本率的差异检验，是临床试验数据分析的重要组成部分。两样本率比较在医学研究中广泛应用，如比较两种治疗方法的成功率、两种药物的有效率、两组患者的不良反应发生率等。

## 模块功能

### 1. 两样本率比较
- 计算两个样本率的差异及其统计显著性
- 使用正态近似法进行检验

### 2. 比率差计算
- 计算两个样本率的差值及置信区间
- 提供精确的区间估计

### 3. 比值比计算
- 计算相对风险比及其置信区间
- 评估相对风险倍数

### 4. Fisher精确检验
- 对于小样本提供精确检验方法
- 适用于四格表数据

### 5. 效应量计算
- 计算差异的效应量指标
- 评估差异的实际意义

## 临床应用价值

- **疗效比较**：比较不同治疗方法的疗效差异
- **药物评价**：比较不同药物的有效性和安全性
- **风险评估**：比较不同暴露因素的风险差异
- **质量改进**：比较不同干预措施的效果

## 统计方法选择指南

### 1. 正态近似法适用条件
- 两个样本量都足够大
- 每组np ≥ 5 且 n(1-p) ≥ 5

### 2. Fisher精确检验适用条件
- 小样本数据
- 某些格子期望频数小于5

### 3. 临床应用场景
- 比较新旧两种治疗方法的有效率
- 比较不同手术方式的并发症发生率
- 比较不同药物的不良反应率

## 主要函数说明

### 两样本比例检验函数 - two_sample_proportion_test

二项分布两样本比较检验。在临床研究中，这用于比较两组数据中阳性率的差异，以评估其统计学显著性。

**参数**：
- `sample1_num`: 样本1的总数
- `sample1_posi_num`: 样本1的阳性数
- `sample2_num`: 样本2的总数
- `sample2_posi_num`: 样本2的阳性数
- `alpha`: 显著性水平 (默认0.05)

**返回值**：
- `Dict[str, Any]`: 包含检验统计量和结论的字典

### 比率差置信区间函数 - calculate_confidence_interval_for_diff

计算两个比例差值的置信区间。在临床研究中，这用于估计两个样本率差值的置信区间，以提供差异的精确范围估计。

**参数**：
- `p1`: 样本1的比率
- `n1`: 样本1的总数
- `p2`: 样本2的比率
- `n2`: 样本2的总数
- `confidence_level`: 置信水平 (默认0.95)

**返回值**：
- `Dict[str, float]`: 包含置信区间上下限的字典

### 相对风险比函数 - calculate_relative_risk

计算相对风险比及其置信区间。在临床研究中，这用于评估一个组相对于另一个组的风险倍数，是流行病学研究中的重要指标。

**参数**：
- `n1`: 样本1的总数
- `x1`: 样本1的阳性数
- `n2`: 样本2的总数
- `x2`: 样本2的阳性数

**返回值**：
- `Dict[str, float]`: 包含相对风险比及其置信区间的字典

### 结果计算函数 - cal_result_bino4

生成二项分布两样本比较统计分析的完整报告字典。此函数整合了二项分布两样本比较的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的二项分布两样本比较分析结果。

**参数**：
- `param`: BaseParamBino4对象，包含sample1_num, sample1_posi_num, sample2_num, sample2_posi_num

**返回值**：
- 包含二项分布两样本比较统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `test_results`: 检验结果
  - `confidence_interval`: 比率差的置信区间
  - `relative_risk`: 相对风险比
  - `validation`: 有效性验证
  - `interpretation`: 结果的专业解释

## 结果解读注意事项

1. **P值解释**：在零假设成立的前提下，观察到当前或更极端结果的概率
2. **置信区间解释**：比率差或比值比的置信区间
3. **显著性水平**：通常使用α=0.05作为判断标准
4. **临床意义**：统计显著性不等同于临床重要性
5. **混杂因素**：注意潜在混杂因素的影响

## 使用示例

```python
from app.stats.base.base_stats_bino4 import cal_result_bino4
from app.schemas.request_data.base_param import BaseParamBino4

# 创建参数对象
param = BaseParamBino4(
    sample1_num=100,         # 样本1总数
    sample1_posi_num=30,     # 样本1阳性数
    sample2_num=120,         # 样本2总数
    sample2_posi_num=40      # 样本2阳性数
)

# 计算二项分布两样本比较
result = cal_result_bino4(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"检验结果: {result['test_results']}")
print(f"比率差置信区间: {result['confidence_interval']}")
print(f"相对风险比: {result['relative_risk']}")
print(f"有效性验证: {result['validation']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范