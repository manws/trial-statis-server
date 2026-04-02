# 临床单样本Wilcoxon符号秩检验统计分析模块 (Clinical Single-Sample Wilcoxon Signed-Rank Test Statistics Module)

## 模块概述

临床单样本Wilcoxon符号秩检验统计分析模块提供全面的单样本Wilcoxon符号秩检验统计分析功能，用于临床数据中单个样本的中位数与假设值的比较，是一种重要的非参数统计分析方法。单样本Wilcoxon符号秩检验通过比较观测值与假设中位数的差值的秩次，来判断样本中位数是否与假设中位数有显著差异，广泛应用于医学研究中的基线比较、疗效评估、异常值检测等领域。

## 模块功能

### 1. 差值计算
- 计算样本数据与假设中位数的差值
- 处理数据的符号和大小

### 2. 符号秩统计
- 计算正负差值的秩和
- 衡量数据偏离假设值的程度

### 3. 检验统计量
- 计算Wilcoxon检验统计量W
- 衡量样本中位数与假设中位数的差异显著性

### 4. 显著性检验
- 执行正态近似检验判断统计量的显著性
- 提供P值评估统计显著性

### 5. 统计解释
- 提供结果的临床意义解释
- 帮助理解检验结果的实际含义

## 临床应用价值

- **基线比较**：比较患者治疗前指标与正常值的差异
- **疗效评估**：评估治疗后指标与目标值的差异
- **异常值检测**：识别与预期值显著不同的观测值
- **质量控制**：检验检测结果与标准值的一致性

## 统计方法选择指南

### 1. 单样本Wilcoxon检验适用条件
- 数据为连续型变量
- 数据分布不要求正态性
- 观测值独立
- 差值分布对称

### 2. 临床应用场景
- 比较患者治疗前血压与正常血压值的差异
- 评估某种药物治疗后血糖水平与目标值的差异
- 分析实验室检测结果与标准值的偏差
- 研究某种干预措施对生物标志物的影响

## 主要函数说明

### 单样本Wilcoxon符号秩检验函数 - wilcoxon_signed_rank_single_sample_test

单样本Wilcoxon符号秩检验。

**参数**：
- `sample_data`: List[float]，样本数据列表
- `hypothesized_median`: float，假设的中位数，默认为0.0

**返回值**：
- 包含检验统计量和结果的字典

### 结果计算函数 - cal_result_rs_single

生成单样本Wilcoxon符号秩检验统计分析的完整报告字典。此函数整合了单样本Wilcoxon符号秩检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的检验结果。

**参数**：
- `sample_data`: List[float]，样本数据列表
- `hypothesized_median`: float，假设的中位数，默认为0.0

**返回值**：
- 包含单样本Wilcoxon符号秩检验统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `sample_statistics`: 样本统计信息
  - `rank_statistics`: 秩统计信息
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释

## 结果解读注意事项

1. **检验统计量解释**：W值越小，表明观测值与假设中位数的差异越显著
2. **P值解释**：P值小于显著性水平（如0.05）时拒绝原假设
3. **临床意义**：统计显著性不等于临床意义，需结合专业知识进行解释

## 使用示例

```python
from app.stats.rs.rs_stats_single import cal_result_rs_single

# 创建样本数据
sample_data = [120, 125, 130, 128, 122, 135, 127, 124, 126, 132]
hypothesized_median = 125.0

# 执行单样本Wilcoxon符号秩检验
result = cal_result_rs_single(sample_data, hypothesized_median)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"样本统计信息: {result['sample_statistics']}")
print(f"秩统计信息: {result['rank_statistics']}")
print(f"检验统计量: {result['test_statistics']}")
print(f"显著性检验: {result['significance_tests']}")
print(f"统计解释: {result['interpretation']}")
print(f"W统计量: {result['rank_statistics']['test_statistic_W']}")
print(f"P值: {result['test_statistics']['p_value_two_sided']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范