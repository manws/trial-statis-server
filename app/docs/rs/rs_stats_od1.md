# 临床两样本Wilcoxon秩和检验统计分析模块 (Clinical Two-Sample Wilcoxon Rank Sum Test Statistics Module)

## 模块概述

临床两样本Wilcoxon秩和检验统计分析模块提供全面的两样本Wilcoxon秩和检验统计分析功能，用于临床数据中两个独立样本的比较，是一种重要的非参数统计分析方法。两样本Wilcoxon秩和检验（也称为Mann-Whitney U检验）通过比较两个独立样本的秩次，来判断两组之间是否存在显著差异，广泛应用于医学研究中的两组比较、等级资料分析、非正态数据检验等领域。

## 模块功能

### 1. 秩和计算
- 计算两组独立样本的秩和
- 衡量数据的分布差异

### 2. U统计量
- 计算Mann-Whitney U检验统计量
- 衡量两组差异的显著性

### 3. 显著性检验
- 执行正态近似检验判断统计量的显著性
- 提供P值评估统计显著性

### 4. 结修正
- 处理数据中相同值的修正
- 确保检验的准确性

### 5. 统计解释
- 提供结果的临床意义解释
- 帮助理解检验结果的实际含义

## 临床应用价值

- **两组比较**：比较两个治疗组的疗效差异
- **等级资料**：分析有序分类变量的差异
- **非正态数据**：处理不满足正态性假设的数据
- **分布自由**：对数据分布形状无特殊要求

## 统计方法选择指南

### 1. Wilcoxon秩和检验适用条件
- 数据为连续型或有序分类变量
- 两组数据独立
- 分布形状相似（但不要求正态分布）
- 观测值相互独立

### 2. 临床应用场景
- 比较两种不同治疗方案的疗效
- 评估不同药物对某项指标的影响
- 分析不同人群某项指标的分布差异
- 研究不同诊断方法的结果差异

## 主要函数说明

### 两样本Wilcoxon秩和检验函数 - wilcoxon_rank_sum_test_ordinal

两样本Wilcoxon秩和检验（适用于等级资料）。

**参数**：
- `group1_data`: List[float]，第一组样本数据（等级资料）
- `group2_data`: List[float]，第二组样本数据（等级资料）

**返回值**：
- 包含检验统计量和结果的字典

### 结果计算函数 - cal_result_rs_od1

生成两样本Wilcoxon秩和检验统计分析的完整报告字典。此函数整合了两样本Wilcoxon秩和检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的检验结果。

**参数**：
- `group1_data`: List[float]，第一组样本数据
- `group2_data`: List[float]，第二组样本数据

**返回值**：
- 包含两样本Wilcoxon秩和检验统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `sample_statistics`: 样本统计信息
  - `rank_statistics`: 秩统计信息
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释

## 结果解读注意事项

1. **检验统计量解释**：U值越偏离期望值，表明两组差异越显著
2. **P值解释**：P值小于显著性水平（如0.05）时拒绝原假设
3. **临床意义**：统计显著性不等于临床意义，需结合专业知识进行解释

## 使用示例

```python
from app.stats.rs.rs_stats_od1 import cal_result_rs_od1

# 创建两组样本数据
group1_data = [120, 125, 130, 128, 122, 135, 127, 124]
group2_data = [115, 120, 125, 123, 118, 130, 122, 119]

# 执行两样本Wilcoxon秩和检验
result = cal_result_rs_od1(group1_data, group2_data)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"样本统计信息: {result['sample_statistics']}")
print(f"秩统计信息: {result['rank_statistics']}")
print(f"检验统计量: {result['test_statistics']}")
print(f"显著性检验: {result['significance_tests']}")
print(f"统计解释: {result['interpretation']}")
print(f"U统计量: {result['rank_statistics']['test_statistic_U']}")
print(f"P值: {result['test_statistics']['p_value_two_sided']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范