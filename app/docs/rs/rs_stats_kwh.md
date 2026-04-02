# 临床Kruskal-Wallis H检验统计分析模块 (Clinical Kruskal-Wallis H Test Statistics Module)

## 模块概述

临床Kruskal-Wallis H检验统计分析模块提供全面的Kruskal-Wallis H检验统计分析功能，用于临床数据中多个独立样本的比较，是一种重要的非参数统计分析方法。Kruskal-Wallis H检验通过对多个独立样本的秩次进行比较，来判断各组之间是否存在显著差异，广泛应用于医学研究中的多组比较、方差分析不满足正态性假设时的替代方法等领域。

## 模块功能

### 1. 秩和计算
- 计算多个独立样本的秩和
- 衡量数据的分布差异

### 2. 检验统计量
- 计算Kruskal-Wallis检验统计量H
- 衡量组间差异的显著性

### 3. 显著性检验
- 执行卡方近似检验判断统计量的显著性
- 提供P值评估统计显著性

### 4. Nemenyi事后检验
- 当主检验显著时，进行成对比较
- 提供更精细的组间差异分析

### 5. 统计解释
- 提供结果的临床意义解释
- 帮助理解检验结果的实际含义

## 临床应用价值

- **多组比较**：比较多个治疗组的疗效差异
- **非正态数据**：处理不满足正态性假设的连续型数据
- **事后分析**：进行成对组间比较
- **分布自由**：对数据分布形状无特殊要求

## 统计方法选择指南

### 1. Kruskal-Wallis检验适用条件
- 数据为连续型或有序分类变量
- 各组数据独立
- 分布形状相似（但不要求正态分布）
- 至少有两组数据

### 2. 临床应用场景
- 比较三种或以上不同治疗方案的疗效
- 评估多个医院间患者满意度的差异
- 分析不同年龄段某项指标的分布差异
- 研究多种药物对某项生物标志物的影响

## 主要函数说明

### Kruskal-Wallis H检验函数 - kruskal_wallis_h_test_kwh

多样本Kruskal-Wallis H秩和检验。

**参数**：
- `groups_data`: List[List[float]]，多组样本数据列表，每组为一个List[float]

**返回值**：
- 包含检验统计量和结果的字典

### Nemenyi事后检验函数 - nemenyi_post_hoc_test

Nemenyi事后检验（用于Kruskal-Wallis检验后的多重比较）。

**参数**：
- `groups_data`: List[List[float]]，多组样本数据列表
- `alpha`: float，显著性水平，默认0.05

**返回值**：
- 包含事后检验结果的字典

### 结果计算函数 - cal_result_rs_kwh

生成Kruskal-Wallis H检验及Nemenyi事后检验统计分析的完整报告字典。此函数整合了Kruskal-Wallis H检验及Nemenyi事后检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的检验结果。

**参数**：
- `groups_data`: List[List[float]]，多组样本数据列表，每组为一个List[float]

**返回值**：
- 包含Kruskal-Wallis H检验及Nemenyi事后检验统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `group_statistics`: 组统计信息
  - `rank_statistics`: 秩统计信息
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释
  - `nemenyi_results`: Nemenyi事后检验结果（如果主检验显著）

## 结果解读注意事项

1. **检验统计量解释**：H值越大，表明组间差异越显著
2. **P值解释**：P值小于显著性水平（如0.05）时拒绝原假设
3. **临床意义**：统计显著性不等于临床意义，需结合专业知识进行解释

## 使用示例

```python
from app.stats.rs.rs_stats_kwh import cal_result_rs_kwh

# 创建多组样本数据
groups_data = [
    [120, 125, 130, 128, 122],  # 第一组
    [115, 120, 125, 123, 118],  # 第二组
    [130, 135, 140, 138, 132]   # 第三组
]

# 执行Kruskal-Wallis H检验
result = cal_result_rs_kwh(groups_data)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"组统计信息: {result['group_statistics']}")
print(f"秩统计信息: {result['rank_statistics']}")
print(f"检验统计量: {result['test_statistics']}")
print(f"显著性检验: {result['significance_tests']}")
print(f"统计解释: {result['interpretation']}")
print(f"H统计量: {result['test_statistics']['chi_square_value']}")
print(f"P值: {result['test_statistics']['p_value']}")
if result['nemenyi_results']:
    print(f"Nemenyi事后检验: {result['nemenyi_results']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范