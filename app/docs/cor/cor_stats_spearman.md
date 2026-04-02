# 临床Spearman秩相关分析统计分析模块 (Clinical Spearman Rank Correlation Analysis Statistics Module)

## 模块概述

临床Spearman秩相关分析统计分析模块提供全面的Spearman秩相关分析统计分析功能，用于临床数据中两个变量之间的单调关系分析，是一种重要的非参数统计分析方法。Spearman秩相关通过对变量的秩次进行Pearson相关分析，得到介于-1到1之间的相关系数，广泛应用于医学研究中的非正态数据关联性分析、等级资料相关性评估、预后因子筛选等领域。

## 模块功能

### 1. 秩相关系数计算
- 计算Spearman秩相关系数
- 衡量两个变量之间的单调关系强度

### 2. 显著性检验
- 执行检验判断相关系数的显著性
- 提供P值评估统计显著性

### 3. 置信区间
- 计算相关系数的置信区间
- 评估估计的可靠性

### 4. 决定系数
- 计算R²以评估解释力度
- 表示一个变量能被另一个变量解释的变异比例

### 5. 统计解释
- 提供结果的临床意义解释
- 帮助理解相关关系的实际含义

## 临床应用价值

- **非正态数据关联性分析**：评估非正态分布变量之间的关系
- **等级资料相关性评估**：分析等级数据间的关联程度
- **预后因子筛选**：识别与预后相关的等级指标
- **研究假设验证**：验证变量间存在单调关系的假设

## 统计方法选择指南

### 1. Spearman相关适用条件
- 两变量不必满足正态分布
- 两变量之间存在单调关系（不一定线性）
- 数据中可能存在异常值
- 观测值独立
- 变量至少为有序分类或连续变量

### 2. 临床应用场景
- 分析疾病严重程度与疗效评分之间的关系
- 评估患者满意度与治疗时间之间的关联
- 比较病理分级与预后评分的相关性
- 研究症状评分与生活质量之间的关系

## 主要函数说明

### Spearman相关分析函数 - spearman_correlation_from_stats_data

基于List[List[float]]的Spearman秩相关分析。在临床研究中，这用于分析两个变量之间的单调关系，特别适用于非正态分布数据或等级数据的相关性分析。

**参数**：
- `stats_data_list`: List[List[float]]，包含两个变量的数据列表，每个子列表代表一个变量的数据

**返回值**：
- 包含相关统计量和结果的字典

### 结果计算函数 - cal_result_cor_spearman

生成Spearman秩相关分析统计分析的完整报告字典。此函数整合了Spearman秩相关分析的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的相关分析结果。

**参数**：
- `stats_data_list`: List[List[float]]，包含两个变量的数据列表，每个子列表代表一个变量的数据

**返回值**：
- 包含Spearman秩相关分析统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `correlation_results`: 相关分析结果
  - `interpretation`: 统计解释

## 结果解读注意事项

1. **相关系数解释**：接近1表示强正单调关系，接近-1表示强负单调关系，接近0表示无单调关系
2. **P值解释**：在零假设（真实相关系数为0）成立的情况下，观察到当前或更极端结果的概率
3. **决定系数解释**：表示一个变量的秩次能够被另一个变量的秩次解释的变异百分比
4. **临床意义**：相关关系不等于因果关系，需结合专业知识进行解释

## 使用示例

```python
from app.stats.cor.cor_stats_spearman import cal_result_cor_spearman

# 创建数据列表
stats_data_list = [
    [1, 2, 2, 3, 3, 4, 4],      # 疾病分期
    [80, 75, 70, 60, 55, 40, 35]  # 疗效评分
]

# 执行Spearman相关分析
result = cal_result_cor_spearman(stats_data_list)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"相关分析结果: {result['correlation_results']}")
print(f"统计解释: {result['interpretation']}")
print(f"相关系数: {result['correlation_results']['correlation_coefficient']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范