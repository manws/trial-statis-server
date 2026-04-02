# 临床完全随机设计方差分析统计分析模块 (Clinical Completely Randomized Design ANOVA Statistics Module)

## 模块概述

临床完全随机设计方差分析统计分析模块提供全面的完全随机设计方差分析统计分析功能，用于临床数据中多个处理组均值差异的比较，是临床试验数据分析的重要组成部分。完全随机设计方差分析在医学研究中广泛应用，如比较不同治疗方案的效果、不同药物剂量的疗效差异等。

## 模块功能

### 1. 平方和计算
- 计算处理间平方和、误差平方和和总平方和
- 量化不同变异来源的贡献

### 2. 自由度计算
- 计算各变异来源的自由度
- 为均方和F检验提供基础

### 3. 均方计算
- 计算处理间均方和误差均方
- 评估每单位自由度的变异程度

### 4. F统计量计算
- 计算F统计量用于检验处理间差异
- 衡量处理效应的显著性

### 5. 显著性检验
- 执行F检验判断处理间差异是否显著
- 提供不同显著性水平的判断

### 6. 效应量计算
- 计算Eta平方评估处理效应大小
- 评估处理效应的实际意义

### 7. 多重比较
- 为处理组间差异提供配对比较信息
- 识别具体哪些处理组间存在差异

## 临床应用价值

- **治疗方案比较**：比较多种治疗方法的疗效差异
- **药物剂量研究**：评估不同剂量药物的疗效差异
- **疗效评估**：判断不同干预措施的效果差异
- **数据质量评估**：识别处理组间的变异来源

## 统计方法选择指南

### 1. 完全随机设计方差分析适用条件
- 各处理组数据独立
- 各处理组数据服从正态分布
- 各处理组方差齐性
- 至少2个处理组（2组时等价于独立样本t检验）

### 2. 临床应用场景
- 比较不同药物治疗方案的疗效
- 评估不同手术方法的效果差异
- 分析不同护理措施的改善程度

## 主要函数说明

### 完全随机设计方差分析函数 - calculate_crd_anova

执行完全随机设计(CRD)的方差分析。在临床研究中，这用于比较多个处理组的均值是否存在显著差异，以评估不同治疗方案或干预措施的效果差异。

**参数**：
- `data_list`: List[List[float]]，每个子列表代表一个处理组的数据

**返回值**：
- 包含CRD ANOVA统计量和结果的字典

### 结果计算函数 - cal_result_anova_crd

生成完全随机设计方差分析统计分析的完整报告字典。此函数整合了完全随机设计方差分析的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的方差分析结果。

**参数**：
- `param`: AnovaParamCRD对象，包含完全随机设计方差分析所需参数

**返回值**：
- 包含完全随机设计方差分析统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `descriptive_statistics`: 描述性统计信息
  - `anova_statistics`: 方差分析统计量信息
  - `significance_tests`: 显著性检验结果
  - `pairwise_comparisons`: 多重比较结果
  - `interpretation`: 结果的专业解释

## 结果解读注意事项

1. **F统计量解释**：衡量处理间变异与误差变异的比值
2. **P值解释**：在零假设（各处理组均值相等）成立的情况下，观察到当前或更极端结果的概率
3. **显著性水平**：通常使用α=0.05作为判断标准
4. **效应量解释**：评估处理效应的实际意义，而不仅仅是统计显著性
5. **临床意义**：统计学显著性不等同于临床重要性，需结合实际意义解读

## 使用示例

```python
from app.stats.anova.anova_stats_crd import cal_result_anova_crd
from app.schemas.request_data.anova_param import AnovaParamCRD
from app.schemas.request_data.stats_data import StatsData

# 创建示例数据，三个处理组
param = AnovaParamCRD(
    stats_data_list=[
        StatsData(field_name="处理组1", data_list=[10.2, 12.1, 9.8, 11.5, 10.9]),
        StatsData(field_name="处理组2", data_list=[13.2, 14.5, 12.8, 13.9, 14.1]),
        StatsData(field_name="处理组3", data_list=[8.7, 9.2, 8.9, 9.5, 8.3])
    ]
)

# 计算完全随机设计方差分析
result = cal_result_anova_crd(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"描述性统计: {result['descriptive_statistics']}")
print(f"方差分析统计量: {result['anova_statistics']}")
print(f"显著性检验: {result['significance_tests']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范