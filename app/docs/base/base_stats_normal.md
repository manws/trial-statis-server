# 临床正态分布统计分析模块 (Clinical Normal Distribution Statistics Module)

## 模块概述

临床正态分布统计分析模块提供全面的正态分布统计分析功能，用于临床数据中正态分布特征的分析和检验，是临床试验数据分析的重要组成部分。正态分布在医学研究中具有核心地位，因为许多生物医学指标都近似服从正态分布，且许多统计方法都基于正态分布假设，是临床统计分析的基础。

## 模块功能

### 1. 正态性检验
- Shapiro-Wilk检验、Kolmogorov-Smirnov检验等
- 提供多种检验方法

### 2. 参数估计
- 均数和标准差的估计
- 参数的置信区间

### 3. 概率计算
- 给定值的概率密度和累积概率
- 分位数计算

### 4. 分布拟合
- 评估数据与正态分布的拟合程度
- 拟合优度检验

### 5. 参考区间
- 建立正常值的参考范围
- 临床决策支持

## 临床应用价值

- **数据预处理**：判断是否满足正态分布假设
- **方法选择**：指导后续统计分析方法的选择
- **质量控制**：监控医疗指标是否在正常范围内
- **参考区间**：建立正常值的参考范围

## 统计方法选择指南

### 1. Shapiro-Wilk检验适用条件
- 小样本数据（n < 5000）
- 完整数据（无缺失值）

### 2. Kolmogorov-Smirnov检验适用条件
- 大样本数据
- 可用于与其他理论分布比较

### 3. 临床应用场景
- 判断实验室指标是否符合正态分布
- 选择合适的统计检验方法
- 建立生理指标的参考范围

## 主要函数说明

### 正态性检验函数 - normality_test

进行正态性检验。在临床研究中，这用于检验数据是否符合正态分布，以决定后续统计分析方法的选择。

**参数**：
- `data`: List[float]，输入数据列表

**返回值**：
- `Dict[str, Any]`: 包含正态性检验结果的字典

### 正态分布参数估计函数 - estimate_normal_parameters

估计正态分布的参数（均值和标准差）。在临床研究中，这用于估计数据的正态分布参数，以描述数据的中心趋势和变异程度。

**参数**：
- `data`: List[float]，输入数据列表

**返回值**：
- `Dict[str, float]`: 包含参数估计值的字典

### 概率计算函数 - calculate_normal_probabilities

计算给定值的正态分布概率。在临床研究中，这用于计算特定值在正态分布下的概率，以评估其出现的可能性。

**参数**：
- `data`: List[float]，输入数据列表
- `values`: List[float]，要计算概率的值列表（如果为None，则使用数据的均值和极值）

**返回值**：
- `Dict[str, Any]`: 包含概率计算结果的字典

### 结果计算函数 - cal_result_normal

生成正态分布统计分析的完整报告字典。此函数整合了正态分布的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的正态分布分析结果。

**参数**：
- `param`: BaseParamNormal对象，包含stats_data_list

**返回值**：
- 包含正态分布统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `normality_tests`: 正态性检验结果
  - `parameter_estimates`: 参数估计结果
  - `probability_calculations`: 概率计算结果
  - `interpretation`: 结果解释
  - `remark`: 备注信息

## 结果解读注意事项

1. **P值解释**：在正态分布假设成立的前提下，观察到当前或更极端结果的概率
2. **检验功效**：小样本可能无法检测出轻微的偏离
3. **实际意义**：轻微偏离正态分布不一定影响分析结果
4. **可视化**：结合直方图、Q-Q图等图形工具判断
5. **临床意义**：统计显著性不等同于临床重要性

## 使用示例

```
from app.stats.base.base_stats_normal import cal_result_normal
from app.schemas.request_data.base_param import BaseParamNormal
from app.schemas.request_data.stats_data import StatsData

# 创建统计数据对象列表
data_list = [
    StatsData(field_name="身高", data_list=[160, 165, 170, 175, 180, 185, 190]),
    StatsData(field_name="体重", data_list=[60, 65, 70, 75, 80, 85, 90])
]

# 创建参数对象
param = BaseParamNormal(
    stats_data_list=data_list
)

# 计算正态分布统计分析
result = cal_result_normal(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"正态性检验: {result['normality_tests']}")
print(f"参数估计: {result['parameter_estimates']}")
print(f"概率计算: {result['probability_calculations']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范