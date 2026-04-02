# 临床频数分布统计分析模块 (Clinical Frequency Distribution Statistics Module)

## 模块概述

临床频数分布统计分析模块提供全面的频数分布统计分析功能，用于临床数据中数据分布情况的可视化和分析，是临床试验数据分析的基础组成部分。频数分布在医学研究中广泛应用，如分析患者的年龄分布、疾病严重程度分布、实验室指标的分布等，为后续的数据分析和临床决策提供直观的分布信息。

## 模块功能

### 1. 频数统计
- 计算各区间或类别的观测频数
- 提供分组统计

### 2. 频率计算
- 计算各区间或类别的相对频率
- 计算百分比

### 3. 累积频数
- 计算累积频数和累积频率
- 提供累积分布信息

### 4. 直方图
- 生成直方图数据
- 可视化分布形状

### 5. 分布拟合
- 评估数据是否符合特定分布
- 提供拟合优度检验

## 临床应用价值

- **数据探索**：了解数据的分布特征
- **异常值检测**：识别数据中的异常值
- **分组分析**：为后续统计分析提供分组依据
- **报告展示**：为研究结果提供直观展示

## 统计方法选择指南

### 1. 等距分组适用条件
- 数据分布相对均匀
- 数值型变量

### 2. 不等距分组适用条件
- 数据分布极度偏斜
- 需要突出某些特殊区间

### 3. 临床应用场景
- 患者年龄分布分析
- 实验室指标分布分析
- 疾病严重程度分布

## 主要函数说明

### 频数分布计算函数 - calculate_frequency_distribution

计算频数分布。在临床研究中，这用于分析数据在不同区间内的分布情况，以揭示数据的分布特征和规律。

**参数**：
- `data`: List[float]，输入数据列表
- `num_classes`: int，分组数（如果为None则自动计算）

**返回值**：
- `Dict[str, Any]`: 包含频数分布统计量的字典

### 分布拟合评估函数 - assess_distribution_fit

评估数据与特定分布的拟合优度。在临床研究中，这用于评估数据是否符合特定的理论分布，以指导后续的统计分析方法选择。

**参数**：
- `data`: List[float]，输入数据列表
- `distribution`: str，要检验的分布类型

**返回值**：
- `Dict[str, Any]`: 包含拟合优度检验结果的字典

### 结果计算函数 - cal_result_freq

生成频数分布统计分析的完整报告字典。此函数整合了频数分布的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的频数分布分析结果。

**参数**：
- `param`: BaseParamFreq对象，包含stats_data_list

**返回值**：
- 包含频数分布统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `frequency_distributions`: 各数据集的频数分布
  - `distribution_assessments`: 分布拟合评估
  - `interpretation`: 结果解释
  - `remark`: 备注信息

## 结果解读注意事项

1. **频数解释**：各区间或类别的观测个数
2. **频率解释**：各区间或类别的占比
3. **累积频率解释**：累计占比
4. **分布形状**：对称、偏斜等
5. **临床意义**：结合专业背景解释分布

## 使用示例

```
from app.stats.base.base_stats_freq import cal_result_freq
from app.schemas.request_data.base_param import BaseParamFreq
from app.schemas.request_data.stats_data import StatsData

# 创建统计数据对象列表
data_list = [
    StatsData(field_name="年龄", data_list=[25, 30, 35, 40, 45, 50, 55, 60, 65]),
    StatsData(field_name="血压", data_list=[120, 125, 130, 135, 140, 145, 150, 155, 160])
]

# 创建参数对象
param = BaseParamFreq(
    stats_data_list=data_list
)

# 计算频数分布
result = cal_result_freq(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"频数分布: {result['frequency_distributions']}")
print(f"分布拟合评估: {result['distribution_assessments']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范