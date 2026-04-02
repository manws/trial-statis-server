# 临床描述性统计量统计分析模块 (Clinical Descriptive Statistics Module)

## 模块概述

临床描述性统计量统计分析模块提供全面的描述性统计量计算功能，用于临床数据的基本特征描述和总结，是临床试验数据分析的基础组成部分。描述性统计在医学研究中广泛应用，如计算患者的平均年龄、疾病持续时间的中位数、实验室指标的变异程度等，为后续的推断统计分析奠定基础。

## 模块功能

### 1. 集中趋势指标
- 均数、中位数、众数
- 几何均数、调和均数

### 2. 离散趋势指标
- 方差、标准差、极差
- 四分位间距、变异系数

### 3. 形状指标
- 偏度、峰度

### 4. 位置指标
- 分位数、百分位数

### 5. 区间估计
- 均数的置信区间

## 临床应用价值

- **数据探索**：初步了解数据的分布特征
- **质量控制**：检查数据的合理性和完整性
- **报告撰写**：为研究结果提供基本统计信息
- **研究设计**：为样本量估算提供基础数据

## 统计方法选择指南

### 1. 均数适用条件
- 数据呈正态分布
- 无明显异常值
- 连续型变量

### 2. 中位数适用条件
- 数据呈偏态分布
- 存在异常值
- 等级数据

### 3. 临床应用场景
- 描述患者基本特征
- 汇总实验室检查结果
- 报告治疗效果指标

## 主要函数说明

### 描述性统计计算函数 - calculate_descriptive_stats

计算描述性统计量。在临床研究中，这用于计算数据的基本统计特征，以描述数据的集中趋势、离散趋势和分布形状等特征。

**参数**：
- `data`: List[float]，输入数据列表

**返回值**：
- `Dict[str, float]`: 包含各种描述性统计量的字典

### 结果计算函数 - cal_result_des

生成描述性统计量统计分析的完整报告字典。此函数整合了描述性统计量的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的描述性统计分析结果。

**参数**：
- `param`: BaseParamDes对象，包含stats_data_list

**返回值**：
- 包含描述性统计量统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `descriptive_statistics`: 描述性统计信息
  - `summary`: 统计摘要信息
  - `interpretation`: 结果解释
  - `remark`: 备注信息

## 结果解读注意事项

1. **均数解释**：反映数据的平均水平
2. **标准差解释**：反映数据的离散程度
3. **变异系数解释**：消除量纲影响的变异程度指标
4. **置信区间解释**：总体均数的可能范围
5. **临床意义**：结合专业背景解释统计结果

## 使用示例

```python
from app.stats.base.base_stats_des import cal_result_des
from app.schemas.request_data.base_param import BaseParamDes
from app.schemas.request_data.stats_data import StatsData

# 创建统计数据对象列表
data_list = [
    StatsData(field_name="年龄", data_list=[25, 30, 35, 40, 45, 50, 55]),
    StatsData(field_name="体重", data_list=[60, 65, 70, 75, 80, 85, 90])
]

# 创建参数对象
param = BaseParamDes(
    stats_data_list=data_list
)

# 计算描述性统计量
result = cal_result_des(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"描述性统计: {result['descriptive_statistics']}")
print(f"统计摘要: {result['summary']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范