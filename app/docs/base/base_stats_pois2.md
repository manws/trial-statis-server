# 临床泊松分布总体均数区间估计统计分析模块 (Clinical Poisson Distribution Population Mean Interval Estimation Statistics Module)

## 模块概述

临床泊松分布总体均数区间估计统计分析模块提供全面的泊松分布总体均数区间估计统计分析功能，用于临床数据中总体均数的置信区间计算，是临床试验数据分析的重要组成部分。总体均数的区间估计在医学研究中广泛应用，如估计某种罕见疾病的年发病率、某种不良事件的发生率、某种病原体的检出率等的置信区间，为临床决策提供精确的均数范围。

## 模块功能

### 1. 精确区间估计
- 使用泊松分布的精确方法计算置信区间
- 适用于所有样本量

### 2. 正态近似区间
- 使用正态近似法计算置信区间
- 适用于大均数情况

### 3. 精度评估
- 评估估计的精度和可靠性
- 计算所需的观察单位

### 4. 结果解释
- 提供统计和临床意义的解释
- 比较不同方法的结果

## 临床应用价值

- **疾病监测**：估计罕见疾病的发病率置信区间
- **安全评估**：估计不良事件发生率的置信区间
- **质量控制**：估计医疗差错发生率的置信区间
- **资源规划**：估计急诊就诊人数的置信区间

## 统计方法选择指南

### 1. 精确方法适用条件
- 适用于所有样本量
- 特别适用于小样本或均数较小的情况

### 2. 正态近似法适用条件
- 均数足够大（通常λ ≥ 5）
- 一般适用于观察单位足够大的情况

### 3. 临床应用场景
- 罕见疾病发病率的区间估计
- 不良事件发生率的区间估计
- 感染事件发生率的区间估计

## 主要函数说明

### 精确置信区间函数 - exact_confidence_interval

使用精确方法计算泊松分布总体均数的置信区间。在临床研究中，这用于计算总体均数的精确置信区间，特别适用于小样本或均数较小的情况。

**参数**：
- `count`: 观察到的事件数
- `confidence_level`: 置信水平

**返回值**：
- `Dict[str, float]`: 包含置信区间上下限的字典

### 正态近似区间函数 - normal_approximation_interval

使用正态近似法计算泊松分布总体均数的置信区间。在临床研究中，这用于计算总体均数的置信区间，适用于大均数的情况，是常用的区间估计方法之一。

**参数**：
- `count`: 观察到的事件数
- `confidence_level`: 置信水平

**返回值**：
- `Dict[str, float]`: 包含置信区间上下限的字典

### 样本量计算函数 - calculate_sample_size_for_precision

计算达到特定精度所需的观察单位大小。在临床研究设计中，这用于确定获得指定精度的置信区间所需的观察单位，有助于研究设计阶段的规划。

**参数**：
- `confidence_level`: 置信水平
- `precision`: 所需精度（置信区间半宽）
- `lambda_estimate`: 均数预估值

**返回值**：
- `float`: 所需的观察单位大小

### 结果计算函数 - cal_result_pois2

生成泊松分布总体均数区间估计统计分析的完整报告字典。此函数整合了泊松分布总体均数区间估计的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的区间估计分析结果。

**参数**：
- `param`: BaseParamPois2对象，包含confidence_level, total_avg, sample_posi_num

**返回值**：
- 包含泊松分布总体均数区间估计统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `exact_method`: 精确法置信区间
  - `normal_approximation_method`: 正态近似法置信区间
  - `interval_widths`: 各方法的区间宽度
  - `validation`: 有效性验证
  - `required_sample_size_for_precision`: 达到精度所需的样本量
  - `interpretation`: 结果解释

## 结果解读注意事项

1. **置信区间解释**：95%置信区间表示真实总体均数有95%的可能性落在该区间内
2. **区间宽度解释**：区间越窄，估计越精确
3. **临床意义**：置信区间应结合临床实际意义进行解读
4. **样本量影响**：观察单位越大，置信区间越窄
5. **临床决策**：置信区间可用于临床决策的参考

## 使用示例

```python
from app.stats.base.base_stats_pois2 import cal_result_pois2
from app.schemas.request_data.base_param import BaseParamPois2

# 创建参数对象
param = BaseParamPois2(
    confidence_level=0.95,  # 95%置信水平
    total_avg=2.0,          # 总体均数为2.0
    sample_posi_num=15      # 样本阳性数（观察到的事件数）为15
)

# 计算泊松分布总体均数区间估计
result = cal_result_pois2(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"精确法置信区间: {result['exact_method']}")
print(f"正态近似法置信区间: {result['normal_approximation_method']}")
print(f"区间宽度: {result['interval_widths']}")
print(f"有效性验证: {result['validation']}")
print(f"结果解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范