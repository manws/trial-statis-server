# 临床Wilcoxon符号秩检验统计分析模块 (Clinical Wilcoxon Signed-Rank Test Statistics Module)

## 模块概述

临床Wilcoxon符号秩检验统计分析模块提供全面的Wilcoxon符号秩检验统计分析功能，用于临床数据中配对样本的比较，是一种重要的非参数统计分析方法。Wilcoxon符号秩检验通过比较配对数据的差值的秩次，来判断配对样本是否存在显著差异，广泛应用于医学研究中的前后对照实验、配对设计研究、重复测量分析等领域。

## 模块功能

### 1. 差值计算
- 计算配对样本的差值
- 处理数据的符号和大小

### 2. 符号秩统计
- 计算正负差值的秩和
- 衡量配对数据的差异程度

### 3. 检验统计量
- 计算Wilcoxon检验统计量W
- 衡量配对样本的差异显著性

### 4. 显著性检验
- 执行正态近似检验判断统计量的显著性
- 提供P值评估统计显著性

### 5. 统计解释
- 提供结果的临床意义解释
- 帮助理解检验结果的实际含义

## 临床应用价值

- **前后对照实验**：比较患者治疗前后的指标变化
- **配对设计研究**：比较不同处理在相同受试者上的效果
- **重复测量分析**：分析同一受试者多次测量结果的差异
- **交叉设计分析**：比较交叉试验中不同阶段的治疗效果

## 统计方法选择指南

### 1. Wilcoxon符号秩检验适用条件
- 数据为连续型变量
- 数据分布不要求正态性
- 配对观测值相互独立
- 差值分布对称

### 2. 临床应用场景
- 比较患者治疗前后血压的差异
- 评估某种药物治疗前后血糖水平的变化
- 分析手术前后某项指标的改变
- 研究不同治疗方案在同一患者的疗效差异

## 主要函数说明

### Wilcoxon符号秩检验函数 - wilcoxon_signed_rank_test

Wilcoxon符号秩检验。

**参数**：
- `data1`: List[float]，第一组配对数据
- `data2`: List[float]，第二组配对数据

**返回值**：
- 包含检验统计量和结果的字典

### 结果计算函数 - cal_result_rs_paired

生成Wilcoxon符号秩检验统计分析的完整报告字典。此函数整合了Wilcoxon符号秩检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的检验结果。

**参数**：
- `data1`: List[float]，第一组配对数据
- `data2`: List[float]，第二组配对数据

**返回值**：
- 包含Wilcoxon符号秩检验统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `sample_statistics`: 样本统计信息
  - `rank_statistics`: 秩统计信息
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释

## 结果解读注意事项

1. **检验统计量解释**：W值越小，表明配对样本的差异越显著
2. **P值解释**：P值小于显著性水平（如0.05）时拒绝原假设
3. **临床意义**：统计显著性不等于临床意义，需结合专业知识进行解释

## 使用示例

```python
from app.stats.rs.rs_stats_paired import cal_result_rs_paired

# 创建配对数据
data1 = [120, 125, 130, 128, 122, 135, 127, 124, 126, 132]
data2 = [115, 120, 125, 123, 118, 130, 122, 119, 121, 127]

# 执行Wilcoxon符号秩检验
result = cal_result_rs_paired(data1, data2)

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