# 临床Friedman M检验统计分析模块 (Clinical Friedman M Test Statistics Module)

## 模块概述

临床Friedman M检验统计分析模块提供全面的Friedman M检验统计分析功能，用于临床数据中随机区组设计的等级资料比较，是一种重要的非参数统计分析方法。Friedman M检验通过对区组内各处理的秩次进行比较，来判断不同处理之间是否存在显著差异，广泛应用于医学研究中的配伍组设计、重复测量分析、多处理比较等领域。

## 模块功能

### 1. 秩和计算
- 计算各处理的秩和
- 衡量数据的分布差异

### 2. 检验统计量
- 计算Friedman检验统计量M
- 衡量处理间差异的显著性

### 3. 显著性检验
- 执行卡方近似检验判断统计量的显著性
- 提供P值评估统计显著性

### 4. Q事后检验
- 当主检验显著时，进行成对比较
- 提供更精细的处理间差异分析

### 5. 统计解释
- 提供结果的临床意义解释
- 帮助理解检验结果的实际含义

## 临床应用价值

- **配伍组设计**：比较同一区组内不同处理的效果
- **重复测量**：分析同一受试者在不同时间点的差异
- **多处理比较**：评估多种治疗方法的差异
- **非参数分析**：处理不满足正态性假设的数据

## 统计方法选择指南

### 1. Friedman检验适用条件
- 数据为连续型或有序分类变量
- 同一区组内数据相关
- 不同区组间数据独立
- 至少有三种处理

### 2. 临床应用场景
- 比较三种或以上不同治疗方案在同一群组患者中的效果
- 评估同一患者在不同时间点的指标变化
- 分析不同药物对同一组患者的影响
- 研究不同护理方法在相同条件下效果差异

## 主要函数说明

### Friedman M检验函数 - friedman_m_test

Friedman M秩和检验（适用于随机区组设计的等级资料）。

**参数**：
- `block_data`: List[List[float]]，区组数据，每个子列表代表一个区组，包含各处理的观测值

**返回值**：
- 包含检验统计量和结果的字典

### Q事后检验函数 - friedman_q_post_hoc_test

Friedman检验后的Q检验（用于多重比较）。

**参数**：
- `block_data`: List[List[float]]，区组数据
- `alpha`: float，显著性水平，默认0.05

**返回值**：
- 包含事后检验结果的字典

### 结果计算函数 - cal_result_rs_fm

生成Friedman M检验及Q检验统计分析的完整报告字典。此函数整合了Friedman M检验及Q检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的检验结果。

**参数**：
- `block_data`: List[List[float]]，区组数据，每个子列表代表一个区组，包含各处理的观测值

**返回值**：
- 包含Friedman M检验及Q检验统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `treatment_statistics`: 处理统计信息
  - `rank_statistics`: 秩统计信息
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释
  - `q_results`: Q事后检验结果（如果主检验显著）

## 结果解读注意事项

1. **检验统计量解释**：M值越大，表明处理间差异越显著
2. **P值解释**：P值小于显著性水平（如0.05）时拒绝原假设
3. **临床意义**：统计显著性不等于临床意义，需结合专业知识进行解释

## 使用示例

```python
from app.stats.rs.rs_stats_fm import cal_result_rs_fm

# 创建区组数据（例如，3个处理，5个区组）
block_data = [
    [120, 125, 130],  # 第一个区组的三个处理
    [115, 120, 125],  # 第二个区组的三个处理
    [130, 135, 140],  # 第三个区组的三个处理
    [125, 130, 135],  # 第四个区组的三个处理
    [118, 123, 128]   # 第五个区组的三个处理
]

# 执行Friedman M检验
result = cal_result_rs_fm(block_data)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"处理统计信息: {result['treatment_statistics']}")
print(f"秩统计信息: {result['rank_statistics']}")
print(f"检验统计量: {result['test_statistics']}")
print(f"显著性检验: {result['significance_tests']}")
print(f"统计解释: {result['interpretation']}")
print(f"M统计量: {result['test_statistics']['chi_square_value']}")
print(f"P值: {result['test_statistics']['p_value']}")
if result['q_results']:
    print(f"Q事后检验: {result['q_results']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范