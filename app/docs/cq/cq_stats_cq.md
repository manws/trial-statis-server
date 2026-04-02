# 临床Cochran's Q检验统计分析模块 (Clinical Cochran's Q Test Statistics Module)

## 模块概述

临床Cochran's Q检验统计分析模块提供全面的Cochran's Q检验统计分析功能，用于临床数据中多个相关样本的二分类数据检验，是一种重要的非参数统计分析方法。Cochran's Q检验用于判断多个相关样本（如不同治疗方法）的二分类结果（成功/失败）是否存在显著差异，广泛应用于医学研究中的治疗效果比较、诊断方法评价、干预措施评估等领域。

## 模块功能

### 1. Cochran's Q统计量计算
- 计算Cochran's Q检验统计量
- 衡量多个相关样本间的差异程度

### 2. 显著性检验
- 执行卡方检验判断Q统计量的显著性
- 提供P值评估统计显著性

### 3. P值计算
- 计算检验的P值
- 评估结果的统计学意义

### 4. 统计解释
- 提供结果的临床意义解释
- 帮助理解检验结果的实际含义

## 临床应用价值

- **治疗效果比较**：比较多种治疗方法的效果差异
- **诊断方法评价**：评估不同诊断方法的准确性差异
- **干预措施评估**：分析不同干预措施的有效性差异
- **研究假设验证**：验证多个相关样本是否存在显著差异

## 统计方法选择指南

### 1. Cochran's Q检验适用条件
- 数据为二分类（通常编码为0和1）
- 样本为相关样本（如同一组受试者接受多种处理）
- 各处理组样本量相同
- 观测值独立

### 2. 临床应用场景
- 比较三种或以上不同药物治疗同一疾病的疗效
- 评估多个诊断测试的准确性差异
- 分析不同手术方式的治疗效果差异
- 研究多种康复方案的疗效对比

## 主要函数说明

### Cochran's Q检验函数 - cochran_q_test

Cochran's Q检验（适用于二分类数据的多个相关样本检验）。

**参数**：
- `data_matrix`: List[List[int]]，数据矩阵，每一行代表一个受试者，每一列代表一种处理/条件，数据应为二分类（0或1）

**返回值**：
- 包含检验统计量和结果的字典

### 结果计算函数 - cal_result_cq

生成Cochran's Q检验统计分析的完整报告字典。此函数整合了Cochran's Q检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的检验结果。

**参数**：
- `data_matrix`: List[List[int]]，数据矩阵，每一行代表一个受试者，每一列代表一种处理/条件，数据应为二分类（0或1）

**返回值**：
- 包含Cochran's Q检验统计分析指标的字典，包括：
  - `input_parameters`: 输入参数信息
  - `treatment_statistics`: 处理统计信息
  - `subject_statistics`: 受试者统计信息
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释

## 结果解读注意事项

1. **Q统计量解释**：Q值越大，表明各处理组间差异越显著
2. **P值解释**：P值小于显著性水平（如0.05）时拒绝原假设
3. **临床意义**：统计显著性不等于临床意义，需结合专业知识进行解释

## 使用示例

```python
from app.stats.cq.cq_stats_cq import cal_result_cq

# 创建数据矩阵（行表示受试者，列表示处理，值为0或1）
data_matrix = [
    [1, 1, 0],  # 受试者1：处理1成功，处理2成功，处理3失败
    [1, 0, 0],  # 受试者2：处理1成功，处理2失败，处理3失败
    [0, 1, 1],  # 受试者3：处理1失败，处理2成功，处理3成功
    [1, 1, 1],  # 受试者4：处理1成功，处理2成功，处理3成功
    [0, 0, 1]   # 受试者5：处理1失败，处理2失败，处理3成功
]

# 执行Cochran's Q检验
result = cal_result_cq(data_matrix)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"处理统计信息: {result['treatment_statistics']}")
print(f"受试者统计信息: {result['subject_statistics']}")
print(f"检验统计量: {result['test_statistics']}")
print(f"显著性检验: {result['significance_tests']}")
print(f"统计解释: {result['interpretation']}")
print(f"Q值: {result['test_statistics']['q_value']}")
print(f"P值: {result['test_statistics']['p_value_two_sided']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范