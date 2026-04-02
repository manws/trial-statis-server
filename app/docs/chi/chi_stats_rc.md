# 临床RC列联表卡方检验统计分析模块 (Clinical R×C Contingency Table Chi-Square Test Statistics Module)

## 模块概述

临床RC列联表卡方检验统计分析模块提供全面的R×C列联表卡方检验统计分析功能，用于临床数据中多个分类变量间的关联性检验，是临床试验数据分析的重要组成部分。RC列联表卡方检验在医学研究中广泛应用，如比较多个治疗组的分类疗效指标、分析多个分层因素与结局的关系等。

## 模块功能

### 1. R×C列联表卡方检验
- 计算R×C列联表的卡方统计量
- 评估多个分类变量间的关联性

### 2. 期望频数计算
- 计算每个单元格的期望频数
- 评估检验的有效性

### 3. 效应量计算
- 计算Cramér's V系数
- 计算Pearson列联系数

### 4. 结果解释
- 提供统计和临床意义的解释
- 评估关联强度

### 5. 适用性检验
- 验证理论频数是否满足检验要求
- 提供替代方案建议

## 临床应用价值

- **多组比较**：比较多个治疗组的分类疗效指标
- **分层分析**：分析不同分层因素与结局的关系
- **关联评估**：评估多个分类变量间的关联性
- **流行病学研究**：分析多个风险因素与疾病的关系

## 统计方法选择指南

### 1. R×C列联表卡方检验适用条件
- 样本量足够大
- 理论频数≥5的单元格比例≥80%
- 不允许有理论频数<1

### 2. 临床应用场景
- 比较三种或以上治疗方法的疗效
- 分析不同亚组的治疗效果
- 评估多个风险因素与结局的关系

## 主要函数说明

### R×C卡方检验函数 - chi_square_test_rc

执行RC列联表的卡方检验。在临床研究中，这用于比较多个分类变量之间的关联性，以评估其统计学显著性。

**参数**：
- `row_num`: 行数
- `col_num`: 列数
- `data_list`: 数据列表，按行优先顺序排列

**返回值**：
- `Dict[str, Any]`: 检验结果字典

### 完整检验流程函数 - perform_chi_square_rc_test

执行完整的RC列联表卡方检验流程。在临床研究中，这用于提供全面的R×C列联表分析结果，包括卡方检验、效应量计算和关联强度评估。

**参数**：
- `row_num`: 行数
- `col_num`: 列数
- `data_list`: 数据列表，按行优先顺序排列

**返回值**：
- `Dict[str, Any]`: 完整的检验结果字典

### 结果计算函数 - cal_result_chi_rc

生成RC列联表卡方检验统计分析的完整报告字典。此函数整合了RC列联表卡方检验的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的卡方检验结果。报告包括输入参数、观察频数表、期望频数分析、检验统计量、显著性检验和统计解释等信息，便于临床医生和研究人员快速理解卡方检验的特征。

**参数**：
- `param`: ChiSquareParamRC对象，包含row_num, col_num, data_list

**返回值**：
- 包含RC列联表卡方检验统计分析指标的字典，包括：
  - `table_name`: 报告表格名称
  - `input_parameters`: 输入参数信息
  - `contingency_table`: 观察频数表
  - `expected_frequencies`: 期望频数分析
  - `test_statistics`: 检验统计量
  - `significance_tests`: 显著性检验结果
  - `interpretation`: 统计解释
  - `remark`: 备注信息

## 结果解读注意事项

1. **P值解释**：在零假设成立的前提下，观察到当前或更极端结果的概率
2. **检验有效性**：检查理论频数是否满足检验要求
3. **显著性水平**：通常使用α=0.05作为判断标准
4. **临床意义**：统计显著性不等同于临床重要性
5. **关联强度**：关注Cramér's V等效应量指标

## 使用示例

```python
from app.stats.chi.chi_stats_rc import cal_result_chi_rc
from app.schemas.request_data.chi_param import ChiSquareParamRC

# 创建参数对象
param = ChiSquareParamRC(
    row_num=2,                    # 行数
    col_num=3,                    # 列数
    data_list=[20, 15, 10, 25, 20, 15]  # 数据列表，按行优先排列
)

# 计算RC列联表卡方检验
result = cal_result_chi_rc(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"列联表: {result['contingency_table']}")
print(f"检验统计量: {result['test_statistics']}")
print(f"统计解释: {result['interpretation']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范