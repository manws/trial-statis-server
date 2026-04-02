# 生存率比较分析 - 两组比较 (Survival Rate Comparison Analysis - Two Groups)

## 概述
两组生存率比较分析使用Log-rank检验和Breslow检验来比较两组生存曲线的差异。这些方法广泛应用于临床试验中，用于比较不同治疗方案或不同患者群体的生存差异。

## 数学原理
### Log-rank检验
- 检验统计量：χ² = (O₁-E₁)²/V，其中O₁为观察事件数，E₁为期望事件数，V为方差
- 对晚期差异更敏感，适合比例风险假设成立的情况

### Breslow检验（广义Wilcoxon检验）
- 加权版本的检验方法，使用权重等于风险集大小
- 对早期差异更敏感，适合早期差异较大的情况

## 功能特点
1. **Log-rank检验**：对晚期差异敏感的非参数检验方法
2. **Breslow检验**：对早期差异敏感的非参数检验方法（广义Wilcoxon检验）
3. **生存曲线比较**：比较两组在各时间点的生存概率
4. **敏感性分析**：分析早期与晚期差异的敏感性
5. **统计推断**：提供显著性检验和置信区间

## 适用条件
### Log-rank检验
- 生存曲线比例风险假设成立
- 对晚期差异更敏感
- 适合风险比恒定的情况

### Breslow检验
- 对早期差异更敏感
- 适合早期差异较大的情况
- 广义Wilcoxon检验的变体

## 临床应用场景
- 临床试验中的主要终点分析
- 随机对照试验的疗效比较
- 预后因子的单因素分析
- 不同治疗组的生存比较

## 输出指标说明
- **table_name**：报告表格名称，固定为"两组生存率比较分析"
- **input_parameters**：输入参数信息
  - group_a_name: A组名称，字符串类型
  - group_b_name: B组名称，字符串类型
  - a_sample_size: A组样本量，整数类型
  - b_sample_size: B组样本量，整数类型
  - a_events: A组事件数，整数类型
  - b_events: B组事件数，整数类型
  - common_time_points: 公共时间点列表，浮点数列表类型
- **time_point_analysis**：时间点分析
  - times: 时间点列表，浮点数列表类型
  - group_a: A组统计信息列表
    - time: 时间点，浮点数类型
    - at_risk: 风险集中人数，整数类型
    - actual_deaths: 实际死亡数，整数类型
    - theoretical_deaths: 理论死亡数，浮点数类型
  - group_b: B组统计信息列表
    - time: 时间点，浮点数类型
    - at_risk: 风险集中人数，整数类型
    - actual_deaths: 实际死亡数，整数类型
    - theoretical_deaths: 理论死亡数，浮点数类型
- **statistical_tests**：统计检验结果
  - log_rank: Log-rank检验结果
    - chi_square_statistic: 卡方统计量，浮点数类型
    - degrees_of_freedom: 自由度，整数类型
    - p_value: P值，浮点数类型
    - observed_events: 观察事件数列表，浮点数列表类型
    - expected_events: 期望事件数列表，浮点数列表类型
    - method_name: 方法名称，字符串类型
    - interpretation: 方法解释，字符串类型
  - breslow: Breslow检验结果
    - chi_square_statistic: 卡方统计量，浮点数类型
    - degrees_of_freedom: 自由度，整数类型
    - p_value: P值，浮点数类型
    - observed_events: 观察事件数列表，浮点数列表类型
    - expected_events: 期望事件数列表，浮点数列表类型
    - method_name: 方法名称，字符串类型
    - interpretation: 方法解释，字符串类型
- **sensitivity_analysis**：敏感性分析
  - early_difference: 早期差异，整数类型
  - late_difference: 晚期差异，整数类型
  - pattern: 差异模式，字符串类型
  - interpretation: 方法建议，字符串类型
- **interpretation**：结果解释
  - log_rank_description: Log-rank描述，字符串类型
  - breslow_description: Breslow描述，字符串类型
  - sensitivity_interpretation: 敏感性解释，字符串类型
- **remark**：备注信息，字符串类型

## 注意事项
- P值解释：在零假设（两组生存率无差异）成立的情况下，观察到当前或更极端结果的概率
- 统计显著性：通常使用α=0.05作为判断标准
- 临床意义：统计显著性不等于临床重要性
- 检验选择：根据研究目的选择合适的检验方法

## 示例
给定两组患者的生存时间、死亡数和删失数数据，通过Log-rank检验和Breslow检验比较两组的生存率差异，分析不同治疗方案的疗效。