# Kaplan-Meier生存分析 - 基于汇总数据 (Kaplan-Meier Survival Analysis - Based on Summary Data)

## 概述
基于汇总数据的Kaplan-Meier生存分析是一种非参数统计方法，用于估计生存函数和计算生存概率。与基于原始数据的方法相比，此模块适用于只有汇总统计信息（时间点、死亡数、删失数）的场景，特别适合数据隐私保护要求较高的研究。

## 数学原理
Kaplan-Meier估计基于汇总数据的计算公式：
- 生存概率：S(t) = ∏(1 - di/ni)，其中di为时间ti的死亡数，ni为时间ti的风险集大小
- 标准误：SE[S(t)] = S(t)√∑[di/(ni(ni-di))] （Greenwood公式）
- 置信区间：使用log(-log)变换计算95%置信区间

## 功能特点
1. **生存概率估计**：计算各时间点的生存概率
2. **标准误计算**：使用Greenwood公式计算生存概率的标准误
3. **置信区间估计**：计算生存概率的置信区间
4. **中位生存时间**：估计中位生存时间
5. **百分位数生存时间**：计算25%、50%、75%生存时间

## 适用条件
- 具备各时间点的生存时间、死亡数、删失数
- 数据完整性良好，无重要缺失
- 时间点划分合理，不过于稀疏或密集

## 临床应用场景
- 多中心生存数据汇总分析
- 已发表研究的再分析
- 机构内部数据汇总报告
- 隐私受限环境下的数据分析

## 输出指标说明
- **table_name**：报告表格名称，固定为"Kaplan-Meier生存分析报告（汇总数据）"
- **input_parameters**：输入参数信息
  - n_time_points: 时间点数，整数类型
  - total_subjects: 总样本数，整数类型
  - total_events: 总事件数，整数类型
  - total_censored: 总删失数，整数类型
  - time_points: 时间点列表，浮点数列表类型
- **survival_analysis**：生存分析结果
  - survival_probabilities: 生存概率列表，浮点数列表类型
  - standard_errors: 标准误列表，浮点数列表类型
  - confidence_lower: 置信区间下限列表，浮点数列表类型
  - confidence_upper: 置信区间上限列表，浮点数列表类型
  - at_risk_counts: 风险集中人数列表，浮点数列表类型
  - event_counts: 事件数列表，浮点数列表类型
  - censored_counts: 删失数列表，浮点数列表类型
  - median_survival: 中位生存时间，浮点数类型或None
- **summary_statistics**：汇总统计
  - overall_survival_rate: 最终生存率，浮点数类型
  - median_survival_time: 中位生存时间，浮点数类型或None
  - survival_at_percentiles: 百分位数生存时间，字典类型
- **remark**：备注信息，字符串类型

## 注意事项
- 生存概率解释：表示在特定时间点仍存活的概率
- 置信区间解释：提供生存概率估计的不确定性范围
- 中位生存时间：50%的个体存活的时间点
- 百分位数时间：不同生存率对应的时间点
- 数据完整性：汇总数据的完整性影响结果准确性

## 示例
给定某研究的汇总数据（各时间点的生存时间、死亡数、删失数），通过Kaplan-Meier方法估计生存概率，分析该研究的总体生存模式。