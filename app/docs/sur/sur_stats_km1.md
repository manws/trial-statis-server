# Kaplan-Meier生存分析 - 基于原始数据 (Kaplan-Meier Survival Analysis - Based on Raw Data)

## 概述
基于原始数据的Kaplan-Meier生存分析是一种非参数统计方法，用于估计生存函数和计算生存概率。该方法适用于存在删失数据的生存分析，通过将时间轴划分为事件发生的时间点，计算每个时间点的生存概率。

## 数学原理
Kaplan-Meier估计基于以下统计量计算：
- 生存概率：S(t) = ∏(1 - di/ni)，其中di为时间ti的死亡数，ni为时间ti的风险集大小
- 标准误：SE[S(t)] = S(t)√∑[di/(ni(ni-di))] （Greenwood公式）
- 置信区间：使用log(-log)变换计算95%置信区间

## 功能特点
1. **生存概率估计**：计算各时间点的生存概率
2. **标准误计算**：使用Greenwood公式计算生存概率的标准误
3. **置信区间估计**：计算生存概率的置信区间
4. **中位生存时间**：估计中位生存时间
5. **生存曲线绘制**：提供绘制生存曲线所需的数据

## 适用条件
- 数据为时间到事件型数据
- 存在右删失数据
- 事件发生时间准确记录
- 样本量适中至较大

## 临床应用场景
- 评估癌症患者的生存时间
- 比较手术与保守治疗的效果
- 评估药物长期疗效
- 评估医疗器械使用寿命

## 输出指标说明
- **table_name**：报告表格名称，固定为"Kaplan-Meier生存分析报告"
- **input_parameters**：输入参数信息
  - n_groups: 分析组数，整数类型
  - group_names: 组名称列表，字符串列表类型
  - analysis_type: 分析类型，字符串类型
- **km_results**：Kaplan-Meier分析结果
  - group_name: 组名称，字符串类型
  - n_total: 总样本数，整数类型
  - n_events: 事件数，整数类型
  - n_censored: 删失数，整数类型
  - unique_times: 唯一时间点列表，浮点数列表类型
  - survival_probabilities: 生存概率列表，浮点数列表类型
  - standard_errors: 标准误列表，浮点数列表类型
  - confidence_lower: 置信区间下限列表，浮点数列表类型
  - confidence_upper: 置信区间上限列表，浮点数列表类型
  - at_risk_counts: 风险集中人数列表，整数列表类型
  - event_counts: 事件数列表，整数列表类型
  - censored_counts: 删失数列表，整数列表类型
  - median_survival: 中位生存时间，浮点数类型或None
- **summary_statistics**：汇总统计
  - total_groups: 总组数，整数类型
  - group_details: 组详细信息列表
    - group_name: 组名称，字符串类型
    - sample_size: 样本大小，整数类型
    - events: 事件数，整数类型
    - censored: 删失数，整数类型
    - median_survival: 中位生存时间，浮点数类型或None
- **remark**：备注信息，字符串类型

## 注意事项
- 生存概率解释：表示在特定时间点仍存活的概率
- 置信区间解释：提供生存概率估计的不确定性范围
- 中位生存时间：50%的个体存活的时间点
- 曲线下降速度：反映事件发生的速度，越陡峭表示事件发生越快

## 示例
给定一组癌症患者的生存时间数据，通过Kaplan-Meier方法估计生存概率，分析不同治疗方案的生存差异。