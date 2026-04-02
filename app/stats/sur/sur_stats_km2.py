"""
Kaplan-Meier生存分析模块 - 基于汇总数据 (Kaplan-Meier Survival Analysis Module - Based on Summary Data)

本模块提供基于汇总数据的Kaplan-Meier生存分析功能，用于估计生存函数和计算生存概率。
与基于原始数据的方法相比，此模块适用于只有汇总统计信息（时间点、死亡数、删失数）的场景，
特别适合数据隐私保护要求较高的研究。

【模块功能概述】:
1. 生存概率估计：计算各时间点的生存概率
2. 标准误计算：使用Greenwood公式计算生存概率的标准误
3. 置信区间估计：计算生存概率的置信区间
4. 中位生存时间：估计中位生存时间
5. 百分位数生存时间：计算25%、50%、75%生存时间

【临床应用价值】:
- 数据汇总分析：当只有汇总数据时进行生存分析
- 多中心研究：整合来自不同中心的汇总数据
- 隐私保护：在不共享个体数据的情况下进行分析
- 文献综述：基于发表文献的汇总数据进行分析

【统计方法选择指南】:
1. Kaplan-Meier基于汇总数据的适用条件：
   - 具备各时间点的生存时间、死亡数、删失数
   - 数据完整性良好，无重要缺失
   - 时间点划分合理，不过于稀疏或密集

2. 临床应用场景：
   - 多中心生存数据汇总分析
   - 已发表研究的再分析
   - 机构内部数据汇总报告
   - 隐私受限环境下的数据分析

【结果解读注意事项】:
1. 生存概率解释：表示在特定时间点仍存活的概率
2. 置信区间解释：提供生存概率估计的不确定性范围
3. 中位生存时间：50%的个体存活的时间点
4. 百分位数时间：不同生存率对应的时间点
5. 数据完整性：汇总数据的完整性影响结果准确性

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "汇总数据可以做Kaplan-Meier分析吗？"
- "如何解释百分位数生存时间？"
- "汇总数据与原始数据分析的区别？"
- "Kaplan-Meier估计的适用条件是什么？"
- "如何处理汇总数据中的缺失值？"

【相关标准和规范】:
- CONSORT 声明：随机临床试验报告规范
- STROBE 声明：观察性研究报告规范
- ICH E9 指导原则：临床试验的统计学原则
- SAMPL 指南：医学研究报告中的统计方法描述规范

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

import numpy as np
from typing import Dict, List, Any


def kaplan_meier_survival_summary_data(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    基于汇总数据的Kaplan-Meier生存分析
    stats_data_list包含三个StatsData对象：
    - field_name="生存时间"的StatsData：data_list为各时间点列表
    - field_name="死亡数"的StatsData：data_list为各时间点死亡数列表  
    - field_name="删失数"的StatsData：data_list为各时间点删失数列表
    
    参数:
    - stats_data_list: 包含三个StatsData对象的列表
    
    返回:
    - 包含生存分析统计量和结果的字典
    """
    # 参数验证
    if not stats_data_list:
        raise ValueError("Kaplan-Meier分析需要提供生存时间、死亡数和删失数数据")
    
    if len(stats_data_list) != 3:
        raise ValueError("Kaplan-Meier分析需要恰好3个StatsData对象：生存时间、死亡数、删失数")
    
    # 验证数据结构
    time_data = None
    death_data = None
    censored_data = None
    
    for stats_data in stats_data_list:
        if stats_data.get("field_name") == '生存时间':
            time_data = stats_data
        elif stats_data.get("field_name") == '死亡数':
            death_data = stats_data
        elif stats_data.get("field_name") == '删失数':
            censored_data = stats_data
        else:
            raise ValueError(f"未知的数据类型: {stats_data.get('field_name')}，应为'生存时间'、'死亡数'或'删失数'")
    
    if time_data is None:
        raise ValueError("缺少生存时间数据")
    
    if death_data is None:
        raise ValueError("缺少死亡数数据")
    
    if censored_data is None:
        raise ValueError("缺少删失数数据")
    
    # 验证数据长度匹配
    if not (len(time_data["data_list"]) == len(death_data["data_list"]) == len(censored_data["data_list"])):
        raise ValueError(f"数据长度不匹配：生存时间({len(time_data['data_list'])})，死亡数({len(death_data['data_list'])})，删失数({len(censored_data['data_list'])})")
    
    # 转换数据
    times = np.array(time_data["data_list"])
    deaths = np.array(death_data["data_list"]).astype(int)
    censored = np.array(censored_data["data_list"]).astype(int)
    
    # 验证数据合理性
    if np.any(deaths < 0) or np.any(censored < 0):
        raise ValueError("死亡数和删失数必须为非负整数")
    
    # 对数据按时间排序
    sort_indices = np.argsort(times)
    sorted_times = times[sort_indices]
    sorted_deaths = deaths[sort_indices]
    sorted_censored = censored[sort_indices]
    
    n_time_points = len(sorted_times)
    total_events = np.sum(sorted_deaths)
    total_censored = np.sum(sorted_censored)
    total_subjects = total_events + total_censored
    
    # 初始化结果数组
    survival_prob = np.ones(n_time_points)
    standard_error = np.zeros(n_time_points)
    confidence_lower = np.ones(n_time_points)
    confidence_upper = np.ones(n_time_points)
    at_risk_counts = np.zeros(n_time_points)
    event_counts = np.zeros(n_time_points)
    censored_counts = np.zeros(n_time_points)
    
    # 计算累积风险和生存概率
    cumulative_hazard = 0.0
    current_at_risk = total_subjects
    
    for i, time_point in enumerate(sorted_times):
        # 设置当前时间点的统计量
        at_risk_counts[i] = current_at_risk
        event_counts[i] = sorted_deaths[i]
        censored_counts[i] = sorted_censored[i]
        
        # 计算该时间点的死亡概率
        if current_at_risk > 0:
            death_prob = sorted_deaths[i] / current_at_risk
            
            # 更新累积风险
            cumulative_hazard += death_prob
            
            # 计算生存概率
            survival_prob[i] = np.exp(-cumulative_hazard)
            
            # 计算标准误（Greenwood公式）
            if survival_prob[i] > 0:
                variance = 0.0
                for j in range(i + 1):
                    risk_j = at_risk_counts[j]
                    events_j = event_counts[j]
                    if risk_j > 0 and events_j > 0:
                        variance += events_j / (risk_j * (risk_j - events_j))
                
                # 避免负方差
                variance = max(variance, 0)
                standard_error[i] = survival_prob[i] * np.sqrt(variance)
                
                # 计算置信区间（使用log(-log)变换）
                if standard_error[i] > 0:
                    z_alpha = 1.96  # 95%置信水平
                    log_log_survival = np.log(-np.log(survival_prob[i]))
                    se_log_log = standard_error[i] / (survival_prob[i] * np.abs(np.log(survival_prob[i])))
                    
                    lower_log_log = log_log_survival - z_alpha * se_log_log
                    upper_log_log = log_log_survival + z_alpha * se_log_log
                    
                    confidence_lower[i] = np.exp(-np.exp(lower_log_log))
                    confidence_upper[i] = np.exp(-np.exp(upper_log_log))
                else:
                    confidence_lower[i] = survival_prob[i]
                    confidence_upper[i] = survival_prob[i]
            else:
                standard_error[i] = 0.0
                confidence_lower[i] = 0.0
                confidence_upper[i] = 0.0
        else:
            survival_prob[i] = 0.0
            standard_error[i] = 0.0
            confidence_lower[i] = 0.0
            confidence_upper[i] = 0.0
        
        # 更新风险集大小
        current_at_risk -= (sorted_deaths[i] + sorted_censored[i])
    
    # 计算中位生存时间
    median_survival = _calculate_median_survival(sorted_times, survival_prob)
    
    return {
        "input_parameters": {
            "n_time_points": int(n_time_points),
            "total_subjects": int(total_subjects),
            "total_events": int(total_events),
            "total_censored": int(total_censored),
            "time_points": sorted_times.tolist()
        },
        "survival_analysis": {
            "survival_probabilities": survival_prob.tolist(),
            "standard_errors": standard_error.tolist(),
            "confidence_lower": confidence_lower.tolist(),
            "confidence_upper": confidence_upper.tolist(),
            "at_risk_counts": at_risk_counts.tolist(),
            "event_counts": event_counts.tolist(),
            "censored_counts": censored_counts.tolist(),
            "median_survival": float(median_survival) if median_survival is not None else None
        },
        "summary_statistics": {
            "overall_survival_rate": float(survival_prob[-1]) if n_time_points > 0 else 0.0,
            "median_survival_time": float(median_survival) if median_survival is not None else None,
            "survival_at_percentiles": _calculate_survival_percentiles(sorted_times, survival_prob)
        }
    }


def _calculate_median_survival(times: np.ndarray, survival_probs: np.ndarray) -> float:
    """
    计算中位生存时间
    
    参数:
    - times: 时间点数组
    - survival_probs: 对应的生存概率数组
    
    返回:
    - 中位生存时间，如果无法计算则返回None
    """
    if len(survival_probs) == 0:
        return None
    
    # 找到生存概率首次降到0.5以下的时间点
    for i, prob in enumerate(survival_probs):
        if prob <= 0.5:
            if i == 0:
                return float(times[0])
            else:
                # 线性插值计算精确的中位生存时间
                prev_time = times[i-1]
                curr_time = times[i]
                prev_prob = survival_probs[i-1]
                curr_prob = survival_probs[i]
                
                if prev_prob > curr_prob:  # 确保概率是递减的
                    fraction = (0.5 - curr_prob) / (prev_prob - curr_prob)
                    median_time = curr_time + fraction * (prev_time - curr_time)
                    return float(median_time)
                else:
                    return float(curr_time)
    
    # 如果生存概率从未降到0.5以下，返回最后一个时间点
    return float(times[-1]) if len(times) > 0 else None


def _calculate_survival_percentiles(times: np.ndarray, survival_probs: np.ndarray) -> Dict:
    """
    计算不同百分位数的生存时间
    
    参数:
    - times: 时间点数组
    - survival_probs: 对应的生存概率数组
    
    返回:
    - 包含各百分位数生存时间的字典
    """
    percentiles = [25, 50, 75]  # 25%, 50%(中位数), 75%
    result = {}
    
    for percentile in percentiles:
        target_prob = percentile / 100.0
        survival_time = None
        
        for i, prob in enumerate(survival_probs):
            if prob <= target_prob:
                if i == 0:
                    survival_time = float(times[0])
                else:
                    # 线性插值
                    prev_time = times[i-1]
                    curr_time = times[i]
                    prev_prob = survival_probs[i-1]
                    curr_prob = survival_probs[i]
                    
                    if prev_prob > curr_prob:
                        fraction = (target_prob - curr_prob) / (prev_prob - curr_prob)
                        survival_time = curr_time + fraction * (prev_time - curr_time)
                        survival_time = float(survival_time)
                    else:
                        survival_time = float(curr_time)
                break
        
        if survival_time is None:
            survival_time = float(times[-1]) if len(times) > 0 else None
            
        result[f"{percentile}%"] = survival_time
    
    return result


def perform_kaplan_meier_survival_summary(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    根据给定参数执行基于汇总数据的Kaplan-Meier生存分析
    
    参数:
    - stats_data_list: 包含生存分析所需参数的列表
    
    返回:
    - 包含Kaplan-Meier生存分析完整结果的字典
    """
    # 验证参数
    if not stats_data_list:
        raise ValueError("生存分析需要提供生存时间、死亡数和删失数数据")
    
    if len(stats_data_list) != 3:
        raise ValueError("生存分析需要恰好3个StatsData对象")
    
    # 执行Kaplan-Meier生存分析
    results = kaplan_meier_survival_summary_data(stats_data_list)
    
    return results


def cal_result_sur_km2(stats_data_list: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    生成基于汇总数据的Kaplan-Meier生存分析的完整报告字典
    
    此函数整合了基于汇总数据的Kaplan-Meier生存分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的生存分析结果。报告包括输入参数、
    生存分析结果和汇总统计等信息，
    便于临床医生和研究人员快速理解生存分析的特征。
    
    Args:
        stats_data_list: 包含生存分析所需参数的列表，必须包含3个字典元素，按顺序为：
            - 第一个字典: field_name='生存时间'，data_list为各时间点列表，浮点数列表类型
            - 第二个字典: field_name='死亡数'，data_list为各时间点死亡数列表，整数列表类型
            - 第三个字典: field_name='删失数'，data_list为各时间点删失数列表，整数列表类型
    
    Returns:
        Dict[str, Any]: 包含生存分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"Kaplan-Meier生存分析报告（汇总数据）"
            - input_parameters: 输入参数信息
                - n_time_points: 时间点数，整数类型
                - total_subjects: 总样本数，整数类型
                - total_events: 总事件数，整数类型
                - total_censored: 总删失数，整数类型
                - time_points: 时间点列表，浮点数列表类型
            - survival_analysis: 生存分析结果
                - survival_probabilities: 生存概率列表，浮点数列表类型
                - standard_errors: 标准误列表，浮点数列表类型
                - confidence_lower: 置信区间下限列表，浮点数列表类型
                - confidence_upper: 置信区间上限列表，浮点数列表类型
                - at_risk_counts: 风险集中人数列表，浮点数列表类型
                - event_counts: 事件数列表，浮点数列表类型
                - censored_counts: 删失数列表，浮点数列表类型
                - median_survival: 中位生存时间，浮点数类型或None
            - summary_statistics: 汇总统计
                - overall_survival_rate: 最终生存率，浮点数类型
                - median_survival_time: 中位生存时间，浮点数类型或None
                - survival_at_percentiles: 百分位数生存时间，字典类型
            - remark: 备注信息，字符串类型
    """
    # 执行Kaplan-Meier生存分析
    results = perform_kaplan_meier_survival_summary(stats_data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Kaplan-Meier生存分析报告（汇总数据）",
        "input_parameters": results["input_parameters"],
        "survival_analysis": results["survival_analysis"],
        "summary_statistics": results["summary_statistics"],
        "remark": f"时间点数: {results['input_parameters']['n_time_points']}, 总样本数: {results['input_parameters']['total_subjects']}"
    }
    
    return result_dict