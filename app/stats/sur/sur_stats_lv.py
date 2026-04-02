"""
生存率比较分析模块 - 两组比较 (Survival Rate Comparison Module - Two Groups)

本模块提供两组生存率比较分析功能，使用Log-rank检验和Breslow检验来比较两组生存曲线的差异。
这些方法广泛应用于临床试验中，用于比较不同治疗方案或不同患者群体的生存差异。

【模块功能概述】:
1. Log-rank检验：对晚期差异敏感的非参数检验方法
2. Breslow检验：对早期差异敏感的非参数检验方法（广义Wilcoxon检验）
3. 生存曲线比较：比较两组在各时间点的生存概率
4. 敏感性分析：分析早期与晚期差异的敏感性
5. 统计推断：提供显著性检验和置信区间

【临床应用价值】:
- 治疗效果比较：比较不同治疗方法的生存差异
- 预后因素分析：评估不同预后因素对生存的影响
- 临床试验：评估新药或新疗法的有效性
- 疗效评价：客观评价治疗措施的长期效果

【统计方法选择指南】:
1. Log-rank检验适用条件：
   - 生存曲线比例风险假设成立
   - 对晚期差异更敏感
   - 适合风险比恒定的情况

2. Breslow检验适用条件：
   - 对早期差异更敏感
   - 适合早期差异较大的情况
   - 广义Wilcoxon检验的变体

3. 临床应用场景：
   - 临床试验中的主要终点分析
   - 随机对照试验的疗效比较
   - 预后因子的单因素分析
   - 不同治疗组的生存比较

【结果解读注意事项】:
1. P值解释：在零假设（两组生存率无差异）成立的情况下，观察到当前或更极端结果的概率
2. 统计显著性：通常使用α=0.05作为判断标准
3. 临床意义：统计显著性不等于临床重要性
4. 检验选择：根据研究目的选择合适的检验方法

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "如何选择Log-rank还是Breslow检验？"
- "生存率比较结果如何解释？"
- "什么是比例风险假设？"
- "生存分析有哪些注意事项？"
- "如何判断早期或晚期差异？"

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
from scipy import stats
from typing import Dict, List, Any
from app.schemas.request_data.sur_param import SurParamLV


def survival_rate_comparison_two_groups(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    两组生存率比较分析
    stats_data_list包含6个StatsData对象，依次为：
    A组生存时间, A组死亡数, A组删失数, B组生存时间, B组死亡数, B组删失数
    
    参数:
    - stats_data_list: 包含6个StatsData对象的列表
    
    返回:
    - 包含生存率比较统计量和结果的字典
    """
    # 参数验证
    if not stats_data_list:
        raise ValueError("生存率比较分析需要提供完整的数据")
    
    if len(stats_data_list) != 6:
        raise ValueError("生存率比较分析需要恰好6个StatsData对象")
    
    # 按顺序提取数据
    a_times = np.array(stats_data_list[0]["data_list"])  # A组生存时间
    a_deaths = np.array(stats_data_list[1]["data_list"]).astype(int)  # A组死亡数
    a_censored = np.array(stats_data_list[2]["data_list"]).astype(int)  # A组删失数
    b_times = np.array(stats_data_list[3]["data_list"])  # B组生存时间
    b_deaths = np.array(stats_data_list[4]["data_list"]).astype(int)  # B组死亡数
    b_censored = np.array(stats_data_list[5]["data_list"]).astype(int)  # B组删失数
    
    # 验证数据长度一致性
    if not (len(a_times) == len(a_deaths) == len(a_censored)):
        raise ValueError("A组数据长度不一致")
    
    if not (len(b_times) == len(b_deaths) == len(b_censored)):
        raise ValueError("B组数据长度不一致")
    
    # 验证数据合理性
    if np.any(a_deaths < 0) or np.any(a_censored < 0) or np.any(b_deaths < 0) or np.any(b_censored < 0):
        raise ValueError("死亡数和删失数必须为非负整数")
    
    # 计算样本量（期初例数）
    a_sample_size = int(np.sum(a_deaths) + np.sum(a_censored))
    b_sample_size = int(np.sum(b_deaths) + np.sum(b_censored))
    
    # 获取共同时间点
    common_times = _get_common_time_points_for_comparison(a_times, b_times)
    
    # 计算各组在公共时间点的统计量
    a_statistics = _calculate_group_statistics_at_times(a_times, a_deaths, a_censored, common_times)
    b_statistics = _calculate_group_statistics_at_times(b_times, b_deaths, b_censored, common_times)
    
    # 执行Log-rank检验
    log_rank_result = _perform_log_rank_test_two_groups(
        a_times, a_deaths, a_censored,
        b_times, b_deaths, b_censored
    )
    
    # 执行Breslow检验
    breslow_result = _perform_breslow_test_two_groups(
        a_times, a_deaths, a_censored,
        b_times, b_deaths, b_censored
    )
    
    # 判断早期和后期敏感性
    early_late_sensitivity = _analyze_early_late_sensitivity(
        a_statistics, b_statistics, common_times
    )
    
    return {
        "input_parameters": {
            "group_a_name": stats_data_list[0].get("field_name") or "A组",
            "group_b_name": stats_data_list[3].get("field_name") or "B组",
            "a_sample_size": a_sample_size,
            "b_sample_size": b_sample_size,
            "a_events": int(np.sum(a_deaths)),
            "b_events": int(np.sum(b_deaths)),
            "common_time_points": common_times.tolist()
        },
        "time_point_analysis": {
            "times": common_times.tolist(),
            "group_a": a_statistics,
            "group_b": b_statistics
        },
        "statistical_tests": {
            "log_rank": log_rank_result,
            "breslow": breslow_result
        },
        "sensitivity_analysis": early_late_sensitivity,
        "interpretation": {
            "log_rank_description": "Log-rank检验对晚期差异敏感，是最常用的生存曲线比较方法",
            "breslow_description": "Breslow检验（广义Wilcoxon）对早期差异更敏感",
            "sensitivity_interpretation": early_late_sensitivity["interpretation"]
        }
    }


def _expand_summary_data(times: np.ndarray, deaths: np.ndarray, censored: np.ndarray) -> np.ndarray:
    """将汇总数据展开为完整生存数据"""
    expanded_data = []
    
    for i in range(len(times)):
        # 添加死亡事件
        expanded_data.extend([times[i]] * deaths[i])
        # 添加删失事件（用负数表示删失）
        expanded_data.extend([-times[i]] * censored[i])
    
    return np.array(expanded_data)


def _get_common_time_points_for_comparison(a_times: np.ndarray, b_times: np.ndarray) -> np.ndarray:
    """获取两组的公共时间点"""
    # 合并所有时间点并去重排序
    all_times = np.concatenate([a_times, b_times])
    unique_times = np.unique(all_times)
    return np.sort(unique_times)


def _calculate_group_statistics_at_times(times: np.ndarray, deaths: np.ndarray, 
                                       censored: np.ndarray, analysis_times: np.ndarray) -> List[Dict]:
    """计算组在指定时间点的统计量"""
    statistics = []
    
    # 按时间排序
    sort_indices = np.argsort(times)
    sorted_times = times[sort_indices]
    sorted_deaths = deaths[sort_indices]
    sorted_censored = censored[sort_indices]
    
    # 计算期初总人数
    total_subjects = int(np.sum(deaths) + np.sum(censored))
    cumulative_entered = total_subjects
    
    time_index = 0
    
    for analysis_time in analysis_times:
        # 计算期初例数（生存时间 >= 分析时间的个体数）
        at_risk = 0
        period_deaths = 0
        period_censored = 0
        
        # 计算到当前分析时间点为止累计发生的事件数
        cumulative_deaths = 0
        cumulative_censored = 0
        
        # 累计到当前分析时间点
        temp_index = 0
        while temp_index < len(sorted_times) and sorted_times[temp_index] < analysis_time:
            cumulative_deaths += sorted_deaths[temp_index]
            cumulative_censored += sorted_censored[temp_index]
            temp_index += 1
        
        # 期初例数 = 总样本数 - 已发生的事件数
        at_risk = total_subjects - cumulative_deaths - cumulative_censored
        
        # 计算当前时间点的实际事件数
        if temp_index < len(sorted_times) and sorted_times[temp_index] == analysis_time:
            period_deaths = sorted_deaths[temp_index]
            period_censored = sorted_censored[temp_index]
        
        # 计算理论死亡数（基于合并的风险）
        remaining_deaths = np.sum(sorted_deaths[sorted_times >= analysis_time])
        remaining_censored = np.sum(sorted_censored[sorted_times >= analysis_time])
        total_remaining = remaining_deaths + remaining_censored
        
        if total_remaining > 0 and at_risk > 0:
            theoretical_deaths = at_risk * remaining_deaths / total_remaining
        else:
            theoretical_deaths = 0.0
        
        statistics.append({
            "time": float(analysis_time),
            "at_risk": int(at_risk),
            "actual_deaths": int(period_deaths),
            "theoretical_deaths": float(theoretical_deaths)
        })
    
    return statistics


def _perform_log_rank_test_two_groups(a_times: np.ndarray, a_deaths: np.ndarray, a_censored: np.ndarray,
                                    b_times: np.ndarray, b_deaths: np.ndarray, b_censored: np.ndarray) -> Dict:
    """执行两组Log-rank检验"""
    # 合并所有时间点
    all_times = np.concatenate([a_times, b_times])
    unique_times = np.unique(all_times)
    
    observed_a = 0.0
    observed_b = 0.0
    expected_a = 0.0
    expected_b = 0.0
    
    for time_point in unique_times:
        # 计算各组在该时间点的风险集
        a_at_risk = np.sum(a_times >= time_point)
        b_at_risk = np.sum(b_times >= time_point)
        total_at_risk = a_at_risk + b_at_risk
        
        if total_at_risk == 0:
            continue
        
        # 计算该时间点的事件数
        a_events = np.sum(a_deaths[a_times == time_point])
        b_events = np.sum(b_deaths[b_times == time_point])
        total_events = a_events + b_events
        
        if total_events == 0:
            continue
        
        # 计算期望事件数
        expected_a += total_events * a_at_risk / total_at_risk
        expected_b += total_events * b_at_risk / total_at_risk
        
        observed_a += a_events
        observed_b += b_events
    
    # 计算检验统计量
    o_minus_e_a = observed_a - expected_a
    o_minus_e_b = observed_b - expected_b
    
    # 方差计算
    variance = 0.0
    for time_point in unique_times:
        a_at_risk = np.sum(a_times >= time_point)
        b_at_risk = np.sum(b_times >= time_point)
        total_at_risk = a_at_risk + b_at_risk
        
        if total_at_risk <= 1:
            continue
        
        total_events = (np.sum(a_deaths[a_times == time_point]) + 
                       np.sum(b_deaths[b_times == time_point]))
        
        if total_events == 0:
            continue
        
        # 计算该时间点对方差的贡献
        term = (total_events * a_at_risk * b_at_risk * (total_at_risk - total_events)) / (total_at_risk**2 * (total_at_risk - 1))
        variance += term
    
    # 计算卡方统计量
    if variance > 0:
        chi_square = (o_minus_e_a**2) / variance
        degrees_of_freedom = 1
        p_value = 1 - stats.chi2.cdf(chi_square, degrees_of_freedom)
    else:
        chi_square = 0.0
        degrees_of_freedom = 1
        p_value = 1.0
    
    return {
        "chi_square_statistic": float(chi_square),
        "degrees_of_freedom": int(degrees_of_freedom),
        "p_value": float(p_value),
        "observed_events": [float(observed_a), float(observed_b)],
        "expected_events": [float(expected_a), float(expected_b)],
        "method_name": "Log-rank检验",
        "interpretation": "Log-rank检验对晚期差异敏感"
    }


def _perform_breslow_test_two_groups(a_times: np.ndarray, a_deaths: np.ndarray, a_censored: np.ndarray,
                                   b_times: np.ndarray, b_deaths: np.ndarray, b_censored: np.ndarray) -> Dict:
    """执行两组Breslow检验（广义Wilcoxon检验）"""
    # 合并所有时间点
    all_times = np.concatenate([a_times, b_times])
    unique_times = np.unique(all_times)
    
    observed_a = 0.0
    observed_b = 0.0
    expected_a = 0.0
    expected_b = 0.0
    
    for time_point in unique_times:
        # 计算各组在该时间点的风险集
        a_at_risk = np.sum(a_times >= time_point)
        b_at_risk = np.sum(b_times >= time_point)
        total_at_risk = a_at_risk + b_at_risk
        
        if total_at_risk == 0:
            continue
        
        # 计算该时间点的事件数
        a_events = np.sum(a_deaths[a_times == time_point])
        b_events = np.sum(b_deaths[b_times == time_point])
        total_events = a_events + b_events
        
        if total_events == 0:
            continue
        
        # Breslow检验使用权重等于风险集大小
        weight = total_at_risk
        
        # 计算加权期望事件数
        expected_a += weight * total_events * a_at_risk / total_at_risk
        expected_b += weight * total_events * b_at_risk / total_at_risk
        
        observed_a += weight * a_events
        observed_b += weight * b_events
    
    # 计算检验统计量
    o_minus_e_a = observed_a - expected_a
    o_minus_e_b = observed_b - expected_b
    
    # 方差计算（加权版本）
    variance = 0.0
    for time_point in unique_times:
        a_at_risk = np.sum(a_times >= time_point)
        b_at_risk = np.sum(b_times >= time_point)
        total_at_risk = a_at_risk + b_at_risk
        
        if total_at_risk <= 1:
            continue
        
        a_events = np.sum(a_deaths[a_times == time_point])
        b_events = np.sum(b_deaths[b_times == time_point])
        total_events = a_events + b_events
        
        if total_events == 0:
            continue
        
        weight = total_at_risk
        
        # 计算该时间点对方差的贡献
        term = weight**2 * (total_events * a_at_risk * b_at_risk * (total_at_risk - total_events)) / (total_at_risk**2 * (total_at_risk - 1))
        variance += term
    
    # 计算卡方统计量
    if variance > 0:
        chi_square = (o_minus_e_a**2) / variance
        degrees_of_freedom = 1
        p_value = 1 - stats.chi2.cdf(chi_square, degrees_of_freedom)
    else:
        chi_square = 0.0
        degrees_of_freedom = 1
        p_value = 1.0
    
    return {
        "chi_square_statistic": float(chi_square),
        "degrees_of_freedom": int(degrees_of_freedom),
        "p_value": float(p_value),
        "observed_events": [float(observed_a), float(observed_b)],
        "expected_events": [float(expected_a), float(expected_b)],
        "method_name": "Breslow检验（广义Wilcoxon）",
        "interpretation": "Breslow检验对早期差异更敏感"
    }


def _analyze_early_late_sensitivity(a_stats: List[Dict], b_stats: List[Dict], times: np.ndarray) -> Dict:
    """分析早期和后期差异的敏感性"""
    if len(times) < 4:
        return {
            "early_difference": "数据点不足，无法分析",
            "late_difference": "数据点不足，无法分析",
            "pattern": "数据点不足，无法分析",
            "interpretation": "时间点数量不足，建议收集更多随访数据"
        }
    
    # 分割时间点为早期和晚期
    mid_point = len(times) // 2
    early_times = times[:mid_point]
    late_times = times[mid_point:]
    
    # 计算早期差异
    early_a_deaths = sum(stat["actual_deaths"] for stat in a_stats[:mid_point])
    early_b_deaths = sum(stat["actual_deaths"] for stat in b_stats[:mid_point])
    early_difference = abs(early_a_deaths - early_b_deaths)
    
    # 计算晚期差异
    late_a_deaths = sum(stat["actual_deaths"] for stat in a_stats[mid_point:])
    late_b_deaths = sum(stat["actual_deaths"] for stat in b_stats[mid_point:])
    late_difference = abs(late_a_deaths - late_b_deaths)
    
    # 判断敏感性
    if early_difference > late_difference:
        pattern = "早期差异更明显"
        interpretation = "Breslow检验可能更适合，因为它对早期差异更敏感"
    elif late_difference > early_difference:
        pattern = "晚期差异更明显"
        interpretation = "Log-rank检验可能更适合，因为它对晚期差异更敏感"
    else:
        pattern = "早期晚期差异相当"
        interpretation = "两种检验方法结果可能相似"
    
    return {
        "early_difference": int(early_difference),
        "late_difference": int(late_difference),
        "pattern": pattern,
        "interpretation": interpretation
    }


def cal_result_sur_lv(param: SurParamLV) -> Dict[str, Any]:
    """
    生成两组生存率比较分析的完整报告字典
    
    此函数整合了两组生存率比较分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的生存率比较结果。报告包括输入参数、
    时间点分析、统计检验结果和敏感性分析等信息，
    便于临床医生和研究人员快速理解生存率比较的特征。
    
    Args:
        stats_data_list: 包含生存分析所需参数的列表，必须包含6个字典元素，按顺序为：
            - 第1个字典: field_name='A组生存时间'，data_list为A组生存时间列表，浮点数列表类型
            - 第2个字典: field_name='A组死亡数'，data_list为A组死亡数列表，整数列表类型
            - 第3个字典: field_name='A组删失数'，data_list为A组删失数列表，整数列表类型
            - 第4个字典: field_name='B组生存时间'，data_list为B组生存时间列表，浮点数列表类型
            - 第5个字典: field_name='B组死亡数'，data_list为B组死亡数列表，整数列表类型
            - 第6个字典: field_name='B组删失数'，data_list为B组删失数列表，整数列表类型
    
    Returns:
        Dict[str, Any]: 包含生存率比较指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"两组生存率比较分析"
            - input_parameters: 输入参数信息
                - group_a_name: A组名称，字符串类型
                - group_b_name: B组名称，字符串类型
                - a_sample_size: A组样本量，整数类型
                - b_sample_size: B组样本量，整数类型
                - a_events: A组事件数，整数类型
                - b_events: B组事件数，整数类型
                - common_time_points: 公共时间点列表，浮点数列表类型
            - time_point_analysis: 时间点分析
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
            - statistical_tests: 统计检验结果
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
            - sensitivity_analysis: 敏感性分析
                - early_difference: 早期差异，整数类型
                - late_difference: 晚期差异，整数类型
                - pattern: 差异模式，字符串类型
                - interpretation: 方法建议，字符串类型
            - interpretation: 结果解释
                - log_rank_description: Log-rank描述，字符串类型
                - breslow_description: Breslow描述，字符串类型
                - sensitivity_interpretation: 敏感性解释，字符串类型
            - remark: 备注信息，字符串类型
    """
    # 从参数对象解构
    stats_data_list = [item.model_dump() for item in param.stats_data_list]

    # 执行两组生存率比较分析
    results = survival_rate_comparison_two_groups(stats_data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "两组生存率比较分析",
        "input_parameters": results["input_parameters"],
        "time_point_analysis": results["time_point_analysis"],
        "statistical_tests": results["statistical_tests"],
        "sensitivity_analysis": results["sensitivity_analysis"],
        "interpretation": results["interpretation"],
        "remark": f"A组: {results['input_parameters']['group_a_name']}, B组: {results['input_parameters']['group_b_name']}"
    }
    
    return result_dict