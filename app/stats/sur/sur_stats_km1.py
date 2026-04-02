"""
Kaplan-Meier生存分析模块 - 基于原始数据 (Kaplan-Meier Survival Analysis Module - Based on Raw Data)

本模块提供基于原始数据的Kaplan-Meier生存分析功能，用于估计生存函数和计算生存概率。
Kaplan-Meier估计是一种非参数统计方法，广泛应用于医学研究、工程可靠性等领域，
特别适用于存在删失数据的情况。

【模块功能概述】:
1. 生存概率估计：计算各时间点的生存概率
2. 标准误计算：使用Greenwood公式计算生存概率的标准误
3. 置信区间估计：计算生存概率的置信区间
4. 中位生存时间：估计中位生存时间
5. 生存曲线绘制：提供绘制生存曲线所需的数据

【临床应用价值】:
- 生存分析：评估患者在特定时间点的生存概率
- 治疗效果比较：比较不同治疗方法的生存差异
- 预后评估：预测疾病进展或复发的可能性
- 临床试验：评估新药或新疗法的有效性

【统计方法选择指南】:
1. Kaplan-Meier适用条件：
   - 数据为时间到事件型数据
   - 存在右删失数据
   - 事件发生时间准确记录
   - 样本量适中至较大

2. 临床应用场景：
   - 评估癌症患者的生存时间
   - 比较手术与保守治疗的效果
   - 评估药物长期疗效
   - 评估医疗器械使用寿命

【结果解读注意事项】:
1. 生存概率解释：表示在特定时间点仍存活的概率
2. 置信区间解释：提供生存概率估计的不确定性范围
3. 中位生存时间：50%的个体存活的时间点
4. 曲线下降速度：反映事件发生的速度，越陡峭表示事件发生越快

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用Kaplan-Meier分析吗？"
- "如何解释生存曲线？"
- "什么是中位生存时间？"
- "Kaplan-Meier估计的前提条件是什么？"
- "如何比较两条生存曲线？"

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
from app.schemas.request_data.sur_param import SurParamKM1


def kaplan_meier_survival_raw_data(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    基于原始数据的Kaplan-Meier生存分析
    每个StatsData代表一个组的数据：
    - field_name: 组名
    - data_list: 该组所有个体的生存时间列表
    
    参数:
    - stats_data_list: StatsData列表，每个代表一个组的生存时间数据
    
    返回:
    - 包含生存分析统计量和结果的字典
    """
    # 参数验证
    if not stats_data_list:
        raise ValueError("Kaplan-Meier分析至少需要1组数据")
    
    # 处理多组数据
    groups_data = []
    group_names = []
    
    for i, stats_data in enumerate(stats_data_list):
        survival_times = np.array(stats_data["data_list"])
        if len(survival_times) == 0:
            raise ValueError(f"第{i+1}组数据为空")
        
        # 所有数据都视为事件发生（无删失）
        outcomes = np.ones(len(survival_times))  # 全部标记为事件发生
        
        groups_data.append((survival_times, outcomes))
        group_names.append(stats_data.get("field_name") or f"组{i+1}")
    
    # 计算每个组的Kaplan-Meier生存曲线
    km_results = []
    
    for i, ((times, outcomes), name) in enumerate(zip(groups_data, group_names)):
        # 对数据按时间排序
        sort_indices = np.argsort(times)
        sorted_times = times[sort_indices]
        sorted_outcomes = outcomes[sort_indices]
        
        n_total = len(sorted_times)
        n_events = np.sum(sorted_outcomes)
        n_censored = n_total - n_events
        
        # 计算风险集大小和事件数
        unique_times, inverse_indices = np.unique(sorted_times, return_inverse=True)
        n_unique = len(unique_times)
        
        # 初始化结果数组
        survival_prob = np.ones(n_unique)
        standard_error = np.zeros(n_unique)
        confidence_lower = np.ones(n_unique)
        confidence_upper = np.ones(n_unique)
        at_risk_counts = np.zeros(n_unique)
        event_counts = np.zeros(n_unique)
        censored_counts = np.zeros(n_unique)
        
        # 计算生存概率（标准 Kaplan-Meier 乘积极限估计）
        cumulative_survival = 1.0
        
        for j, time_point in enumerate(unique_times):
            # 计算在该时间点的风险集大小
            at_risk = np.sum(sorted_times >= time_point)
            at_risk_counts[j] = at_risk
            
            # 计算该时间点的事件数
            events_at_time = np.sum((sorted_times == time_point) & (sorted_outcomes == 1))
            event_counts[j] = events_at_time
            
            # 计算该时间点的删失数
            censored_at_time = np.sum((sorted_times == time_point) & (sorted_outcomes == 0))
            censored_counts[j] = censored_at_time
            
            # 标准 KM 乘积极限: S(t) = S(t-1) * (1 - d/n)
            if at_risk > 0:
                cumulative_survival *= (1 - events_at_time / at_risk)
                survival_prob[j] = cumulative_survival
                
                # 计算标准误（Greenwood公式）
                if survival_prob[j] > 0:
                    variance = 0.0
                    for k in range(j + 1):
                        risk_k = at_risk_counts[k]
                        events_k = event_counts[k]
                        if risk_k > 0 and events_k > 0:
                            variance += events_k / (risk_k * (risk_k - events_k))
                    
                    # 避免负方差
                    variance = max(variance, 0)
                    standard_error[j] = survival_prob[j] * np.sqrt(variance)
                    
                    # 计算置信区间（使用log(-log)变换）
                    if standard_error[j] > 0:
                        z_alpha = 1.96  # 95%置信水平
                        log_log_survival = np.log(-np.log(survival_prob[j]))
                        se_log_log = standard_error[j] / (survival_prob[j] * np.abs(np.log(survival_prob[j])))
                        
                        lower_log_log = log_log_survival - z_alpha * se_log_log
                        upper_log_log = log_log_survival + z_alpha * se_log_log
                        
                        confidence_lower[j] = np.exp(-np.exp(lower_log_log))
                        confidence_upper[j] = np.exp(-np.exp(upper_log_log))
                    else:
                        confidence_lower[j] = survival_prob[j]
                        confidence_upper[j] = survival_prob[j]
                else:
                    standard_error[j] = 0.0
                    confidence_lower[j] = 0.0
                    confidence_upper[j] = 0.0
            else:
                survival_prob[j] = 0.0
                standard_error[j] = 0.0
                confidence_lower[j] = 0.0
                confidence_upper[j] = 0.0
        
        # 计算中位生存时间
        median_survival = _calculate_median_survival(unique_times, survival_prob)
        
        km_results.append({
            "group_name": name,
            "n_total": int(n_total),
            "n_events": int(n_events),
            "n_censored": int(n_censored),
            "unique_times": unique_times.tolist(),
            "survival_probabilities": survival_prob.tolist(),
            "standard_errors": standard_error.tolist(),
            "confidence_lower": confidence_lower.tolist(),
            "confidence_upper": confidence_upper.tolist(),
            "at_risk_counts": at_risk_counts.tolist(),
            "event_counts": event_counts.tolist(),
            "censored_counts": censored_counts.tolist(),
            "median_survival": float(median_survival) if median_survival is not None else None
        })
    
    return {
        "input_parameters": {
            "n_groups": len(groups_data),
            "group_names": group_names,
            "analysis_type": "Kaplan-Meier生存分析（原始数据）"
        },
        "km_results": km_results,
        "summary_statistics": {
            "total_groups": len(groups_data),
            "group_details": [
                {
                    "group_name": result["group_name"],
                    "sample_size": result["n_total"],
                    "events": result["n_events"],
                    "censored": result["n_censored"],
                    "median_survival": result["median_survival"]
                }
                for result in km_results
            ]
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


def perform_kaplan_meier_survival_raw(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    根据给定参数执行Kaplan-Meier生存分析
    
    参数:
    - stats_data_list: 包含生存分析所需参数的列表
    
    返回:
    - 包含Kaplan-Meier生存分析完整结果的字典
    """
    # 验证参数
    if not stats_data_list:
        raise ValueError("生存分析至少需要1组数据")
    
    # 执行Kaplan-Meier生存分析
    results = kaplan_meier_survival_raw_data(stats_data_list)
    
    return results


def cal_result_sur_km1(param: SurParamKM1) -> Dict[str, Any]:
    """
    生成Kaplan-Meier生存分析的完整报告字典
    
    此函数整合了Kaplan-Meier生存分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的生存分析结果。报告包括输入参数、
    生存分析结果和汇总统计等信息，
    便于临床医生和研究人员快速理解生存分析的特征。
    
    Args:
        stats_data_list: 包含生存分析所需参数的列表，每个元素是一个字典，包含：
            - field_name: 组名，字符串类型，可选，若不提供则默认为"组i"（i为组索引）
            - data_list: 该组所有个体的生存时间列表，浮点数列表类型
    
    Returns:
        Dict[str, Any]: 包含生存分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"Kaplan-Meier生存分析报告"
            - input_parameters: 输入参数信息
                - n_groups: 分析组数，整数类型
                - group_names: 组名称列表，字符串列表类型
                - analysis_type: 分析类型，字符串类型
            - km_results: Kaplan-Meier分析结果
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
            - summary_statistics: 汇总统计
                - total_groups: 总组数，整数类型
                - group_details: 组详细信息列表
                    - group_name: 组名称，字符串类型
                    - sample_size: 样本大小，整数类型
                    - events: 事件数，整数类型
                    - censored: 删失数，整数类型
                    - median_survival: 中位生存时间，浮点数类型或None
            - remark: 备注信息，字符串类型
    """
    # 执行Kaplan-Meier生存分析
    # 从参数对象解构
    stats_data_list = [item.model_dump() for item in param.stats_data_list]

    results = perform_kaplan_meier_survival_raw(stats_data_list)
    
    # 构建备注信息
    group_names = results['input_parameters'].get('group_names', [])
    group_info = ', '.join(['{name}(n={n})'.format(name=name, n=len(sd['data_list'])) for name, sd in zip(group_names, stats_data_list)])

    # 构建结果字典
    result_dict = {
        "table_name": "Kaplan-Meier生存分析报告",
        "input_parameters": results["input_parameters"],
        "km_results": results["km_results"],
        "summary_statistics": results["summary_statistics"],
        "remark": f"分析组数: {results['input_parameters'].get('n_groups', len(stats_data_list))}, 组别: {group_info}"
    }
    
    return result_dict