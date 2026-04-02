"""
寿命表法生存分析模块 (Life Table Survival Analysis Module)

本模块提供基于寿命表法的生存分析功能，用于估计生存函数和计算生存概率。
寿命表法是一种经典的生存分析方法，特别适用于大样本数据，通过将时间轴划分为区间，
计算每个时间区间内的生存概率。

【模块功能概述】:
1. 生存概率估计：计算各时间区间的生存概率
2. 死亡概率计算：计算各时间区间的死亡概率
3. 风险集估计：估计各时间区间的平均风险集大小
4. 中位生存时间：估计中位生存时间
5. 生存标准误：计算生存概率的标准误

【临床应用价值】:
- 人群寿命研究：分析大规模人群的生存模式
- 长期随访研究：处理长期随访的生存数据
- 分组比较：比较不同组别的生存模式
- 保险精算：计算生命表用于保险精算

【统计方法选择指南】:
1. 寿命表法适用条件：
   - 大样本数据（通常n>30）
   - 数据按时间区间分组
   - 适用于分组数据的生存分析
   - 适合处理间隔删失数据

2. 临床应用场景：
   - 人口普查生存分析
   - 长期疾病预后研究
   - 保险公司的生命表编制
   - 医疗政策制定参考

【结果解读注意事项】:
1. 生存概率解释：表示在特定时间区间开始时尚存活的概率
2. 死亡概率解释：表示在特定时间区间内发生死亡的概率
3. 中位生存时间：50%的个体存活的时间区间
4. 风险集大小：每个时间区间内的平均风险人数

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "寿命表法与Kaplan-Meier法有何区别？"
- "如何解释寿命表中的各项指标？"
- "寿命表法适用于什么类型的数据？"
- "如何选择合适的区间长度？"
- "如何比较不同组的寿命表？"

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
from app.schemas.request_data.sur_param import SurParamLT


def life_table_survival_analysis(stats_data_list: List[Dict[str, Any]], time_intervals: List[str]) -> Dict:
    """
    基于汇总数据的寿命表法生存分析
    stats_data_list包含两个StatsData对象：
    - field_name='期内死亡数'的StatsData: data_list是各时间区间的死亡数列表
    - field_name='期内删失数'的StatsData: data_list是各时间区间的删失数列表
    time_intervals是时间区间标识列表
    
    参数:
    - stats_data_list: 包含生存分析所需参数的列表
    - time_intervals: 时间区间标识列表
    
    返回:
    - 包含寿命表分析统计量和结果的字典
    """
    # 参数验证
    if not stats_data_list:
        raise ValueError("寿命表分析至少需要2组数据：期内死亡数和期内删失数")
    
    if len(stats_data_list) != 2:
        raise ValueError("寿命表分析需要恰好2组数据：期内死亡数和期内删失数")
    
    if not time_intervals:
        raise ValueError("时间区间列表不能为空")
    
    # 验证数据结构
    death_data = None
    censored_data = None
    
    for stats_data in stats_data_list:
        if stats_data.get("field_name") == '期内死亡数':
            death_data = stats_data
        elif stats_data.get("field_name") == '期内删失数':
            censored_data = stats_data
        else:
            raise ValueError(f"未知的数据类型: {stats_data.get('field_name')}，应为'期内死亡数'或'期内删失数'")
    
    if death_data is None:
        raise ValueError("缺少期内死亡数数据")
    
    if censored_data is None:
        raise ValueError("缺少期内删失数数据")
    
    # 获取时间区间列表
    n_intervals = len(time_intervals)
    
    # 验证数据长度
    if len(death_data["data_list"]) != n_intervals:
        raise ValueError(f"期内死亡数数据长度({len(death_data['data_list'])})与时间区间数({n_intervals})不匹配")
    
    if len(censored_data["data_list"]) != n_intervals:
        raise ValueError(f"期内删失数数据长度({len(censored_data['data_list'])})与时间区间数({n_intervals})不匹配")
    
    # 转换数据为整数数组
    deaths = np.array(death_data["data_list"]).astype(int)
    censored = np.array(censored_data["data_list"]).astype(int)
    
    # 验证数据合理性
    if np.any(deaths < 0) or np.any(censored < 0):
        raise ValueError("死亡数和删失数必须为非负整数")
    
    # 计算统计量
    total_deaths = np.sum(deaths)
    total_censored = np.sum(censored)
    total_subjects = total_deaths + total_censored
    
    # 初始化寿命表数组
    interval_starts = []  # 区间起始时间
    interval_ends = []    # 区间结束时间
    entered = []          # 期初人数
    died = []             # 区间内死亡数
    censored_in_interval = []  # 区间内删失数
    at_risk_midpoint = [] # 平均风险集大小
    death_probability = [] # 死亡概率
    survival_probability = [] # 生存概率
    survival_se = []      # 生存标准误
    
    # 计算每个区间的生命表统计量
    current_at_risk = total_subjects
    
    for i in range(n_intervals):
        # 基本统计量
        interval_died = deaths[i]
        interval_censored = censored[i]
        
        # 期初人数
        entered.append(current_at_risk)
        
        # 区间内死亡数和删失数
        died.append(interval_died)
        censored_in_interval.append(interval_censored)
        
        # 平均风险集大小（期初人数 - 删失数/2）
        avg_at_risk = current_at_risk - interval_censored / 2
        at_risk_midpoint.append(avg_at_risk)
        
        # 死亡概率
        if avg_at_risk > 0:
            prob_death = interval_died / avg_at_risk
            death_probability.append(prob_death)
        else:
            prob_death = 0
            death_probability.append(0)
        
        # 生存概率（累积生存概率）
        if i == 0:
            surv_prob = 1 - prob_death
        else:
            surv_prob = survival_probability[i-1] * (1 - prob_death)
        survival_probability.append(surv_prob)
        
        # 生存标准误（使用 Greenwood 公式的近似）
        if surv_prob > 0 and avg_at_risk > 0 and interval_died > 0:
            variance = (interval_died / (avg_at_risk * (avg_at_risk - interval_died))) if avg_at_risk > interval_died else 0
            se = surv_prob * np.sqrt(variance)
        else:
            se = 0
        survival_se.append(se)
        
        # 更新下一个区间的期初人数
        current_at_risk -= (interval_died + interval_censored)
    
    # 计算中位生存时间
    median_survival = _calculate_life_table_median(time_intervals, survival_probability)
    
    return {
        "input_parameters": {
            "n_intervals": int(n_intervals),
            "time_intervals": time_intervals,
            "total_subjects": int(total_subjects),
            "total_deaths": int(total_deaths),
            "total_censored": int(total_censored)
        },
        "life_table": {
            "interval_starts": interval_starts,
            "interval_ends": interval_ends,
            "entered": entered,
            "died": died,
            "censored": censored_in_interval,
            "at_risk_midpoint": at_risk_midpoint,
            "death_probability": death_probability,
            "survival_probability": survival_probability,
            "survival_se": survival_se
        },
        "summary_statistics": {
            "median_survival": float(median_survival) if median_survival is not None else None,
            "final_survival_rate": float(survival_probability[-1]) if survival_probability else 0.0
        }
    }


def _calculate_life_table_median(time_intervals: List[str], survival_probs: List[float]) -> float:
    """
    计算寿命表法的中位生存时间
    
    参数:
    - time_intervals: 时间区间标识列表
    - survival_probs: 对应的生存概率列表
    
    返回:
    - 中位生存时间估算值，如果无法计算则返回None
    """
    if len(survival_probs) == 0:
        return None
    
    # 找到生存概率首次降到0.5以下的区间
    for i, prob in enumerate(survival_probs):
        if prob <= 0.5:
            # 简单估算：返回该区间的中点时间
            try:
                interval_desc = time_intervals[i]
                # 尝试解析区间描述，例如"[0,1)" 或 "0-1年"
                if '-' in interval_desc:
                    parts = interval_desc.replace('[','').replace(']','').replace('(','').replace(')','').split('-')
                    if len(parts) == 2:
                        start_time = float(parts[0])
                        end_time = float(parts[1])
                        median_time = (start_time + end_time) / 2
                        return median_time
                elif ',' in interval_desc:
                    parts = interval_desc.replace('[','').replace(']','').replace('(','').replace(')','').split(',')
                    if len(parts) == 2:
                        start_time = float(parts[0])
                        end_time = float(parts[1])
                        median_time = (start_time + end_time) / 2
                        return median_time
            except:
                # 如果无法解析区间描述，返回区间的序号作为粗略估计
                return float(i + 1)
            
            # 默认返回区间的序号
            return float(i + 1)
    
    # 如果生存概率从未降到0.5以下，返回最后一个区间
    return float(len(time_intervals))


def perform_life_table_survival(stats_data_list: List[Dict[str, Any]], time_intervals: List[str]) -> Dict:
    """
    根据给定参数执行寿命表法生存分析
    
    参数:
    - stats_data_list: 包含寿命表分析所需参数的列表
    - time_intervals: 时间区间标识列表
    
    返回:
    - 包含寿命表分析完整结果的字典
    """
    # 验证参数
    if not stats_data_list:
        raise ValueError("寿命表分析至少需要期内死亡数和期内删失数数据")
    
    if len(stats_data_list) != 2:
        raise ValueError("寿命表分析需要恰好2组数据")
    
    # 执行寿命表分析
    results = life_table_survival_analysis(stats_data_list, time_intervals)
    
    return results


def cal_result_sur_lt(param: SurParamLT) -> Dict[str, Any]:
    """
    生成寿命表法生存分析的完整报告字典
    
    此函数整合了寿命表法生存分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的生存分析结果。报告包括输入参数、
    寿命表结果和汇总统计等信息，
    便于临床医生和研究人员快速理解寿命表分析的特征。
    
    Args:
        stats_data_list: 包含生存分析所需参数的列表，必须包含2个字典元素，按顺序为：
            - 第一个字典: field_name='期内死亡数'，data_list为各时间区间的死亡数列表，整数列表类型
            - 第二个字典: field_name='期内删失数'，data_list为各时间区间的删失数列表，整数列表类型
        time_intervals: 时间区间标识列表，字符串列表类型，例如['[0,1)', '[1,2)', ...]
    
    Returns:
        Dict[str, Any]: 包含生存分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"寿命表法生存分析报告"
            - input_parameters: 输入参数信息
                - n_intervals: 区间数，整数类型
                - time_intervals: 时间区间列表，字符串列表类型
                - total_subjects: 总样本数，整数类型
                - total_deaths: 总死亡数，整数类型
                - total_censored: 总删失数，整数类型
            - life_table: 寿命表结果
                - entered: 期初人数列表，整数列表类型
                - died: 死亡数列表，整数列表类型
                - censored: 删失数列表，整数列表类型
                - at_risk_midpoint: 平均风险集大小列表，浮点数列表类型
                - death_probability: 死亡概率列表，浮点数列表类型
                - survival_probability: 生存概率列表，浮点数列表类型
                - survival_se: 生存标准误列表，浮点数列表类型
            - summary_statistics: 汇总统计
                - median_survival: 中位生存时间，浮点数类型或None
                - final_survival_rate: 最终生存率，浮点数类型
            - remark: 备注信息，字符串类型
    """
    # 执行寿命表分析
    # 从参数对象解构
    stats_data_list = [item.model_dump() for item in param.stats_data_list]
    time_intervals = param.stats_name.name_list

    results = perform_life_table_survival(stats_data_list, time_intervals)
    
    # 构建结果字典
    result_dict = {
        "table_name": "寿命表法生存分析报告",
        "input_parameters": results["input_parameters"],
        "life_table": results["life_table"],
        "summary_statistics": results["summary_statistics"],
        "remark": f"区间数: {results['input_parameters']['n_intervals']}, 总样本数: {results['input_parameters']['total_subjects']}"
    }
    
    return result_dict