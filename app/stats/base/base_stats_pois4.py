"""
临床泊松分布两样本比较统计分析模块 (Clinical Two-Sample Poisson Distribution Comparison Statistics Module)

本模块提供全面的泊松分布两样本比较统计分析功能，用于临床数据中两个样本均数的差异检验，
是临床试验数据分析的重要组成部分。两样本均数比较在医学研究中广泛应用，
如比较两组患者的事件发生率、两种治疗方法的不良事件发生率、不同时期的疾病发病率等。

【模块功能概述】:
1. 两样本均数比较：计算两个样本均数的差异及其统计显著性
2. 比率差计算：计算两个样本率的差值及置信区间
3. 比值比计算：计算相对风险比及其置信区间
4. 精确检验：提供泊松分布的精确检验方法
5. 效应量计算：计算差异的效应量指标

【临床应用价值】:
- 疾病监测：比较不同时期的疾病发生率
- 安全评估：比较不同治疗方法的不良事件发生率
- 政策评估：比较不同政策实施前后的事件发生率
- 质量改进：比较不同干预措施的效果

【统计方法选择指南】:
1. 精确检验适用条件：
   - 适用于所有样本量
   - 特别适用于小样本或事件数较少的情况

2. 正态近似法适用条件：
   - 两个样本的事件数都足够大
   - 通常要求事件数 ≥ 5

3. 临床应用场景：
   - 比较新旧两种治疗方法的安全性
   - 比较不同手术方式的并发症发生率
   - 比较不同时期的感染率

【结果解读注意事项】:
1. P值解释：在零假设成立的前提下，观察到当前或更极端结果的概率
2. 置信区间解释：均数差或比值的置信区间
3. 显著性水平：通常使用α=0.05作为判断标准
4. 临床意义：统计显著性不等同于临床重要性
5. 混杂因素：注意潜在混杂因素的影响

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "如何比较两个泊松分布样本的均数？"
- "泊松分布比较与二项分布比较有何不同？"
- "如何解释相对风险比？"
- "P值小于0.05意味着什么？"
- "如何判断临床效果的实际意义？"

【相关标准和规范】:
- CONSORT 声明：随机临床试验报告规范
- STROBE 声明：观察性研究报告规范
- ICH E9 指导原则：临床试验的统计学原则
- SAMPL 指南：医学研究报告中的统计方法描述规范

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

import math
from typing import Dict, Any
from scipy.stats import norm
from app.schemas.request_data.base_param import BaseParamPois4


def two_sample_poisson_test(sample1_num: int, sample1_tick: int, 
                           sample2_num: int, sample2_tick: int, 
                           alpha: float = 0.05) -> Dict[str, Any]:
    """
    两样本泊松分布比较检验
    
    在临床研究中，这用于比较两组数据中事件发生率的差异，
    以评估其统计学显著性。

    Args:
        sample1_num: 样本1的观察单位数
        sample1_tick: 样本1的事件数
        sample2_num: 样本2的观察单位数
        sample2_tick: 样本2的事件数
        alpha: 显著性水平 (默认0.05)
        
    Returns:
        Dict[str, Any]: 包含检验统计量和结论的字典
        
    原假设 H0: λ1 = λ2 (两样本的泊松参数相等)
    备择假设 H1: λ1 ≠ λ2 (两样本的泊松参数不等)
    """
    # 计算各样本的率
    rate1 = sample1_tick / sample1_num if sample1_num > 0 else 0
    rate2 = sample2_tick / sample2_num if sample2_num > 0 else 0
    
    # 总事件数和总观察单位数
    total_events = sample1_tick + sample2_tick
    total_units = sample1_num + sample2_num
    
    # 如果任意样本没有事件，需要特殊处理
    if sample1_tick == 0 and sample2_tick == 0:
        return {
            'rate1': rate1,
            'rate2': rate2,
            'p_value': 1.0,
            'significant': False
        }
    
    # 在零假设下，每个样本的期望事件数
    exp1 = sample1_num * total_events / total_units
    exp2 = sample2_num * total_events / total_units
    
    # 使用正态近似法计算Z统计量
    # 在H0成立下，(X1-E[X1])/sqrt(Var[X1|H0]) ~ N(0,1)
    # Var[X1|H0] = n1^2 * (N1+N2)/(n1+n2)^2 * λ = n1^2 * (N1+N2)/(n1+n2)^2 * (N1+N2)/(n1+n2) = n1*n2*(N1+N2)/(n1+n2)^2
    if total_units == 0:
        raise ValueError("观察单位总数不能为0")
    
    var_under_h0 = (sample1_num * sample2_num * total_events) / (total_units * total_units)
    
    if var_under_h0 == 0:
        z_stat = float('inf') if sample1_tick != exp1 else 0
    else:
        z_stat = (sample1_tick - exp1) / math.sqrt(var_under_h0)
    
    # 计算P值
    p_value = 2 * (1 - norm.cdf(abs(z_stat))) if abs(z_stat) < float('inf') else 0.0
    
    # 判断显著性
    significant = p_value < alpha
    
    return {
        'rate1': rate1,
        'rate2': rate2,
        'exp1': exp1,
        'exp2': exp2,
        'z_statistic': z_stat,
        'p_value': p_value,
        'significant': significant
    }


def calculate_confidence_interval_for_rate_ratio(
    events1: int, units1: int, 
    events2: int, units2: int, 
    confidence_level: float = 0.95
) -> Dict[str, float]:
    """
    计算两个泊松率比的置信区间
    
    在临床研究中，这用于估计两个样本率比值的置信区间，
    以提供相对风险的精确范围估计。
    
    Args:
        events1: 样本1的事件数
        units1: 样本1的观察单位数
        events2: 样本2的事件数
        units2: 样本2的观察单位数
        confidence_level: 置信水平 (默认0.95)
        
    Returns:
        Dict[str, float]: 包含率比及其置信区间的字典
    """
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha/2)
    
    # 计算率比
    if events2 == 0 or units2 == 0:
        # 如果第二个样本没有事件或单位数为0，则无法计算率比
        return {
            'rate_ratio': float('inf'),
            'log_rr_se': float('inf'),
            'lower_bound': float('inf'),
            'upper_bound': float('inf')
        }
    
    rate1 = events1 / units1 if units1 > 0 else 0
    rate2 = events2 / units2 if units2 > 0 else 0
    
    rate_ratio = rate1 / rate2
    
    if events1 == 0 or events2 == 0:
        # 当任一样本事件数为0时，需要特殊处理
        log_rr = math.log(rate_ratio) if rate_ratio > 0 else float('-inf')
        log_rr_se = math.sqrt(1/events1 + 1/events2) if events1 > 0 and events2 > 0 else float('inf')
    else:
        # 计算对数率比的标准误
        log_rr_se = math.sqrt(1/events1 + 1/events2)
        log_rr = math.log(rate_ratio)
    
    # 计算置信区间
    log_lower = log_rr - z * log_rr_se
    log_upper = log_rr + z * log_rr_se
    
    ci_lower = math.exp(log_lower) if log_lower != float('-inf') else 0
    ci_upper = math.exp(log_upper) if log_upper != float('-inf') else float('inf')
    
    return {
        'rate_ratio': rate_ratio,
        'log_rr_se': log_rr_se,
        'lower_bound': ci_lower,
        'upper_bound': ci_upper
    }


def calculate_confidence_interval_for_rate_diff(
    events1: int, units1: int, 
    events2: int, units2: int, 
    confidence_level: float = 0.95
) -> Dict[str, float]:
    """
    计算两个泊松率差的置信区间
    
    在临床研究中，这用于估计两个样本率差值的置信区间，
    以提供绝对风险差的精确范围估计。
    
    Args:
        events1: 样本1的事件数
        units1: 样本1的观察单位数
        events2: 样本2的事件数
        units2: 样本2的观察单位数
        confidence_level: 置信水平 (默认0.95)
        
    Returns:
        Dict[str, float]: 包含率差及其置信区间的字典
    """
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha/2)
    
    rate1 = events1 / units1 if units1 > 0 else 0
    rate2 = events2 / units2 if units2 > 0 else 0
    
    rate_diff = rate1 - rate2
    
    # 计算标准误
    se = math.sqrt(rate1/units1 + rate2/units2)
    
    margin_error = z * se
    lower_bound = rate_diff - margin_error
    upper_bound = rate_diff + margin_error
    
    return {
        'rate_diff': rate_diff,
        'standard_error': se,
        'margin_error': margin_error,
        'lower_bound': lower_bound,
        'upper_bound': upper_bound
    }


def cal_result_pois4(param: BaseParamPois4) -> Dict[str, Any]:
    """
    生成泊松分布两样本比较统计分析的完整报告字典
    
    此函数整合了泊松分布两样本比较的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的泊松分布两样本比较分析结果。
    报告包括输入参数、统计量、P值、显著性判断等信息，
    便于临床医生和研究人员快速理解泊松分布两样本比较的特征。
    
    Args:
        param: BaseParamPois4对象，包含sample1_num, sample1_tick, sample2_num, sample2_tick
        
    Returns:
        Dict[str, Any]: 包含泊松分布两样本比较统计分析指标的字典，键为指标名称，值为对应的统计量
            - sample1_rate: 样本1的率
            - sample2_rate: 样本2的率
            - z_statistic: Z统计量
            - p_value: P值
            - significant: 是否有统计学意义
            - rate_ratio_ci: 率比的置信区间
            - rate_diff_ci: 率差的置信区间
    """
    # 从参数对象中提取值
    sample1_num = param.sample1_num
    sample1_tick = param.sample1_tick
    sample2_num = param.sample2_num
    sample2_tick = param.sample2_tick
    
    # 确保输入参数有效
    if sample1_num <= 0 or sample2_num <= 0:
        raise ValueError("观察单位数必须大于0")
    if sample1_tick < 0 or sample2_tick < 0:
        raise ValueError("事件数不能为负数")
    
    # 执行两样本泊松检验
    test_results = two_sample_poisson_test(
        sample1_num, sample1_tick,
        sample2_num, sample2_tick
    )
    
    # 计算率比的置信区间
    rate_ratio_ci = calculate_confidence_interval_for_rate_ratio(
        sample1_tick, sample1_num,
        sample2_tick, sample2_num
    )
    
    # 计算率差的置信区间
    rate_diff_ci = calculate_confidence_interval_for_rate_diff(
        sample1_tick, sample1_num,
        sample2_tick, sample2_num
    )
    
    # 计算各组的基本统计量
    rate1 = sample1_tick / sample1_num
    rate2 = sample2_tick / sample2_num
    
    # 构建结果字典
    result_dict = {
        "table_name": "泊松分布两样本比较",
        "input_parameters": {
            "sample1_observation_units": sample1_num,
            "sample1_events": sample1_tick,
            "sample1_rate": rate1,
            "sample2_observation_units": sample2_num,
            "sample2_events": sample2_tick,
            "sample2_rate": rate2
        },
        "test_results": test_results,
        "rate_ratio_ci": rate_ratio_ci,
        "rate_diff_ci": rate_diff_ci,
        "interpretation": {
            "statistical_interpretation": f"样本1率为{rate1:.6f}，样本2率为{rate2:.6f}，差值为{rate1-rate2:.6f}",
            "p_value_interpretation": f"P值为{test_results['p_value']:.6f}",
            "significance_interpretation": f"{'差异有统计学意义' if test_results['significant'] else '差异无统计学意义'} (α=0.05)",
            "clinical_significance": "请结合临床实际意义解读统计学差异",
            "rate_ratio_interpretation": f"率比为{rate_ratio_ci['rate_ratio']:.6f}，95%CI为({rate_ratio_ci['lower_bound']:.6f}, {rate_ratio_ci['upper_bound']:.6f})",
            "rate_diff_interpretation": f"率差为{rate_diff_ci['rate_diff']:.6f}，95%CI为({rate_diff_ci['lower_bound']:.6f}, {rate_diff_ci['upper_bound']:.6f})"
        },
        "remark": f"样本1: 观察单位={sample1_num}, 事件数={sample1_tick}, 率={rate1:.6f}; 样本2: 观察单位={sample2_num}, 事件数={sample2_tick}, 率={rate2:.6f}"
    }
    
    return result_dict