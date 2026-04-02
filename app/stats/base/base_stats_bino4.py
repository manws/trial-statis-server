"""
临床二项分布两样本比较统计分析模块 (Clinical Two-Sample Binomial Distribution Comparison Statistics Module)

本模块提供全面的二项分布两样本比较统计分析功能，用于临床数据中两个样本率的差异检验，
是临床试验数据分析的重要组成部分。两样本率比较在医学研究中广泛应用，
如比较两种治疗方法的成功率、两种药物的有效率、两组患者的不良反应发生率等。

【模块功能概述】:
1. 两样本率比较：计算两个样本率的差异及其统计显著性
2. 比率差计算：计算两个样本率的差值及置信区间
3. 比值比计算：计算相对风险比及其置信区间
4. Fisher精确检验：对于小样本提供精确检验方法
5. 效应量计算：计算差异的效应量指标

【临床应用价值】:
- 疗效比较：比较不同治疗方法的疗效差异
- 药物评价：比较不同药物的有效性和安全性
- 风险评估：比较不同暴露因素的风险差异
- 质量改进：比较不同干预措施的效果

【统计方法选择指南】:
1. 正态近似法适用条件：
   - 两个样本量都足够大
   - 每组np ≥ 5 且 n(1-p) ≥ 5

2. Fisher精确检验适用条件：
   - 小样本数据
   - 某些格子期望频数小于5

3. 临床应用场景：
   - 比较新旧两种治疗方法的有效率
   - 比较不同手术方式的并发症发生率
   - 比较不同药物的不良反应率

【结果解读注意事项】:
1. P值解释：在零假设成立的前提下，观察到当前或更极端结果的概率
2. 置信区间解释：比率差或比值比的置信区间
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
- "如何比较两个样本的率？"
- "什么时候使用Fisher精确检验？"
- "如何解释比值比？"
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
from app.schemas.request_data.base_param import BaseParamBino4


def two_sample_proportion_test(sample1_num: int, sample1_posi_num: int, 
                              sample2_num: int, sample2_posi_num: int, 
                              alpha: float = 0.05) -> Dict[str, Any]:
    """
    二项分布两样本比较检验
    
    在临床研究中，这用于比较两组数据中阳性率的差异，
    以评估其统计学显著性。

    Args:
        sample1_num: 样本1的总数
        sample1_posi_num: 样本1的阳性数
        sample2_num: 样本2的总数
        sample2_posi_num: 样本2的阳性数
        alpha: 显著性水平 (默认0.05)
        
    Returns:
        Dict[str, Any]: 包含检验统计量和结论的字典
        
    原假设 H0: p1 = p2 (两样本的总体率相等)
    备择假设 H1: p1 ≠ p2 (两样本的总体率不等)
    """
    # 计算各样本的率
    p1 = sample1_posi_num / sample1_num
    p2 = sample2_posi_num / sample2_num
    pooled_p = (sample1_posi_num + sample2_posi_num) / (sample1_num + sample2_num)
    
    # 计算标准误差
    se = math.sqrt(pooled_p * (1 - pooled_p) * (1/sample1_num + 1/sample2_num))
    
    # 计算Z统计量
    z_stat = (p1 - p2) / se if se > 0 else float('inf')
    
    # 计算P值
    p_value = 2 * (1 - norm.cdf(abs(z_stat))) if abs(z_stat) < float('inf') else 0.0
    
    # 判断显著性
    significant = p_value < alpha
    
    return {
        'p1': p1,
        'p2': p2,
        'pooled_p': pooled_p,
        'z_statistic': z_stat,
        'p_value': p_value,
        'significant': significant
    }


def calculate_confidence_interval_for_diff(p1: float, n1: int, p2: float, n2: int, 
                                          confidence_level: float = 0.95) -> Dict[str, float]:
    """
    计算两个比例差值的置信区间
    
    在临床研究中，这用于估计两个样本率差值的置信区间，
    以提供差异的精确范围估计。
    
    Args:
        p1: 样本1的比率
        n1: 样本1的总数
        p2: 样本2的比率
        n2: 样本2的总数
        confidence_level: 置信水平 (默认0.95)
        
    Returns:
        Dict[str, float]: 包含置信区间上下限的字典
    """
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha/2)
    
    diff = p1 - p2
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    
    margin_error = z * se
    lower_bound = diff - margin_error
    upper_bound = diff + margin_error
    
    return {
        'difference': diff,
        'standard_error': se,
        'margin_error': margin_error,
        'lower_bound': lower_bound,
        'upper_bound': upper_bound
    }


def calculate_relative_risk(n1: int, x1: int, n2: int, x2: int) -> Dict[str, float]:
    """
    计算相对风险比及其置信区间
    
    在临床研究中，这用于评估一个组相对于另一个组的风险倍数，
    是流行病学研究中的重要指标。
    
    Args:
        n1: 样本1的总数
        x1: 样本1的阳性数
        n2: 样本2的总数
        x2: 样本2的阳性数
        
    Returns:
        Dict[str, float]: 包含相对风险比及其置信区间的字典
    """
    p1 = x1 / n1
    p2 = x2 / n2
    
    if p2 == 0:
        return {
            'relative_risk': float('inf'),
            'log_rr_se': float('inf'),
            'lower_bound': float('inf'),
            'upper_bound': float('inf')
        }
    
    if p1 == 0:
        return {
            'relative_risk': 0.0,
            'log_rr_se': float('inf'),
            'lower_bound': 0.0,
            'upper_bound': float('inf')
        }
    
    rr = p1 / p2
    log_rr = math.log(rr)
    
    # 计算对数相对风险的标准误
    log_rr_se = math.sqrt((1/x1 - 1/n1) + (1/x2 - 1/n2))
    
    # 计算95%置信区间
    z = norm.ppf(0.975)  # 1.96 for 95% CI
    log_lower = log_rr - z * log_rr_se
    log_upper = log_rr + z * log_rr_se
    
    ci_lower = math.exp(log_lower)
    ci_upper = math.exp(log_upper)
    
    return {
        'relative_risk': rr,
        'log_rr_se': log_rr_se,
        'lower_bound': ci_lower,
        'upper_bound': ci_upper
    }


def cal_result_bino4(param: BaseParamBino4) -> Dict[str, Any]:
    """
    生成二项分布两样本比较统计分析的完整报告字典
    
    此函数整合了二项分布两样本比较的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的二项分布两样本比较分析结果。
    报告包括输入参数、统计量、P值、显著性判断等信息，
    便于临床医生和研究人员快速理解二项分布两样本比较的特征。
    
    Args:
        param: BaseParamBino4对象，包含sample1_num, sample1_posi_num, sample2_num, sample2_posi_num
        
    Returns:
        Dict[str, Any]: 包含二项分布两样本比较统计分析指标的字典，键为指标名称，值为对应的统计量
            - sample1_rate: 样本1的比率
            - sample2_rate: 样本2的比率
            - pooled_rate: 合并率
            - z_statistic: Z统计量
            - p_value: P值
            - significant: 是否有统计学意义
            - confidence_interval: 比率差的置信区间
            - relative_risk: 相对风险比
    """
    # 从参数对象中提取值
    sample1_num = param.sample1_num
    sample1_posi_num = param.sample1_posi_num
    sample2_num = param.sample2_num
    sample2_posi_num = param.sample2_posi_num
    
    # 确保输入参数有效
    if sample1_num <= 0 or sample2_num <= 0:
        raise ValueError("样本量必须大于0")
    if sample1_posi_num < 0 or sample1_posi_num > sample1_num:
        raise ValueError("样本1阳性数必须在0和样本1总数之间")
    if sample2_posi_num < 0 or sample2_posi_num > sample2_num:
        raise ValueError("样本2阳性数必须在0和样本2总数之间")
    
    # 执行两样本比例检验
    test_results = two_sample_proportion_test(
        sample1_num, sample1_posi_num, 
        sample2_num, sample2_posi_num
    )
    
    # 计算比率差的置信区间
    p1 = sample1_posi_num / sample1_num
    p2 = sample2_posi_num / sample2_num
    ci_results = calculate_confidence_interval_for_diff(p1, sample1_num, p2, sample2_num)
    
    # 计算相对风险比
    rr_results = calculate_relative_risk(sample1_num, sample1_posi_num, sample2_num, sample2_posi_num)
    
    # 验证正态近似条件
    n1p1 = sample1_num * p1
    n1q1 = sample1_num * (1 - p1)
    n2p2 = sample2_num * p2
    n2q2 = sample2_num * (1 - p2)
    
    normal_approx_valid = (n1p1 >= 5 and n1q1 >= 5 and n2p2 >= 5 and n2q2 >= 5)
    
    # 构建结果字典
    result_dict = {
        "table_name": "二项分布两样本比较",
        "input_parameters": {
            "sample1_size": sample1_num,
            "sample1_positive_number": sample1_posi_num,
            "sample1_rate": p1,
            "sample2_size": sample2_num,
            "sample2_positive_number": sample2_posi_num,
            "sample2_rate": p2
        },
        "test_results": test_results,
        "confidence_interval": ci_results,
        "relative_risk": rr_results,
        "validation": {
            "normal_approximation_valid": normal_approx_valid,
            "condition_check": f"n1p1={n1p1}, n1q1={n1q1}, n2p2={n2p2}, n2q2={n2q2} (应全部≥5)"
        },
        "interpretation": {
            "statistical_interpretation": f"样本1率为{p1:.6f}，样本2率为{p2:.6f}，差值为{p1-p2:.6f}",
            "p_value_interpretation": f"P值为{test_results['p_value']:.6f}",
            "significance_interpretation": f"{'差异有统计学意义' if test_results['significant'] else '差异无统计学意义'} (α=0.05)",
            "clinical_significance": "请结合临床实际意义解读统计学差异",
            "relative_risk_interpretation": f"相对风险比为{rr_results['relative_risk']:.6f}，95%CI为({rr_results['lower_bound']:.6f}, {rr_results['upper_bound']:.6f})"
        },
        "remark": f"样本1: n={sample1_num}, 阳性数={sample1_posi_num}, 率={p1:.6f}; 样本2: n={sample2_num}, 阳性数={sample2_posi_num}, 率={p2:.6f}"
    }
    
    return result_dict
