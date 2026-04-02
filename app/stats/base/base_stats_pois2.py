"""
临床泊松分布总体均数区间估计统计分析模块 (Clinical Poisson Distribution Population Mean Interval Estimation Statistics Module)

本模块提供全面的泊松分布总体均数区间估计统计分析功能，用于临床数据中总体均数的置信区间计算，
是临床试验数据分析的重要组成部分。总体均数的区间估计在医学研究中广泛应用，
如估计某种罕见疾病的年发病率、某种不良事件的发生率、某种病原体的检出率等的置信区间，
为临床决策提供精确的均数范围。

【模块功能概述】:
1. 精确区间估计：使用泊松分布的精确方法计算置信区间
2. 正态近似区间：使用正态近似法计算置信区间
3. 精度评估：评估估计的精度和可靠性
4. 结果解释：提供统计和临床意义的解释

【临床应用价值】:
- 疾病监测：估计罕见疾病的发病率置信区间
- 安全评估：估计不良事件发生率的置信区间
- 质量控制：估计医疗差错发生率的置信区间
- 资源规划：估计急诊就诊人数的置信区间

【统计方法选择指南】:
1. 精确方法适用条件：
   - 适用于所有样本量
   - 特别适用于小样本或均数较小的情况

2. 正态近似法适用条件：
   - 均数足够大（通常λ ≥ 5）
   - 一般适用于观察单位足够大的情况

3. 临床应用场景：
   - 罕见疾病发病率的区间估计
   - 不良事件发生率的区间估计
   - 感染事件发生率的区间估计

【结果解读注意事项】:
1. 置信区间解释：95%置信区间表示真实总体均数有95%的可能性落在该区间内
2. 区间宽度解释：区间越窄，估计越精确
3. 临床意义：置信区间应结合临床实际意义进行解读
4. 样本量影响：观察单位越大，置信区间越窄
5. 临床决策：置信区间可用于临床决策的参考

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用哪种区间估计方法？"
- "如何解释总体均数的置信区间？"
- "精确法和正态近似法有什么区别？"
- "置信区间的宽度意味着什么？"
- "如何基于置信区间做出临床决策？"

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
from scipy.stats import chi2
from app.schemas.request_data.base_param import BaseParamPois2


def exact_confidence_interval(count: int, confidence_level: float) -> Dict[str, float]:
    """
    使用精确方法计算泊松分布总体均数的置信区间
    
    在临床研究中，这用于计算总体均数的精确置信区间，
    特别适用于小样本或均数较小的情况。
    
    Args:
        count: 观察到的事件数
        confidence_level: 置信水平
        
    Returns:
        Dict[str, float]: 包含置信区间上下限的字典
    """
    alpha = 1 - confidence_level
    
    # 下限：使用χ²分布
    lower_bound = chi2.ppf(alpha/2, 2*count) / 2 if count > 0 else 0
    
    # 上限：使用χ²分布
    upper_bound = chi2.ppf(1 - alpha/2, 2*count + 2) / 2
    
    return {
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'point_estimate': count
    }


def normal_approximation_interval(count: int, confidence_level: float) -> Dict[str, float]:
    """
    使用正态近似法计算泊松分布总体均数的置信区间
    
    在临床研究中，这用于计算总体均数的置信区间，适用于大均数的情况，
    是常用的区间估计方法之一。
    
    Args:
        count: 观察到的事件数
        confidence_level: 置信水平
        
    Returns:
        Dict[str, float]: 包含置信区间上下限的字典
    """
    from scipy.stats import norm
    
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha/2)
    
    # 泊松分布的性质：均数 = 方差
    se = math.sqrt(count)
    margin_error = z * se
    
    lower_bound = count - margin_error
    upper_bound = count + margin_error
    
    # 确保下限非负
    lower_bound = max(0, lower_bound)
    
    return {
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'point_estimate': count,
        'standard_error': se
    }


def calculate_sample_size_for_precision(confidence_level: float, precision: float, lambda_estimate: float) -> float:
    """
    计算达到特定精度所需的观察单位大小
    
    在临床研究设计中，这用于确定获得指定精度的置信区间所需的观察单位，
    有助于研究设计阶段的规划。
    
    Args:
        confidence_level: 置信水平
        precision: 所需精度（置信区间半宽）
        lambda_estimate: 均数预估值
        
    Returns:
        float: 所需的观察单位大小
    """
    from scipy.stats import norm
    
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha/2)
    
    # 精度公式：precision = z * sqrt(lambda / n)
    # 解得：n = (z^2 * lambda) / precision^2
    n = (z ** 2 * lambda_estimate) / (precision ** 2)
    
    return n


def cal_result_pois2(param: BaseParamPois2) -> Dict[str, Any]:
    """
    生成泊松分布总体均数区间估计统计分析的完整报告字典
    
    此函数整合了泊松分布总体均数区间估计的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的区间估计分析结果。
    报告包括精确法和正态近似法的置信区间、均数估计、区间宽度等信息，
    便于临床医生和研究人员快速理解区间估计的特点。
    
    Args:
        param: BaseParamPois2对象，包含confidence_level, total_avg, sample_posi_num
        
    Returns:
        Dict[str, Any]: 包含泊松分布总体均数区间估计统计分析指标的字典，键为指标名称，值为对应的统计量
            - sample_positive_number: 样本阳性数（观察到的事件数）
            - confidence_level: 置信水平
            - exact_ci: 精确法置信区间
            - normal_ci: 正态近似法置信区间
            - interval_widths: 各方法的区间宽度
    """
    # 从参数对象中提取值
    confidence_level = param.confidence_level
    total_avg = param.total_avg
    sample_posi_num = param.sample_posi_num
    
    # 确保输入参数有效
    if confidence_level <= 0 or confidence_level >= 1:
        raise ValueError("置信水平必须在0和1之间")
    if sample_posi_num < 0:
        raise ValueError("样本阳性数（事件数）不能为负数")
    
    # 计算精确法置信区间
    exact_ci = exact_confidence_interval(sample_posi_num, confidence_level)
    
    # 计算正态近似法置信区间
    normal_ci = normal_approximation_interval(sample_posi_num, confidence_level)
    
    # 计算区间宽度
    exact_width = exact_ci['upper_bound'] - exact_ci['lower_bound']
    normal_width = normal_ci['upper_bound'] - normal_ci['lower_bound']
    
    # 检查正态近似条件
    normal_approx_valid = sample_posi_num >= 5
    
    # 计算达到特定精度所需的观察单位
    desired_precision = (exact_ci['upper_bound'] - exact_ci['lower_bound']) / 2
    try:
        required_sample_size = calculate_sample_size_for_precision(
            confidence_level, desired_precision, sample_posi_num
        )
    except:
        required_sample_size = None
    
    # 构建结果字典
    result_dict = {
        "table_name": "泊松分布总体均数区间估计",
        "input_parameters": {
            "observed_event_count": sample_posi_num,
            "confidence_level": confidence_level,
            "total_average": total_avg
        },
        "exact_method": exact_ci,
        "normal_approximation_method": normal_ci,
        "interval_widths": {
            "exact_method": exact_width,
            "normal_approximation": normal_width
        },
        "validation": {
            "normal_approximation_valid": normal_approx_valid,
            "normal_approx_condition": f"观察事件数={sample_posi_num}{'≥5，正态近似有效' if normal_approx_valid else '<5，建议使用精确法'}"
        },
        "required_sample_size_for_precision": required_sample_size,
        "interpretation": {
            "exact_method_interpretation": f"精确法95%置信区间为[{exact_ci['lower_bound']:.6f}, {exact_ci['upper_bound']:.6f}]",
            "normal_method_interpretation": f"正态近似法95%置信区间为[{normal_ci['lower_bound']:.6f}, {normal_ci['upper_bound']:.6f}]",
            "recommended_method": "精确法" if not normal_approx_valid else "两种方法均可，精确法更可靠",
            "clinical_significance": "置信区间提供了总体均数的可能范围，有助于临床决策"
        },
        "remark": f"观察事件数={sample_posi_num}, 置信水平={confidence_level}, 总体均数={total_avg}"
    }
    
    return result_dict