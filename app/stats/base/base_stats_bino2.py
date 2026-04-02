"""
临床总体率的区间估计统计分析模块 (Clinical Overall Rate Interval Estimation Statistics Module)

本模块提供全面的总体率区间估计统计分析功能，用于临床数据中总体率的置信区间计算，
是临床试验数据分析的重要组成部分。总体率的区间估计在医学研究中广泛应用，
如计算某种疾病的发病率、某种治疗的成功率、某种不良反应的发生率等的置信区间，
为临床决策提供精确的概率范围。

【模块功能概述】:
1. 区间估计：计算总体率的置信区间（Wilson得分法、正态近似法等）
2. 方法比较：提供多种区间估计方法的结果对比
3. 精度评估：评估估计的精度和可靠性
4. 结果解释：提供统计和临床意义的解释

【临床应用价值】:
- 疗效评估：评估治疗成功率的置信区间
- 疾病发生率：估计特定人群中疾病发生率的置信区间
- 诊断准确性：估计诊断测试敏感性/特异性的置信区间
- 临床决策：为临床决策提供概率范围而非单一估计值

【统计方法选择指南】:
1. Wilson得分法适用条件：
   - 适用于小样本或比例接近0或1的情况
   - 不受正态近似的限制
   - 在边界情况下表现更好

2. 正态近似法适用条件：
   - 样本量足够大
   - np ≥ 5 且 n(1-p) ≥ 5
   - 比例远离0或1

3. 临床应用场景：
   - 治疗成功率的置信区间估计
   - 疾病发生率的区间估计
   - 不良反应发生率的区间估计

【结果解读注意事项】:
1. 置信区间解释：95%置信区间表示真实总体率有95%的可能性落在该区间内
2. 区间宽度解释：区间越窄，估计越精确
3. 临床意义：置信区间应结合临床实际意义进行解读
4. 样本量影响：样本量越大，置信区间越窄
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
- "如何解释总体率的置信区间？"
- "Wilson得分法和正态近似法有什么区别？"
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
from scipy.stats import norm
from app.schemas.request_data.base_param import BaseParamBino2


def wilson_score_interval(n: int, x: int, confidence_level: float) -> Dict[str, float]:
    """
    计算Wilson得分置信区间
    
    在临床研究中，这用于计算总体率的置信区间，特别适用于小样本或比例接近0或1的情况，
    以提供更可靠的区间估计。
    
    Args:
        n: 样本量
        x: 样本阳性数
        confidence_level: 置信水平
        
    Returns:
        Dict[str, float]: 包含置信区间上下限的字典
    """
    p = x / n
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha / 2)
    z2 = z * z
    
    # Wilson得分区间的计算公式
    center = (x + z2 / 2) / (n + z2)
    denominator = 1 + z2 / n
    se = z * math.sqrt((p * (1 - p) / n) + (z2 / (4 * n * n)))
    
    lower_bound = (center - se / denominator)
    upper_bound = (center + se / denominator)
    
    return {
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'center': center
    }


def normal_approximation_interval(n: int, x: int, confidence_level: float) -> Dict[str, float]:
    """
    计算正态近似置信区间
    
    在临床研究中，这用于计算总体率的置信区间，适用于大样本且比例远离0或1的情况，
    是最常用的区间估计方法之一。
    
    Args:
        n: 样本量
        x: 样本阳性数
        confidence_level: 置信水平
        
    Returns:
        Dict[str, float]: 包含置信区间上下限的字典
    """
    p = x / n
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha / 2)
    
    se = math.sqrt(p * (1 - p) / n)
    margin_error = z * se
    
    lower_bound = p - margin_error
    upper_bound = p + margin_error
    
    # 确保区间在[0,1]范围内
    lower_bound = max(0, lower_bound)
    upper_bound = min(1, upper_bound)
    
    return {
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'center': p
    }


def agresti_coull_interval(n: int, x: int, confidence_level: float) -> Dict[str, float]:
    """
    计算Agresti-Coull置信区间
    
    在临床研究中，这用于计算总体率的置信区间，是对正态近似法的一种改进，
    通过增加2个成功和2个失败来提高估计精度。
    
    Args:
        n: 样本量
        x: 样本阳性数
        confidence_level: 置信水平
        
    Returns:
        Dict[str, float]: 包含置信区间上下限的字典
    """
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha / 2)
    
    # 增广样本量和增广比例
    n_tilde = n + z * z
    p_tilde = (x + z * z / 2) / n_tilde
    
    se = math.sqrt(p_tilde * (1 - p_tilde) / n_tilde)
    margin_error = z * se
    
    lower_bound = p_tilde - margin_error
    upper_bound = p_tilde + margin_error
    
    return {
        'lower_bound': lower_bound,
        'upper_bound': upper_bound,
        'center': p_tilde
    }


def calculate_sample_size_for_ci(confidence_level: float, precision: float, p_estimate: float = 0.5) -> int:
    """
    计算达到特定精度所需的样本量
    
    在临床研究设计中，这用于确定获得指定精度的置信区间所需的样本量，
    有助于研究设计阶段的样本量估算。
    
    Args:
        confidence_level: 置信水平
        precision: 所需精度（置信区间半宽）
        p_estimate: 比例预估值（如果未知，默认为0.5以获得最大样本量）
        
    Returns:
        int: 所需样本量
    """
    alpha = 1 - confidence_level
    z = norm.ppf(1 - alpha / 2)
    
    n = (z ** 2 * p_estimate * (1 - p_estimate)) / (precision ** 2)
    
    return math.ceil(n)


def cal_result_bino2(param: BaseParamBino2) -> Dict[str, Any]:
    """
    生成总体率区间估计统计分析的完整报告字典
    
    此函数整合了总体率区间估计的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的区间估计分析结果。
    报告包括多种方法的置信区间、样本率、区间宽度等信息，
    便于临床医生和研究人员快速理解区间估计的特点。
    
    Args:
        param: BaseParamBino2对象，包含confidence_level, sample_size, sample_posi_num
        
    Returns:
        Dict[str, Any]: 包含总体率区间估计统计分析指标的字典，键为指标名称，值为对应的统计量
            - sample_size: 样本量
            - sample_positive_number: 样本阳性数
            - sample_rate: 样本率
            - confidence_level: 置信水平
            - wilson_ci: Wilson得分法置信区间
            - normal_ci: 正态近似法置信区间
            - agresti_coull_ci: Agresti-Coull法置信区间
            - interval_widths: 各方法的区间宽度
    """
    # 从参数对象中提取值
    confidence_level = param.confidence_level
    sample_size = param.sample_size
    sample_posi_num = param.sample_posi_num
    
    # 确保输入参数有效
    if sample_size <= 0:
        raise ValueError("样本量必须大于0")
    if sample_posi_num < 0 or sample_posi_num > sample_size:
        raise ValueError("样本阳性数必须在0和样本量之间")
    if confidence_level <= 0 or confidence_level >= 1:
        raise ValueError("置信水平必须在0和1之间")
    
    # 计算样本率
    sample_rate = sample_posi_num / sample_size
    
    # 计算各种方法的置信区间
    wilson_ci = wilson_score_interval(sample_size, sample_posi_num, confidence_level)
    normal_ci = normal_approximation_interval(sample_size, sample_posi_num, confidence_level)
    agresti_coull_ci = agresti_coull_interval(sample_size, sample_posi_num, confidence_level)
    
    # 计算区间宽度
    wilson_width = wilson_ci['upper_bound'] - wilson_ci['lower_bound']
    normal_width = normal_ci['upper_bound'] - normal_ci['lower_bound']
    agresti_coull_width = agresti_coull_ci['upper_bound'] - agresti_coull_ci['lower_bound']
    
    # 检查是否满足正态近似条件
    np_condition = sample_size * sample_rate >= 5
    nq_condition = sample_size * (1 - sample_rate) >= 5
    normal_approx_valid = np_condition and nq_condition
    
    # 构建结果字典
    result_dict = {
        "table_name": "总体率区间估计",
        "sample_size": sample_size,
        "sample_positive_number": sample_posi_num,
        "sample_rate": sample_rate,
        "confidence_level": confidence_level,
        "wilson_ci": wilson_ci,
        "normal_ci": normal_ci,
        "agresti_coull_ci": agresti_coull_ci,
        "interval_widths": {
            "wilson_method": wilson_width,
            "normal_approximation": normal_width,
            "agresti_coull_method": agresti_coull_width
        },
        "normal_approximation_valid": normal_approx_valid,
        "interpretation": {
            "normal_approx_condition": f"np={sample_size * sample_rate:.2f}>=5 and n(1-p)={sample_size * (1 - sample_rate):.2f}>=5",
            "validity": "正态近似法适用于当前数据" if normal_approx_valid else "正态近似法可能不适用于当前数据，请优先参考Wilson得分法结果",
            "method_comparison": "通常情况下，Wilson得分法在各种条件下都表现良好，是首选方法"
        },
        "remark": f"样本量={sample_size}, 样本阳性数={sample_posi_num}, 样本率={sample_rate:.6f}, 置信水平={confidence_level}"
    }
    
    return result_dict
