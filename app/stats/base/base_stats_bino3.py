"""
临床样本率与总体率比较统计分析模块 (Clinical Sample Rate vs Population Rate Comparison Statistics Module)

本模块提供全面的样本率与总体率比较统计分析功能，用于临床数据中样本率与已知总体率的差异检验，
是临床试验数据分析的重要组成部分。样本率与总体率比较在医学研究中广泛应用，
如评估某种新治疗方法的成功率是否显著高于历史平均水平、某种疫苗接种率是否达到目标要求等。

【模块功能概述】:
1. 精确检验：使用二项分布精确计算P值
2. 正态近似检验：使用正态近似法计算Z统计量和P值
3. 单侧双侧检验：支持左侧、右侧和双侧检验
4. 效应量计算：计算差异的效应量
5. 检验效能：计算检验的效能

【临床应用价值】:
- 疗效评估：评估新疗法是否显著优于既往标准
- 质量控制：评估医疗指标是否达到既定标准
- 政策评估：评估公共卫生干预措施是否达到预期目标
- 安全监控：监测不良事件发生率是否超过安全阈值

【统计方法选择指南】:
1. 精确检验适用条件：
   - 任何样本量下都适用
   - 特别适用于小样本或比例接近0或1的情况

2. 正态近似法适用条件：
   - 样本量较大
   - np0 ≥ 5 且 n(1-p0) ≥ 5（其中p0为总体率）

3. 临床应用场景：
   - 新疗法成功率是否显著高于历史水平
   - 疫苗接种率是否达到目标要求
   - 院内感染率是否控制在标准范围内

【结果解读注意事项】:
1. P值解释：在零假设成立的前提下，观察到当前或更极端结果的概率
2. 显著性水平：通常使用α=0.05作为判断标准
3. 单双侧选择：根据研究目的选择合适的检验方向
4. 临床意义：统计显著性不等同于临床重要性
5. 检验效能：注意小样本可能导致的低效能问题

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "如何比较样本率与已知总体率？"
- "二项检验和正态近似法有什么区别？"
- "如何选择单侧还是双侧检验？"
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
from scipy.stats import norm, binom
from app.schemas.request_data.base_param import BaseParamBino3


def exact_test_binomial(x: int, n: int, p0: float, alternative: str = 'two-sided') -> Dict[str, float]:
    """
    二项分布精确检验
    
    在临床研究中，这用于检验样本阳性数与理论阳性数之间的差异，
    通过计算二项分布的确切概率来判断差异的统计学意义。
    
    Args:
        x: 样本阳性数
        n: 样本量
        p0: 总体阳性概率（零假设下的概率）
        alternative: 备择假设类型 ('two-sided', 'less', 'greater')
        
    Returns:
        Dict[str, float]: 包含检验统计量和P值的字典
    """
    # 计算观测值的概率
    obs_prob = binom.pmf(x, n, p0)
    
    if alternative == 'two-sided':
        # 找到所有等于或低于观测概率的x值，并计算其概率之和
        p_value = 0.0
        # 计算期望值附近的概率
        expected_x = n * p0
        
        # 计算所有pmf值小于等于obs_prob的x值
        for i in range(n + 1):
            if binom.pmf(i, n, p0) <= obs_prob + 1e-10:  # 加一个小的容差
                p_value += binom.pmf(i, n, p0)
                
        # 确保p值不超过1
        p_value = min(p_value, 1.0)
    elif alternative == 'less':
        # 左侧检验：P(X ≤ x)
        p_value = binom.cdf(x, n, p0)
    else:  # greater
        # 右侧检验：P(X ≥ x)
        p_value = 1 - binom.cdf(x - 1, n, p0)
    
    return {
        'p_value': p_value,
        'observed_prob': x/n,
        'theoretical_prob': p0
    }


def normal_approximation_test(x: int, n: int, p0: float, alternative: str = 'two-sided') -> Dict[str, float]:
    """
    正态近似法检验
    
    在临床研究中，这用于检验样本率与总体率之间的差异，
    通过正态近似计算Z统计量来判断差异的统计学意义。
    注意：此方法要求样本量足够大，且满足np0 ≥ 5和n(1-p0) ≥ 5。
    
    Args:
        x: 样本阳性数
        n: 样本量
        p0: 总体阳性概率（零假设下的概率）
        alternative: 备择假设类型 ('two-sided', 'less', 'greater')
        
    Returns:
        Dict[str, float]: 包含检验统计量和P值的字典
    """
    # 检查正态近似的适用条件
    if n * p0 < 5 or n * (1 - p0) < 5:
        print(f"警告：正态近似条件可能不满足 (np0={n*p0}, n(1-p0)={n*(1-p0)})，建议使用精确检验")
    
    sample_p = x / n
    
    # 计算标准误
    se = math.sqrt(p0 * (1 - p0) / n)
    
    # 计算Z统计量
    if se == 0:
        # 如果标准误为0，说明p0是0或1，无法进行检验
        raise ValueError("总体率不能为0或1，无法进行正态近似检验")
    
    z_stat = (sample_p - p0) / se
    
    # 计算P值
    if alternative == 'two-sided':
        p_value = 2 * (1 - norm.cdf(abs(z_stat)))
    elif alternative == 'less':
        p_value = norm.cdf(z_stat)
    else:  # greater
        p_value = 1 - norm.cdf(z_stat)
    
    return {
        'z_statistic': z_stat,
        'p_value': p_value,
        'observed_prob': sample_p,
        'theoretical_prob': p0,
        'se': se
    }


def calculate_power(x: int, n: int, p0: float, alpha: float = 0.05) -> float:
    """
    计算检验效能
    
    在临床研究中，这用于评估当前样本量下检测特定差异的能力，
    有助于研究设计和结果解释。
    
    Args:
        x: 样本阳性数
        n: 样本量
        p0: 总体阳性概率
        alpha: 显著性水平
        
    Returns:
        float: 检验效能值
    """
    # 这是一个简化的功效计算，实际功效分析会更复杂
    observed_p = x / n
    effect_size = abs(observed_p - p0) / math.sqrt(p0 * (1 - p0) / n) if p0 * (1 - p0) / n > 0 else 0
    
    # 简化的功效估计
    # 实际的功效分析会考虑更多因素，这里仅作示意
    power = norm.cdf(effect_size - norm.ppf(1 - alpha/2)) if effect_size > 0 else 0.5
    
    return power


def cal_result_bino3(param: BaseParamBino3) -> Dict[str, Any]:
    """
    生成样本率与总体率比较统计分析的完整报告字典
    
    此函数整合了样本率与总体率比较的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的比较分析结果。
    报告包括精确检验、正态近似检验、P值、显著性判断等信息，
    便于临床医生和研究人员快速理解比较分析的特征。
    
    Args:
        param: BaseParamBino3对象，包含total_posi_p, sample_size, sample_posi_num
        
    Returns:
        Dict[str, Any]: 包含样本率与总体率比较统计分析指标的字典，键为指标名称，值为对应的统计量
            - total_positive_probability: 总体阳性概率
            - sample_size: 样本量
            - sample_positive_number: 样本阳性数
            - sample_rate: 样本率
            - exact_test: 精确检验结果
            - normal_approximation_test: 正态近似检验结果
            - power: 检验效能
            - interpretation: 结果的专业解释
    """
    # 从参数对象中提取值
    total_posi_p = param.total_posi_p
    sample_size = param.sample_size
    sample_posi_num = param.sample_posi_num
    
    # 确保输入参数有效
    if sample_size <= 0:
        raise ValueError("样本量必须大于0")
    if sample_posi_num < 0 or sample_posi_num > sample_size:
        raise ValueError("样本阳性数必须在0和样本量之间")
    if total_posi_p < 0 or total_posi_p > 1:
        raise ValueError("总体阳性概率必须在0和1之间")
    
    # 计算样本率
    sample_rate = sample_posi_num / sample_size
    
    # 执行精确检验
    exact_results = exact_test_binomial(sample_posi_num, sample_size, total_posi_p)
    
    # 执行正态近似检验
    normal_results = {}
    try:
        normal_results = normal_approximation_test(sample_posi_num, sample_size, total_posi_p)
    except ValueError as e:
        normal_results = {'error': str(e)}
    
    # 计算检验效能
    try:
        power = calculate_power(sample_posi_num, sample_size, total_posi_p)
    except:
        power = None
    
    # 检查正态近似条件
    normal_approx_valid = sample_size * total_posi_p >= 5 and sample_size * (1 - total_posi_p) >= 5
    
    # 构建结果字典
    result_dict = {
        "table_name": "样本率与总体率比较",
        "input_parameters": {
            "total_positive_probability": total_posi_p,
            "sample_size": sample_size,
            "sample_positive_number": sample_posi_num,
            "sample_rate": sample_rate
        },
        "exact_test": exact_results,
        "normal_approximation_test": normal_results,
        "power": power,
        "validation": {
            "normal_approximation_valid": normal_approx_valid,
            "condition_check": f"np0={sample_size * total_posi_p}, n(1-p0)={sample_size * (1 - total_posi_p)}"
        },
        "interpretation": {
            "statistical_interpretation": f"样本率为{sample_rate:.6f}，总体率为{total_posi_p}",
            "exact_test_p_value": f"P值为{exact_results['p_value']:.6f}",
            "normal_test_p_value": f"正态近似P值为{normal_results.get('p_value', 'N/A'):.6f}" if 'p_value' in normal_results else "正态近似不可用",
            "significance_05": f"{'差异有统计学意义' if exact_results['p_value'] < 0.05 else '差异无统计学意义'} (α=0.05)",
            "clinical_significance": "请结合临床实际意义解读统计学差异"
        },
        "remark": f"总体率={total_posi_p}, 样本量={sample_size}, 样本阳性数={sample_posi_num}, 样本率={sample_rate:.6f}"
    }
    
    return result_dict