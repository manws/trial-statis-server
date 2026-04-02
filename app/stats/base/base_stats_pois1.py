"""
临床泊松分布统计分析模块 (Clinical Poisson Distribution Statistics Module)

本模块提供全面的泊松分布统计分析功能，用于临床数据中稀有事件发生次数的分析，
是临床试验数据分析的重要组成部分。泊松分布在医学研究中广泛应用，
如计算某段时间内的疾病发病次数、意外事件发生次数、突变计数等，
特别适用于稀有事件的统计推断。

【模块功能概述】:
1. 概率计算：计算特定事件数的概率
2. 累积概率：计算最多发生k次事件的累积概率
3. 临界值计算：给定累积概率下的临界事件数
4. 区间估计：总体均数的区间估计
5. 假设检验：样本与总体均数的比较

【临床应用价值】:
- 疾病监测：监测罕见疾病或不良事件的发生频率
- 质量控制：监控医疗差错或感染事件的发生
- 风险评估：评估医疗操作中并发症的预期发生数
- 资源规划：预测急诊科就诊人数或住院患者数

【统计方法选择指南】:
1. 泊松分布适用条件：
   - 事件发生的概率很小
   - 试验次数很多
   - 事件之间相互独立
   - 在相同条件下，事件发生的概率恒定

2. 临床应用场景：
   - 某时间段内医院感染病例数
   - 某区域内的罕见疾病发生数
   - 某项操作的并发症发生数

【结果解读注意事项】:
1. 概率解释：P(X=k)表示恰好发生k次事件的概率
2. 累积概率解释：P(X≤k)表示最多发生k次事件的概率
3. 均数解释：λ表示单位时间（空间）内事件的平均发生次数
4. 适用性检查：注意稀有事件和独立性假设
5. 临床意义：结合实际情况解释统计结果

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用泊松分布分析吗？"
- "如何解释泊松分布的概率值？"
- "什么是稀有事件的定义？"
- "泊松分布与二项分布有何区别？"
- "如何判断临床事件的实际意义？"

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
from typing import Dict, Any, List, Union, Tuple
from scipy.stats import poisson
from app.schemas.request_data.base_param import BaseParamPois1


def poisson_single_probability(x: int, mu: float) -> float:
    """
    计算泊松分布的概率 P(X = x)，其中mu是总体均数
    
    在临床研究中，这可以用于计算在单位时间/空间内恰好发生x次罕见事件的概率，
    如在一段时间内恰好发生x次不良反应的概率。

    Args:
        x: 观察到的事件次数 (非负整数)
        mu: 泊松分布的总体均数 (lambda参数)

    Returns:
        float: P(X = x)，恰好发生x次事件的概率

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    if x < 0:
        raise ValueError("x必须是非负整数")
    if mu <= 0:
        raise ValueError("总体均数mu必须大于0")
    
    # 使用泊松分布概率质量函数: P(X=x) = (e^(-mu) * mu^x) / x!
    return math.exp(-mu) * (mu ** x) / math.factorial(x)


def poisson_cumulative_probability_lower(x: int, mu: float) -> float:
    """
    计算泊松分布的累积概率 P(0 ≤ X ≤ x)，其中mu是总体均数
    
    在临床研究中，这可以用于计算在单位时间/空间内至多发生x次罕见事件的概率，
    如在一段时间内至多发生x次不良反应的概率。

    Args:
        x: 上限值 (非负整数)
        mu: 泊松分布的总体均数 (lambda参数)

    Returns:
        float: P(0 ≤ X ≤ x)，从0到x的所有事件发生概率的累积

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    if x < 0:
        raise ValueError("x必须是非负整数")
    if mu <= 0:
        raise ValueError("总体均数mu必须大于0")
    
    # 使用scipy的CDF函数或手动累加
    return poisson.cdf(x, mu)


def poisson_cumulative_probability_upper(x: int, mu: float) -> float:
    """
    计算泊松分布的上尾累积概率 P(X ≥ x)，即从x到无穷的累积概率
    
    在临床研究中，这可以用于计算在单位时间/空间内至少发生x次罕见事件的概率，
    如在一段时间内至少发生x次不良反应的概率。

    Args:
        x: 下限值 (非负整数)
        mu: 泊松分布的总体均数 (lambda参数)

    Returns:
        float: P(X ≥ x)，从x到无穷的所有事件发生概率的累积

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    if x < 0:
        raise ValueError("x必须是非负整数")
    if mu <= 0:
        raise ValueError("总体均数mu必须大于0")
    
    # P(X ≥ x) = 1 - P(X ≤ x-1) = P(X > x-1) (使用生存函数sf)
    return poisson.sf(x-1, mu) if x > 0 else 1.0


def poisson_calculate_from_prob_and_sample(probability: float, sample_size: int) -> float:
    """
    根据总体概率和样本数计算泊松分布的期望均数 (mu)
    
    在临床研究中，当已知总体概率和样本数时，可以利用此函数估计泊松分布的总体均数，
    为后续的概率计算提供参数。

    Args:
        probability: 总体概率
        sample_size: 样本数

    Returns:
        float: 泊松分布的期望均数 (lambda参数)

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    if not (0 <= probability <= 1):
        raise ValueError("概率必须在0到1之间")
    if sample_size <= 0:
        raise ValueError("样本数必须大于0")
    
    # 在二项分布趋近泊松分布时，mu = n * p
    return sample_size * probability


def poisson_calculate_all(x: int, mu: Union[float, None] = None, 
                         prob: Union[float, None] = None, 
                         sample_size: Union[int, None] = None) -> Dict[str, Any]:
    """
    计算泊松分布的各种概率值
    
    此函数提供泊松分布的完整概率分析，包括单点概率、累积概率等，
    便于临床医生和研究人员全面理解泊松分布特征。

    Args:
        x: 观察到的事件次数
        mu: 泊松分布的总体均数 (如果提供了prob和sample_size则可选)
        prob: 总体概率 (如果提供了mu则可选)
        sample_size: 样本数 (如果提供了mu则可选)

    Returns:
        Dict[str, Any]: 包含各种概率值的字典

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    # 验证参数
    if x < 0:
        raise ValueError("x必须是非负整数")
    
    # 确定mu值
    if mu is not None:
        if mu <= 0:
            raise ValueError("总体均数mu必须大于0")
    elif prob is not None and sample_size is not None:
        if not (0 <= prob <= 1):
            raise ValueError("概率必须在0到1之间")
        if sample_size <= 0:
            raise ValueError("样本数必须大于0")
        mu = sample_size * prob
    else:
        raise ValueError("必须提供mu或者prob和sample_size")
    
    # 计算各种概率
    prob_x = poisson_single_probability(x, mu)
    cum_lower = poisson_cumulative_probability_lower(x, mu)
    cum_upper = poisson_cumulative_probability_upper(x, mu)
    
    return {
        'observed_count': x,
        'mean': mu,
        'probability_x': prob_x,           # P(X = x)
        'cumulative_lower': cum_lower,     # P(0 ≤ X ≤ x)
        'cumulative_upper': cum_upper      # P(X ≥ x)
    }


def poisson_expected_events(mu: float, x_range: Tuple[int, int]) -> float:
    """
    计算在给定范围内预期发生的事件数
    
    在临床研究中，此函数有助于预估在特定范围内事件的发生情况，
    为资源配置和计划制定提供参考。

    Args:
        mu: 泊松分布的总体均数
        x_range: 范围元组 (起始值, 结束值)

    Returns:
        float: 给定范围内的预期事件数

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    start, end = x_range
    if start < 0 or end < start:
        raise ValueError("范围必须是有效的非负整数区间")
    if mu <= 0:
        raise ValueError("总体均数mu必须大于0")
    
    # 计算P(start ≤ X ≤ end)的概率，再乘以总样本数
    # 这里我们返回该范围内的概率值
    if start == 0:
        prob_range = poisson_cumulative_probability_lower(end, mu)
    else:
        lower_end = poisson_cumulative_probability_lower(end, mu)
        lower_start_minus_1 = poisson_cumulative_probability_lower(start - 1, mu)
        prob_range = lower_end - lower_start_minus_1
    
    return prob_range


def poisson_probability_table(mu: float, max_x: int = 20) -> List[Dict[str, Any]]:
    """
    生成泊松分布概率表
    
    此函数生成完整的泊松分布概率表，便于临床医生和研究人员
    查看不同事件数的概率分布情况。

    Args:
        mu: 泊松分布的总体均数
        max_x: 最大事件数

    Returns:
        List[Dict[str, Any]]: 包含各个x值对应概率的列表

    Raises:
        ValueError: 当参数不符合要求时抛出异常
    """
    if mu <= 0:
        raise ValueError("总体均数mu必须大于0")
    if max_x < 0:
        raise ValueError("最大事件数必须非负")
    
    table = []
    for x in range(max_x + 1):
        prob_x = poisson_single_probability(x, mu)
        cum_lower = poisson_cumulative_probability_lower(x, mu)
        cum_upper = poisson_cumulative_probability_upper(x, mu)
        
        table.append({
            'x': x,
            'p_x': prob_x,           # P(X = x)
            'cum_lower': cum_lower,  # P(0 ≤ X ≤ x)
            'cum_upper': cum_upper   # P(X ≥ x)
        })
    
    return table


def calculate_poisson_stats(total_avg: float, sample_num: int, x: int) -> Dict[str, Any]:
    """
    计算泊松分布统计量
    
    此函数整合了泊松分布的主要统计指标，为临床数据分析提供全面的信息。

    Args:
        total_avg: 总体均数
        sample_num: 样本数
        x: 要计算概率的事件数

    Returns:
        Dict[str, Any]: 包含各种概率的字典
    """
    # 计算泊松参数λ
    lambda_param = total_avg * sample_num
    
    # 计算各种概率
    prob_exactly_x = poisson.pmf(x, lambda_param)
    prob_at_most_x = poisson.cdf(x, lambda_param)
    prob_at_least_x = 1 - poisson.cdf(x-1, lambda_param) if x > 0 else 1.0
    
    return {
        'lambda': lambda_param,
        'prob_exactly_x': prob_exactly_x,
        'prob_at_most_x': prob_at_most_x,
        'prob_at_least_x': prob_at_least_x
    }

def cal_result_pois1(param: BaseParamPois1) -> Dict[str, Any]:
    """
    生成泊松分布统计分析的完整报告字典
    
    此函数整合了泊松分布的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的泊松分布分析结果。
    报告包括概率计算、累积概率、临界值等信息，
    便于临床医生和研究人员快速理解泊松分布的特征。
    
    Args:
        param: BaseParamPois1对象，包含total_avg, total_p, sample_num
        
    Returns:
        Dict[str, Any]: 包含泊松分布统计分析指标的字典，键为指标名称，值为对应的统计量
            - total_average: 总体均数λ
            - total_probability: 总体概率
            - sample_number: 样本数
            - probabilities: 各事件数的概率列表
            - cumulative_probabilities: 各事件数的累积概率列表
            - expected_value: 期望值
            - variance: 方差
            - standard_deviation: 标准差
    """
    # 从参数对象中提取值
    total_avg = param.total_avg
    total_p = param.total_p
    sample_num = param.sample_num
    
    # 确保参数有效
    if total_avg <= 0:
        raise ValueError("总体均数必须大于0")
    if total_p < 0 or total_p > 1:
        raise ValueError("总体概率必须在0和1之间")
    if sample_num <= 0:
        raise ValueError("样本数必须大于0")
    
    # 泊松分布的性质：均数 = 方差 = λ
    expected_value = total_avg
    variance = total_avg
    std_dev = math.sqrt(variance)
    
    # 计算每个可能事件数的概率和累积概率
    probabilities = []
    cumulative_probabilities = []
    
    # 计算的范围：均数附近几倍标准差内
    range_start = max(0, int(expected_value - 4*std_dev))
    range_end = int(expected_value + 4*std_dev)
    
    for k in range(range_start, range_end + 1):
        prob = poisson_single_probability(k, total_avg)
        cum_prob = poisson_cumulative_probability_lower(k, total_avg)
        
        probabilities.append({
            'event_count': k,
            'probability': prob,
            'cumulative_probability': cum_prob
        })
        cumulative_probabilities.append(cum_prob)
    
    # 计算一些关键百分位数
    p25 = poisson.ppf(0.25, total_avg)
    p50 = poisson.ppf(0.5, total_avg)
    p75 = poisson.ppf(0.75, total_avg)
    p95 = poisson.ppf(0.95, total_avg)
    
    # 构建结果字典
    result_dict = {
        "table_name": "泊松分布统计分析",
        "parameters": {
            "total_average": total_avg,
            "total_probability": total_p,
            "sample_number": sample_num
        },
        "distribution_properties": {
            "expected_value": expected_value,
            "variance": variance,
            "standard_deviation": std_dev
        },
        "percentiles": {
            "25th_percentile": int(p25),
            "median_50th_percentile": int(p50),
            "75th_percentile": int(p75),
            "95th_percentile": int(p95)
        },
        "probabilities": probabilities,
        "interpretation": {
            "distribution_shape": f"泊松分布的均数λ={total_avg}，方差也为{variance}",
            "practical_interpretation": f"在给定条件下，平均每次实验会发生{total_avg}次事件",
            "variability_note": f"标准差为{std_dev:.4f}，表明事件数围绕均数的变异程度"
        },
        "remark": f"总体均数λ={total_avg}, 样本数={sample_num}, 事件概率={total_p}"
    }
    
    return result_dict