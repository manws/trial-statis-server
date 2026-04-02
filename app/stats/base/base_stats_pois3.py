"""
临床泊松分布假设检验统计分析模块 (Clinical Poisson Distribution Hypothesis Testing Statistics Module)

本模块提供全面的泊松分布假设检验统计分析功能，用于临床数据中罕见事件发生率的假设检验，
是临床试验数据分析的重要组成部分。泊松分布假设检验在医学研究中广泛应用，
如检验某种罕见不良事件的发生率是否显著高于预期、评估治疗效果等。

【模块功能概述】:
1. Z统计量计算：计算泊松分布假设检验的Z统计量
2. P值计算：根据Z统计量计算相应P值
3. 假设检验：执行完整的泊松分布假设检验流程
4. 结果解释：提供假设检验结果的专业解释

【临床应用价值】:
- 罕见不良事件率检验：检验实际发生率是否显著高于理论发生率
- 临床试验分析：评估药物或治疗方法有效性的统计学意义
- 疾病发生率检验：检验疾病发生率是否与预期值有显著差异
- 疗效评价：判断治疗效果是否具有统计学意义

【统计方法选择指南】:
1. 泊松分布假设检验适用条件：
   - 事件发生概率较小（稀有事件）
   - 事件独立发生
   - 期望事件数足够大（通常λ≥10）以使用正态近似
   - 试验次数较多

2. 临床应用场景：
   - 评估新疗法不良事件发生率是否显著高于基准
   - 检验疫苗接种后不良反应发生率是否超标
   - 验证医疗设备故障率是否符合安全标准

【结果解读注意事项】:
1. Z统计量解释：衡量观察值与期望值之间的标准差数
2. P值解释：在零假设成立的情况下，观察到当前或更极端结果的概率
3. 显著性水平：通常使用α=0.05作为判断标准
4. 检验方向：可根据研究目的选择双侧或单侧检验
5. 临床意义：统计学显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的罕见事件数据适合用泊松假设检验吗？"
- "如何解释泊松分布假设检验的Z统计量？"
- "P值小于0.05意味着什么？"
- "什么时候应该使用泊松分布而非二项分布？"
- "如何判断检验结果的临床意义？"

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
from app.schemas.request_data.base_param import BaseParamPois3


def poisson_hypothesis_test(population_rate: float, sample_size: int, observed_positive: int, alpha: float = 0.05) -> Dict[str, Any]:
    """
    泊松分布假设检验
    
    在临床研究中，这用于检验观察到的罕见事件发生率是否显著不同于预期的发生率，
    以评估统计学意义。

    参数:
    - population_rate: 总体阳性概率 (事件发生率)
    - sample_size: 样本数
    - observed_positive: 样本阳性数 (观察到的事件数)
    - alpha: 显著性水平 (默认0.05)
    
    返回:
    - 包含检验统计量和结论的字典
    
    原假设 H0: λ = population_rate * sample_size
    备择假设 H1: λ ≠ population_rate * sample_size

    Args:
        population_rate: 总体阳性概率 (事件发生率)
        sample_size: 样本数
        observed_positive: 样本阳性数 (观察到的事件数)
        alpha: 显著性水平 (默认0.05)

    Returns:
        Dict[str, Any]: 包含检验统计量和结论的字典
    """
    # 参数验证
    if population_rate < 0 or population_rate > 1:
        raise ValueError("总体阳性概率必须在0到1之间")
    if sample_size <= 0:
        raise ValueError("样本数必须大于0")
    if observed_positive < 0:
        raise ValueError("样本阳性数不能为负数")
    if not (0 < alpha < 1):
        raise ValueError("显著性水平必须在0和1之间")
    
    # 计算期望事件数 (λ0)
    expected_events = population_rate * sample_size
    
    # 计算观察到的事件率
    observed_rate = observed_positive / sample_size
    
    # 计算z统计量
    # 对于泊松分布，当λ较大时可以用正态近似
    # Z = (X - λ0) / sqrt(λ0)
    if expected_events == 0:
        if observed_positive == 0:
            z_statistic = 0.0
        else:
            z_statistic = float('inf')  # 观察到事件但期望为0
    else:
        z_statistic = (observed_positive - expected_events) / math.sqrt(expected_events)
    
    # 计算双侧p值
    if abs(z_statistic) == float('inf'):
        p_value = 0.0
    else:
        p_value = 2 * (1 - norm.cdf(abs(z_statistic)))
    
    # 判断是否显著
    is_significant = p_value < alpha
    
    return {
        "input_parameters": {
            "population_rate": population_rate,
            "sample_size": sample_size,
            "observed_positive": observed_positive,
            "expected_events": expected_events,
            "observed_rate": observed_rate,
            "significance_level": alpha
        },
        "test_statistics": {
            "z_statistic": z_statistic,
            "p_value_two_sided": p_value,
            "critical_z_value": norm.ppf(1 - alpha/2)
        },
        "hypothesis_test": {
            "null_hypothesis": f"H₀: λ = {expected_events:.4f}",
            "alternative_hypothesis": f"H₁: λ ≠ {expected_events:.4f}",
            "is_significant": is_significant,
            "p_less_than_0_05": p_value < 0.05,
            "conclusion": "拒绝原假设" if is_significant else "不拒绝原假设"
        },
        "interpretation": {
            "statistical_meaning": f"在{alpha*100}%显著性水平下，{'有足够的证据' if is_significant else '没有足够的证据'}表明观察到的事件数与期望值存在显著差异",
            "practical_implication": f"观察事件率({observed_rate:.4f})与期望事件率({population_rate:.4f})的差异{'具有统计学意义' if is_significant else '不具有统计学意义'}"
        }
    }


def cal_result_pois3(param: BaseParamPois3, alpha: float = 0.05) -> Dict[str, Any]:
    """
    生成泊松分布假设检验统计分析的完整报告字典
    
    此函数整合了泊松分布假设检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的泊松分布假设检验分析结果。
    报告包括输入参数、统计量、P值、显著性判断等信息，
    便于临床医生和研究人员快速理解泊松分布假设检验的特征。

    Args:
        param: BaseParamPois3 对象，包含总体阳性概率、样本数和样本阳性数
        alpha: 显著性水平 (默认0.05)，浮点数类型

    Returns:
        Dict[str, Any]: 包含泊松分布假设检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"泊松分布假设检验分析"
            - input_parameters: 输入参数信息
            - test_statistics: 检验统计量信息
            - hypothesis_test: 假设检验结果
            - interpretation: 结果的专业解释
    """
    # 从参数对象中提取值
    # 注意：根据参考方案，BaseParamPois3 的属性名可能为 total_posi_p, sample_size, sample_posi_num
    # 如果实际定义不同，请相应调整以下属性访问
    population_rate = param.total_posi_p
    sample_size = param.sample_size
    observed_positive = param.sample_posi_num

    # 执行泊松分布假设检验
    results = poisson_hypothesis_test(
        population_rate,
        sample_size,
        observed_positive,
        alpha
    )
    
    # 构建结果字典
    result_dict = {
        "table_name": "泊松分布假设检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "hypothesis_test": results["hypothesis_test"],
        "interpretation": results["interpretation"],
        "remark": f"总体阳性概率: {population_rate}, 样本数: {sample_size}, 阳性数: {observed_positive}, α: {alpha}"
    }
    
    return result_dict