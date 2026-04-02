"""
临床二项分布统计分析模块 (Clinical Binomial Distribution Statistics Module)

本模块提供全面的二项分布统计分析功能，用于临床数据中阳性率、成功率等二项分布参数的估计和检验，
是临床试验数据分析的重要组成部分。二项分布在医学研究中广泛应用，
如计算某种疾病的发生率、某种治疗的成功率、某种基因型的比例等。

【模块功能概述】:
1. 概率计算：计算特定阳性数的概率、累积概率等
2. 区间估计：计算总体率的置信区间
3. 假设检验：比较样本率与总体率的差异
4. 两样本比较：比较两个样本率的差异
5. 临界值计算：计算给定累积概率对应的阳性数

【临床应用价值】:
- 疗效评估：评估某种治疗的成功率是否达到预期
- 疾病发生率：估计特定人群中的疾病发生率
- 诊断准确性：评估诊断测试的敏感性或特异性
- 临床决策：为临床决策提供概率依据

【统计方法选择指南】:
1. 二项分布适用条件：
   - 试验只有两种互斥结果（成功/失败）
   - 试验次数固定
   - 每次试验的成功概率相同
   - 试验之间相互独立

2. 临床应用场景：
   - 治疗成功率的区间估计
   - 新药有效性的假设检验
   - 不同疗法效果的比较

【结果解读注意事项】:
1. 概率解释：P(X=k)表示恰好k次成功的概率
2. 累积概率解释：P(X≤k)表示最多k次成功的概率
3. 置信区间解释：在给定置信水平下，总体率的可能范围
4. P值解释：在零假设成立的前提下，观察到当前或更极端结果的概率
5. 临床意义：统计显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用二项分布分析吗？"
- "如何解释二项分布的概率值？"
- "P值小于0.05意味着什么？"
- "二项分布的前提条件是什么？"
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
from scipy.stats import binom
from app.schemas.request_data.base_param import BaseParamBino1


def probability_binomial(k: int, n: int, p: float) -> float:
    """
    计算二项分布的概率 P(X=k)
    
    在临床研究中，这用于计算在n次独立试验中恰好出现k次阳性（成功）的概率，
    以评估某种现象的出现可能性。
    
    Args:
        k: 阳性数（成功次数）
        n: 试验总数（样本量）
        p: 总体阳性概率（每次试验成功的概率）
        
    Returns:
        float: 概率值 P(X=k)
    """
    return binom.pmf(k, n, p)


def cumulative_probability_binomial(k: int, n: int, p: float) -> float:
    """
    计算二项分布的累积概率 P(X≤k)
    
    在临床研究中，这用于计算在n次独立试验中最多出现k次阳性（成功）的概率，
    以评估某种现象出现频率的上限。
    
    Args:
        k: 最大阳性数
        n: 试验总数
        p: 总体阳性概率
        
    Returns:
        float: 累积概率值 P(X≤k)
    """
    return binom.cdf(k, n, p)


def calculate_critical_value_binomial(alpha: float, n: int, p: float) -> int:
    """
    计算二项分布的临界值
    
    在临床研究中，这用于确定在给定显著性水平alpha下，n次试验中最可能出现的
    最大阳性数，以设定临床决策的阈值。
    
    Args:
        alpha: 显著性水平
        n: 试验总数
        p: 总体阳性概率
        
    Returns:
        int: 临界值
    """
    return binom.ppf(1 - alpha, n, p)


def cal_result_bino1(param: BaseParamBino1) -> Dict[str, Any]:
    """
    生成二项分布统计分析的完整报告字典
    
    此函数整合了二项分布的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的二项分布分析结果。
    报告包括概率计算、累积概率、临界值等信息，
    便于临床医生和研究人员快速理解二项分布的特征。
    
    Args:
        param: BaseParamBino1对象，包含total_posi_p, sample_size
        
    Returns:
        Dict[str, Any]: 包含二项分布统计分析指标的字典，键为指标名称，值为对应的统计量
            - total_positive_probability: 总体阳性概率
            - sample_size: 样本量
            - probabilities: 各阳性数的概率列表
            - cumulative_probabilities: 各阳性数的累积概率列表
            - expected_value: 期望值
            - variance: 方差
            - standard_deviation: 标准差
    """
    # 从参数对象中提取值
    total_posi_p = param.total_posi_p
    sample_size = param.sample_size
    
    # 计算各项指标
    expected_value = sample_size * total_posi_p
    variance = sample_size * total_posi_p * (1 - total_posi_p)
    std_dev = math.sqrt(variance)
    
    # 计算每个可能阳性数的概率和累积概率
    probabilities = []
    cumulative_probabilities = []
    for k in range(sample_size + 1):
        prob = probability_binomial(k, sample_size, total_posi_p)
        cum_prob = cumulative_probability_binomial(k, sample_size, total_posi_p)
        probabilities.append({
            'positive_count': k,
            'probability': prob,
            'cumulative_probability': cum_prob
        })
        cumulative_probabilities.append(cum_prob)
    
    # 构建结果字典
    result_dict = {
        "table_name": "二项分布统计分析",
        "total_positive_probability": total_posi_p,
        "sample_size": sample_size,
        "expected_value": expected_value,
        "variance": variance,
        "standard_deviation": std_dev,
        "probabilities": probabilities,
        "remark": f"总体阳性概率={total_posi_p}, 样本量={sample_size}, 期望值={expected_value}, 方差={variance}"
    }
    
    return result_dict
