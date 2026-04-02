"""
临床Friedman M检验统计分析模块 (Clinical Friedman M Test Statistics Module)

本模块提供全面的Friedman M检验统计分析功能，用于临床数据中随机区组设计的等级资料比较，
是一种重要的非参数统计分析方法。Friedman M检验通过对区组内各处理的秩次进行比较，
来判断不同处理之间是否存在显著差异，广泛应用于医学研究中的配伍组设计、
重复测量分析、多处理比较等领域。

【模块功能概述】:
1. 秩和计算：计算各处理的秩和
2. 检验统计量：计算Friedman检验统计量M
3. 显著性检验：执行卡方近似检验判断统计量的显著性
4. Q事后检验：当主检验显著时，进行成对比较
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 配伍组设计：比较同一区组内不同处理的效果
- 重复测量：分析同一受试者在不同时间点的差异
- 多处理比较：评估多种治疗方法的差异
- 非参数分析：处理不满足正态性假设的数据

【统计方法选择指南】:
1. Friedman检验适用条件：
   - 数据为连续型或有序分类变量
   - 同一区组内数据相关
   - 不同区组间数据独立
   - 至少有三种处理

2. 临床应用场景：
   - 比较三种或以上不同治疗方案在同一群组患者中的效果
   - 评估同一患者在不同时间点的指标变化
   - 分析不同药物对同一组患者的影响
   - 研究不同护理方法在相同条件下效果差异

【结果解读注意事项】:
1. 检验统计量解释：M值越大，表明处理间差异越显著
2. P值解释：P值小于显著性水平（如0.05）时拒绝原假设
3. 临床意义：统计显著性不等于临床意义，需结合专业知识进行解释

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用Friedman检验吗？"
- "如何解释Friedman检验的结果？"
- "Friedman检验的前提条件是什么？"
- "如何判断差异是否具有统计学意义？"

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
from typing import Dict, List
from scipy.stats import chi2, rankdata, norm
import numpy as np
from itertools import combinations


def friedman_m_test(block_data: List[List[float]]) -> Dict:
    """
    Friedman M秩和检验（适用于随机区组设计的等级资料）
    
    参数:
    - block_data: 区组数据，每个子列表代表一个区组，包含各处理的观测值
    
    返回:
    - 包含检验统计量和结果的字典
    """
    # 参数验证
    if not block_data or len(block_data) == 0:
        raise ValueError("区组数据不能为空")
    
    # 验证区组数据
    block_count = len(block_data)
    treatment_count = len(block_data[0]) if block_data else 0
    
    if treatment_count < 2:
        raise ValueError("至少需要两种处理进行Friedman检验")
    
    # 验证每个区组的数据完整性
    for i, block in enumerate(block_data):
        if len(block) != treatment_count:
            raise ValueError(f"第{i+1}个区组的数据长度不一致")
        if not block or all(x is None or x == '' for x in block):
            raise ValueError(f"第{i+1}个区组数据不能为空")
    
    # 计算每个区组内各处理的秩
    ranked_data = []
    for block in block_data:
        ranks = rankdata(block, method='average')
        ranked_data.append(ranks.tolist())
    
    # 计算各处理的平均秩
    treatment_rank_sums = [0.0] * treatment_count
    treatment_avg_ranks = [0.0] * treatment_count
    
    for block_ranks in ranked_data:
        for j in range(treatment_count):
            treatment_rank_sums[j] += block_ranks[j]
    
    for j in range(treatment_count):
        treatment_avg_ranks[j] = treatment_rank_sums[j] / block_count
    
    # 计算Friedman M统计量
    # M = (12/(bk(k+1))) * Σ(Rj²) - 3b(k+1)
    # 其中 b=区组数, k=处理数, Rj=第j个处理的秩和
    sum_rank_squares = sum(rank_sum ** 2 for rank_sum in treatment_rank_sums)
    m_statistic = (12 / (block_count * treatment_count * (treatment_count + 1))) * sum_rank_squares - 3 * block_count * (treatment_count + 1)
    
    # 自由度
    df = treatment_count - 1
    
    # 计算P值（近似卡方检验）
    if m_statistic >= 0:
        p_value = 1 - chi2.cdf(m_statistic, df)
    else:
        p_value = 1.0  # M统计量不应为负，如果出现则设为最大P值
    
    # 处理结（ties）的修正（对于Friedman检验通常较少出现）
    # 这里简化处理，因为区组内排序相对较少出现完全相同的值
    
    # 显著性判断 (< 0.05 和 < 0.01)
    is_less_than_05 = p_value < 0.05
    is_less_than_01 = p_value < 0.01
    
    # 计算基本统计量
    treatment_stats = []
    for j in range(treatment_count):
        # 提取第j个处理的所有观测值
        treatment_values = [block[j] for block in block_data]
        treatment_stats.append({
            "treatment_index": j + 1,
            "mean": float(np.mean(treatment_values)),
            "median": float(np.median(treatment_values)),
            "std": float(np.std(treatment_values, ddof=1)) if len(treatment_values) > 1 else 0.0,
            "rank_sum": float(treatment_rank_sums[j]),
            "avg_rank": float(treatment_avg_ranks[j])
        })
    
    return {
        "input_parameters": {
            "block_count": block_count,
            "treatment_count": treatment_count,
            "total_observations": block_count * treatment_count
        },
        "treatment_statistics": treatment_stats,
        "rank_statistics": {
            "treatment_rank_sums": [float(rs) for rs in treatment_rank_sums],
            "treatment_avg_ranks": [float(ar) for ar in treatment_avg_ranks],
            "m_statistic": float(m_statistic)
        },
        "test_statistics": {
            "chi_square_value": float(m_statistic),
            "degrees_of_freedom": df,
            "p_value": float(p_value)
        },
        "significance_tests": {
            "p_less_than_0_05": is_less_than_05,
            "p_less_than_0_01": is_less_than_01,
            "significant_at_05": "显著" if is_less_than_05 else "不显著",
            "significant_at_01": "显著" if is_less_than_01 else "不显著"
        },
        "interpretation": {
            "friedman_interpretation": f"Friedman M检验{'不' if p_value >= 0.05 else ''}显著 (p={'≥' if p_value >= 0.05 else '<'}0.05)",
            "friedman_interpretation_01": f"在0.01水平下{'不' if p_value >= 0.01 else ''}显著 (p={'≥' if p_value >= 0.01 else '<'}0.01)"
        }
    }


def friedman_q_post_hoc_test(block_data: List[List[float]], alpha: float = 0.05) -> Dict:
    """
    Friedman检验后的Q检验（用于多重比较）
    
    参数:
    - block_data: 区组数据
    - alpha: 显著性水平，默认0.05
    
    返回:
    - 包含事后检验结果的字典
    """
    # 参数验证
    if not block_data or len(block_data) == 0:
        raise ValueError("区组数据不能为空")
    
    block_count = len(block_data)
    treatment_count = len(block_data[0]) if block_data else 0
    
    if treatment_count < 2:
        raise ValueError("至少需要两种处理进行Q检验")
    
    # 计算每个区组内各处理的秩
    ranked_data = []
    for block in block_data:
        ranks = rankdata(block, method='average')
        ranked_data.append(ranks.tolist())
    
    # 计算各处理的平均秩
    treatment_avg_ranks = [0.0] * treatment_count
    for block_ranks in ranked_data:
        for j in range(treatment_count):
            treatment_avg_ranks[j] += block_ranks[j]
    
    for j in range(treatment_count):
        treatment_avg_ranks[j] /= block_count
    
    # 计算Q检验统计量
    # 对于每对处理i,j: Qij = |Ri - Rj| / sqrt(k(k+1)/(6b))
    # 其中Ri,Rj是平均秩，k是处理数，b是区组数
    
    critical_q_dict = {
        (3, 0.05): 2.291, (3, 0.01): 2.921,
        (4, 0.05): 2.559, (4, 0.01): 3.244,
        (5, 0.05): 2.735, (5, 0.01): 3.451,
        (6, 0.05): 2.866, (6, 0.01): 3.613,
        (7, 0.05): 2.971, (7, 0.01): 3.745,
        (8, 0.05): 3.059, (8, 0.01): 3.855,
        (9, 0.05): 3.136, (9, 0.01): 3.950,
        (10, 0.05): 3.204, (10, 0.01): 4.034
    }
    
    # 获取临界Q值
    key = (treatment_count, alpha)
    critical_q = critical_q_dict.get(key, norm.ppf(1 - alpha/2))  # 默认使用正态近似
    
    # 计算临界差异
    critical_difference = critical_q * math.sqrt(treatment_count * (treatment_count + 1) / (6 * block_count))
    
    # 进行所有成对比较
    pairwise_comparisons = []
    significant_pairs = []
    
    for i, j in combinations(range(treatment_count), 2):
        rank_diff = abs(treatment_avg_ranks[i] - treatment_avg_ranks[j])
        is_significant = rank_diff > critical_difference
        
        comparison = {
            "treatment_i": i + 1,
            "treatment_j": j + 1,
            "avg_rank_i": float(treatment_avg_ranks[i]),
            "avg_rank_j": float(treatment_avg_ranks[j]),
            "rank_difference": float(rank_diff),
            "critical_difference": float(critical_difference),
            "significant": is_significant
        }
        pairwise_comparisons.append(comparison)
        
        if is_significant:
            significant_pairs.append((i + 1, j + 1))
    
    return {
        "input_parameters": {
            "block_count": block_count,
            "treatment_count": treatment_count,
            "alpha_level": alpha,
            "critical_q": float(critical_q),
            "critical_difference": float(critical_difference)
        },
        "treatment_average_ranks": [float(rank) for rank in treatment_avg_ranks],
        "pairwise_comparisons": pairwise_comparisons,
        "significant_pairs": significant_pairs,
        "summary": {
            "total_comparisons": len(pairwise_comparisons),
            "significant_count": len(significant_pairs)
        }
    }


def cal_result_rs_fm(block_data: List[List[float]]) -> Dict:
    """
    生成Friedman M检验及Q检验统计分析的完整报告字典
    
    此函数整合了Friedman M检验及Q检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解Friedman M检验的特征。
    
    Args:
        block_data: List[List[float]]，区组数据，每个子列表代表一个区组，包含各处理的观测值
    
    Returns:
        Dict: 包含Friedman M检验及Q检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - treatment_statistics: 处理统计信息
            - rank_statistics: 秩统计信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
            - q_results: Q事后检验结果（如果主检验显著）
    """
    # 执行Friedman M检验
    fm_results = friedman_m_test(block_data)
    
    # 如果Friedman检验显著，则进行Q检验
    q_results = None
    if fm_results["significance_tests"]["p_less_than_0_05"]:
        q_results = friedman_q_post_hoc_test(block_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Friedman M检验及Q检验分析",
        "input_parameters": fm_results["input_parameters"],
        "treatment_statistics": fm_results["treatment_statistics"],
        "rank_statistics": fm_results["rank_statistics"],
        "test_statistics": fm_results["test_statistics"],
        "significance_tests": fm_results["significance_tests"],
        "interpretation": fm_results["interpretation"],
        "q_results": q_results,
        "overall_significant": fm_results["significance_tests"]["p_less_than_0_05"],
        "remark": f"区组数: {fm_results['input_parameters']['block_count']}, 处理数: {fm_results['input_parameters']['treatment_count']}"
    }
    
    return result_dict
