"""
临床两样本Wilcoxon秩和检验统计分析模块 (Clinical Two-Sample Wilcoxon Rank Sum Test Statistics Module)

本模块提供全面的两样本Wilcoxon秩和检验统计分析功能，用于临床数据中两个独立样本的比较，
是一种重要的非参数统计分析方法。两样本Wilcoxon秩和检验（也称为Mann-Whitney U检验）
通过比较两个独立样本的秩次，来判断两组之间是否存在显著差异，广泛应用于医学研究中的
两组比较、等级资料分析、非正态数据检验等领域。

【模块功能概述】:
1. 秩和计算：计算两组独立样本的秩和
2. U统计量：计算Mann-Whitney U检验统计量
3. 显著性检验：执行正态近似检验判断统计量的显著性
4. 结修正：处理数据中相同值的修正
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 两组比较：比较两个治疗组的疗效差异
- 等级资料：分析有序分类变量的差异
- 非正态数据：处理不满足正态性假设的数据
- 分布自由：对数据分布形状无特殊要求

【统计方法选择指南】:
1. Wilcoxon秩和检验适用条件：
   - 数据为连续型或有序分类变量
   - 两组数据独立
   - 分布形状相似（但不要求正态分布）
   - 观测值相互独立

2. 临床应用场景：
   - 比较两种不同治疗方案的疗效
   - 评估不同药物对某项指标的影响
   - 分析不同人群某项指标的分布差异
   - 研究不同诊断方法的结果差异

【结果解读注意事项】:
1. 检验统计量解释：U值越偏离期望值，表明两组差异越显著
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
- "我的数据适合用两样本Wilcoxon秩和检验吗？"
- "如何解释两样本Wilcoxon秩和检验的结果？"
- "两样本Wilcoxon秩和检验的前提条件是什么？"
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
from scipy.stats import norm, rankdata
from app.schemas.request_data.rs_param import RSParamOD1
import numpy as np


def wilcoxon_rank_sum_test_ordinal(group1_data: List[float], group2_data: List[float]) -> Dict:
    """
    两样本Wilcoxon秩和检验（适用于等级资料）
    
    参数:
    - group1_data: 第一组样本数据（等级资料）
    - group2_data: 第二组样本数据（等级资料）
    
    返回:
    - 包含检验统计量和结果的字典
    """
    # 参数验证
    if not group1_data or len(group1_data) == 0:
        raise ValueError("第一组数据不能为空")
    
    if not group2_data or len(group2_data) == 0:
        raise ValueError("第二组数据不能为空")
    
    n1 = len(group1_data)
    n2 = len(group2_data)
    total_n = n1 + n2
    
    # 合并所有数据并计算秩
    combined_data = group1_data + group2_data
    ranks = rankdata(combined_data, method='average')
    
    # 分离两组的秩
    group1_ranks = ranks[:n1]
    group2_ranks = ranks[n1:]
    
    # 计算秩和
    R1 = sum(group1_ranks)
    R2 = sum(group2_ranks)
    
    # 计算U统计量（Mann-Whitney U检验等价于Wilcoxon秩和检验）
    # U1 = n1*n2 + n1*(n1+1)/2 - R1
    # U2 = n1*n2 + n2*(n2+1)/2 - R2
    U1 = n1 * n2 + n1 * (n1 + 1) / 2 - R1
    U2 = n1 * n2 + n2 * (n2 + 1) / 2 - R2
    
    # 使用较小的U值作为检验统计量
    U_statistic = min(U1, U2)
    
    # 计算期望值和方差（无结时）
    E_U = n1 * n2 / 2
    Var_U = n1 * n2 * (n1 + n2 + 1) / 12
    
    # 处理等级资料中的结（相同等级值）
    # 计算结的修正因子
    unique_values, counts = np.unique(combined_data, return_counts=True)
    tie_correction = sum([count**3 - count for count in counts if count > 1])
    if tie_correction > 0:
        Var_U -= tie_correction * n1 * n2 / (12 * (total_n - 1))
    
    # 计算标准化检验统计量Z（近似Z检验）
    if Var_U > 0:
        Z_statistic = (U_statistic - E_U) / math.sqrt(Var_U)
    else:
        Z_statistic = 0.0
    
    # 计算双侧P值
    # P值 = 2 × P(Z ≥ |z|)
    p_value_two_sided = 2 * (1 - norm.cdf(abs(Z_statistic)))
    
    # 显著性判断 (< 0.05 和 < 0.01)
    is_less_than_05 = p_value_two_sided < 0.05
    is_less_than_01 = p_value_two_sided < 0.01
    
    return {
        "input_parameters": {
            "group1_size": n1,
            "group2_size": n2,
            "total_samples": total_n
        },
        "sample_statistics": {
            "group1_mean": float(np.mean(group1_data)),
            "group1_median": float(np.median(group1_data)),
            "group1_std": float(np.std(group1_data, ddof=1)),
            "group2_mean": float(np.mean(group2_data)),
            "group2_median": float(np.median(group2_data)),
            "group2_std": float(np.std(group2_data, ddof=1))
        },
        "rank_statistics": {
            "group1_ranks_sum": float(R1),
            "group2_ranks_sum": float(R2),
            "u_statistic_1": float(U1),
            "u_statistic_2": float(U2),
            "test_statistic_U": float(U_statistic),
            "expected_U": float(E_U),
            "variance_U": float(Var_U),
            "tie_correction": float(tie_correction) if tie_correction > 0 else 0.0
        },
        "test_statistics": {
            "z_value": float(Z_statistic),
            "p_value_two_sided": float(p_value_two_sided)
        },
        "significance_tests": {
            "p_less_than_0_05": is_less_than_05,
            "p_less_than_0_01": is_less_than_01,
            "significant_at_05": "显著" if is_less_than_05 else "不显著",
            "significant_at_01": "显著" if is_less_than_01 else "不显著"
        },
        "interpretation": {
            "wilcoxon_interpretation": f"两样本Wilcoxon秩和检验{'不' if p_value_two_sided >= 0.05 else ''}显著 (p={'≥' if p_value_two_sided >= 0.05 else '<'}0.05)",
            "wilcoxon_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided >= 0.01 else ''}显著 (p={'≥' if p_value_two_sided >= 0.01 else '<'}0.01)"
        }
    }


def cal_result_rs_od1(param: RSParamOD1) -> Dict:
    """
    生成两样本Wilcoxon秩和检验统计分析的完整报告字典
    
    此函数整合了两样本Wilcoxon秩和检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解两样本Wilcoxon秩和检验的特征。
    
    Args:
        group1_data: List[float]，第一组样本数据
        group2_data: List[float]，第二组样本数据
    
    Returns:
        Dict: 包含两样本Wilcoxon秩和检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - sample_statistics: 样本统计信息
            - rank_statistics: 秩统计信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 参数验证
    if not group1_data or not group2_data:
        raise ValueError("两组数据都不能为空")
    
    # 执行Wilcoxon秩和检验
    # 从参数对象解构
    group1_data = param.stats_data_list[0].data_list
    group2_data = param.stats_data_list[1].data_list

    results = wilcoxon_rank_sum_test_ordinal(group1_data, group2_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "两样本Wilcoxon秩和检验分析（等级资料）",
        "input_parameters": results["input_parameters"],
        "sample_statistics": results["sample_statistics"],
        "rank_statistics": results["rank_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"第一组样本量：{results['input_parameters']['group1_size']}, 第二组样本量：{results['input_parameters']['group2_size']}"
    }
    
    return result_dict