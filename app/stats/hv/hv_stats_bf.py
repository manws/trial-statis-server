"""
临床Brown-Forsythe方差齐性检验统计分析模块 (Clinical Brown-Forsythe Variance Homogeneity Test Statistics Module)

本模块提供全面的Brown-Forsythe方差齐性检验统计分析功能，用于临床数据中多组方差齐性的检验，
是检验多组数据方差是否相等的经典统计方法。Brown-Forsythe检验是Levene检验的一个变种，
使用整体中位数而非各组中位数，对数据分布形状不敏感，比Bartlett检验更稳健。

【模块功能概述】:
1. Brown-Forsythe统计量计算：计算多组数据方差齐性的BF统计量
2. 绝对偏差计算：计算各组数据相对于整体中位数的绝对偏差
3. 显著性检验：执行F检验判断方差是否齐性
4. 方差分解：计算组间、组内和总平方和
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 方差齐性检验：验证多组数据方差是否齐性
- 统计方法选择：为后续统计分析方法的选择提供依据
- 数据质量评估：评估多组间变异程度的差异
- 假设检验验证：验证方差分析等方法的前提条件

【统计方法选择指南】:
1. Brown-Forsythe检验适用条件：
   - 各组数据相互独立
   - 适用于2组及以上数据的方差比较
   - 对数据分布形状不敏感，比Bartlett检验更稳健
   - 适用于偏态分布或含有异常值的数据
   - 相比Levene检验，使用整体中位数而非组内中位数

2. 临床应用场景：
   - 比较多个治疗组某项指标的变异性差异
   - 验证方差分析前的方差齐性假设
   - 检验不同处理组间变异程度的差异（特别是当数据不满足正态性时）

【结果解读注意事项】:
1. F统计量解释：统计量越大表示方差差异越显著
2. P值解释：在零假设（各组方差相等）成立的情况下，观察到当前或更极端结果的概率
3. 显著性水平：通常使用α=0.05作为判断方差齐性的标准
4. 临床意义：方差齐性是方差分析等统计方法的重要前提条件

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的多组数据适合用Brown-Forsythe检验比较方差吗？"
- "如何解释Brown-Forsythe检验的结果？"
- "P值小于0.05意味着什么？"
- "Brown-Forsythe检验的前提条件是什么？"
- "Brown-Forsythe检验和Levene检验有何区别？"

【相关标准和规范】:
- CONSORT 声明：随机临床试验报告规范
- STROBE 声明：观察性研究报告规范
- ICH E9 指导原则：临床试验的统计学原则
- SAMPL 指南：医学研究报告中的统计方法描述规范

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

import numpy as np
from scipy import stats
from typing import Dict, List, Any

def calculate_bf_test_variance_homogeneity(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    执行Brown-Forsythe方差齐性检验
    
    在临床研究中，这用于检验多组数据的方差是否齐性，
    即比较多个组数据的变异程度是否存在显著差异。
    Brown-Forsythe检验是Levene检验的一个变种，使用整体中位数而非各组中位数，
    对数据分布形状不敏感，比Bartlett检验更稳健。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据
        
    Returns:
        Dict[str, Any]: 包含Brown-Forsythe检验统计量和结果的字典
        
    Raises:
        ValueError: 当数据组数少于2组或数据无效时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) < 2:
        raise ValueError("Brown-Forsythe检验至少需要2组数据")
    
    # 提取数据和基本信息
    groups_data = []
    groups_info = []
    
    for i, data in enumerate(data_list):
        data_array = np.array(data)
        n = len(data_array)
        var = np.var(data_array, ddof=1)  # 样本方差
        mean = np.mean(data_array)
        median = np.median(data_array)
        std = np.std(data_array, ddof=1)
        
        if n < 2:
            raise ValueError(f"第{i+1}组数据样本量不足，至少需要2个观测值")
        
        groups_data.append(data_array)
        groups_info.append({
            "name": f"第{i+1}组",
            "sample_size": n,
            "mean": float(mean),
            "median": float(median),
            "std": float(std),
            "variance": float(var)
        })
    
    k = len(groups_data)  # 组数
    total_n = sum(info["sample_size"] for info in groups_info)  # 总样本量
    
    # 计算所有数据的中位数（Brown-Forsythe使用整体中位数而非均值）
    all_data = np.concatenate(groups_data)
    overall_median = np.median(all_data)
    
    # 计算每组相对于整体中位数的绝对偏差
    abs_deviations = []
    for data in groups_data:
        deviations = np.abs(data - overall_median)
        abs_deviations.append(deviations)
    
    # 计算每组绝对偏差的均值
    group_means_abs_dev = []
    for deviations in abs_deviations:
        group_means_abs_dev.append(np.mean(deviations))
    
    # 执行单因素方差分析（ANOVA）on绝对偏差
    # 计算总平方和
    all_abs_deviations = np.concatenate(abs_deviations)
    grand_mean_abs_dev = np.mean(all_abs_deviations)
    
    # 总平方和 (SST)
    sst = np.sum((all_abs_deviations - grand_mean_abs_dev) ** 2)
    
    # 组间平方和 (SSB)
    ssb = 0
    for i, (deviations, group_info) in enumerate(zip(abs_deviations, groups_info)):
        group_mean = group_means_abs_dev[i]
        n_i = group_info["sample_size"]
        ssb += n_i * (group_mean - grand_mean_abs_dev) ** 2
    
    # 组内平方和 (SSW)
    ssw = sst - ssb
    
    # 自由度
    df_between = k - 1
    df_within = total_n - k
    df_total = total_n - 1
    
    # 均方
    msb = ssb / df_between if df_between > 0 else 0
    msw = ssw / df_within if df_within > 0 else 0
    
    # F统计量
    f_statistic = msb / msw if msw > 0 else 0
    
    # 计算P值
    if msw > 0 and df_between > 0 and df_within > 0:
        p_value = 1 - stats.f.cdf(f_statistic, df_between, df_within)
    else:
        p_value = 1.0
    
    # 确保P值在合理范围内
    p_value = max(min(p_value, 1.0), 0.0)
    
    # 显著性判断
    significant_01 = p_value < 0.01
    significant_05 = p_value < 0.05
    significant_10 = p_value < 0.10
    
    return {
        "input_parameters": {
            "group_count": int(k),
            "total_sample_size": int(total_n),
            "overall_median": float(overall_median),
            "groups_info": groups_info
        },
        "test_statistics": {
            "f_statistic": float(f_statistic),
            "ss_between": float(ssb),
            "ss_within": float(ssw),
            "ss_total": float(sst),
            "ms_between": float(msb),
            "ms_within": float(msw),
            "degrees_of_freedom_between": int(df_between),
            "degrees_of_freedom_within": int(df_within),
            "grand_mean_abs_dev": float(grand_mean_abs_dev)
        },
        "significance_tests": {
            "p_value": float(p_value),
            "significant_at_01": bool(significant_01),
            "significant_at_05": bool(significant_05),
            "significant_at_10": bool(significant_10),
            "p_value_lt_01": bool(p_value < 0.01),
            "p_value_lt_05": bool(p_value < 0.05),
            "p_value_lt_10": bool(p_value < 0.10)
        },
        "interpretation": {
            "statistic_interpretation": f"F统计量为{f_statistic:.4f}，组间自由度为{df_between}，组内自由度为{df_within}",
            "p_value_interpretation": f"P值为{p_value:.6f}",
            "significance_interpretation": _get_significance_interpretation(p_value),
            "homogeneity_assessment": _get_homogeneity_assessment(p_value),
            "assumption_note": "注意：Brown-Forsythe检验对数据的分布形状不敏感，比Bartlett检验更稳健"
        }
    }


def _get_significance_interpretation(p_value: float) -> str:
    """获取显著性解释"""
    if p_value < 0.01:
        return f"P值={p_value:.6f} < 0.01，方差差异极显著"
    elif p_value < 0.05:
        return f"P值={p_value:.6f} < 0.05，方差差异显著"
    elif p_value < 0.10:
        return f"P值={p_value:.6f} < 0.10，方差差异边缘显著"
    else:
        return f"P值={p_value:.6f} ≥ 0.10，方差差异不显著"


def _get_homogeneity_assessment(p_value: float) -> str:
    """获取方差齐性评估"""
    if p_value < 0.05:
        return "方差不齐（拒绝方差齐性假设）"
    else:
        return "方差齐性（不拒绝方差齐性假设）"


def cal_result_hv_bf(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    生成Brown-Forsythe方差齐性检验统计分析的完整报告字典
    
    此函数整合了Brown-Forsythe方差齐性检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的Brown-Forsythe检验结果。报告包括输入参数、
    检验统计量、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解Brown-Forsythe检验的特征。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据
    
    Returns:
        Dict[str, Any]: 包含Brown-Forsythe方差齐性检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 执行Brown-Forsythe检验
    results = calculate_bf_test_variance_homogeneity(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Brown-Forsythe方差齐性检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"组数: {results['input_parameters']['group_count']}, 总样本量: {results['input_parameters']['total_sample_size']}, 整体中位数: {results['input_parameters']['overall_median']:.4f}"
    }
    
    return result_dict