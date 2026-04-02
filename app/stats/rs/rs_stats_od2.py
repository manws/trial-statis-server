"""
临床多样本 Kruskal-Wallis H 检验统计分析模块 (Clinical Multi-Sample Kruskal-Wallis H Test Statistics Module)

本模块提供全面的多样本 Kruskal-Wallis H 检验统计分析功能，用于临床数据中多个独立样本的比较，
是一种重要的非参数统计分析方法。多样本 Kruskal-Wallis H 检验通过比较多个独立样本的秩次，
来判断各组之间是否存在显著差异，广泛应用于医学研究中的多组比较、
等级资料分析、非正态数据检验等领域。

【模块功能概述】:
1. 秩和计算：计算多个独立样本的秩和
2. 检验统计量：计算 Kruskal-Wallis 检验统计量 H
3. 显著性检验：执行卡方近似检验判断统计量的显著性
4. 结修正：处理数据中相同值的修正
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 多组比较：比较多个治疗组的疗效差异
- 等级资料：分析有序分类变量的差异
- 非正态数据：处理不满足正态性假设的数据
- 分布自由：对数据分布形状无特殊要求

【统计方法选择指南】:
1. Kruskal-Wallis 检验适用条件：
   - 数据为连续型或有序分类变量
   - 各组数据独立
   - 分布形状相似（但不要求正态分布）
   - 至少有两组数据

2. 临床应用场景：
   - 比较三种或以上不同治疗方案的疗效
   - 评估多个医院间患者满意度的差异
   - 分析不同年龄段某项指标的分布差异
   - 研究多种药物对某项生物标志物的影响

【结果解读注意事项】:
1. 检验统计量解释：H 值越大，表明组间差异越显著
2. P 值解释：P 值小于显著性水平（如 0.05）时拒绝原假设
3. 临床意义：统计显著性不等于临床意义，需结合专业知识进行解释

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用 Kruskal-Wallis 检验吗？"
- "如何解释 Kruskal-Wallis 检验的结果？"
- "Kruskal-Wallis 检验的前提条件是什么？"
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
from typing import Dict, List, Any
from scipy.stats import chi2, rankdata
import numpy as np
from app.schemas.request_data.rs_param import RSParamOD2


def kruskal_wallis_h_test(groups_data: List[List[float]]) -> Dict:
    """
    多样本Kruskal-Wallis H秩和检验（适用于等级资料）
    
    参数:
    - groups_data: 多组样本数据列表，每组为一个List[float]
    
    返回:
    - 包含检验统计量和结果的字典
    """
    # 参数验证
    if not groups_data or len(groups_data) < 2:
        raise ValueError("至少需要两组数据进行Kruskal-Wallis检验")
    
    # 验证每组数据
    for i, group_data in enumerate(groups_data):
        if not group_data or len(group_data) == 0:
            raise ValueError(f"第{i+1}组数据不能为空")
    
    k = len(groups_data)  # 组数
    n_total = sum(len(group) for group in groups_data)  # 总样本量
    
    # 合并所有数据并计算秩
    combined_data = []
    group_indices = []  # 记录每个数据点属于哪一组
    
    for group_idx, group_data in enumerate(groups_data):
        combined_data.extend(group_data)
        group_indices.extend([group_idx] * len(group_data))
    
    ranks = rankdata(combined_data, method='average')
    
    # 计算每组的秩和
    group_rank_sums = [0.0] * k
    group_sizes = [0] * k
    
    for i, group_idx in enumerate(group_indices):
        group_rank_sums[group_idx] += ranks[i]
        group_sizes[group_idx] += 1
    
    # 计算Kruskal-Wallis H统计量
    # H = (12/(N(N+1))) * Σ(Ri²/ni) - 3(N+1)
    # 其中 Ri 是第i组的秩和，ni 是第i组的样本量，N 是总样本量
    
    h_statistic = 0.0
    for i in range(k):
        if group_sizes[i] > 0:
            h_statistic += (group_rank_sums[i] ** 2) / group_sizes[i]
    
    h_statistic = (12 / (n_total * (n_total + 1))) * h_statistic - 3 * (n_total + 1)
    
    # 自由度
    df = k - 1
    
    # 计算P值（近似卡方检验）
    if h_statistic >= 0:
        p_value = 1 - chi2.cdf(h_statistic, df)
    else:
        p_value = 1.0  # H统计量不应为负，如果出现则设为最大P值
    
    # 处理结（ties）的修正
    # 计算结修正因子
    unique_values, counts = np.unique(combined_data, return_counts=True)
    tie_correction = sum([count**3 - count for count in counts if count > 1])
    
    if tie_correction > 0 and n_total > 1:
        # 应用结修正
        correction_factor = 1 - tie_correction / (n_total**3 - n_total)
        if correction_factor > 0:
            corrected_h = h_statistic / correction_factor
            corrected_p_value = 1 - chi2.cdf(corrected_h, df)
        else:
            corrected_h = h_statistic
            corrected_p_value = p_value
    else:
        corrected_h = h_statistic
        corrected_p_value = p_value
    
    # 显著性判断 (< 0.05 和 < 0.01)
    is_less_than_05 = corrected_p_value < 0.05
    is_less_than_01 = corrected_p_value < 0.01
    
    # 计算每组的基本统计量
    group_stats = []
    for i, group_data in enumerate(groups_data):
        group_stats.append({
            "group_index": i + 1,
            "sample_size": len(group_data),
            "mean": float(np.mean(group_data)),
            "median": float(np.median(group_data)),
            "std": float(np.std(group_data, ddof=1)) if len(group_data) > 1 else 0.0,
            "rank_sum": float(group_rank_sums[i])
        })
    
    return {
        "input_parameters": {
            "number_of_groups": k,
            "total_samples": n_total,
            "group_sizes": group_sizes
        },
        "group_statistics": group_stats,
        "rank_statistics": {
            "group_rank_sums": [float(rs) for rs in group_rank_sums],
            "h_statistic_original": float(h_statistic),
            "h_statistic_corrected": float(corrected_h),
            "tie_correction": float(tie_correction) if tie_correction > 0 else 0.0
        },
        "test_statistics": {
            "chi_square_value": float(corrected_h),
            "degrees_of_freedom": df,
            "p_value": float(corrected_p_value)
        },
        "significance_tests": {
            "p_less_than_0_05": is_less_than_05,
            "p_less_than_0_01": is_less_than_01,
            "significant_at_05": "显著" if is_less_than_05 else "不显著",
            "significant_at_01": "显著" if is_less_than_01 else "不显著"
        },
        "interpretation": {
            "kruskal_wallis_interpretation": f"Kruskal-Wallis H检验{'不' if corrected_p_value >= 0.05 else ''}显著 (p={'≥' if corrected_p_value >= 0.05 else '<'}0.05)",
            "kruskal_wallis_interpretation_01": f"在0.01水平下{'不' if corrected_p_value >= 0.01 else ''}显著 (p={'≥' if corrected_p_value >= 0.01 else '<'}0.01)"
        }
    }


def perform_kruskal_wallis_h_test(groups_data: List[List[float]]) -> Dict:
    """
    根据给定参数执行多样本Kruskal-Wallis H检验
    
    参数:
    - groups_data: 多组样本数据列表，每组为一个List[float]
    
    返回:
    - 包含Kruskal-Wallis H检验完整结果的字典
    """
    # 验证参数
    if not groups_data or len(groups_data) < 2:
        raise ValueError("至少需要两组数据进行Kruskal-Wallis检验")
    
    # 验证每组数据
    for i, group_data in enumerate(groups_data):
        if not group_data or len(group_data) == 0:
            raise ValueError(f"第{i+1}组数据不能为空")
    
    # 执行Kruskal-Wallis H检验
    results = kruskal_wallis_h_test(groups_data)
    
    return results


def cal_result_rs_od2(param: RSParamOD2) -> Dict[str, Any]:
    """
    生成多样本Kruskal-Wallis H检验分析的完整报告字典
    
    此函数整合了多样本Kruskal-Wallis H检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解Kruskal-Wallis H检验的特征。
    
    Args:
        groups_data: List[List[float]]，多组样本数据列表，每组为一个List[float]
    
    Returns:
        Dict[str, Any]: 包含多样本Kruskal-Wallis H检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"多样本Kruskal-Wallis H检验分析"
            - input_parameters: 输入参数信息
            - group_statistics: 组统计信息
            - rank_statistics: 秩统计信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 执行Kruskal-Wallis H检验
    # 从参数对象解构
    groups_data = [item.data_list for item in param.stats_data_list]

    results = perform_kruskal_wallis_h_test(groups_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "多样本Kruskal-Wallis H检验分析",
        "input_parameters": results["input_parameters"],
        "group_statistics": results["group_statistics"],
        "rank_statistics": results["rank_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"组数: {results['input_parameters']['number_of_groups']}, 总样本量: {results['input_parameters']['total_samples']}"
    }
    
    return result_dict