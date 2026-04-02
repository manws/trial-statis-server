"""
临床Hartley方差齐性检验统计分析模块 (Clinical Hartley Variance Homogeneity Test Statistics Module)

本模块提供全面的Hartley方差齐性检验统计分析功能，用于临床数据中多组方差齐性的检验，
是检验多组数据方差是否相等的经典统计方法。Hartley检验特别适用于等样本量的多组数据，
通过比较最大方差与最小方差的比值来判断方差齐性。

【模块功能概述】:
1. Hartley统计量计算：计算多组数据方差齐性的Hartley统计量
2. 方差比计算：计算最大方差与最小方差的比值
3. 显著性检验：执行检验判断方差是否齐性
4. 等样本量验证：验证所有组样本量是否相等
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 方差齐性检验：验证多组数据方差是否齐性
- 统计方法选择：为后续统计分析方法的选择提供依据
- 数据质量评估：评估多组间变异程度的差异
- 假设检验验证：验证方差分析等方法的前提条件

【统计方法选择指南】:
1. Hartley检验适用条件：
   - 各组数据相互独立
   - 所有组的样本量必须相等
   - 适用于3组及以上数据的方差比较
   - 对数据正态性假设敏感
   - 相比其他方法，计算简单直观

2. 临床应用场景：
   - 比较多个治疗组某项指标的变异性差异（等样本量情况）
   - 验证方差分析前的方差齐性假设（等样本量情况）
   - 检验不同处理组间变异程度的差异（等样本量情况）

【结果解读注意事项】:
1. Hartley统计量解释：统计量越接近1表示方差越齐，值越大表示方差差异越显著
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
- "我的多组等样本量数据适合用Hartley检验比较方差吗？"
- "如何解释Hartley检验的结果？"
- "P值小于0.05意味着什么？"
- "Hartley检验的前提条件是什么？"
- "Hartley检验和其他方差齐性检验有何区别？"

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

def calculate_hartley_test_variance_homogeneity(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    执行Hartley方差齐性检验
    
    在临床研究中，这用于检验多组数据的方差是否齐性，
    即比较多个组数据的变异程度是否存在显著差异。
    Hartley检验特别适用于等样本量的多组数据，
    通过比较最大方差与最小方差的比值来判断方差齐性。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据，所有组的样本量必须相等
        
    Returns:
        Dict[str, Any]: 包含Hartley检验统计量和结果的字典
        
    Raises:
        ValueError: 当数据组数少于2组、样本量不相等或数据无效时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) < 2:
        raise ValueError("Hartley检验至少需要2组数据")
    
    # 提取数据和基本信息
    groups_data = []
    groups_info = []
    
    for i, data in enumerate(data_list):
        data_array = np.array(data)
        n = len(data_array)
        var = np.var(data_array, ddof=1)  # 样本方差
        mean = np.mean(data_array)
        std = np.std(data_array, ddof=1)
        
        if n < 2:
            raise ValueError(f"第{i+1}组数据样本量不足，至少需要2个观测值")
        
        groups_data.append(data_array)
        groups_info.append({
            "name": f"第{i+1}组",
            "sample_size": n,
            "mean": float(mean),
            "std": float(std),
            "variance": float(var)
        })
    
    k = len(groups_data)  # 组数
    
    # 检查样本量是否相等（Hartley检验要求等样本量）
    sample_sizes = [info["sample_size"] for info in groups_info]
    if len(set(sample_sizes)) > 1:
        raise ValueError("Hartley检验要求所有组的样本量相等")
    
    n_per_group = sample_sizes[0]  # 每组样本量
    
    # 提取各组方差
    variances = [info["variance"] for info in groups_info]
    
    # 计算最大和最小方差
    max_variance = max(variances)
    min_variance = min(variances)
    
    # Hartley统计量：F_max = max(s²) / min(s²)
    if min_variance > 0:
        hartley_statistic = max_variance / min_variance
    else:
        raise ValueError("存在方差为0的情况，无法计算Hartley统计量")
    
    # 计算临界值和P值
    # Hartley检验使用F-max分布
    df = n_per_group - 1  # 每组的自由度
    
    # 计算P值（基于F-max分布的近似）
    # 这里使用保守的近似方法
    if df > 0:
        # 使用F分布近似计算P值
        # F-max统计量近似服从F(df, df)分布
        p_value = 1 - stats.f.cdf(hartley_statistic, df, df)
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
            "sample_size_per_group": int(n_per_group),
            "total_sample_size": int(k * n_per_group),
            "groups_info": groups_info,
            "variances": variances,
            "max_variance": float(max_variance),
            "min_variance": float(min_variance)
        },
        "test_statistics": {
            "hartley_statistic": float(hartley_statistic),
            "degrees_of_freedom": int(df),
            "variance_ratio": float(hartley_statistic)
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
            "statistic_interpretation": f"Hartley统计量为{hartley_statistic:.4f}，自由度为{df}",
            "p_value_interpretation": f"P值为{p_value:.6f}",
            "significance_interpretation": _get_significance_interpretation(p_value),
            "homogeneity_assessment": _get_homogeneity_assessment(p_value),
            "assumption_note": "注意：Hartley检验要求所有组样本量相等，且对正态性假设敏感"
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


def cal_result_hv_hartley(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    生成Hartley方差齐性检验统计分析的完整报告字典
    
    此函数整合了Hartley方差齐性检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的Hartley检验结果。报告包括输入参数、
    检验统计量、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解Hartley检验的特征。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据，所有组的样本量必须相等
    
    Returns:
        Dict[str, Any]: 包含Hartley方差齐性检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 执行Hartley检验
    results = calculate_hartley_test_variance_homogeneity(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Hartley方差齐性检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"组数: {results['input_parameters']['group_count']}, 每组样本量: {results['input_parameters']['sample_size_per_group']}"
    }
    
    return result_dict