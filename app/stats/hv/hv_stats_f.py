"""
临床方差齐性F检验统计分析模块 (Clinical Variance Homogeneity F-Test Statistics Module)

本模块提供全面的方差齐性F检验统计分析功能，用于临床数据中方差齐性的检验，
是比较两组数据方差是否相等的经典统计方法。F检验在医学研究中广泛应用，
如比较两组患者的方差是否存在显著差异、验证方差齐性假设等。

【模块功能概述】:
1. F统计量计算：计算两组数据方差比的F统计量
2. 自由度计算：计算分子和分母的自由度
3. 显著性检验：执行双侧F检验判断方差是否齐性
4. 方差比计算：计算两组方差的比值
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 方差齐性检验：验证不同组别间方差是否齐性
- 统计方法选择：为后续统计分析方法的选择提供依据
- 数据质量评估：评估不同组别间变异程度的差异
- 假设检验验证：验证统计方法应用的前提条件

【统计方法选择指南】:
1. F检验适用条件：
   - 两组数据均服从正态分布
   - 数据相互独立
   - 适用于小样本情况下的方差比较

2. 临床应用场景：
   - 比较两组患者某项指标的变异性差异
   - 验证方差分析前的方差齐性假设
   - 检验不同处理组间变异程度的差异

【结果解读注意事项】:
1. F统计量解释：衡量两组方差的比值，F值越接近1表示方差越齐
2. P值解释：在零假设（两组方差相等）成立的情况下，观察到当前或更极端结果的概率
3. 显著性水平：通常使用α=0.05作为判断方差齐性的标准
4. 临床意义：方差齐性是许多统计方法的前提条件，不满足时需选用其他方法

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的两组数据适合用F检验比较方差吗？"
- "如何解释F检验的结果？"
- "P值小于0.05意味着什么？"
- "F检验的前提条件是什么？"
- "如何判断方差是否齐性？"

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
from app.schemas.request_data.hv_param import HVParamF

def calculate_f_test_variance_homogeneity(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    执行方差齐性F检验（方差比检验）
    
    在临床研究中，这用于检验两组数据的方差是否齐性，
    即比较两组数据的变异程度是否存在显著差异。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据
        
    Returns:
        Dict[str, Any]: 包含F检验统计量和结果的字典
        
    Raises:
        ValueError: 当数据组数不是2组或数据无效时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) != 2:
        raise ValueError("F检验仅适用于两组数据的方差比较")
    
    # 提取数据
    group1_data = np.array(data_list[0])
    group2_data = np.array(data_list[1])
    
    # 基本统计量
    n1 = len(group1_data)
    n2 = len(group2_data)
    var1 = np.var(group1_data, ddof=1)  # 样本方差
    var2 = np.var(group2_data, ddof=1)  # 样本方差
    mean1 = np.mean(group1_data)
    mean2 = np.mean(group2_data)
    std1 = np.std(group1_data, ddof=1)
    std2 = np.std(group2_data, ddof=1)
    
    # 自由度
    df1 = n1 - 1
    df2 = n2 - 1
    
    # 方差为零时的特殊处理
    if var1 == 0 and var2 == 0:
        # 两组方差都为零，完全齐性
        f_statistic = 1.0
        numerator_df = df1
        denominator_df = df2
        larger_var_group = "第一组"
        smaller_var_group = "第二组"
    elif var2 == 0 and var1 > 0:
        raise ValueError("第二组方差为0（所有值相同），无法计算F统计量")
    elif var1 == 0 and var2 > 0:
        raise ValueError("第一组方差为0（所有值相同），无法计算F统计量")
    elif var1 >= var2:
        f_statistic = var1 / var2
        numerator_df = df1
        denominator_df = df2
        larger_var_group = "第一组"
        smaller_var_group = "第二组"
    else:
        f_statistic = var2 / var1
        numerator_df = df2
        denominator_df = df1
        larger_var_group = "第二组"
        smaller_var_group = "第一组"
    
    # 计算P值（双侧检验）
    # F检验通常是双侧的，但我们使用scipy的单侧检验然后乘以2
    p_value = 2 * (1 - stats.f.cdf(f_statistic, numerator_df, denominator_df))
    # 确保P值不超过1
    p_value = min(p_value, 1.0)
    
    # 显著性判断（转换为Python原生布尔值）
    significant_01 = bool(p_value < 0.01)
    significant_05 = bool(p_value < 0.05)
    significant_10 = bool(p_value < 0.10)
    
    # 方差比
    variance_ratio = var1 / var2 if var1 >= var2 else var2 / var1
    
    return {
        "input_parameters": {
            "group1_sample_size": int(n1),
            "group2_sample_size": int(n2),
            "group1_mean": float(mean1),
            "group2_mean": float(mean2),
            "group1_std": float(std1),
            "group2_std": float(std2),
            "group1_variance": float(var1),
            "group2_variance": float(var2)
        },
        "test_statistics": {
            "f_statistic": float(f_statistic),
            "numerator_df": int(numerator_df),
            "denominator_df": int(denominator_df),
            "variance_ratio": float(variance_ratio),
            "larger_variance_group": larger_var_group,
            "smaller_variance_group": smaller_var_group
        },
        "significance_tests": {
            "p_value": float(p_value),
            "significant_at_01": significant_01,
            "significant_at_05": significant_05,
            "significant_at_10": significant_10,
            "p_value_lt_01": bool(p_value < 0.01),
            "p_value_lt_05": bool(p_value < 0.05),
            "p_value_lt_10": bool(p_value < 0.10)
        },
        "interpretation": {
            "f_statistic_interpretation": f"F统计量为{f_statistic:.4f}，自由度为({numerator_df}, {denominator_df})",
            "p_value_interpretation": f"双侧P值为{p_value:.6f}",
            "variance_comparison": f"{larger_var_group}的方差({max(var1, var2):.4f})是{smaller_var_group}方差({min(var1, var2):.4f})的{variance_ratio:.4f}倍",
            "significance_interpretation": _get_significance_interpretation(p_value),
            "homogeneity_assessment": _get_homogeneity_assessment(p_value)
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


def cal_result_hv_f(param: HVParamF) -> Dict[str, Any]:
    """
    生成方差齐性F检验统计分析的完整报告字典
    
    此函数整合了方差齐性F检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的F检验结果。报告包括输入参数、
    检验统计量、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解F检验的特征。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据
    
    Returns:
        Dict[str, Any]: 包含方差齐性F检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 从参数对象解构
    data_list = [item.data_list for item in param.stats_data_list]

    # 执行F检验
    results = calculate_f_test_variance_homogeneity(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "方差齐性F检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"两组数据方差比较: 第一组样本量={results['input_parameters']['group1_sample_size']}, 第二组样本量={results['input_parameters']['group2_sample_size']}"
    }
    
    return result_dict