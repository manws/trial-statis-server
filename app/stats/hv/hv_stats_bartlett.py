"""
临床Bartlett方差齐性检验统计分析模块 (Clinical Bartlett Variance Homogeneity Test Statistics Module)

本模块提供全面的Bartlett方差齐性检验统计分析功能，用于临床数据中多组方差齐性的检验，
是检验多组数据方差是否相等的经典统计方法。Bartlett检验在医学研究中广泛应用，
如比较多个治疗组方差的齐性、验证方差分析的前提条件等。

【模块功能概述】:
1. Bartlett统计量计算：计算多组数据方差齐性的Bartlett统计量
2. 修正统计量计算：应用修正因子调整统计量
3. 显著性检验：执行卡方检验判断方差是否齐性
4. 合并方差计算：计算各组的合并方差
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 方差齐性检验：验证多组数据方差是否齐性
- 统计方法选择：为后续统计分析方法的选择提供依据
- 数据质量评估：评估多组间变异程度的差异
- 假设检验验证：验证方差分析等方法的前提条件

【统计方法选择指南】:
1. Bartlett检验适用条件：
   - 各组数据均服从正态分布
   - 数据相互独立
   - 适用于3组及以上数据的方差比较
   - 对正态性假设较为敏感

2. 临床应用场景：
   - 比较多个治疗组某项指标的变异性差异
   - 验证方差分析前的方差齐性假设
   - 检验不同处理组间变异程度的差异

【结果解读注意事项】:
1. Bartlett统计量解释：统计量越大表示方差差异越显著
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
- "我的多组数据适合用Bartlett检验比较方差吗？"
- "如何解释Bartlett检验的结果？"
- "P值小于0.05意味着什么？"
- "Bartlett检验的前提条件是什么？"
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
from app.schemas.request_data.hv_param import HVParamBartlett

def calculate_bartlett_test_variance_homogeneity(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    执行Bartlett方差齐性检验
    
    在临床研究中，这用于检验多组数据的方差是否齐性，
    即比较多个组数据的变异程度是否存在显著差异。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据
        
    Returns:
        Dict[str, Any]: 包含Bartlett检验统计量和结果的字典
        
    Raises:
        ValueError: 当数据组数少于2组或数据无效时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) < 2:
        raise ValueError("Bartlett检验至少需要2组数据")
    
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
    total_n = sum(info["sample_size"] for info in groups_info)  # 总样本量
    
    # 计算合并方差
    pooled_variance = 0
    denominator = 0
    for info in groups_info:
        pooled_variance += (info["sample_size"] - 1) * info["variance"]
        denominator += (info["sample_size"] - 1)
    
    if denominator == 0:
        raise ValueError("无法计算合并方差")
    
    pooled_variance /= denominator
    
    # 计算Bartlett检验统计量
    # Bartlett统计量: T = (N-k)*ln(Sp²) - Σ[(ni-1)*ln(si²)]
    # 其中N是总样本量，k是组数，Sp²是合并方差，si²是第i组方差
    
    N = total_n
    ln_pooled_var = np.log(pooled_variance)
    
    # 计算第一项
    first_term = (N - k) * ln_pooled_var
    
    # 计算第二项
    second_term = 0
    for info in groups_info:
        ni = info["sample_size"]
        si_squared = info["variance"]
        if si_squared <= 0:
            raise ValueError(f"第{groups_info.index(info)+1}组方差必须大于0")
        second_term += (ni - 1) * np.log(si_squared)
    
    # Bartlett统计量
    bartlett_statistic = first_term - second_term
    
    # 计算修正因子C
    # C = 1 + (1/(3*(k-1))) * (Σ(1/(ni-1)) - 1/(N-k))
    sum_reciprocal_df = sum(1/(info["sample_size"] - 1) for info in groups_info)
    c_factor = 1 + (1/(3*(k-1))) * (sum_reciprocal_df - 1/(N-k))
    
    # 修正后的统计量
    corrected_statistic = bartlett_statistic / c_factor
    
    # 自由度
    df = k - 1
    
    # 计算P值（基于卡方分布）
    p_value = 1 - stats.chi2.cdf(corrected_statistic, df)
    
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
            "groups_info": groups_info
        },
        "test_statistics": {
            "bartlett_statistic": float(bartlett_statistic),
            "corrected_statistic": float(corrected_statistic),
            "chi_square_statistic": float(corrected_statistic),
            "degrees_of_freedom": int(df),
            "c_factor": float(c_factor),
            "pooled_variance": float(pooled_variance)
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
            "statistic_interpretation": f"Bartlett统计量为{corrected_statistic:.4f}，自由度为{df}",
            "p_value_interpretation": f"P值为{p_value:.6f}",
            "significance_interpretation": _get_significance_interpretation(p_value),
            "homogeneity_assessment": _get_homogeneity_assessment(p_value),
            "assumption_note": "注意：Bartlett检验对数据的正态性假设较为敏感"
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


def cal_result_hv_bartlett(param: HVParamBartlett) -> Dict[str, Any]:
    """
    生成Bartlett方差齐性检验统计分析的完整报告字典
    
    此函数整合了Bartlett方差齐性检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的Bartlett检验结果。报告包括输入参数、
    检验统计量、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解Bartlett检验的特征。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一组数据
    
    Returns:
        Dict[str, Any]: 包含Bartlett方差齐性检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 执行Bartlett检验
    # 从参数对象解构
    data_list = [item.data_list for item in param.stats_data_list]

    results = calculate_bartlett_test_variance_homogeneity(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Bartlett方差齐性检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"组数: {results['input_parameters']['group_count']}, 总样本量: {results['input_parameters']['total_sample_size']}"
    }
    
    return result_dict