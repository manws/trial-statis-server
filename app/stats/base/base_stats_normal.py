"""
临床正态性检验统计分析模块 (Clinical Normality Test Statistics Module)

本模块提供全面的正态性检验统计分析功能，使用矩法检验数据是否符合正态分布，
是临床数据分析的重要组成部分。正态性检验在医学研究中广泛应用，
如检验生理指标、实验室数据等是否符合正态分布假设。

【模块功能概述】:
1. 描述性统计量计算：计算样本量、均值、标准差等基本统计量
2. 矩法正态性检验：计算偏度和峰度系数及其显著性检验
3. 正态性判断：根据偏度和峰度检验结果判断数据是否符合正态分布
4. 统计量标准化：提供标准化的正态性检验报告

【临床应用价值】:
- 数据分布检验：检验生理指标是否符合正态分布假设
- 统计方法选择：为后续参数检验提供分布假设依据
- 实验数据质量评估：识别偏离正态分布的异常数据
- 统计推断基础：为基于正态分布的统计推理提供前提

【统计方法选择指南】:
1. 矩法正态性检验适用条件：
   - 样本量适中至较大（通常n≥20）
   - 连续型数据
   - 数据无严重离群值

2. 临床应用场景：
   - 检验血压、血糖等生理指标的分布
   - 验证实验室检测结果的正态性
   - 评估医学影像数据的分布特征

【结果解读注意事项】:
1. 偏度系数：衡量分布的对称性，正态分布偏度约为0
2. 峰度系数：衡量分布的尖锐程度，正态分布峰度约为0（超额峰度）
3. P值解释：P值>0.05表示数据符合正态分布
4. 样本量影响：大样本容易产生显著性，需结合效应量判断
5. 临床意义：统计学正态性不等同于临床合理性，需结合专业知识

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的临床数据适合用矩法正态性检验吗？"
- "如何解释偏度和峰度系数？"
- "P值小于0.05意味着什么？"
- "正态性检验的前提条件是什么？"
- "如何判断数据是否符合正态分布？"

【相关标准和规范】:
- CONSORT 声明：随机临床试验报告规范
- STROBE 声明：观察性研究报告规范
- ICH E9 指导原则：临床试验的统计学原则
- SAMPL 指南：医学研究报告中的统计方法描述规范

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from typing import List, Dict, Any
import math
from scipy.stats import norm
from app.schemas.request_data.base_param import BaseParamNormal
from app.schemas.request_data.stats_data import StatsData


def sample_size(data: List[float]) -> int:
    """计算List[float]的样本量"""
    return len(data)


def mean(data: List[float]) -> float:
    """计算List[float]的平均数"""
    if not data:
        raise ValueError("Cannot compute mean of an empty list")
    return sum(data) / len(data)


def std_dev(data: List[float]) -> float:
    """计算List[float]的标准差"""
    if len(data) < 2:
        raise ValueError("Standard deviation requires at least two values")
    
    mean_val = mean(data)
    variance = sum((x - mean_val) ** 2 for x in data) / (len(data) - 1)
    return math.sqrt(variance)


def skewness(data: List[float]) -> float:
    """
    计算List[float]的偏度系数
    
    偏度是衡量数据分布对称性的统计量。在临床数据分析中，
    偏度可以帮助我们了解数据分布是否对称，这对于选择合适的统计方法至关重要。
    
    Args:
        data: 待分析的浮点数列表

    Returns:
        float: 偏度系数，正态分布的偏度约为0

    Raises:
        ValueError: 当数据长度不足或标准差为0时抛出异常
    """
    if len(data) < 3:
        raise ValueError("Skewness requires at least three values")
    
    n = len(data)
    mean_val = mean(data)
    std = std_dev(data)
    
    if std == 0:
        raise ValueError("Cannot compute skewness when standard deviation is zero")
    
    skew_sum = sum(((x - mean_val) / std) ** 3 for x in data)
    # 使用修正系数调整偏度
    corrected_skewness = (n / ((n - 1) * (n - 2))) * skew_sum
    return corrected_skewness


def kurtosis(data: List[float]) -> float:
    """
    计算List[float]的峰度系数
    
    峰度是衡量数据分布尖锐程度的统计量。在临床数据分析中，
    峰度可以帮助我们了解数据分布相对于正态分布的尖锐程度，
    这对于评估数据的离散程度和异常值的潜在存在非常重要。
    
    Args:
        data: 待分析的浮点数列表

    Returns:
        float: 峰度系数（超额峰度），正态分布的峰度约为0

    Raises:
        ValueError: 当数据长度不足或标准差为0时抛出异常
    """
    if len(data) < 4:
        raise ValueError("Kurtosis requires at least four values")
    
    n = len(data)
    mean_val = mean(data)
    std = std_dev(data)
    
    if std == 0:
        raise ValueError("Cannot compute kurtosis when standard deviation is zero")
    
    kurt_sum = sum(((x - mean_val) / std) ** 4 for x in data)
    # 计算超额峰度 (Fisher's definition)，减去3使其正态分布的峰度为0
    excess_kurtosis = (n * (n + 1) * kurt_sum - 3 * (n - 1) ** 2) / ((n - 1) * (n - 2) * (n - 3))
    return excess_kurtosis


def skewness_z_score(data: List[float]) -> float:
    """
    计算偏度系数的z值
    
    z值用于检验偏度系数是否显著不同于0，是正态性检验的重要组成部分。
    在临床研究中，偏度z值可以帮助判断数据分布是否显著偏离对称性。
    
    Args:
        data: 待分析的浮点数列表

    Returns:
        float: 偏度系数的z值

    Raises:
        ValueError: 当样本量不足时抛出异常
    """
    n = len(data)
    if n < 3:
        raise ValueError("Sample size must be at least 3 for skewness test")
    
    skew = skewness(data)
    # 偏度的标准误差
    se_skewness = math.sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))
    
    return skew / se_skewness


def skewness_p_value(data: List[float]) -> float:
    """
    计算偏度系数的p值
    
    p值用于判断偏度系数是否显著不同于0，是正态性检验的关键指标。
    在临床数据分析中，偏度p值帮助判断数据分布是否显著偏离对称性。
    
    Args:
        data: 待分析的浮点数列表

    Returns:
        float: 偏度系数的双侧p值
    """
    z_score = skewness_z_score(data)
    # 使用标准正态分布计算双侧p值
    p_value = 2 * (1 - norm.cdf(abs(z_score)))
    return p_value


def kurtosis_z_score(data: List[float]) -> float:
    """
    计算峰度系数的z值
    
    z值用于检验峰度系数是否显著不同于0，是正态性检验的重要组成部分。
    在临床研究中，峰度z值可以帮助判断数据分布是否显著偏离正态分布的尖锐程度。
    
    Args:
        data: 待分析的浮点数列表

    Returns:
        float: 峰度系数的z值

    Raises:
        ValueError: 当样本量不足时抛出异常
    """
    n = len(data)
    if n < 4:
        raise ValueError("Sample size must be at least 4 for kurtosis test")
    
    kurt = kurtosis(data)
    # 峰度的标准误差
    se_kurtosis = math.sqrt((24 * n * (n - 1) ** 2) / ((n - 2) * (n - 3) * (n + 3) * (n + 5)))
    
    return kurt / se_kurtosis


def kurtosis_p_value(data: List[float]) -> float:
    """
    计算峰度系数的p值
    
    p值用于判断峰度系数是否显著不同于0，是正态性检验的关键指标。
    在临床数据分析中，峰度p值帮助判断数据分布是否显著偏离正态分布的尖锐程度。
    
    Args:
        data: 待分析的浮点数列表

    Returns:
        float: 峰度系数的双侧p值
    """
    z_score = kurtosis_z_score(data)
    # 使用标准正态分布计算双侧p值
    p_value = 2 * (1 - norm.cdf(abs(z_score)))
    return p_value


def cal_result_normal(data: List[float]) -> Dict[str, Any]:
    """
    矩法正态性检验统计分析的完整报告字典
    
    此函数整合了矩法正态性检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的正态性检验分析结果。
    检验假设 H0: γ1=0 且 γ2=0 (总体服从正态分布)
    
    在临床数据分析中，正态性检验是许多统计方法应用的前提条件，
    本函数提供了完整的偏度和峰度检验，帮助研究人员判断数据是否符合正态分布假设。
    
    Args:
        data: 待检验的浮点数列表

    Returns:
        Dict[str, Any]: 包含正态性检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - normality_params: 正态分布参数信息
            - normality_test: 矩法正态检验结果
            - inference: 推断检验结果
    """
    if not data or len(data) < 4:
        raise ValueError("Sample size must be at least 4 for normality test")
    
    # 计算各项统计量
    sample_size_val = sample_size(data)
    mean_val = mean(data)
    std_dev_val = std_dev(data)
    skewness_val = skewness(data)
    kurtosis_val = kurtosis(data)
    skewness_z = skewness_z_score(data)
    skewness_p = skewness_p_value(data)
    kurtosis_z = kurtosis_z_score(data)
    kurtosis_p = kurtosis_p_value(data)
    
    # 确定显著性
    skew_p_range = "不显著" if skewness_p > 0.05 else "显著"
    kurt_p_range = "不显著" if kurtosis_p > 0.05 else "显著"
    
    # α=0.05 和 α=0.1 水平下的判断
    alpha_05 = "显著" if abs(skewness_z) > 1.96 or abs(kurtosis_z) > 1.96 else "不显著"
    alpha_01 = "显著" if abs(skewness_z) > 2.576 or abs(kurtosis_z) > 2.576 else "不显著"
    
    # 构建结果字典
    result_dict = {
        "table_name": "矩法正态性检验",
        "normality_params": {
            "sample_size": sample_size_val,
            "mean": round(mean_val, 4),
            "std_dev": round(std_dev_val, 4),
            "skewness": round(skewness_val, 4),
            "kurtosis": round(kurtosis_val, 4)
        },
        "normality_test": {
            "skewness_z": round(skewness_z, 4),
            "skewness_p": round(skewness_p, 4),
            "skewness_p_range": skew_p_range,
            "kurtosis_z": round(kurtosis_z, 4),
            "kurtosis_p": round(kurtosis_p, 4),
            "kurtosis_p_range": kurt_p_range
        },
        "inference": {
            "alpha_05": alpha_05,
            "alpha_01": alpha_01
        },
        "interpretation": {
            "skewness_interpretation": f"偏度系数为{round(skewness_val, 4)}，表明数据分布{'右偏' if skewness_val > 0 else '左偏' if skewness_val < 0 else '对称'}",
            "kurtosis_interpretation": f"峰度系数为{round(kurtosis_val, 4)}，表明数据分布比正态分布{'更尖锐' if kurtosis_val > 0 else '更平坦' if kurtosis_val < 0 else '尖锐程度相似'}",
            "normality_interpretation": f"在0.05显著性水平下，{'数据不服从正态分布' if alpha_05 == '显著' else '数据服从正态分布'}"
        },
        "remark": f"样本量: {sample_size_val}, 用于检验数据是否服从正态分布"
    }
    
    return result_dict