"""
临床单样本Wilcoxon符号秩检验统计分析模块 (Clinical Single-Sample Wilcoxon Signed-Rank Test Statistics Module)

本模块提供全面的单样本Wilcoxon符号秩检验统计分析功能，用于临床数据中单个样本的中位数与假设值的比较，
是一种重要的非参数统计分析方法。单样本Wilcoxon符号秩检验通过比较观测值与假设中位数的差值的秩次，
来判断样本中位数是否与假设中位数有显著差异，广泛应用于医学研究中的基线比较、
疗效评估、异常值检测等领域。

【模块功能概述】:
1. 差值计算：计算样本数据与假设中位数的差值
2. 符号秩统计：计算正负差值的秩和
3. 检验统计量：计算Wilcoxon检验统计量W
4. 显著性检验：执行正态近似检验判断统计量的显著性
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 基线比较：比较患者治疗前指标与正常值的差异
- 疗效评估：评估治疗后指标与目标值的差异
- 异常值检测：识别与预期值显著不同的观测值
- 质量控制：检验检测结果与标准值的一致性

【统计方法选择指南】:
1. 单样本Wilcoxon检验适用条件：
   - 数据为连续型变量
   - 数据分布不要求正态性
   - 观测值独立
   - 差值分布对称

2. 临床应用场景：
   - 比较患者治疗前血压与正常血压值的差异
   - 评估某种药物治疗后血糖水平与目标值的差异
   - 分析实验室检测结果与标准值的偏差
   - 研究某种干预措施对生物标志物的影响

【结果解读注意事项】:
1. 检验统计量解释：W值越小，表明观测值与假设中位数的差异越显著
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
- "我的数据适合用单样本Wilcoxon符号秩检验吗？"
- "如何解释单样本Wilcoxon符号秩检验的结果？"
- "单样本Wilcoxon符号秩检验的前提条件是什么？"
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
from typing import Dict, List, Tuple
from scipy.stats import norm, rankdata
import numpy as np


def wilcoxon_signed_rank_single_sample_test(sample_data: List[float], hypothesized_median: float = 0.0) -> Dict:
    """
    单样本Wilcoxon符号秩检验
    
    参数:
    - sample_data: 样本数据列表
    - hypothesized_median: 假设的中位数，默认为0.0
    
    返回:
    - 包含检验统计量和结果的字典
    """
    # 参数验证
    if not sample_data or len(sample_data) == 0:
        raise ValueError("样本数据不能为空")
    
    if len(sample_data) < 2:
        raise ValueError("样本量至少需要2个观测值")
    
    # 计算观测值与假设中位数的差值
    differences = [x - hypothesized_median for x in sample_data]
    
    # 移除差值为0的数据点
    non_zero_diffs = [(diff, i) for i, diff in enumerate(differences) if diff != 0]
    
    if len(non_zero_diffs) == 0:
        raise ValueError("所有差值都为0，无法进行检验")
    
    n = len(non_zero_diffs)
    
    # 提取非零差值
    diffs_only = [diff for diff, _ in non_zero_diffs]
    
    # 计算绝对值和符号
    abs_diffs = [abs(diff) for diff in diffs_only]
    signs = [1 if diff > 0 else -1 for diff in diffs_only]
    
    # 计算秩
    ranks = rankdata(abs_diffs, method='average')
    
    # 计算正秩和和负秩和
    positive_ranks = [ranks[i] for i in range(n) if signs[i] > 0]
    negative_ranks = [ranks[i] for i in range(n) if signs[i] < 0]
    
    W_plus = sum(positive_ranks)
    W_minus = sum(negative_ranks)
    
    # 使用较小的秩和作为检验统计量
    W_statistic = min(W_plus, W_minus)
    
    # 计算期望值和方差（无结时）
    E_W = n * (n + 1) / 4
    Var_W = n * (n + 1) * (2 * n + 1) / 24
    
    # 如果存在结（相同绝对值），需要调整方差
    # 计算结的修正因子
    unique_abs_diffs, counts = np.unique(abs_diffs, return_counts=True)
    tie_correction = sum([count**3 - count for count in counts if count > 1])
    if tie_correction > 0:
        Var_W -= tie_correction / 48
    
    # 计算标准化检验统计量Z
    if Var_W > 0:
        Z_statistic = (W_statistic - E_W) / math.sqrt(Var_W)
    else:
        Z_statistic = 0.0
    
    # 计算双侧P值
    # P值 = 2 × P(Z ≥ |z|)
    p_value_two_sided = 2 * (1 - norm.cdf(abs(Z_statistic)))
    
    # 显著性判断
    is_less_than_05 = p_value_two_sided < 0.05
    is_less_than_01 = p_value_two_sided < 0.01
    
    # 计算样本基本统计量
    sample_mean = float(np.mean(sample_data))
    sample_median = float(np.median(sample_data))
    sample_std = float(np.std(sample_data, ddof=1))
    
    return {
        "input_parameters": {
            "sample_size": len(sample_data),
            "hypothesized_median": hypothesized_median,
            "non_zero_differences": n,
            "zero_differences": len(sample_data) - n
        },
        "sample_statistics": {
            "sample_mean": sample_mean,
            "sample_median": sample_median,
            "sample_std": sample_std
        },
        "rank_statistics": {
            "positive_ranks_sum": float(W_plus),
            "negative_ranks_sum": float(W_minus),
            "test_statistic_W": float(W_statistic),
            "expected_W": float(E_W),
            "variance_W": float(Var_W)
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
            "wilcoxon_interpretation": f"单样本Wilcoxon符号秩检验{'不' if p_value_two_sided >= 0.05 else ''}显著 (p={'≥' if p_value_two_sided >= 0.05 else '<'}0.05)",
            "wilcoxon_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided >= 0.01 else ''}显著 (p={'≥' if p_value_two_sided >= 0.01 else '<'}0.01)"
        }
    }


def cal_result_rs_single(sample_data: List[float], hypothesized_median: float = 0.0) -> Dict:
    """
    生成单样本Wilcoxon符号秩检验统计分析的完整报告字典
    
    此函数整合了单样本Wilcoxon符号秩检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解单样本Wilcoxon符号秩检验的特征。
    
    Args:
        sample_data: List[float]，样本数据列表
        hypothesized_median: float，假设的中位数，默认为0.0
    
    Returns:
        Dict: 包含单样本Wilcoxon符号秩检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - sample_statistics: 样本统计信息
            - rank_statistics: 秩统计信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 执行单样本Wilcoxon符号秩检验
    results = wilcoxon_signed_rank_single_sample_test(sample_data, hypothesized_median)
    
    # 构建结果字典
    result_dict = {
        "table_name": "单样本Wilcoxon符号秩检验分析",
        "input_parameters": results["input_parameters"],
        "sample_statistics": results["sample_statistics"],
        "rank_statistics": results["rank_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {len(sample_data)}, 假设中位数: {hypothesized_median}"
    }
    
    return result_dict
