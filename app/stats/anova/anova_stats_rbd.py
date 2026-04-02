"""
临床随机区组设计方差分析统计分析模块 (Clinical Randomized Block Design ANOVA Statistics Module)

本模块提供全面的随机区组设计方差分析统计分析功能，用于临床数据中多个处理组在控制区组效应后的均值差异比较，
是临床试验数据分析的重要组成部分。随机区组设计方差分析在医学研究中广泛应用，
如比较不同治疗方案的效果（控制个体差异）、不同药物剂量的疗效差异（控制病情差异）等。

【模块功能概述】:
1. 平方和计算：计算处理间平方和、区组间平方和、误差平方和和总平方和
2. 自由度计算：计算各变异来源的自由度
3. 均方计算：计算处理间均方、区组间均方和误差均方
4. F统计量计算：分别计算处理和区组的F统计量
5. 显著性检验：执行F检验判断处理和区组效应是否显著
6. 效应量计算：计算Eta平方评估处理和区组效应大小
7. 事后比较：为处理组间差异提供LSD-t检验等事后比较

【临床应用价值】:
- 治疗方案比较：比较多种治疗方法的疗效差异，控制个体差异
- 药物剂量研究：评估不同剂量药物的疗效差异，控制基线差异
- 交叉设计分析：分析交叉试验中的处理效应，控制个体和时期效应
- 数据质量评估：识别处理组和区组的变异来源

【统计方法选择指南】:
1. 随机区组设计方差分析适用条件：
   - 各处理组数据独立
   - 各处理组数据服从正态分布
   - 各处理组方差齐性
   - 每个区组内各处理的观测值数量相等
   - 不存在处理与区组的交互效应

2. 临床应用场景：
   - 比较不同药物治疗方案的疗效（控制个体差异）
   - 评估不同手术方法的效果差异（控制病情差异）
   - 分析不同护理措施的改善程度（控制基线差异）

【结果解读注意事项】:
1. 处理F统计量解释：衡量处理间变异与误差变异的比值
2. 区组F统计量解释：衡量区组间变异与误差变异的比值
3. P值解释：在零假设（处理/区组均值相等）成立的情况下，观察到当前或更极端结果的概率
4. 显著性水平：通常使用α=0.05作为判断标准
5. 效应量解释：评估处理和区组效应的实际意义，而不仅仅是统计显著性
6. 临床意义：统计学显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的多组临床数据适合用随机区组设计方差分析吗？"
- "如何解释随机区组设计方差分析的F统计量？"
- "P值小于0.05意味着什么？"
- "随机区组设计方差分析的前提条件是什么？"
- "如何判断处理效应的临床意义？"

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
from app.schemas.request_data.anova_param import AnovaParamRBD


def calculate_rbd_anova(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    执行随机区组设计(RBD)的方差分析
    
    在临床研究中，这用于比较多个处理组在控制区组效应后的均值是否存在显著差异，
    以评估不同治疗方案或干预措施的效果差异，同时控制可能的混杂因素（区组效应）。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一个区组的数据，每个区组内包含各处理的观测值
        
    Returns:
        Dict[str, Any]: 包含RBD ANOVA统计量和结果的字典
        
    Raises:
        ValueError: 当区组数不足或数据结构不一致时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) < 2:
        raise ValueError("随机区组设计方差分析至少需要2个区组")
    
    # 提取数据和基本信息
    blocks_data = []
    block_names = []
    treatment_count = len(data_list[0])  # 处理数
    
    # 验证数据结构一致性
    for i, data in enumerate(data_list):
        data_array = np.array(data)
        
        if len(data_array) != treatment_count:
            raise ValueError(f"第{i+1}个区组的数据长度({len(data_array)})与其他区组不一致({treatment_count})")
        
        if len(data_array) < 2:
            raise ValueError(f"每个区组至少需要2个处理观测值")
        
        blocks_data.append(data_array)
        block_names.append(f"区组_{i+1}")
    
    b = len(blocks_data)  # 区组数
    t = treatment_count   # 处理数
    total_n = b * t       # 总观测数
    
    # 重构数据为处理×区组的矩阵形式
    data_matrix = np.array(blocks_data)  # shape: (b, t)
    
    # 计算各种均值
    grand_mean = np.mean(data_matrix)
    treatment_means = np.mean(data_matrix, axis=0)  # 每个处理的均值
    block_means = np.mean(data_matrix, axis=1)      # 每个区组的均值
    
    # 计算总平方和 (SST)
    ss_total = np.sum((data_matrix - grand_mean) ** 2)
    
    # 计算处理间平方和 (SSTreatments)
    ss_treatments = b * np.sum((treatment_means - grand_mean) ** 2)
    
    # 计算区组间平方和 (SSBlocks)
    ss_blocks = t * np.sum((block_means - grand_mean) ** 2)
    
    # 计算误差平方和 (SSE)
    ss_error = ss_total - ss_treatments - ss_blocks
    
    # 自由度计算
    df_treatments = t - 1
    df_blocks = b - 1
    df_error = (t - 1) * (b - 1)
    df_total = total_n - 1
    
    # 均方计算
    ms_treatments = ss_treatments / df_treatments if df_treatments > 0 else 0
    ms_blocks = ss_blocks / df_blocks if df_blocks > 0 else 0
    ms_error = ss_error / df_error if df_error > 0 else 0
    
    # F统计量
    f_treatments = ms_treatments / ms_error if ms_error > 0 else 0
    f_blocks = ms_blocks / ms_error if ms_error > 0 else 0
    
    # 计算P值
    p_treatments = 1 - stats.f.cdf(f_treatments, df_treatments, df_error) if ms_error > 0 else 1.0
    p_blocks = 1 - stats.f.cdf(f_blocks, df_blocks, df_error) if ms_error > 0 else 1.0
    
    # 确保P值在合理范围内
    p_treatments = max(min(p_treatments, 1.0), 0.0)
    p_blocks = max(min(p_blocks, 1.0), 0.0)
    
    # 显著性判断
    sig_treat_01 = p_treatments < 0.01
    sig_treat_05 = p_treatments < 0.05
    sig_treat_10 = p_treatments < 0.10
    
    sig_block_01 = p_blocks < 0.01
    sig_block_05 = p_blocks < 0.05
    sig_block_10 = p_blocks < 0.10
    
    # 计算效应量
    eta_squared_treatments = ss_treatments / ss_total if ss_total > 0 else 0
    eta_squared_blocks = ss_blocks / ss_total if ss_total > 0 else 0
    
    # 事后比较（LSD-t检验）
    post_hoc_comparisons = []
    if t >= 3 and ms_error > 0 and df_error > 0:
        try:
            # LSD临界值
            lsd_critical = stats.t.ppf(0.975, df_error) * np.sqrt(2 * ms_error / b)
            
            for i in range(t):
                for j in range(i + 1, t):
                    mean_diff = treatment_means[i] - treatment_means[j]
                    se_diff = np.sqrt(2 * ms_error / b)
                    t_stat = abs(mean_diff) / se_diff if se_diff > 0 else 0
                    p_value = 2 * (1 - stats.t.cdf(t_stat, df_error)) if t_stat > 0 else 1.0
                    
                    post_hoc_comparisons.append({
                        "comparison": f"处理{i+1} vs 处理{j+1}",
                        "mean_difference": float(mean_diff),
                        "standard_error": float(se_diff),
                        "t_statistic": float(t_stat),
                        "p_value": float(p_value),
                        "significant": p_value < 0.05,
                        "lsd_significant": abs(mean_diff) > lsd_critical if lsd_critical > 0 else False
                    })
        except Exception as e:
            # 如果事后比较计算出错，返回空列表但不影响主分析结果
            post_hoc_comparisons = []
    
    return {
        "input_parameters": {
            "block_count": int(b),
            "treatment_count": int(t),
            "total_sample_size": int(total_n),
            "block_names": block_names
        },
        "descriptive_statistics": {
            "treatment_means": [float(m) for m in treatment_means],
            "block_means": [float(m) for m in block_means],
            "grand_mean": float(grand_mean),
            "treatment_stds": [float(np.std(data_matrix[:, i], ddof=1)) for i in range(t)],
            "block_stds": [float(np.std(data_matrix[i, :], ddof=1)) for i in range(b)]
        },
        "anova_statistics": {
            "ss_treatments": float(ss_treatments),
            "ss_blocks": float(ss_blocks),
            "ss_error": float(ss_error),
            "ss_total": float(ss_total),
            "ms_treatments": float(ms_treatments),
            "ms_blocks": float(ms_blocks),
            "ms_error": float(ms_error),
            "f_treatments": float(f_treatments),
            "f_blocks": float(f_blocks),
            "degrees_of_freedom_treatments": int(df_treatments),
            "degrees_of_freedom_blocks": int(df_blocks),
            "degrees_of_freedom_error": int(df_error),
            "eta_squared_treatments": float(eta_squared_treatments),
            "eta_squared_blocks": float(eta_squared_blocks)
        },
        "significance_tests": {
            "p_value_treatments": float(p_treatments),
            "p_value_blocks": float(p_blocks),
            "significant_treatments_at_01": bool(sig_treat_01),
            "significant_treatments_at_05": bool(sig_treat_05),
            "significant_treatments_at_10": bool(sig_treat_10),
            "significant_blocks_at_01": bool(sig_block_01),
            "significant_blocks_at_05": bool(sig_block_05),
            "significant_blocks_at_10": bool(sig_block_10)
        },
        "post_hoc_comparisons": post_hoc_comparisons,
        "interpretation": {
            "treatment_interpretation": _get_treatment_significance_interpretation(p_treatments),
            "block_interpretation": _get_block_significance_interpretation(p_blocks),
            "effect_size_treatment": _get_effect_size_interpretation(eta_squared_treatments),
            "effect_size_block": _get_effect_size_interpretation(eta_squared_blocks),
            "assumption_note": "注意：随机区组设计方差分析假设数据满足正态性、方差齐性和各区组内处理的独立性"
        }
    }


def _get_treatment_significance_interpretation(p_value: float) -> str:
    """获取处理效应显著性解释"""
    if p_value < 0.01:
        return f"处理间P值={p_value:.6f} < 0.01，处理效应极显著"
    elif p_value < 0.05:
        return f"处理间P值={p_value:.6f} < 0.05，处理效应显著"
    elif p_value < 0.10:
        return f"处理间P值={p_value:.6f} < 0.10，处理效应边缘显著"
    else:
        return f"处理间P值={p_value:.6f} ≥ 0.10，处理效应不显著"


def _get_block_significance_interpretation(p_value: float) -> str:
    """获取区组效应显著性解释"""
    if p_value < 0.01:
        return f"区组间P值={p_value:.6f} < 0.01，区组效应极显著"
    elif p_value < 0.05:
        return f"区组间P值={p_value:.6f} < 0.05，区组效应显著"
    elif p_value < 0.10:
        return f"区组间P值={p_value:.6f} < 0.10，区组效应边缘显著"
    else:
        return f"区组间P值={p_value:.6f} ≥ 0.10，区组效应不显著"


def _get_effect_size_interpretation(eta_squared: float) -> str:
    """获取效应量解释"""
    if eta_squared >= 0.14:
        return f"η²={eta_squared:.4f}，大效应"
    elif eta_squared >= 0.06:
        return f"η²={eta_squared:.4f}，中等效应"
    elif eta_squared >= 0.01:
        return f"η²={eta_squared:.4f}，小效应"
    else:
        return f"η²={eta_squared:.4f}，效应量很小"


def cal_result_anova_rbd(param: AnovaParamRBD) -> Dict[str, Any]:
    """
    生成随机区组设计方差分析统计分析的完整报告字典
    
    此函数整合了随机区组设计方差分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的方差分析结果。报告包括输入参数、
    描述性统计、方差分析统计量、显著性检验结果等信息，
    便于临床医生和研究人员快速理解方差分析的特征。
    
    Args:
        param: AnovaParamRBD对象，包含随机区组设计方差分析所需参数
        
    Returns:
        Dict[str, Any]: 包含随机区组设计方差分析统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - descriptive_statistics: 描述性统计信息
            - anova_statistics: 方差分析统计量信息
            - significance_tests: 显著性检验结果
            - post_hoc_comparisons: 事后比较结果
            - interpretation: 结果的专业解释
    """
    # 提取数据列表
    data_list = [data.data_list for data in param.stats_data_list]
    
    # 执行随机区组设计方差分析
    results = calculate_rbd_anova(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "随机区组设计方差分析",
        "input_parameters": results["input_parameters"],
        "descriptive_statistics": results["descriptive_statistics"],
        "anova_statistics": results["anova_statistics"],
        "significance_tests": results["significance_tests"],
        "post_hoc_comparisons": results["post_hoc_comparisons"],
        "interpretation": results["interpretation"],
        "remark": f"区组数: {results['input_parameters']['block_count']}, 处理数: {results['input_parameters']['treatment_count']}, 总样本量: {results['input_parameters']['total_sample_size']}"
    }
    
    return result_dict


# 导出主要函数
__all__ = ['cal_result_anova_rbd', 'calculate_rbd_anova']