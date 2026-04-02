"""
临床完全随机设计方差分析统计分析模块 (Clinical Completely Randomized Design ANOVA Statistics Module)

本模块提供全面的完全随机设计方差分析统计分析功能，用于临床数据中多个处理组均值差异的比较，
是临床试验数据分析的重要组成部分。完全随机设计方差分析在医学研究中广泛应用，
如比较不同治疗方案的效果、不同药物剂量的疗效差异等。

【模块功能概述】:
1. 平方和计算：计算处理间平方和、误差平方和和总平方和
2. 自由度计算：计算各变异来源的自由度
3. 均方计算：计算处理间均方和误差均方
4. F统计量计算：计算F统计量用于检验处理间差异
5. 显著性检验：执行F检验判断处理间差异是否显著
6. 效应量计算：计算Eta平方评估处理效应大小
7. 多重比较：为处理组间差异提供配对比较信息

【临床应用价值】:
- 治疗方案比较：比较多种治疗方法的疗效差异
- 药物剂量研究：评估不同剂量药物的疗效差异
- 疗效评估：判断不同干预措施的效果差异
- 数据质量评估：识别处理组间的变异来源

【统计方法选择指南】:
1. 完全随机设计方差分析适用条件：
   - 各处理组数据独立
   - 各处理组数据服从正态分布
   - 各处理组方差齐性
   - 至少3个处理组（2组时等价于独立样本t检验）

2. 临床应用场景：
   - 比较不同药物治疗方案的疗效
   - 评估不同手术方法的效果差异
   - 分析不同护理措施的改善程度

【结果解读注意事项】:
1. F统计量解释：衡量处理间变异与误差变异的比值
2. P值解释：在零假设（各处理组均值相等）成立的情况下，观察到当前或更极端结果的概率
3. 显著性水平：通常使用α=0.05作为判断标准
4. 效应量解释：评估处理效应的实际意义，而不仅仅是统计显著性
5. 临床意义：统计学显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的多组临床数据适合用完全随机设计方差分析吗？"
- "如何解释方差分析的F统计量？"
- "P值小于0.05意味着什么？"
- "方差分析的前提条件是什么？"
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
from app.schemas.request_data.anova_param import AnovaParamCRD


def calculate_crd_anova(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    执行完全随机设计(CRD)的方差分析
    
    在临床研究中，这用于比较多个处理组的均值是否存在显著差异，
    以评估不同治疗方案或干预措施的效果差异。
    
    Args:
        data_list: List[List[float]]，每个子列表代表一个处理组的数据
        
    Returns:
        Dict[str, Any]: 包含CRD ANOVA统计量和结果的字典
        
    Raises:
        ValueError: 当处理组数不足或某个处理组样本量不足时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) < 2:
        raise ValueError("完全随机设计方差分析至少需要2个处理组")
    
    # 提取数据和基本信息
    treatments_data = []
    treatments_info = []
    
    for i, data in enumerate(data_list):
        data_array = np.array(data)
        n = len(data_array)
        mean = np.mean(data_array)
        std = np.std(data_array, ddof=1)
        var = np.var(data_array, ddof=1)
        
        if n < 2:
            raise ValueError(f"第{i+1}个处理组样本量不足，至少需要2个观测值")
        
        treatments_data.append(data_array)
        treatments_info.append({
            "group_name": f"处理组_{i+1}",
            "sample_size": n,
            "mean": float(mean),
            "std": float(std),
            "variance": float(var)
        })
    
    t = len(treatments_data)  # 处理数组
    total_n = sum(info["sample_size"] for info in treatments_info)  # 总样本量
    
    # 计算总均值
    all_data = np.concatenate(treatments_data)
    grand_mean = np.mean(all_data)
    
    # 计算处理间平方和 (SSTreatments)
    ss_treatments = 0
    for info, data in zip(treatments_info, treatments_data):
        treatment_mean = info["mean"]
        n_i = info["sample_size"]
        ss_treatments += n_i * (treatment_mean - grand_mean) ** 2
    
    # 计算误差平方和 (SSE)
    ss_error = 0
    for info, data in zip(treatments_info, treatments_data):
        treatment_mean = info["mean"]
        ss_error += np.sum((data - treatment_mean) ** 2)
    
    # 计算总平方和 (SSTotal)
    ss_total = ss_treatments + ss_error
    
    # 自由度计算
    df_treatments = t - 1
    df_error = total_n - t
    df_total = total_n - 1
    
    # 均方计算
    ms_treatments = ss_treatments / df_treatments if df_treatments > 0 else 0
    ms_error = ss_error / df_error if df_error > 0 else 0
    
    # F统计量
    f_statistic = ms_treatments / ms_error if ms_error > 0 else 0
    
    # 计算P值
    if ms_error > 0 and df_treatments > 0 and df_error > 0:
        p_value = 1 - stats.f.cdf(f_statistic, df_treatments, df_error)
    else:
        p_value = 1.0
    
    # 确保P值在合理范围内
    p_value = max(min(p_value, 1.0), 0.0)
    
    # 显著性判断
    significant_01 = p_value < 0.01
    significant_05 = p_value < 0.05
    significant_10 = p_value < 0.10
    
    # 计算效应量 (Eta-squared)
    eta_squared = ss_treatments / ss_total if ss_total > 0 else 0
    
    # 计算各处理组间的多重比较准备（仅当处理组≥3时）
    pairwise_comparisons = []
    if t >= 3:
        for i in range(t):
            for j in range(i + 1, t):
                mean_diff = treatments_info[i]["mean"] - treatments_info[j]["mean"]
                pooled_std = np.sqrt(ms_error * (1/treatments_info[i]["sample_size"] + 
                                               1/treatments_info[j]["sample_size"]))
                pairwise_comparisons.append({
                    "comparison": f"{treatments_info[i]['group_name']} vs {treatments_info[j]['group_name']}",
                    "mean_difference": float(mean_diff),
                    "standard_error": float(pooled_std)
                })
    
    return {
        "input_parameters": {
            "treatment_count": int(t),
            "total_sample_size": int(total_n),
            "treatment_names": [info["group_name"] for info in treatments_info],
            "sample_sizes": [info["sample_size"] for info in treatments_info]
        },
        "descriptive_statistics": {
            "treatment_means": [info["mean"] for info in treatments_info],
            "treatment_stds": [info["std"] for info in treatments_info],
            "treatment_vars": [info["variance"] for info in treatments_info],
            "grand_mean": float(grand_mean)
        },
        "anova_statistics": {
            "ss_treatments": float(ss_treatments),
            "ss_error": float(ss_error),
            "ss_total": float(ss_total),
            "ms_treatments": float(ms_treatments),
            "ms_error": float(ms_error),
            "f_statistic": float(f_statistic),
            "degrees_of_freedom_treatments": int(df_treatments),
            "degrees_of_freedom_error": int(df_error),
            "eta_squared": float(eta_squared)
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
        "pairwise_comparisons": pairwise_comparisons,
        "interpretation": {
            "statistic_interpretation": f"F统计量为{f_statistic:.4f}，处理间自由度为{df_treatments}，误差自由度为{df_error}",
            "p_value_interpretation": f"P值为{p_value:.6f}",
            "significance_interpretation": _get_significance_interpretation(p_value),
            "effect_size_interpretation": _get_effect_size_interpretation(eta_squared),
            "assumption_note": "注意：完全随机设计方差分析假设数据满足正态性、方差齐性和独立性"
        }
    }


def _get_significance_interpretation(p_value: float) -> str:
    """获取显著性解释"""
    if p_value < 0.01:
        return f"P值={p_value:.6f} < 0.01，处理间差异极显著"
    elif p_value < 0.05:
        return f"P值={p_value:.6f} < 0.05，处理间差异显著"
    elif p_value < 0.10:
        return f"P值={p_value:.6f} < 0.10，处理间差异边缘显著"
    else:
        return f"P值={p_value:.6f} ≥ 0.10，处理间差异不显著"


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


def cal_result_anova_crd(param: AnovaParamCRD) -> Dict[str, Any]:
    """
    生成完全随机设计方差分析统计分析的完整报告字典
    
    此函数整合了完全随机设计方差分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的方差分析结果。报告包括输入参数、
    描述性统计、方差分析统计量、显著性检验结果等信息，
    便于临床医生和研究人员快速理解方差分析的特征。
    
    Args:
        param: AnovaParamCRD对象，包含完全随机设计方差分析所需参数
        
    Returns:
        Dict[str, Any]: 包含完全随机设计方差分析统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - descriptive_statistics: 描述性统计信息
            - anova_statistics: 方差分析统计量信息
            - significance_tests: 显著性检验结果
            - pairwise_comparisons: 多重比较结果
            - interpretation: 结果的专业解释
    """
    # 提取数据列表
    data_list = [data.data_list for data in param.stats_data_list]
    
    # 执行完全随机设计方差分析
    results = calculate_crd_anova(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "完全随机设计方差分析",
        "input_parameters": results["input_parameters"],
        "descriptive_statistics": results["descriptive_statistics"],
        "anova_statistics": results["anova_statistics"],
        "significance_tests": results["significance_tests"],
        "pairwise_comparisons": results["pairwise_comparisons"],
        "interpretation": results["interpretation"],
        "remark": f"处理数组: {results['input_parameters']['treatment_count']}, 总样本量: {results['input_parameters']['total_sample_size']}"
    }
    
    return result_dict