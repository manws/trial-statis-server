"""
临床RR列联表卡方检验统计分析模块 (Clinical RR Contingency Table Chi-Square Test Statistics Module)

本模块提供全面的RR列联表卡方检验统计分析功能，用于临床数据中两个分类变量之间的关联性检验，
特别是当行数和列数大于2的情况。RR列联表卡方检验在医学研究中广泛应用，
如比较不同治疗方案在多个分类结果上的差异、评估多个暴露因素与疾病的关系等。

【模块功能概述】:
1. 卡方统计量计算：计算Pearson卡方统计量和P值
2. 期望频数计算：计算列联表中各格子的期望频数
3. 自由度计算：计算卡方分布的自由度
4. 条件验证：验证期望频数条件是否满足
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 多分类疗效评价：比较不同治疗方案在多个分类结果上的差异
- 多因素分析：评估多个暴露因素与疾病的关系
- 诊断试验：评价多分类诊断试验的准确性
- 预后分析：评估多个预后因子与多分类结局的关系

【统计方法选择指南】:
1. RR列联表卡方检验适用条件：
   - 总样本量较大（一般>40）
   - 期望频数小于5的格子不超过总格子数的20%
   - 所有格子的期望频数均大于1
   - 当不满足条件时，考虑使用Fisher精确检验

2. 临床应用场景：
   - 比较多个治疗组在多个分类结局上的差异
   - 评估多个暴露因素与疾病发生的关系
   - 比较多分类诊断方法的准确性

【结果解读注意事项】:
1. 卡方统计量解释：衡量实际频数与理论频数的偏离程度
2. P值解释：在零假设（两变量独立）成立的情况下，观察到当前或更极端结果的概率
3. 自由度解释：df=(行数-1)×(列数-1)，决定了卡方分布的形状
4. 期望频数检查：确保检验结果的可靠性
5. 临床意义：统计学显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的R×C列联表适合用卡方检验吗？"
- "如何解释RR列联表卡方检验的结果？"
- "P值小于0.05意味着什么？"
- "RR列联表卡方检验的前提条件是什么？"
- "如何判断结果的临床意义？"

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
import numpy as np
from scipy.stats import chi2
from app.schemas.request_data.chi_param import ChiSquareParamRR


def chi_square_test_rr(n: int, data_list: List[int]) -> Dict[str, Any]:
    """
    执行RR列联表的卡方检验
    
    Args:
        n: 分类数
        data_list: 数据列表，长度应为n*n，按行优先顺序排列
    
    Returns:
        Dict[str, Any]: 检验结果字典
    """
    # 验证输入数据
    if len(data_list) != n * n:
        raise ValueError(f"数据列表长度({len(data_list)})必须等于n×n({n}×{n})")
    
    if any(x < 0 for x in data_list):
        raise ValueError("数据列表中的所有值必须为非负数")
    
    # 将数据列表转换为二维矩阵
    obs = [[data_list[i * n + j] for j in range(n)] for i in range(n)]
    
    # 计算行合计和列合计
    row_totals = [sum(obs[i][j] for j in range(n)) for i in range(n)]
    col_totals = [sum(obs[i][j] for i in range(n)) for j in range(n)]
    grand_total = sum(row_totals)
    
    # 如果总和为0，无法进行检验
    if grand_total == 0:
        raise ValueError("数据总和不能为0")
    
    # 计算期望频数
    exp = [[(row_totals[i] * col_totals[j]) / grand_total for j in range(n)] for i in range(n)]
    
    # 计算卡方统计量
    chi_square = 0.0
    for i in range(n):
        for j in range(n):
            if exp[i][j] > 0:
                chi_square += (obs[i][j] - exp[i][j])**2 / exp[i][j]
    
    # 计算自由度
    df = (n - 1) * (n - 1)
    
    # 计算P值
    p_value = 1 - chi2.cdf(chi_square, df)
    
    # 计算Cramér's V (列联表相关系数)
    cramers_v = math.sqrt(chi_square / (grand_total * (n-1))) if n-1 > 0 else 0.0
    
    # 计算Pearson相关比
    pearson_contingency = math.sqrt(chi_square / (chi_square + grand_total)) if chi_square + grand_total > 0 else 0.0
    
    return {
        "chi_square": chi_square,
        "degrees_of_freedom": df,
        "p_value": p_value,
        "cramers_v": cramers_v,
        "pearson_contingency": pearson_contingency,
        "observed_table": obs,
        "expected_table": exp,
        "row_totals": row_totals,
        "column_totals": col_totals,
        "grand_total": grand_total
    }


def perform_chi_square_rr_test(n: int, data_list: List[int]) -> Dict[str, Any]:
    """
    执行完整的RR列联表卡方检验流程
    
    Args:
        n: 分类数
        data_list: 数据列表，长度应为n*n，按行优先顺序排列
    
    Returns:
        Dict[str, Any]: 完整的检验结果字典
    """
    # 执行RR列联表卡方检验
    results = chi_square_test_rr(n, data_list)
    
    # 统计解释
    interpretation = {
        "chi_square_value": f"卡方值为 {results['chi_square']:.4f}",
        "degrees_of_freedom": f"自由度为 {(n - 1) * (n - 1)}",
        "p_value_interpretation": f"P值为 {results['p_value']:.6f}",
        "significance_95": "有显著性差异" if results["p_value"] < 0.05 else "无显著性差异",
        "cramers_v_interpretation": f"Cramér's V为 {results['cramers_v']:.4f}",
        "association_strength": (
            "强关联" if results['cramers_v'] >= 0.5 else
            "中等关联" if results['cramers_v'] >= 0.3 else
            "弱关联" if results['cramers_v'] >= 0.1 else
            "几乎无关联"
        )
    }
    
    # 显著性检验
    significance_tests = {
        "p_value": results["p_value"],
        "significant_at_0.05": results["p_value"] < 0.05,
        "significant_at_0.01": results["p_value"] < 0.01,
        "interpretation_05": "在0.05水平上显著" if results["p_value"] < 0.05 else "在0.05水平上不显著",
        "interpretation_01": "在0.01水平上显著" if results["p_value"] < 0.01 else "在0.01水平上不显著"
    }
    
    # 期望频数分析
    exp_analysis = {
        "expected_frequencies_table": results["expected_table"],
        "min_expected_frequency": min(min(row) for row in results["expected_table"]),
        "max_expected_frequency": max(max(row) for row in results["expected_table"]),
        "validity_check": "满足要求" if all(exp_val >= 5 for row in results["expected_table"] for exp_val in row) else "不满足要求"
    }
    
    return {
        "input_parameters": {
            "classification_number": n,
            "data_list": data_list,
            "total_cells": len(data_list)
        },
        "contingency_table": {
            "observed_frequencies": results["observed_table"],
            "row_totals": results["row_totals"],
            "column_totals": results["column_totals"],
            "grand_total": results["grand_total"]
        },
        "expected_frequencies": exp_analysis,
        "test_statistics": {
            "chi_square": results["chi_square"],
            "degrees_of_freedom": results["degrees_of_freedom"],
            "p_value": results["p_value"],
            "effect_sizes": {
                "cramers_v": results["cramers_v"],
                "pearson_contingency": results["pearson_contingency"]
            }
        },
        "significance_tests": significance_tests,
        "interpretation": interpretation
    }


def cal_result_chi_rr(param: ChiSquareParamRR) -> Dict[str, Any]:
    """
    生成RR列联表卡方检验统计分析的完整报告字典
    
    此函数整合了RR列联表卡方检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的卡方检验结果。报告包括输入参数、
    观察频数表、期望频数分析、检验统计量、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解卡方检验的特征。
    
    Args:
        param: ChiSquareParamRR对象，包含n, data_list参数
    
    Returns:
        Dict[str, Any]: 包含RR列联表卡方检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"RR列联表卡方检验分析"
            - input_parameters: 输入参数信息
            - contingency_table: 观察频数表
            - expected_frequencies: 期望频数分析
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 从参数对象中提取值
    n = param.n
    data_list = param.data_list
    
    # 执行RR列联表卡方检验
    results = perform_chi_square_rr_test(n, data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "RR列联表卡方检验分析",
        "input_parameters": results["input_parameters"],
        "contingency_table": results["contingency_table"],
        "expected_frequencies": results["expected_frequencies"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"分类数: {n}, 数据点数: {len(data_list)}"
    }
    
    return result_dict