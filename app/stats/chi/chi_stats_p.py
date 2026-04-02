"""
临床卡方P值计算统计分析模块 (Clinical Chi-square P-value Calculation Statistics Module)

本模块提供全面的卡方P值计算统计分析功能，用于临床数据中卡方统计量的P值计算，
是临床试验数据分析的重要组成部分。卡方P值计算在医学研究中广泛应用，
如评估分类数据间的关联性、检验独立性等。

【模块功能概述】:
1. P值计算：根据卡方统计量和自由度计算P值
2. 统计解释：提供统计显著性解释
3. 结果汇总：生成完整的分析报告

【临床应用价值】:
- 关联性分析：评估分类变量间的关联性
- 独立性检验：检验两个分类变量是否独立
- 拟合优度：检验观察频数与理论频数的拟合程度

【统计方法选择指南】:
1. 卡方检验适用条件：
   - 样本量足够大
   - 理论频数不太小
   - 分类变量间相互独立

2. 临床应用场景：
   - 比较两个分类变量的分布
   - 检验暴露因素与疾病的关系
   - 分析治疗效果与分类指标的关联

【结果解读注意事项】:
1. P值解释：P值越小表示证据越强，拒绝原假设的可能性越大
2. 显著性水平：通常使用α=0.05作为判断标准
3. 临床意义：统计显著性不等同于临床重要性
4. 检验前提：注意卡方检验的应用条件
5. 多重比较：当进行多次检验时需要注意调整显著性水平

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的分类数据适合用卡方检验吗？"
- "如何解释卡方检验的P值？"
- "卡方检验的前提条件是什么？"
- "卡方检验与Fisher精确检验有什么区别？"
- "如何判断卡方检验结果的临床意义？"

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
from scipy.stats import chi2
import numpy as np
from app.schemas.request_data.chi_param import ChiSquareParamP


def chi_square_p_value(chi_square: float, df: int) -> Dict[str, Any]:
    """
    根据卡方值和自由度计算P值
    
    在临床研究中，这用于根据已知的卡方统计量和自由度计算相应的P值，
    以判断统计结果的显著性。
    
    Args:
        chi_square: 卡方统计量值
        df: 自由度
        
    Returns:
        Dict[str, Any]: 包含P值计算结果的字典
        
    Raises:
        ValueError: 当卡方值为负数或自由度非正时抛出异常
    """
    # 参数验证
    if chi_square < 0:
        raise ValueError("卡方值必须为非负数")
    
    if df <= 0:
        raise ValueError("自由度必须为正整数")
    
    # 计算P值：使用卡方分布的生存函数（1-CDF）
    p_value = chi2.sf(chi_square, df)
    
    # 计算其他相关信息
    # 临界值（α=0.05）
    critical_value_05 = chi2.ppf(0.95, df)
    # 临界值（α=0.01）
    critical_value_01 = chi2.ppf(0.99, df)
    
    return {
        "chi_square": chi_square,
        "degrees_of_freedom": df,
        "p_value": p_value,
        "critical_value_alpha_05": critical_value_05,
        "critical_value_alpha_01": critical_value_01,
        "significant_05": chi_square > critical_value_05,
        "significant_01": chi_square > critical_value_01
    }


def perform_chi_square_p_test(chi_square: float, df: int) -> Dict[str, Any]:
    """
    执行卡方P值计算的完整分析
    
    在临床研究中，这用于执行完整的卡方P值计算过程，
    以提供全面的统计分析结果，便于临床医生和研究人员进行全面评估。
    
    Args:
        chi_square: 卡方统计量值
        df: 自由度
    
    Returns:
        Dict[str, Any]: 包含完整分析结果的字典
    """
    # 计算P值
    results = chi_square_p_value(chi_square, df)
    
    # 显著性检验
    significance_tests = {
        "p_value": results["p_value"],
        "significant_at_0.05": results["significant_05"],
        "significant_at_0.01": results["significant_01"],
        "interpretation_05": "在0.05水平上显著" if results["significant_05"] else "在0.05水平上不显著",
        "interpretation_01": "在0.01水平上显著" if results["significant_01"] else "在0.01水平上不显著"
    }
    
    # 统计解释
    interpretation = {
        "chi_square_value": f"卡方值为 {chi_square}",
        "degrees_of_freedom": f"自由度为 {df}",
        "p_value_interpretation": f"P值为 {results['p_value']:.6f}",
        "clinical_significance": "具有统计学意义" if results["p_value"] < 0.05 else "不具有统计学意义",
        "evidence_against_null": "证据支持原假设" if results["p_value"] >= 0.05 else "证据反对原假设"
    }
    
    return {
        "input_parameters": {
            "chi_square_value": chi_square,
            "degrees_of_freedom": df
        },
        "test_results": results,
        "significance_tests": significance_tests,
        "interpretation": interpretation
    }


def cal_result_chi_p(param: ChiSquareParamP) -> Dict[str, Any]:
    """
    生成卡方P值计算统计分析的完整报告字典
    
    此函数整合了卡方P值计算的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的卡方P值计算结果。报告包括输入参数、
    检验结果、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解卡方P值计算的特征。
    
    Args:
        param: ChiSquareParamP对象，包含chi, df参数
    
    Returns:
        Dict[str, Any]: 包含卡方P值计算统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"卡方P值计算分析"
            - input_parameters: 输入参数信息
            - test_results: 检验结果
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 从参数对象中提取值
    chi_square = param.chi
    df = param.df
    
    # 执行卡方P值计算
    results = perform_chi_square_p_test(chi_square, df)
    
    # 构建结果字典
    result_dict = {
        "table_name": "卡方P值计算分析",
        "input_parameters": results["input_parameters"],
        "test_results": results["test_results"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"卡方值: {chi_square}, 自由度: {df}"
    }
    
    return result_dict