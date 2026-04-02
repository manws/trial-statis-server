from typing import Dict, List, Any
import math
from scipy.stats import chi2
from app.schemas.request_data.chi_param import ChiSquareParamRC


def chi_square_test_rc(row_num: int, col_num: int, data_list: List[int]) -> Dict[str, Any]:
    """
    执行RC列联表的卡方检验
    
    Args:
        row_num: 行数
        col_num: 列数
        data_list: 数据列表，按行优先顺序排列
    
    Returns:
        Dict[str, Any]: 检验结果字典
    """
    # 验证输入数据
    if len(data_list) != row_num * col_num:
        raise ValueError(f"数据列表长度({len(data_list)})必须等于行数×列数({row_num}×{col_num})")
    
    if any(x < 0 for x in data_list):
        raise ValueError("数据列表中的所有值必须为非负数")
    
    # 将数据列表转换为二维矩阵
    obs = [[data_list[i * col_num + j] for j in range(col_num)] for i in range(row_num)]
    
    # 计算行合计和列合计
    row_totals = [sum(obs[i][j] for j in range(col_num)) for i in range(row_num)]
    col_totals = [sum(obs[i][j] for i in range(row_num)) for j in range(col_num)]
    grand_total = sum(row_totals)
    
    # 如果总和为0，无法进行检验
    if grand_total == 0:
        raise ValueError("数据总和不能为0")
    
    # 计算期望频数
    exp = [[(row_totals[i] * col_totals[j]) / grand_total for j in range(col_num)] for i in range(row_num)]
    
    # 计算卡方统计量
    chi_square = 0.0
    for i in range(row_num):
        for j in range(col_num):
            if exp[i][j] > 0:
                chi_square += (obs[i][j] - exp[i][j])**2 / exp[i][j]
    
    # 计算自由度
    df = (row_num - 1) * (col_num - 1)
    
    # 计算P值
    p_value = 1 - chi2.cdf(chi_square, df)
    
    # 计算Cramér's V (列联表相关系数)
    cramers_v = math.sqrt(chi_square / (grand_total * min(row_num-1, col_num-1))) if min(row_num-1, col_num-1) > 0 else 0.0
    
    # 计算Phi系数（仅适用于2×2表）
    phi = math.sqrt(chi_square / grand_total) if row_num == 2 and col_num == 2 else None
    
    # 计算Pearson相关比
    pearson_contingency = math.sqrt(chi_square / (chi_square + grand_total)) if chi_square + grand_total > 0 else 0.0
    
    return {
        "chi_square": chi_square,
        "degrees_of_freedom": df,
        "p_value": p_value,
        "cramers_v": cramers_v,
        "phi_coefficient": phi,
        "pearson_contingency": pearson_contingency,
        "observed_table": obs,
        "expected_table": exp,
        "row_totals": row_totals,
        "column_totals": col_totals,
        "grand_total": grand_total
    }


def perform_chi_square_rc_test(row_num: int, col_num: int, data_list: List[int]) -> Dict[str, Any]:
    """
    执行完整的RC列联表卡方检验流程
    
    Args:
        row_num: 行数
        col_num: 列数
        data_list: 数据列表，按行优先顺序排列
    
    Returns:
        Dict[str, Any]: 完整的检验结果字典
    """
    # 执行RC列联表卡方检验
    results = chi_square_test_rc(row_num, col_num, data_list)
    
    # 统计解释
    interpretation = {
        "chi_square_value": f"卡方值为 {results['chi_square']:.4f}",
        "degrees_of_freedom": f"自由度为 {(row_num - 1) * (col_num - 1)}",
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
            "row_number": row_num,
            "column_number": col_num,
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
                "phi_coefficient": results["phi_coefficient"],
                "pearson_contingency": results["pearson_contingency"]
            }
        },
        "significance_tests": significance_tests,
        "interpretation": interpretation
    }


def cal_result_chi_rc(param: ChiSquareParamRC) -> Dict[str, Any]:
    """
    生成RC列联表卡方检验统计分析的完整报告字典
    
    此函数整合了RC列联表卡方检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的卡方检验结果。报告包括输入参数、
    观察频数表、期望频数分析、检验统计量、显著性检验和统计解释等信息，
    便于临床医生和研究人员快速理解卡方检验的特征。
    
    Args:
        param: ChiSquareParamRC对象，包含row_num, col_num, data_list参数
    
    Returns:
        Dict[str, Any]: 包含RC列联表卡方检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"RC列联表卡方检验分析"
            - input_parameters: 输入参数信息
            - contingency_table: 观察频数表
            - expected_frequencies: 期望频数分析
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 从参数对象中提取值
    row_num = param.row_num
    col_num = param.col_num
    data_list = param.data_list
    
    # 执行RC列联表卡方检验
    results = perform_chi_square_rc_test(row_num, col_num, data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "RC列联表卡方检验分析",
        "input_parameters": results["input_parameters"],
        "contingency_table": results["contingency_table"],
        "expected_frequencies": results["expected_frequencies"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"行列数: {row_num}×{col_num}, 数据点数: {len(data_list)}"
    }
    
    return result_dict