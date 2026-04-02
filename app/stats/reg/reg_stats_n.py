"""
多元线性回归统计分析模块 (Multiple Linear Regression Statistics Module)

本模块提供全面的多元线性回归统计分析功能，用于探索一个因变量与多个自变量之间的线性关系，
是医学研究和数据分析中最重要的统计方法之一。多元线性回归在临床研究中广泛应用，
如分析多个因素对疗效的影响、构建预后模型等。

【模块功能概述】:
1. 回归系数计算：计算截距和各变量的回归系数
2. 模型拟合评估：计算决定系数R²、调整R²等指标
3. 参数显著性检验：执行t检验判断各自变量系数是否显著
4. 模型整体检验：执行F检验判断模型整体显著性
5. 残差分析：计算预测值和残差，用于模型诊断
6. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 多因素分析：同时考虑多个影响因素
- 预测建模：基于多个变量预测目标变量
- 因素筛选：识别重要的预测变量
- 统计控制：控制混杂因素的影响

【统计方法选择指南】:
1. 多元线性回归适用条件：
   - 因变量为连续变量
   - 自变量与因变量之间存在线性关系
   - 残差服从正态分布
   - 残差方差齐性
   - 观测值相互独立
   - 自变量间不存在多重共线性

2. 临床应用场景：
   - 分析多个因素对疗效的综合影响
   - 构建疾病预后预测模型
   - 研究多个生物标志物与疾病的关联
   - 控制混杂因素的影响

【结果解读注意事项】:
1. 回归系数解释：在其他变量不变的情况下，该变量每增加一个单位，因变量的平均变化量
2. R²解释：决定系数表示模型解释的变异百分比
3. P值解释：在零假设成立的情况下，观察到当前或更极端结果的概率
4. 临床意义：统计显著性不等于临床重要性，需结合专业知识判断

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用多元线性回归分析吗？"
- "如何解释各自变量的回归系数？"
- "调整R²与R²的区别是什么？"
- "多元回归分析的前提条件是什么？"
- "如何判断多重共线性问题？"

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


def calculate_multiple_linear_regression(data_list: List[List[float]]) -> Dict:
    """
    基于数据列表的多元线性回归分析
    最后一项作为Y变量，其余作为X变量
    
    在临床研究中，这用于探索多个自变量与一个因变量之间的线性关系，
    如分析多个因素对疗效的综合影响、构建预后模型等。
    通过最小二乘法估计回归系数，评估模型拟合优度。
    
    Args:
        data_list: 数据列表，最后一项为Y变量，其余为X变量
        
    Returns:
        Dict[str, Any]: 包含回归统计量和结果的字典
        
    Raises:
        ValueError: 当数据无效或不满足回归前提时抛出异常
    """
    # 参数验证
    if not data_list or len(data_list) < 2:
        raise ValueError("多元线性回归至少需要2个变量（1个因变量和至少1个自变量）")
    
    # 提取变量数据
    n_variables = len(data_list)
    y_data = data_list[-1]  # 最后一项为Y变量
    x_data_list = data_list[:-1]  # 其余为X变量
    
    # 验证数据长度一致性
    n = len(y_data)
    if n < 3:
        raise ValueError("数据点至少需要3个")
    
    for i, x_data in enumerate(x_data_list):
        if len(x_data) != n:
            raise ValueError(f"第{i+1}个X变量数据长度与Y变量不一致")
    
    k = len(x_data_list)  # 自变量个数
    if n <= k:  # 修改这里：样本量必须严格大于自变量个数
        raise ValueError(f"样本量({n})必须大于自变量个数({k})")
    
    # 转换为numpy数组
    y = np.array(y_data)
    X = np.column_stack([np.array(x_data) for x_data in x_data_list])
    
    # 添加截距项
    X_with_intercept = np.column_stack([np.ones(n), X])
    
    # 计算回归系数 (β = (X'X)^(-1)X'y)
    try:
        beta = np.linalg.solve(X_with_intercept.T @ X_with_intercept, X_with_intercept.T @ y)
    except np.linalg.LinAlgError:
        raise ValueError("设计矩阵奇异，无法计算回归系数")
    
    # 计算预测值和残差
    y_pred = X_with_intercept @ beta
    residuals = y - y_pred
    
    # 计算各种平方和
    ss_res = np.sum(residuals ** 2)  # 残差平方和
    ss_tot = np.sum((y - np.mean(y)) ** 2)  # 总平方和
    
    # 计算决定系数R²
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
    
    # 调整R²
    adjusted_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - k - 1)
    
    # 计算标准误差
    mse = ss_res / (n - k - 1)  # 均方误差
    rmse = np.sqrt(mse)         # 均方根误差
    
    # 计算协方差矩阵和标准误
    cov_matrix = mse * np.linalg.inv(X_with_intercept.T @ X_with_intercept)
    se_beta = np.sqrt(np.diag(cov_matrix))
    
    # 计算t统计量和p值
    t_stats = beta / se_beta
    p_values = 2 * (1 - stats.t.cdf(np.abs(t_stats), n - k - 1))
    
    # 确保p值在合理范围内
    p_values = np.maximum(np.minimum(p_values, 1.0), 0.0)
    
    # 计算F统计量（模型整体显著性）
    f_statistic = (r_squared / k) / ((1 - r_squared) / (n - k - 1)) if r_squared != 1 else 0
    f_p_value = 1 - stats.f.cdf(f_statistic, k, n - k - 1) if r_squared != 1 else 1.0
    f_p_value = max(min(f_p_value, 1.0), 0.0)
    
    # 计算各变量的基本统计量
    variable_stats = []
    for i, data in enumerate(data_list):
        data_array = np.array(data)
        variable_stats.append({
            "name": f"Variable_{i+1}",
            "mean": float(np.mean(data_array)),
            "std": float(np.std(data_array, ddof=1)),
            "min": float(np.min(data_array)),
            "max": float(np.max(data_array))
        })
    
    # 构造回归系数信息
    coefficients_info = []
    coefficient_names = ["截距"] + [f"X{i+1}" for i in range(k)]
    
    for i, (name, coef, se, t_stat, p_val) in enumerate(zip(coefficient_names, beta, se_beta, t_stats, p_values)):
        coefficients_info.append({
            "name": name,
            "coefficient": float(coef),
            "standard_error": float(se),
            "t_statistic": float(t_stat),
            "p_value": float(p_val),
            "significant_05": bool(p_val < 0.05),
            "significant_01": bool(p_val < 0.01)
        })
    
    return {
        "input_parameters": {
            "sample_size": int(n),
            "number_of_predictors": int(k),
            "variable_statistics": variable_stats
        },
        "regression_coefficients": coefficients_info,
        "model_statistics": {
            "r_squared": float(r_squared),
            "adjusted_r_squared": float(adjusted_r_squared),
            "rmse": float(rmse),
            "f_statistic": float(f_statistic),
            "f_p_value": float(f_p_value),
            "degrees_of_freedom_model": int(k),
            "degrees_of_freedom_residual": int(n - k - 1),
            "degrees_of_freedom_total": int(n - 1)
        },
        "predictions": {
            "predicted_values": [float(val) for val in y_pred],
            "residuals": [float(res) for res in residuals]
        },
        "interpretation": {
            "model_equation": _generate_model_equation(coefficients_info),
            "r_squared_interpretation": f"模型解释了{r_squared*100:.2f}%的因变量变异",
            "model_significance": "模型显著" if f_p_value < 0.05 else "模型不显著",
            "assumption_note": "注意：多元线性回归假设数据满足线性关系、独立性、正态性、方差齐性和无多重共线性"
        }
    }


def _generate_model_equation(coefficients_info: List[Dict]) -> str:
    """生成回归方程字符串"""
    equation_parts = []
    for i, coef_info in enumerate(coefficients_info):
        coef = coef_info['coefficient']
        name = coef_info['name']
        
        if i == 0:  # 截距
            if coef != 0:
                equation_parts.append(f"{coef:.4f}")
        else:  # 自变量系数
            if coef != 0:
                sign = "+" if coef > 0 else "-"
                abs_coef = abs(coef)
                equation_parts.append(f"{sign} {abs_coef:.4f}{name}")
    
    if not equation_parts:
        return "Ŷ = 0"
    
    # 处理第一个项的符号
    first_part = equation_parts[0]
    if first_part.startswith("+ "):
        first_part = first_part[2:]
    elif first_part.startswith("- "):
        first_part = "-" + first_part[2:]
    
    equation = "Ŷ = " + first_part
    for part in equation_parts[1:]:
        equation += " " + part
    
    return equation


def cal_result_reg_n(data_list: List[List[float]]) -> Dict[str, Any]:
    """
    生成多元线性回归分析的完整报告字典
    
    此函数整合了多元线性回归的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的回归分析结果。报告包括输入参数、
    回归系数、模型统计量、预测值和统计解释等信息，
    便于临床医生和研究人员快速理解回归分析的特征。
    
    Args:
        data_list: 数据列表，最后一项为Y变量，其余为X变量
    
    Returns:
        Dict[str, Any]: 包含多元线性回归分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - regression_coefficients: 回归系数
            - model_statistics: 模型统计量
            - predictions: 预测值和残差
            - interpretation: 统计解释
    """
    # 执行多元线性回归分析
    results = calculate_multiple_linear_regression(data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "多元线性回归分析",
        "input_parameters": results["input_parameters"],
        "regression_coefficients": results["regression_coefficients"],
        "model_statistics": results["model_statistics"],
        "predictions": results["predictions"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}, 自变量数: {results['input_parameters']['number_of_predictors']}"
    }
    
    return result_dict