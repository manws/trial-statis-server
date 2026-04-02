"""
一元线性回归统计分析模块 (Simple Linear Regression Statistics Module)

本模块提供全面的一元线性回归统计分析功能，用于探索两个连续变量之间的线性关系，
是医学研究和数据分析中最常用的统计方法之一。一元线性回归在临床研究中广泛应用，
如分析药物剂量与疗效的关系、生理指标间的关联等。

【模块功能概述】:
1. 回归系数计算：计算截距和斜率参数
2. 模型拟合评估：计算决定系数R²和调整R²
3. 参数显著性检验：执行t检验判断回归系数是否显著
4. 模型整体检验：执行F检验判断模型整体显著性
5. 残差分析：计算预测值和残差，用于模型诊断
6. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 关联性分析：探索两个变量之间的线性关系
- 预测建模：基于一个变量预测另一个变量
- 因果关系初探：分析变量间的依赖关系
- 统计控制：控制混杂因素的影响

【统计方法选择指南】:
1. 一元线性回归适用条件：
   - 两个变量均为连续变量
   - 变量间存在线性关系
   - 残差服从正态分布
   - 残差方差齐性
   - 观测值相互独立

2. 临床应用场景：
   - 分析药物剂量与疗效的关系
   - 研究生理指标间的关联
   - 预测患者预后指标
   - 探索生物标志物与疾病的关系

【结果解读注意事项】:
1. 回归系数解释：斜率表示X每增加一个单位，Y的平均变化量
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
- "我的数据适合用一元线性回归分析吗？"
- "如何解释回归系数的意义？"
- "R²值为0.5意味着什么？"
- "回归分析的前提条件是什么？"
- "如何判断回归模型的拟合效果？"

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


def calculate_simple_linear_regression(x_data: List[float], y_data: List[float]) -> Dict:
    """
    执行一元线性回归分析
    
    在临床研究中，这用于探索两个连续变量之间的线性关系，
    如药物剂量与疗效的关系、生理指标间的关联等。
    通过最小二乘法估计回归系数，评估模型拟合优度。
    
    Args:
        x_data: 自变量X的数据列表
        y_data: 因变量Y的数据列表
        
    Returns:
        Dict[str, Any]: 包含回归统计量和结果的字典
        
    Raises:
        ValueError: 当数据无效或不满足回归前提时抛出异常
    """
    # 参数验证
    if not x_data or not y_data:
        raise ValueError("X和Y变量不能为空")
    
    if len(x_data) != len(y_data):
        raise ValueError("X和Y变量的数据长度必须相等")
    
    if len(x_data) < 3:
        raise ValueError("数据点至少需要3个")
    
    n = len(x_data)
    
    # 转换为numpy数组
    x = np.array(x_data)
    y = np.array(y_data)
    
    # 计算基本统计量
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_std = np.std(x, ddof=1)
    y_std = np.std(y, ddof=1)
    
    # 计算回归系数
    # 斜率: β1 = Σ[(xi-x̄)(yi-ȳ)] / Σ[(xi-x̄)²]
    # 截距: β0 = ȳ - β1*x̄
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    
    if denominator == 0:
        raise ValueError("X变量方差为0，无法计算回归系数")
    
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    
    # 计算预测值和残差
    y_pred = intercept + slope * x
    residuals = y - y_pred
    
    # 计算残差平方和
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - y_mean) ** 2)
    
    # 计算决定系数R²
    r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0
    
    # 计算标准误差
    mse = ss_res / (n - 2)  # 均方误差
    rmse = np.sqrt(mse)     # 均方根误差
    
    # 计算回归系数的标准误
    slope_se = np.sqrt(mse / denominator)
    intercept_se = np.sqrt(mse * (1/n + x_mean**2/denominator))
    
    # 计算t统计量和p值
    slope_t = slope / slope_se if slope_se != 0 else 0
    intercept_t = intercept / intercept_se if intercept_se != 0 else 0
    
    slope_p = 2 * (1 - stats.t.cdf(abs(slope_t), n - 2)) if slope_se != 0 else 1.0
    intercept_p = 2 * (1 - stats.t.cdf(abs(intercept_t), n - 2)) if intercept_se != 0 else 1.0
    
    # 计算F统计量
    f_statistic = (r_squared / 1) / ((1 - r_squared) / (n - 2)) if r_squared != 1 else 0
    f_p_value = 1 - stats.f.cdf(f_statistic, 1, n - 2) if r_squared != 1 else 1.0
    
    # 确保p值在合理范围内
    slope_p = max(min(slope_p, 1.0), 0.0)
    intercept_p = max(min(intercept_p, 1.0), 0.0)
    f_p_value = max(min(f_p_value, 1.0), 0.0)
    
    return {
        "input_parameters": {
            "sample_size": int(n),
            "x_mean": float(x_mean),
            "y_mean": float(y_mean),
            "x_std": float(x_std),
            "y_std": float(y_std)
        },
        "regression_coefficients": {
            "intercept": float(intercept),
            "slope": float(slope),
            "intercept_se": float(intercept_se),
            "slope_se": float(slope_se),
            "intercept_t": float(intercept_t),
            "slope_t": float(slope_t),
            "intercept_p": float(intercept_p),
            "slope_p": float(slope_p)
        },
        "model_statistics": {
            "r_squared": float(r_squared),
            "adjusted_r_squared": float(1 - (1 - r_squared) * (n - 1) / (n - 2)),
            "rmse": float(rmse),
            "f_statistic": float(f_statistic),
            "f_p_value": float(f_p_value),
            "degrees_of_freedom_model": 1,
            "degrees_of_freedom_residual": int(n - 2)
        },
        "predictions": {
            "predicted_values": [float(val) for val in y_pred],
            "residuals": [float(res) for res in residuals]
        },
        "interpretation": {
            "regression_equation": f"Ŷ = {intercept:.4f} + {slope:.4f}X",
            "slope_interpretation": f"X每增加1个单位，Y平均变化{slope:.4f}个单位",
            "r_squared_interpretation": f"模型解释了{r_squared*100:.2f}%的Y变量变异",
            "model_significance": "模型显著" if f_p_value < 0.05 else "模型不显著",
            "assumption_note": "注意：线性回归假设数据满足线性关系、独立性、正态性和方差齐性"
        }
    }


def cal_result_reg_1(x_data: List[float], y_data: List[float]) -> Dict[str, Any]:
    """
    生成一元线性回归分析的完整报告字典
    
    此函数整合了一元线性回归的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的回归分析结果。报告包括输入参数、
    回归系数、模型统计量、预测值和统计解释等信息，
    便于临床医生和研究人员快速理解回归分析的特征。
    
    Args:
        x_data: 自变量X的数据列表
        y_data: 因变量Y的数据列表
    
    Returns:
        Dict[str, Any]: 包含一元线性回归分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - regression_coefficients: 回归系数
            - model_statistics: 模型统计量
            - predictions: 预测值和残差
            - interpretation: 统计解释
    """
    # 执行一元线性回归分析
    results = calculate_simple_linear_regression(x_data, y_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "一元线性回归分析",
        "input_parameters": results["input_parameters"],
        "regression_coefficients": results["regression_coefficients"],
        "model_statistics": results["model_statistics"],
        "predictions": results["predictions"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}"
    }
    
    return result_dict