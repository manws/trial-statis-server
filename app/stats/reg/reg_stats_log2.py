"""
二元逻辑回归统计分析模块 (Binary Logistic Regression Statistics Module)

本模块提供全面的二元逻辑回归统计分析功能，用于分析一个二分类因变量与一个或多个自变量之间的关系，
是医学研究和流行病学中最重要的统计方法之一。二元逻辑回归在临床研究中广泛应用，
如预测疾病发生概率、分析危险因素等。

【模块功能概述】:
1. 参数估计：使用最大似然法估计回归系数
2. 模型拟合评估：计算对数似然值、AIC、BIC等指标
3. 参数显著性检验：执行Wald检验判断各自变量系数是否显著
4. 模型整体检验：计算偏差统计量和McFadden's R²
5. 分类性能评估：计算准确率、敏感性、特异性等指标
6. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 疾病预测：预测疾病发生的概率
- 危险因素分析：识别疾病相关的危险因素
- 治疗效果评估：评估治疗成功的概率
- 分类建模：构建二分类预测模型

【统计方法选择指南】:
1. 二元逻辑回归适用条件：
   - 因变量为二分类变量（0/1）
   - 自变量可以是连续变量或分类变量
   - 观测值相互独立
   - 不存在完全分离现象
   - 线性关系假设（在对数几率尺度上）

2. 临床应用场景：
   - 预测疾病发生风险
   - 分析治疗成功的概率
   - 识别疾病危险因素
   - 构建诊断预测模型

【结果解读注意事项】:
1. 回归系数解释：表示自变量每增加一个单位，因变量对数几率的变化量
2. 优势比解释：表示自变量每增加一个单位，因变量发生与不发生比值的变化倍数
3. P值解释：在零假设成立的情况下，观察到当前或更极端结果的概率
4. 临床意义：关注效应大小的临床重要性，而不仅仅是统计显著性

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用二元逻辑回归分析吗？"
- "如何解释优势比的意义？"
- "McFadden's R²与线性回归的R²有何区别？"
- "逻辑回归分析的前提条件是什么？"
- "如何评估逻辑回归模型的性能？"

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
from scipy.optimize import minimize
from typing import Dict, List, Any
from app.schemas.request_data.reg_param import RegParamLog2

def calculate_binary_logistic_regression(x_data_list: List[List[float]], y_data: List[float]) -> Dict:
    """
    基于数据列表的二元逻辑回归分析
    y_data作为因变量Y（0/1二分类），x_data_list作为自变量X
    
    在临床研究中，这用于分析多个因素对二分类结果的影响，
    如预测疾病发生概率、分析治疗成功与否的因素等。
    通过最大似然法估计回归系数，评估模型拟合效果。
    
    Args:
        x_data_list: 自变量X的数据列表
        y_data: 因变量Y的数据列表（0/1二分类）
        
    Returns:
        Dict[str, Any]: 包含逻辑回归统计量和结果的字典
        
    Raises:
        ValueError: 当数据无效或不满足回归前提时抛出异常
    """
    # 参数验证
    if not x_data_list or not y_data:
        raise ValueError("自变量和因变量不能为空")
    
    if len(x_data_list) < 1:
        raise ValueError("至少需要1个自变量")
    
    # 验证数据长度一致性
    n = len(y_data)
    if n < 3:
        raise ValueError("数据点至少需要3个")
    
    for i, x_data in enumerate(x_data_list):
        if len(x_data) != n:
            raise ValueError(f"第{i+1}个X变量数据长度与Y变量不一致")
    
    # 验证Y变量为二分类数据
    unique_y = set(y_data)
    if not unique_y.issubset({0, 1}):
        raise ValueError("因变量Y必须为0/1二分类数据")
    
    if len(unique_y) < 2:
        raise ValueError("因变量Y必须包含0和1两个类别")
    
    k = len(x_data_list)  # 自变量个数
    if n <= k:
        raise ValueError(f"样本量({n})必须大于自变量个数({k})")
    
    # 转换为numpy数组
    y = np.array(y_data)
    X = np.column_stack([np.array(x_data) for x_data in x_data_list])
    
    # 添加截距项
    X_with_intercept = np.column_stack([np.ones(n), X])
    
    # 定义负对数似然函数
    def neg_log_likelihood(beta):
        linear_predictor = X_with_intercept @ beta
        # 防止数值溢出
        linear_predictor = np.clip(linear_predictor, -500, 500)
        mu = 1 / (1 + np.exp(-linear_predictor))
        # 添加小的epsilon防止log(0)
        epsilon = 1e-15
        mu = np.clip(mu, epsilon, 1 - epsilon)
        return -np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu))
    
    # 定义梯度函数
    def gradient(beta):
        linear_predictor = X_with_intercept @ beta
        linear_predictor = np.clip(linear_predictor, -500, 500)
        mu = 1 / (1 + np.exp(-linear_predictor))
        return -X_with_intercept.T @ (y - mu)
    
    # 使用牛顿-拉夫森方法或BFGS优化
    # 初始值设为0
    initial_beta = np.zeros(k + 1)
    
    # 优化求解
    try:
        result = minimize(
            neg_log_likelihood,
            initial_beta,
            method='BFGS',
            jac=gradient,
            options={'maxiter': 1000, 'disp': False}
        )
        
        if not result.success:
            raise ValueError(f"优化失败: {result.message}")
            
        beta = result.x
    except Exception as e:
        raise ValueError(f"参数估计失败: {str(e)}")
    
    # 计算预测概率和残差
    linear_predictor = X_with_intercept @ beta
    linear_predictor = np.clip(linear_predictor, -500, 500)
    predicted_prob = 1 / (1 + np.exp(-linear_predictor))
    residuals = y - predicted_prob
    
    # 计算对数似然值
    log_likelihood = -neg_log_likelihood(beta)
    
    # 计算空模型的对数似然（只有截距项）
    null_beta = np.array([np.log(np.mean(y) / (1 - np.mean(y))) if np.mean(y) != 0 and np.mean(y) != 1 else 0])
    null_X = np.ones((n, 1))
    null_linear = null_X @ null_beta
    null_linear = np.clip(null_linear, -500, 500)
    null_mu = 1 / (1 + np.exp(-null_linear))
    null_mu = np.clip(null_mu, 1e-15, 1 - 1e-15)
    null_log_likelihood = np.sum(y * np.log(null_mu) + (1 - y) * np.log(1 - null_mu))
    
    # 计算模型统计量
    deviance = -2 * (null_log_likelihood - log_likelihood)
    mcfadden_r2 = 1 - (log_likelihood / null_log_likelihood) if null_log_likelihood != 0 else 0
    
    # 计算协方差矩阵和标准误
    # 使用Hessian矩阵的逆作为协方差矩阵
    try:
        # 数值计算Hessian矩阵
        hessian = np.zeros((k + 1, k + 1))
        eps = 1e-6
        
        for i in range(k + 1):
            beta_plus = beta.copy()
            beta_minus = beta.copy()
            beta_plus[i] += eps
            beta_minus[i] -= eps
            
            grad_plus = gradient(beta_plus)
            grad_minus = gradient(beta_minus)
            hessian[:, i] = (grad_plus - grad_minus) / (2 * eps)
        
        # 确保Hessian矩阵对称
        hessian = (hessian + hessian.T) / 2
        
        # 计算协方差矩阵
        cov_matrix = np.linalg.inv(hessian)
        se_beta = np.sqrt(np.diag(cov_matrix))
        
        # 计算Wald统计量和p值
        wald_stats = beta / se_beta
        p_values = 2 * (1 - stats.norm.cdf(np.abs(wald_stats)))
        
        # 确保p值在合理范围内
        p_values = np.maximum(np.minimum(p_values, 1.0), 0.0)
        
    except np.linalg.LinAlgError:
        # 如果Hessian矩阵奇异，使用简化方法
        se_beta = np.ones_like(beta) * np.inf
        wald_stats = np.zeros_like(beta)
        p_values = np.ones_like(beta)
    
    # 计算分类性能指标
    # 预测类别（阈值=0.5）
    predicted_class = (predicted_prob > 0.5).astype(int)
    accuracy = np.mean(predicted_class == y)
    
    # 计算混淆矩阵
    tp = np.sum((y == 1) & (predicted_class == 1))  # 真正例
    tn = np.sum((y == 0) & (predicted_class == 0))  # 真负例
    fp = np.sum((y == 0) & (predicted_class == 1))  # 假正例
    fn = np.sum((y == 1) & (predicted_class == 0))  # 假负例
    
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0  # 敏感性/召回率
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0  # 特异性
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0    # 精确率
    
    # 构造回归系数信息
    coefficients_info = []
    coefficient_names = ["截距"] + [f"X{i+1}" for i in range(k)]
    
    for i, (name, coef, se, wald_stat, p_val) in enumerate(zip(coefficient_names, beta, se_beta, wald_stats, p_values)):
        coefficients_info.append({
            "name": name,
            "coefficient": float(coef),
            "standard_error": float(se),
            "wald_statistic": float(wald_stat),
            "p_value": float(p_val),
            "odds_ratio": float(np.exp(coef)),
            "significant_05": bool(p_val < 0.05),
            "significant_01": bool(p_val < 0.01)
        })
    
    return {
        "input_parameters": {
            "sample_size": int(n),
            "number_of_predictors": int(k)
        },
        "regression_coefficients": coefficients_info,
        "model_statistics": {
            "log_likelihood": float(log_likelihood),
            "null_log_likelihood": float(null_log_likelihood),
            "deviance": float(deviance),
            "mcfadden_r2": float(mcfadden_r2),
            "aic": float(-2 * log_likelihood + 2 * (k + 1)),
            "bic": float(-2 * log_likelihood + (k + 1) * np.log(n))
        },
        "classification_metrics": {
            "accuracy": float(accuracy),
            "sensitivity": float(sensitivity),
            "specificity": float(specificity),
            "precision": float(precision),
            "true_positives": int(tp),
            "true_negatives": int(tn),
            "false_positives": int(fp),
            "false_negatives": int(fn)
        },
        "predictions": {
            "predicted_probabilities": [float(prob) for prob in predicted_prob],
            "predicted_classes": [int(cls) for cls in predicted_class],
            "residuals": [float(res) for res in residuals]
        },
        "interpretation": {
            "model_equation": _generate_logit_equation(coefficients_info),
            "mcfadden_r2_interpretation": f"McFadden's R² = {mcfadden_r2:.4f}，解释了{mcfadden_r2*100:.2f}%的模型拟合改善",
            "model_significance": "模型显著" if deviance > stats.chi2.ppf(0.95, k) else "模型不显著",
            "assumption_note": "注意：逻辑回归假设数据满足独立性、线性关系(对数几率尺度)和无完全分离"
        }
    }


def _generate_logit_equation(coefficients_info: List[Dict]) -> str:
    """生成逻辑回归方程字符串"""
    # 构建线性预测器部分
    linear_parts = []
    for i, coef_info in enumerate(coefficients_info):
        coef = coef_info['coefficient']
        name = coef_info['name']
        
        if i == 0:  # 截距
            if coef != 0:
                linear_parts.append(f"{coef:.4f}")
        else:  # 自变量系数
            if coef != 0:
                sign = "+" if coef > 0 else "-"
                abs_coef = abs(coef)
                linear_parts.append(f"{sign} {abs_coef:.4f}{name}")
    
    if not linear_parts:
        linear_eq = "0"
    else:
        # 处理第一个项的符号
        first_part = linear_parts[0]
        if first_part.startswith("+ "):
            first_part = first_part[2:]
        elif first_part.startswith("- "):
            first_part = "-" + first_part[2:]
        
        linear_eq = first_part
        for part in linear_parts[1:]:
            linear_eq += " " + part
    
    return f"logit(P) = {linear_eq}"


def cal_result_reg_log2(param: RegParamLog2) -> Dict[str, Any]:
    """
    生成二元逻辑回归分析的完整报告字典
    
    此函数整合了二元逻辑回归的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的回归分析结果。报告包括输入参数、
    回归系数、模型统计量、分类性能指标和统计解释等信息，
    便于临床医生和研究人员快速理解回归分析的特征。
    
    Args:
        x_data_list: 自变量X的数据列表
        y_data: 因变量Y的数据列表（0/1二分类）
    
    Returns:
        Dict[str, Any]: 包含二元逻辑回归分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - regression_coefficients: 回归系数
            - model_statistics: 模型统计量
            - classification_metrics: 分类性能指标
            - predictions: 预测概率和类别
            - interpretation: 统计解释
    """
    # 执行二元逻辑回归分析
    # 从参数对象解构: 最后一个是 y_data，其余是 x_data_list
    all_data = [item.data_list for item in param.stats_data_list]
    x_data_list = all_data[:-1]
    y_data = all_data[-1]

    results = calculate_binary_logistic_regression(x_data_list, y_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "二元逻辑回归分析",
        "input_parameters": results["input_parameters"],
        "regression_coefficients": results["regression_coefficients"],
        "model_statistics": results["model_statistics"],
        "classification_metrics": results["classification_metrics"],
        "predictions": results["predictions"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}, 自变量数: {results['input_parameters']['number_of_predictors']}"
    }
    
    return result_dict