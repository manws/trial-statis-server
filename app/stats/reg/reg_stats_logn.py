"""
多项逻辑回归统计分析模块 (Multinomial Logistic Regression Statistics Module)

本模块提供全面的多项逻辑回归统计分析功能，用于分析一个多项分类因变量与一个或多个自变量之间的关系，
是医学研究和社会科学中重要的统计方法之一。多项逻辑回归在临床研究中广泛应用，
如分析多个类别结果的概率、预测多分类结果等。

【模块功能概述】:
1. 参数估计：使用最大似然法估计各对比组的回归系数
2. 模型拟合评估：计算对数似然值、AIC、BIC等指标
3. 参数显著性检验：执行Wald检验判断各自变量系数是否显著
4. 模型整体检验：计算偏差统计量和McFadden's R²
5. 分类性能评估：计算准确率和混淆矩阵
6. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 多类别预测：预测多类别结果的发生概率
- 疾病分型分析：分析不同疾病类型的危险因素
- 治疗选择评估：评估影响治疗选择的因素
- 分类建模：构建多分类预测模型

【统计方法选择指南】:
1. 多项逻辑回归适用条件：
   - 因变量为多项分类变量（≥3类）
   - 自变量可以是连续变量或分类变量
   - 观测值相互独立
   - 满足比例优势假设（无完全分离现象）
   - 独立性无关选项假设（IIA）

2. 临床应用场景：
   - 预测疾病亚型发生概率
   - 分析治疗选择的影响因素
   - 识别疾病严重程度的危险因素
   - 构建多分类诊断模型

【结果解读注意事项】:
1. 回归系数解释：表示自变量每增加一个单位，相对基准类别的对数优势变化量
2. 优势比解释：表示自变量每增加一个单位，相对基准类别的优势比变化倍数
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
- "我的数据适合用多项逻辑回归分析吗？"
- "如何解释相对基准类别的优势比？"
- "McFadden's R²与线性回归的R²有何区别？"
- "多项逻辑回归分析的前提条件是什么？"
- "如何评估多项逻辑回归模型的性能？"

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
from app.schemas.request_data.reg_param import RegParamLogN


def calculate_multinomial_logistic_regression(x_data_list: List[List[float]], y_data: List[float]) -> Dict:
    """
    基于数据列表的多项逻辑回归分析
    y_data作为因变量Y（多分类），x_data_list作为自变量X
    
    在临床研究中，这用于分析多个因素对多项分类结果的影响，
    如预测疾病亚型、分析治疗选择等。
    通过比较各组与基准组的相对概率，评估模型拟合效果。
    
    Args:
        x_data_list: 自变量X的数据列表
        y_data: 因变量Y的数据列表（多分类）
        
    Returns:
        Dict[str, Any]: 包含多项逻辑回归统计量和结果的字典
        
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
    
    # 验证Y变量为分类数据
    unique_y = sorted(list(set(y_data)))
    if len(unique_y) < 2:
        raise ValueError("因变量Y必须包含至少2个不同的类别")
    
    if len(unique_y) == 2:
        raise ValueError("因变量只有2个类别，请使用二元逻辑回归")
    
    # 检查Y变量是否为整数类别
    if not all(isinstance(y, (int, float)) and y == int(y) for y in unique_y):
        raise ValueError("因变量Y必须为整数类别编码")
    
    k = len(x_data_list)  # 自变量个数
    c = len(unique_y)     # 类别数
    
    if n <= k * (c - 1):
        raise ValueError(f"样本量({n})相对于参数数量({k*(c-1)})过小")
    
    # 转换为numpy数组
    y = np.array(y_data)
    X = np.column_stack([np.array(x_data) for x_data in x_data_list])
    
    # 添加截距项
    X_with_intercept = np.column_stack([np.ones(n), X])
    p = X_with_intercept.shape[1]  # 包含截距的参数个数
    
    # 将Y转换为one-hot编码
    y_onehot = np.zeros((n, c))
    for i, category in enumerate(unique_y):
        y_onehot[:, i] = (y == category).astype(int)
    
    # 选择基准类别（通常是第一个类别）
    base_category = unique_y[0]
    comparison_categories = unique_y[1:]  # 与其他类别的比较
    
    # 为每个非基准类别拟合一个二元逻辑回归
    category_models = {}
    
    for i, category in enumerate(comparison_categories):
        # 创建二分类变量：当前类别vs基准类别
        binary_y = (y == category).astype(int)
        
        # 定义负对数似然函数
        def neg_log_likelihood(beta):
            linear_predictor = X_with_intercept @ beta
            # 防止数值溢出
            linear_predictor = np.clip(linear_predictor, -500, 500)
            mu = 1 / (1 + np.exp(-linear_predictor))
            # 添加小的epsilon防止log(0)
            epsilon = 1e-15
            mu = np.clip(mu, epsilon, 1 - epsilon)
            return -np.sum(binary_y * np.log(mu) + (1 - binary_y) * np.log(1 - mu))
        
        # 定义梯度函数
        def gradient(beta):
            linear_predictor = X_with_intercept @ beta
            linear_predictor = np.clip(linear_predictor, -500, 500)
            mu = 1 / (1 + np.exp(-linear_predictor))
            return -X_with_intercept.T @ (binary_y - mu)
        
        # 初始值设为0
        initial_beta = np.zeros(p)
        
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
                raise ValueError(f"类别{category}优化失败: {result.message}")
                
            beta = result.x
            
            # 计算标准误（简化方法）
            try:
                # 数值计算Hessian矩阵
                hessian = np.zeros((p, p))
                eps = 1e-6
                
                for j in range(p):
                    beta_plus = beta.copy()
                    beta_minus = beta.copy()
                    beta_plus[j] += eps
                    beta_minus[j] -= eps
                    
                    grad_plus = gradient(beta_plus)
                    grad_minus = gradient(beta_minus)
                    hessian[:, j] = (grad_plus - grad_minus) / (2 * eps)
                
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
            
            category_models[category] = {
                'coefficients': beta,
                'standard_errors': se_beta,
                'wald_statistics': wald_stats,
                'p_values': p_values,
                'log_likelihood': -neg_log_likelihood(beta)
            }
            
        except Exception as e:
            raise ValueError(f"类别{category}参数估计失败: {str(e)}")
    
    # 计算整体模型统计量
    # 空模型对数似然（多数类预测）
    null_predictions = np.bincount(y.astype(int))
    null_category = np.argmax(null_predictions)
    null_prob = np.max(null_predictions) / n
    null_log_likelihood = n * (null_prob * np.log(null_prob) + (1 - null_prob) * np.log(1 - null_prob))
    
    # 计算模型对数似然
    total_log_likelihood = sum(model['log_likelihood'] for model in category_models.values())
    
    # 计算模型统计量
    deviance = -2 * (null_log_likelihood - total_log_likelihood)
    mcfadden_r2 = 1 - (total_log_likelihood / null_log_likelihood) if null_log_likelihood != 0 else 0
    
    # 计算预测概率和分类
    predicted_probs = np.zeros((n, c))
    predicted_probs[:, 0] = 1  # 基准类别的概率分子为1
    
    # 计算每个非基准类别的概率分子
    for i, category in enumerate(comparison_categories):
        beta = category_models[category]['coefficients']
        linear_predictor = X_with_intercept @ beta
        linear_predictor = np.clip(linear_predictor, -500, 500)
        predicted_probs[:, i + 1] = np.exp(linear_predictor)
    
    # 标准化概率（除以总和）
    prob_sums = np.sum(predicted_probs, axis=1, keepdims=True)
    predicted_probs = predicted_probs / prob_sums
    
    # 预测类别
    predicted_classes = np.argmax(predicted_probs, axis=1)
    predicted_categories = [unique_y[i] for i in predicted_classes]
    
    # 计算分类性能
    accuracy = np.mean(np.array(predicted_categories) == y)
    
    # 构造回归系数信息
    coefficients_info = []
    
    # 基准类别（通常设为0）
    baseline_info = {
        "category": str(base_category),
        "category_name": f"类别{base_category}(基准)",
        "coefficients": [{"name": "截距", "coefficient": 0.0, "standard_error": 0.0, 
                         "wald_statistic": 0.0, "p_value": 1.0, "odds_ratio": 1.0,
                         "significant_05": False, "significant_01": False}]
    }
    coefficients_info.append(baseline_info)
    
    # 其他类别
    for i, category in enumerate(comparison_categories):
        beta = category_models[category]['coefficients']
        se = category_models[category]['standard_errors']
        wald_stats = category_models[category]['wald_statistics']
        p_vals = category_models[category]['p_values']
        
        category_coef_info = []
        coefficient_names = ["截距"] + [f"X{j+1}" for j in range(k)]
        
        for j, (name, coef, se_val, wald_stat, p_val) in enumerate(zip(coefficient_names, beta, se, wald_stats, p_vals)):
            # 确保p_val是数值类型，避免比较出错
            p_val_float = float(p_val)
            category_coef_info.append({
                "name": name,
                "coefficient": float(coef),
                "standard_error": float(se_val),
                "wald_statistic": float(wald_stat),
                "p_value": p_val_float,
                "odds_ratio": float(np.exp(coef)) if j > 0 else float(np.exp(coef)),
                "significant_05": bool(p_val_float < 0.05),
                "significant_01": bool(p_val_float < 0.01)
            })
        
        coefficients_info.append({
            "category": str(category),
            "category_name": f"类别{category} vs 类别{base_category}",
            "coefficients": category_coef_info
        })
    
    return {
        "input_parameters": {
            "sample_size": int(n),
            "number_of_predictors": int(k),
            "number_of_categories": int(c),
            "categories": [str(cat) for cat in unique_y],
            "base_category": str(base_category),
            "comparison_categories": [str(cat) for cat in comparison_categories]
        },
        "regression_coefficients": coefficients_info,
        "model_statistics": {
            "total_log_likelihood": float(total_log_likelihood),
            "null_log_likelihood": float(null_log_likelihood),
            "deviance": float(deviance),
            "mcfadden_r2": float(mcfadden_r2),
            "aic": float(-2 * total_log_likelihood + 2 * k * (c - 1)),
            "bic": float(-2 * total_log_likelihood + k * (c - 1) * np.log(n))
        },
        "classification_metrics": {
            "accuracy": float(accuracy),
            "confusion_matrix": _compute_confusion_matrix(y, predicted_classes, unique_y)
        },
        "predictions": {
            "predicted_probabilities": predicted_probs.tolist(),
            "predicted_categories": predicted_categories,
            "predicted_class_indices": predicted_classes.tolist()
        },
        "interpretation": {
            "model_description": f"多项逻辑回归模型，{c}个类别，以类别{base_category}为基准",
            "mcfadden_r2_interpretation": f"McFadden's R² = {mcfadden_r2:.4f}，解释了{mcfadden_r2*100:.2f}%的模型拟合改善",
            "model_significance": "模型显著" if deviance > stats.chi2.ppf(0.95, k * (c - 1)) else "模型不显著",
            "assumption_note": "注意：多项逻辑回归假设数据满足独立性、广义线性关系和无完全分离"
        }
    }


def _compute_confusion_matrix(true_labels, predicted_classes, categories):
    """计算混淆矩阵"""
    c = len(categories)
    confusion_matrix = np.zeros((c, c), dtype=int)
    
    for true_cat_idx, pred_cat_idx in zip(true_labels.astype(int), predicted_classes):
        true_idx = categories.index(true_cat_idx)
        pred_idx = pred_cat_idx
        confusion_matrix[true_idx, pred_idx] += 1
    
    return confusion_matrix.tolist()


def cal_result_reg_logn(param: RegParamLogN) -> Dict[str, Any]:
    """
    生成多项逻辑回归分析的完整报告字典
    
    此函数整合了多项逻辑回归的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的回归分析结果。报告包括输入参数、
    回归系数、模型统计量、分类性能指标和统计解释等信息，
    便于临床医生和研究人员快速理解回归分析的特征。
    
    Args:
        x_data_list: 自变量X的数据列表
        y_data: 因变量Y的数据列表（多分类）
    
    Returns:
        Dict[str, Any]: 包含多项逻辑回归分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - regression_coefficients: 回归系数（按类别分组）
            - model_statistics: 模型统计量
            - classification_metrics: 分类性能指标
            - predictions: 预测概率和类别
            - interpretation: 统计解释
    """
    # 执行多项逻辑回归分析
    # 从参数对象解构: 最后一个是 y_data，其余是 x_data_list
    all_data = [item.data_list for item in param.stats_data_list]
    x_data_list = all_data[:-1]
    y_data = all_data[-1]

    results = calculate_multinomial_logistic_regression(x_data_list, y_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "多项逻辑回归分析",
        "input_parameters": results["input_parameters"],
        "regression_coefficients": results["regression_coefficients"],
        "model_statistics": results["model_statistics"],
        "classification_metrics": results["classification_metrics"],
        "predictions": results["predictions"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}, 类别数: {results['input_parameters']['number_of_categories']}"
    }
    
    return result_dict