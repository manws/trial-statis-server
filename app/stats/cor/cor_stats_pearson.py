"""
临床Pearson相关分析统计分析模块 (Clinical Pearson Correlation Analysis Statistics Module)

本模块提供全面的Pearson相关分析统计分析功能，用于临床数据中两个连续变量之间的线性关系分析，
是一种重要的双变量统计分析方法。Pearson相关分析通过对变量间的协方差与各自标准差的比值计算，
得到介于-1到1之间的相关系数，广泛应用于医学研究中的变量关联性分析、指标相关性评估、
预后因子筛选等领域。

【模块功能概述】:
1. 相关系数计算：计算Pearson积矩相关系数
2. 显著性检验：执行t检验判断相关系数的显著性
3. 置信区间：计算相关系数的置信区间
4. 决定系数：计算R²以评估解释力度
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 变量关联性分析：评估两个连续变量之间的线性关系
- 指标相关性评估：分析临床指标间的关联程度
- 预后因子筛选：识别与预后相关的指标
- 研究假设验证：验证变量间存在线性关系的假设

【统计方法选择指南】:
1. Pearson相关适用条件：
   - 两变量均为连续变量
   - 两变量大致呈线性关系
   - 两变量联合正态分布
   - 数据中无明显异常值
   - 观测值独立

2. 临床应用场景：
   - 分析年龄与血压之间的关系
   - 评估药物剂量与疗效之间的关联
   - 比较体重与血糖水平的相关性
   - 研究生物标志物与疾病严重程度的关系

【结果解读注意事项】:
1. 相关系数解释：接近1表示强正相关，接近-1表示强负相关，接近0表示无线性关系
2. P值解释：在零假设（真实相关系数为0）成立的情况下，观察到当前或更极端结果的概率
3. 决定系数解释：表示一个变量能够被另一个变量线性解释的变异百分比
4. 临床意义：相关关系不等于因果关系，需结合专业知识进行解释

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用Pearson相关分析吗？"
- "如何解释Pearson相关的结果？"
- "相关系数的大小有什么意义？"
- "Pearson相关的前提条件是什么？"
- "如何判断相关是否具有统计学意义？"

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


def pearson_correlation_from_stats_data(stats_data_list: List[List[float]]) -> Dict[str, Any]:
    """
    基于List[List[float]]的Pearson直线相关分析
    
    在临床研究中，这用于分析两个连续变量之间的线性关系，
    如分析年龄与血压之间的关系或药物剂量与疗效之间的关联。
    
    Args:
        stats_data_list: List[List[float]]，包含两个变量的数据列表，每个子列表代表一个变量的数据
    
    Returns:
        Dict[str, Any]: 包含相关统计量和结果的字典
        
    Raises:
        ValueError: 当数据不符合分析条件时抛出异常
    """
    # 参数验证
    if not stats_data_list or len(stats_data_list) != 2:
        raise ValueError("Pearson相关分析需要恰好2组数据")
    
    # 提取X和Y变量数据
    x_data = stats_data_list[0]
    y_data = stats_data_list[1]
    
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
    
    # 计算Pearson相关系数
    # r = Σ[(xi-x̄)(yi-ȳ)] / √[Σ(xi-x̄)² * Σ(yi-ȳ)²]
    numerator = np.sum((x - x_mean) * (y - y_mean))
    sum_x_sq = np.sum((x - x_mean) ** 2)
    sum_y_sq = np.sum((y - y_mean) ** 2)
    denominator = np.sqrt(sum_x_sq * sum_y_sq)
    
    if denominator == 0:
        raise ValueError("变量方差为0，无法计算相关系数")
    
    correlation = numerator / denominator
    
    # 计算t统计量和p值
    # t = r * √[(n-2) / (1-r²)]
    if abs(correlation) == 1:
        t_stat = np.inf if correlation == 1 else -np.inf
        p_value = 0.0
    else:
        t_stat = correlation * np.sqrt((n - 2) / (1 - correlation ** 2))
        # 双侧检验
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), n - 2))
    
    # 计算决定系数
    r_squared = correlation ** 2
    
    # 计算置信区间 (使用Fisher Z变换)
    # Fisher Z = 0.5 * ln[(1+r)/(1-r)]
    if abs(correlation) < 1:
        fisher_z = 0.5 * np.log((1 + correlation) / (1 - correlation))
        fisher_z_se = 1 / np.sqrt(n - 3)
        
        # 95%置信区间
        z_critical = stats.norm.ppf(0.975)  # 1.96
        z_lower = fisher_z - z_critical * fisher_z_se
        z_upper = fisher_z + z_critical * fisher_z_se
        
        # 反变换回相关系数
        ci_lower = (np.exp(2 * z_lower) - 1) / (np.exp(2 * z_lower) + 1)
        ci_upper = (np.exp(2 * z_upper) - 1) / (np.exp(2 * z_upper) + 1)
    else:
        ci_lower = correlation
        ci_upper = correlation
    
    return {
        "input_parameters": {
            "sample_size": int(n),
            "x_variable": "X变量",
            "y_variable": "Y变量",
        },
        "correlation_results": {
            "correlation_coefficient": float(correlation),
            "t_statistic": float(t_stat),
            "degrees_of_freedom": int(n - 2),
            "p_value": float(p_value),
            "r_squared": float(r_squared),
            "confidence_interval": {
                "lower": float(ci_lower),
                "upper": float(ci_upper),
                "confidence_level": 0.95
            }
        },
        "interpretation": {
            "strength": _interpret_correlation_strength(abs(correlation)),
            "direction": "正相关" if correlation > 0 else "负相关" if correlation < 0 else "无线性关系",
            "significance": "显著" if p_value < 0.05 else "不显著",
            "r_squared_interpretation": f"决定系数R² = {r_squared:.4f}，表示X变量解释了Y变量{r_squared*100:.2f}%的变异"
        }
    }


def _interpret_correlation_strength(abs_r: float) -> str:
    """解释相关系数强度"""
    if abs_r >= 0.9:
        return "很强"
    elif abs_r >= 0.7:
        return "强"
    elif abs_r >= 0.5:
        return "中等"
    elif abs_r >= 0.3:
        return "弱"
    else:
        return "很弱"


def cal_result_cor_pearson(param: CorParamPearson) -> Dict[str, Any]:
    """
    生成Pearson相关分析统计分析的完整报告字典
    
    此函数整合了Pearson相关分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的相关分析结果。报告包括输入参数、
    相关结果和统计解释等信息，
    便于临床医生和研究人员快速理解Pearson相关分析的特征。
    
    Args:
        stats_data_list: List[List[float]]，包含两个变量的数据列表，每个子列表代表一个变量的数据
    
    Returns:
        Dict[str, Any]: 包含Pearson相关分析统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - correlation_results: 相关分析结果
            - interpretation: 统计解释
    """
    # 从参数对象解构
    stats_data_list = [data.data_list for data in param.stats_data_list]

    # 执行Pearson相关分析
    results = pearson_correlation_from_stats_data(stats_data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Pearson直线相关分析",
        "input_parameters": results["input_parameters"],
        "correlation_results": results["correlation_results"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}, X变量: {results['input_parameters']['x_variable']}, Y变量: {results['input_parameters']['y_variable']}"
    }
    
    return result_dict