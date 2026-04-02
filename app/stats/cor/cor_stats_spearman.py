"""
临床Spearman秩相关分析统计分析模块 (Clinical Spearman Rank Correlation Analysis Statistics Module)

本模块提供全面的Spearman秩相关分析统计分析功能，用于临床数据中两个变量之间的单调关系分析，
是一种重要的非参数统计分析方法。Spearman秩相关通过对变量的秩次进行Pearson相关分析，
得到介于-1到1之间的相关系数，广泛应用于医学研究中的非正态数据关联性分析、
等级资料相关性评估、预后因子筛选等领域。

【模块功能概述】:
1. 秩相关系数计算：计算Spearman秩相关系数
2. 显著性检验：执行检验判断相关系数的显著性
3. 置信区间：计算相关系数的置信区间
4. 决定系数：计算R²以评估解释力度
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 非正态数据关联性分析：评估非正态分布变量之间的关系
- 等级资料相关性评估：分析等级数据间的关联程度
- 预后因子筛选：识别与预后相关的等级指标
- 研究假设验证：验证变量间存在单调关系的假设

【统计方法选择指南】:
1. Spearman相关适用条件：
   - 两变量不必满足正态分布
   - 两变量之间存在单调关系（不一定线性）
   - 数据中可能存在异常值
   - 观测值独立
   - 变量至少为有序分类或连续变量

2. 临床应用场景：
   - 分析疾病严重程度与疗效评分之间的关系
   - 评估患者满意度与治疗时间之间的关联
   - 比较病理分级与预后评分的相关性
   - 研究症状评分与生活质量之间的关系

【结果解读注意事项】:
1. 相关系数解释：接近1表示强正单调关系，接近-1表示强负单调关系，接近0表示无单调关系
2. P值解释：在零假设（真实相关系数为0）成立的情况下，观察到当前或更极端结果的概率
3. 决定系数解释：表示一个变量的秩次能够被另一个变量的秩次解释的变异百分比
4. 临床意义：相关关系不等于因果关系，需结合专业知识进行解释

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用Spearman相关分析吗？"
- "如何解释Spearman相关的结果？"
- "Spearman和Pearson相关有什么区别？"
- "Spearman相关的前提条件是什么？"
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


def spearman_correlation_from_stats_data(stats_data_list: List[List[float]]) -> Dict[str, Any]:
    """
    基于List[List[float]]的Spearman秩相关分析
    
    在临床研究中，这用于分析两个变量之间的单调关系，
    特别适用于非正态分布数据或等级数据的相关性分析。
    
    Args:
        stats_data_list: List[List[float]]，包含两个变量的数据列表，每个子列表代表一个变量的数据
    
    Returns:
        Dict[str, Any]: 包含相关统计量和结果的字典
        
    Raises:
        ValueError: 当数据不符合分析条件时抛出异常
    """
    # 参数验证
    if not stats_data_list or len(stats_data_list) != 2:
        raise ValueError("Spearman相关分析需要恰好2组数据")
    
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
    
    # 计算Spearman秩相关系数
    # 使用scipy.stats.spearmanr函数
    correlation_result = stats.spearmanr(x, y)
    correlation = correlation_result.correlation
    p_value = correlation_result.pvalue
    
    # 计算t统计量（用于大样本近似）
    # 对于Spearman相关，当n较大时可以用t检验近似
    if abs(correlation) < 1 and n > 10:
        t_stat = correlation * np.sqrt((n - 2) / (1 - correlation ** 2))
        degrees_of_freedom = n - 2
    else:
        t_stat = float('nan')
        degrees_of_freedom = n - 2
    
    # 计算决定系数
    r_squared = correlation ** 2
    
    # 计算置信区间 (使用Fisher Z变换，适用于大样本)
    if abs(correlation) < 1 and n > 10:
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
            "x_variable": "X变量（秩次）",
            "y_variable": "Y变量（秩次）",
        },
        "correlation_results": {
            "correlation_coefficient": float(correlation),
            "t_statistic": float(t_stat),
            "degrees_of_freedom": int(degrees_of_freedom),
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
            "direction": "正相关" if correlation > 0 else "负相关" if correlation < 0 else "无单调关系",
            "significance": "显著" if p_value < 0.05 else "不显著",
            "method_note": "Spearman秩相关衡量变量间的单调关系，适用于非正态分布数据",
            "r_squared_interpretation": f"决定系数R² = {r_squared:.4f}，表示秩次变量解释了{r_squared*100:.2f}%的变异"
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


def cal_result_cor_spearman(stats_data_list: List[List[float]]) -> Dict[str, Any]:
    """
    生成Spearman秩相关分析统计分析的完整报告字典
    
    此函数整合了Spearman秩相关分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的相关分析结果。报告包括输入参数、
    相关结果和统计解释等信息，
    便于临床医生和研究人员快速理解Spearman秩相关分析的特征。
    
    Args:
        stats_data_list: List[List[float]]，包含两个变量的数据列表，每个子列表代表一个变量的数据
    
    Returns:
        Dict[str, Any]: 包含Spearman秩相关分析统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - correlation_results: 相关分析结果
            - interpretation: 统计解释
    """
    # 执行Spearman相关分析
    results = spearman_correlation_from_stats_data(stats_data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Spearman秩相关分析",
        "input_parameters": results["input_parameters"],
        "correlation_results": results["correlation_results"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}, X变量: {results['input_parameters']['x_variable']}, Y变量: {results['input_parameters']['y_variable']}"
    }
    
    return result_dict