"""
独立样本t检验统计分析模块 (Independent Samples T-Test Statistics Module)

本模块提供独立样本t检验功能，用于检验两个独立样本的均数是否有显著性差异。
该检验适用于两组数据相互独立，且来自正态分布总体，方差齐性的情况。

【模块功能概述】:
1. 方差齐性假设下的t检验：计算独立样本t检验的统计量
2. P值计算：基于t分布计算双侧P值
3. 显著性检验：判断统计结果的显著性
4. 结果解释：提供统计结果的专业解释

【临床应用价值】:
- 比较两种治疗方法的效果差异：如比较两组患者的治疗后效果
- 比较不同人群的生理指标：如比较男女的血压水平
- 比较实验组和对照组的差异：如比较用药组和安慰剂组
- 评估不同诊断方法的一致性：如比较两种检测方法的结果

【统计方法选择指南】:
1. 独立样本t检验适用条件：
   - 两组数据相互独立
   - 数据来自正态分布总体
   - 两总体方差相等（方差齐性）
   - 样本量较小（通常n<30）
   - 数据独立

2. 临床应用场景：
   - 比较两种治疗方案的疗效差异
   - 比较不同性别、年龄组的生理指标
   - 比较实验组与对照组的结局指标
   - 验证干预措施的有效性

【结果解读注意事项】:
1. P值解释：在零假设（两组均数无差异）成立的情况下，观察到当前或更极端结果的概率
2. 显著性水平：通常使用α=0.05作为判断标准
3. 统计显著性不代表临床重要性
4. 需注意t检验对正态分布和方差齐性的假设

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用独立样本t检验吗？"
- "如何解释独立样本t检验的P值？"
- "t值的含义是什么？"
- "独立样本t检验的前提条件是什么？"
- "如何判断统计结果的显著性？"

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
from scipy.stats import t
import numpy as np


def independent_t_test_from_stats_data(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    基于StatsData列表的独立资料t检验
    
    参数:
    - stats_data_list: 包含两组数据的列表，每组数据是一个字典，包含field_name和data_list字段
    
    返回:
    - 包含独立t检验统计量和结果的字典
    """
    # 参数验证
    if not stats_data_list or len(stats_data_list) != 2:
        raise ValueError("独立资料t检验需要恰好2组数据")
    
    # 提取两组数据
    group1_data = stats_data_list[0]["data_list"]
    group2_data = stats_data_list[1]["data_list"]
    
    if len(group1_data) < 2 or len(group2_data) < 2:
        raise ValueError("每组数据至少需要2个观测值")
    
    # 计算各组统计量
    n1, n2 = len(group1_data), len(group2_data)
    mean1, mean2 = np.mean(group1_data), np.mean(group2_data)
    std1, std2 = np.std(group1_data, ddof=1), np.std(group2_data, ddof=1)
    
    # 特殊情况处理：如果某组标准差为0
    if std1 == 0 and std2 == 0:
        # 两组都是常数，比较均值差异
        mean_diff = mean1 - mean2
        if mean_diff == 0:
            # 均值相同
            return {
                "input_parameters": {
                    "group1_name": stats_data_list[0]["field_name"],
                    "group2_name": stats_data_list[1]["field_name"],
                    "group1_size": int(n1),
                    "group2_size": int(n2),
                    "group1_data": group1_data,
                    "group2_data": group2_data
                },
                "group_statistics": {
                    "group1_mean": float(mean1),
                    "group1_std": float(std1),
                    "group1_min": float(np.min(group1_data)),
                    "group1_max": float(np.max(group1_data)),
                    "group2_mean": float(mean2),
                    "group2_std": float(std2),
                    "group2_min": float(np.min(group2_data)),
                    "group2_max": float(np.max(group2_data)),
                    "pooled_std": 0.0
                },
                "test_statistics": {
                    "degrees_of_freedom": n1 + n2 - 2,
                    "t_value": 0.0,
                    "standard_error": 0.0,
                    "p_value_two_sided": 1.0
                },
                "significance_tests": {
                    "p_greater_than_0_05": True,
                    "p_greater_than_0_01": True,
                    "significant_at_05": "不显著",
                    "significant_at_01": "不显著"
                },
                "interpretation": {
                    "t_test_interpretation": "独立样本t检验不显著 (p>0.05)",
                    "t_test_interpretation_01": "在0.01水平下不显著 (p>0.01)"
                }
            }
        else:
            # 均值不同但标准差都为0，这在实际中几乎不可能，但理论上t值为无穷大
            return {
                "input_parameters": {
                    "group1_name": stats_data_list[0]["field_name"],
                    "group2_name": stats_data_list[1]["field_name"],
                    "group1_size": int(n1),
                    "group2_size": int(n2),
                    "group1_data": group1_data,
                    "group2_data": group2_data
                },
                "group_statistics": {
                    "group1_mean": float(mean1),
                    "group1_std": float(std1),
                    "group1_min": float(np.min(group1_data)),
                    "group1_max": float(np.max(group1_data)),
                    "group2_mean": float(mean2),
                    "group2_std": float(std2),
                    "group2_min": float(np.min(group2_data)),
                    "group2_max": float(np.max(group2_data)),
                    "pooled_std": 0.0
                },
                "test_statistics": {
                    "degrees_of_freedom": n1 + n2 - 2,
                    "t_value": float('inf') if mean_diff > 0 else float('-inf'),
                    "standard_error": 0.0,
                    "p_value_two_sided": 0.0
                },
                "significance_tests": {
                    "p_greater_than_0_05": False,
                    "p_greater_than_0_01": False,
                    "significant_at_05": "显著",
                    "significant_at_01": "显著"
                },
                "interpretation": {
                    "t_test_interpretation": "独立样本t检验显著 (p≤0.05)",
                    "t_test_interpretation_01": "在0.01水平下显著 (p≤0.01)"
                }
            }
    
    # 如果只有一组标准差为0，这在实际中也很少见，但我们仍需处理
    if std1 <= 0 or std2 <= 0:
        raise ValueError("各组标准差必须大于0")
    
    # 计算合并标准差（假设方差齐性）
    pooled_std = math.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
    
    # 计算标准误
    standard_error = pooled_std * math.sqrt(1/n1 + 1/n2)
    
    # 计算t统计量
    t_statistic = (mean1 - mean2) / standard_error
    
    # 计算自由度
    degrees_of_freedom = n1 + n2 - 2
    
    # 计算双侧p值
    # P值 = 2 × P(T ≥ |t|)
    p_value_two_sided = 2 * (1 - t.cdf(abs(t_statistic), df=degrees_of_freedom))
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "group1_name": stats_data_list[0]["field_name"],
            "group2_name": stats_data_list[1]["field_name"],
            "group1_size": int(n1),
            "group2_size": int(n2),
            "group1_data": group1_data,
            "group2_data": group2_data
        },
        "group_statistics": {
            "group1_mean": float(mean1),
            "group1_std": float(std1),
            "group1_min": float(np.min(group1_data)),
            "group1_max": float(np.max(group1_data)),
            "group2_mean": float(mean2),
            "group2_std": float(std2),
            "group2_min": float(np.min(group2_data)),
            "group2_max": float(np.max(group2_data)),
            "pooled_std": float(pooled_std)
        },
        "test_statistics": {
            "degrees_of_freedom": degrees_of_freedom,
            "t_value": float(t_statistic),
            "standard_error": float(standard_error),
            "p_value_two_sided": float(p_value_two_sided)
        },
        "significance_tests": {
            "p_greater_than_0_05": is_greater_than_05,
            "p_greater_than_0_01": is_greater_than_01,
            "significant_at_05": "不显著" if is_greater_than_05 else "显著",
            "significant_at_01": "不显著" if is_greater_than_01 else "显著"
        },
        "interpretation": {
            "t_test_interpretation": f"独立样本t检验{'不' if p_value_two_sided > 0.05 else ''}显著 (p={'>' if p_value_two_sided > 0.05 else '≤'}0.05)",
            "t_test_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided > 0.01 else ''}显著 (p={'>' if p_value_two_sided > 0.01 else '≤'}0.01)"
        }
    }


def perform_independent_t_test(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    根据给定参数执行独立资料t检验
    
    参数:
    - stats_data_list: 包含两组数据的列表
    
    返回:
    - 包含独立t检验完整结果的字典
    """
    # 验证参数
    if not stats_data_list or len(stats_data_list) != 2:
        raise ValueError("独立资料t检验需要恰好2组数据")
    
    # 执行独立t检验
    results = independent_t_test_from_stats_data(stats_data_list)
    
    return results


def cal_result_t_indep(stats_data_list: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    生成独立资料t检验分析的完整报告字典
    
    此函数整合了独立资料t检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的t检验结果。报告包括输入参数、
    组统计量、检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解t检验的特征。
    
    Args:
        stats_data_list: 包含两组数据的列表，每组数据是一个字典，包含：
            - field_name: 组名，字符串类型
            - data_list: 该组数据列表，浮点数列表类型
    
    Returns:
        Dict[str, Any]: 包含t检验分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"独立资料t检验分析"
            - input_parameters: 输入参数信息
                - group1_name: 组1名称，字符串类型
                - group2_name: 组2名称，字符串类型
                - group1_size: 组1样本量，整数类型
                - group2_size: 组2样本量，整数类型
            - group_statistics: 组统计量
                - group1_mean: 组1均数，浮点数类型
                - group1_std: 组1标准差，浮点数类型
                - group1_min: 组1最小值，浮点数类型
                - group1_max: 组1最大值，浮点数类型
                - group2_mean: 组2均数，浮点数类型
                - group2_std: 组2标准差，浮点数类型
                - group2_min: 组2最小值，浮点数类型
                - group2_max: 组2最大值，浮点数类型
                - pooled_std: 合并标准差，浮点数类型
            - test_statistics: 检验统计量
                - degrees_of_freedom: 自由度，整数类型
                - t_value: t值，浮点数类型
                - standard_error: 标准误，浮点数类型
                - p_value_two_sided: 双侧P值，浮点数类型
            - significance_tests: 显著性检验结果
                - p_greater_than_0_05: P值是否大于0.05，布尔类型
                - p_greater_than_0_01: P值是否大于0.01，布尔类型
                - significant_at_05: 0.05水平显著性，字符串类型
                - significant_at_01: 0.01水平显著性，字符串类型
            - interpretation: 统计解释
                - t_test_interpretation: t检验解释，字符串类型
                - t_test_interpretation_01: 0.01水平解释，字符串类型
            - remark: 备注信息，字符串类型
    """
    # 执行独立t检验
    results = perform_independent_t_test(stats_data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "独立资料t检验分析",
        "input_parameters": results["input_parameters"],
        "group_statistics": results["group_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"组1: {stats_data_list[0]['field_name']} (n={len(stats_data_list[0]['data_list'])}), "
                  f"组2: {stats_data_list[1]['field_name']} (n={len(stats_data_list[1]['data_list'])})"
    }
    
    return result_dict