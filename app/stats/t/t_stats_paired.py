"""
配对样本t检验统计分析模块 (Paired Samples T-Test Statistics Module)

本模块提供配对样本t检验功能，用于检验配对数据的均数是否有显著性差异。
该检验适用于配对设计数据，如同一受试者在不同时间点的测量值、或配对的受试者之间的对比。

【模块功能概述】:
1. 差值计算：计算配对数据的差值
2. t统计量计算：计算配对样本t检验的统计量
3. P值计算：基于t分布计算双侧P值
4. 显著性检验：判断统计结果的显著性
5. 结果解释：提供统计结果的专业解释

【临床应用价值】:
- 比较同一受试者治疗前后的差异：如比较患者用药前后的血压变化
- 比较配对受试者的不同处理效果：如比较同卵双胞胎接受不同处理的反应
- 比较两种测量方法的一致性：如比较两种仪器对同一标本的测量结果
- 比较病例-对照匹配研究：如比较疾病患者与其匹配对照的指标

【统计方法选择指南】:
1. 配对样本t检验适用条件：
   - 数据为配对设计
   - 差值来自正态分布总体
   - 样本量较小（通常n<30）
   - 配对数据独立

2. 临床应用场景：
   - 治疗前后比较：评估干预措施的有效性
   - 交叉试验：比较两种治疗方法的效果
   - 重复测量：同一受试者在不同时间点的比较
   - 方法学比较：比较两种检测方法的一致性

【结果解读注意事项】:
1. P值解释：在零假设（配对均数无差异）成立的情况下，观察到当前或更极端结果的概率
2. 显著性水平：通常使用α=0.05作为判断标准
3. 统计显著性不代表临床重要性
4. 需注意t检验对差值正态分布的假设

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用配对样本t检验吗？"
- "如何解释配对样本t检验的P值？"
- "配对t检验与独立样本t检验的区别是什么？"
- "配对样本t检验的前提条件是什么？"
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


def paired_t_test_from_stats_data(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    基于StatsData列表的配对资料t检验
    
    参数:
    - stats_data_list: 包含两组配对数据的列表，每组数据是一个字典，包含field_name和data_list字段
    
    返回:
    - 包含配对t检验统计量和结果的字典
    """
    # 参数验证
    if not stats_data_list or len(stats_data_list) != 2:
        raise ValueError("配对资料t检验需要恰好2个变量（列）")
    
    # 提取两个变量的数据
    var1_data = stats_data_list[0]["data_list"]
    var2_data = stats_data_list[1]["data_list"]
    
    # 验证数据长度一致
    if len(var1_data) != len(var2_data):
        raise ValueError(f"两个变量的数据长度不一致: {len(var1_data)} vs {len(var2_data)}")
    
    if len(var1_data) < 2:
        raise ValueError("每组数据至少需要2个观测值")
    
    # 计算差值
    differences = np.array(var2_data) - np.array(var1_data)
    
    # 计算差值的统计量
    n_pairs = len(differences)
    mean_diff = np.mean(differences)
    std_diff = np.std(differences, ddof=1)  # 使用样本标准差
    
    # 特殊情况处理：如果所有差值都相同（标准差为0）
    if std_diff == 0:
        if mean_diff == 0:
            # 完全相同的两组数据
            return {
                "input_parameters": {
                    "variable1_name": stats_data_list[0]["field_name"],
                    "variable2_name": stats_data_list[1]["field_name"],
                    "number_of_pairs": int(n_pairs),
                    "variable1_data": var1_data,
                    "variable2_data": var2_data
                },
                "difference_statistics": {
                    "mean_difference": float(mean_diff),
                    "std_difference": float(std_diff),
                    "min_difference": float(np.min(differences)),
                    "max_difference": float(np.max(differences))
                },
                "test_statistics": {
                    "degrees_of_freedom": n_pairs - 1,
                    "t_value": 0.0,
                    "p_value_two_sided": 1.0
                },
                "significance_tests": {
                    "p_greater_than_0_05": True,
                    "p_greater_than_0_01": True,
                    "significant_at_05": "不显著",
                    "significant_at_01": "不显著"
                },
                "interpretation": {
                    "t_test_interpretation": "配对t检验不显著 (p>0.05)",
                    "t_test_interpretation_01": "在0.01水平下不显著 (p>0.01)"
                }
            }
        else:
            # 差值恒定但不为0，这是理论上不可能的情况，但在数值计算中可能出现
            # 在这种情况下，我们可以认为差异非常显著
            return {
                "input_parameters": {
                    "variable1_name": stats_data_list[0]["field_name"],
                    "variable2_name": stats_data_list[1]["field_name"],
                    "number_of_pairs": int(n_pairs),
                    "variable1_data": var1_data,
                    "variable2_data": var2_data
                },
                "difference_statistics": {
                    "mean_difference": float(mean_diff),
                    "std_difference": float(std_diff),
                    "min_difference": float(np.min(differences)),
                    "max_difference": float(np.max(differences))
                },
                "test_statistics": {
                    "degrees_of_freedom": n_pairs - 1,
                    "t_value": float('inf') if mean_diff > 0 else float('-inf'),
                    "p_value_two_sided": 0.0
                },
                "significance_tests": {
                    "p_greater_than_0_05": False,
                    "p_greater_than_0_01": False,
                    "significant_at_05": "显著",
                    "significant_at_01": "显著"
                },
                "interpretation": {
                    "t_test_interpretation": "配对t检验显著 (p≤0.05)",
                    "t_test_interpretation_01": "在0.01水平下显著 (p≤0.01)"
                }
            }
    
    # 计算自由度
    degrees_of_freedom = n_pairs - 1
    
    # 计算t统计量
    # t = mean_diff / (std_diff / sqrt(n_pairs))
    t_statistic = mean_diff / (std_diff / math.sqrt(n_pairs))
    
    # 计算双侧p值
    # P值 = 2 × P(T ≥ |t|)
    p_value_two_sided = 2 * (1 - t.cdf(abs(t_statistic), df=degrees_of_freedom))
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "variable1_name": stats_data_list[0]["field_name"],
            "variable2_name": stats_data_list[1]["field_name"],
            "number_of_pairs": int(n_pairs),
            "variable1_data": var1_data,
            "variable2_data": var2_data
        },
        "difference_statistics": {
            "mean_difference": float(mean_diff),
            "std_difference": float(std_diff),
            "min_difference": float(np.min(differences)),
            "max_difference": float(np.max(differences))
        },
        "test_statistics": {
            "degrees_of_freedom": degrees_of_freedom,
            "t_value": float(t_statistic),
            "p_value_two_sided": float(p_value_two_sided)
        },
        "significance_tests": {
            "p_greater_than_0_05": is_greater_than_05,
            "p_greater_than_0_01": is_greater_than_01,
            "significant_at_05": "不显著" if is_greater_than_05 else "显著",
            "significant_at_01": "不显著" if is_greater_than_01 else "显著"
        },
        "interpretation": {
            "t_test_interpretation": f"配对t检验{'不' if p_value_two_sided > 0.05 else ''}显著 (p={'>' if p_value_two_sided > 0.05 else '≤'}0.05)",
            "t_test_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided > 0.01 else ''}显著 (p={'>' if p_value_two_sided > 0.01 else '≤'}0.01)"
        }
    }


def perform_paired_t_test(stats_data_list: List[Dict[str, Any]]) -> Dict:
    """
    根据给定参数执行配对资料t检验
    
    参数:
    - stats_data_list: 包含两组配对数据的列表
    
    返回:
    - 包含配对t检验完整结果的字典
    """
    # 验证参数
    if not stats_data_list or len(stats_data_list) != 2:
        raise ValueError("配对资料t检验需要恰好2个变量的数据")
    
    # 执行配对t检验
    results = paired_t_test_from_stats_data(stats_data_list)
    
    return results


def cal_result_t_paired(stats_data_list: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    生成配对资料t检验分析的完整报告字典
    
    此函数整合了配对资料t检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的t检验结果。报告包括输入参数、
    差值统计量、检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解t检验的特征。
    
    Args:
        stats_data_list: 包含两组配对数据的列表，每组数据是一个字典，包含：
            - field_name: 变量名，字符串类型
            - data_list: 该变量数据列表，浮点数列表类型
    
    Returns:
        Dict[str, Any]: 包含t检验分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"配对资料t检验分析"
            - input_parameters: 输入参数信息
                - variable1_name: 变量1名称，字符串类型
                - variable2_name: 变量2名称，字符串类型
                - number_of_pairs: 配对数，整数类型
            - difference_statistics: 差值统计量
                - mean_difference: 平均差值，浮点数类型
                - std_difference: 差值标准差，浮点数类型
                - min_difference: 最小差值，浮点数类型
                - max_difference: 最大差值，浮点数类型
            - test_statistics: 检验统计量
                - degrees_of_freedom: 自由度，整数类型
                - t_value: t值，浮点数类型
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
    # 从参数对象解构为底层函数所需的 List[Dict] 格式
    stats_data_list = [item.model_dump() for item in param.stats_data_list]

    # 执行配对t检验
    results = perform_paired_t_test(stats_data_list)
    
    # 构建结果字典
    result_dict = {
        "table_name": "配对资料t检验分析",
        "input_parameters": results["input_parameters"],
        "difference_statistics": results["difference_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"变量1: {stats_data_list[0]['field_name']}, "
                  f"变量2: {stats_data_list[1]['field_name']}, "
                  f"配对数: {len(stats_data_list[0]['data_list'])}"
    }
    
    return result_dict