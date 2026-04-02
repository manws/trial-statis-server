"""
基于原始数据的独立样本Z检验统计分析模块 (Independent Samples Z-Test Statistics Module - Based on Raw Data)

本模块提供基于原始数据的独立样本Z检验功能，用于检验两个独立样本的均数是否有显著性差异。
该检验适用于两组数据相互独立，总体标准差已知且仅有原始数据的情况，模块将自动计算样本统计量。

【模块功能概述】:
1. 样本统计量计算：计算两组样本的均数、标准差等描述性统计量
2. Z统计量计算：计算独立样本Z检验的统计量
3. P值计算：基于标准正态分布计算双侧P值
4. 显著性检验：判断统计结果的显著性
5. 结果解释：提供统计结果的专业解释

【临床应用价值】:
- 比较两种治疗方法的效果差异：如比较两组患者的治疗后效果
- 比较不同人群的生理指标：如比较男女的血压水平
- 比较实验组和对照组的差异：如比较用药组和安慰剂组
- 评估不同诊断方法的一致性：如比较两种检测方法的结果

【统计方法选择指南】:
1. 独立样本Z检验适用条件：
   - 两组数据相互独立
   - 总体标准差已知
   - 样本量较大（通常n≥30）
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
4. 需注意Z检验对总体标准差已知的假设

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用独立样本Z检验吗？"
- "如何解释独立样本Z检验的P值？"
- "Z值的含义是什么？"
- "独立样本Z检验的前提条件是什么？"
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
from scipy.stats import norm
import numpy as np
from app.schemas.request_data.z_param import ZParamIndep2


def independent_sample_z_test_from_raw_data(raw_data_group1: List[float], raw_data_group2: List[float]) -> Dict:
    """
    基于原始数据的独立样本Z检验
    
    参数:
    - raw_data_group1: 第一组原始数据列表
    - raw_data_group2: 第二组原始数据列表
    
    返回:
    - 包含Z检验统计量和结果的字典
    """
    # 参数验证
    if not raw_data_group1 or len(raw_data_group1) == 0:
        raise ValueError("第一组原始数据不能为空")
    
    if not raw_data_group2 or len(raw_data_group2) == 0:
        raise ValueError("第二组原始数据不能为空")
    
    # 计算两组样本统计量
    n1 = len(raw_data_group1)
    n2 = len(raw_data_group2)
    
    mean1 = np.mean(raw_data_group1)
    mean2 = np.mean(raw_data_group2)
    
    std1 = np.std(raw_data_group1, ddof=1)  # 样本标准差
    std2 = np.std(raw_data_group2, ddof=1)  # 样本标准差
    
    # 计算两组均数差
    mean_diff = mean1 - mean2
    
    # 计算标准误
    standard_error = math.sqrt((std1**2 / n1) + (std2**2 / n2))
    
    # 处理标准误为0的情况（两组数据完全相同）
    if standard_error == 0:
        if mean_diff == 0:
            # 两组完全相同
            z_statistic = 0.0
            p_value_two_sided = 1.0
        else:
            # 均数不同但标准差为0（理论上不可能，但防止除零错误）
            z_statistic = float('inf') if mean_diff > 0 else float('-inf')
            p_value_two_sided = 0.0
    else:
        # 计算Z统计量
        z_statistic = mean_diff / standard_error
        
        # 计算双侧P值
        # P值 = 2 × P(Z ≥ |z|)
        p_value_two_sided = 2 * (1 - norm.cdf(abs(z_statistic)))
    
    # 显著性判断 (> 0.05, > 0.01)
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "sample1_n": int(n1),
            "sample1_mean": float(mean1),
            "sample1_std": float(std1),
            "sample2_n": int(n2),
            "sample2_mean": float(mean2),
            "sample2_std": float(std2)
        },
        "sample_statistics": {
            "sample1_min": float(np.min(raw_data_group1)),
            "sample1_max": float(np.max(raw_data_group1)),
            "sample2_min": float(np.min(raw_data_group2)),
            "sample2_max": float(np.max(raw_data_group2))
        },
        "test_statistics": {
            "mean_difference": float(mean_diff),
            "standard_error": float(standard_error),
            "z_value": float(z_statistic),
            "p_value_two_sided": float(p_value_two_sided)
        },
        "significance_tests": {
            "p_greater_than_0_05": is_greater_than_05,
            "p_greater_than_0_01": is_greater_than_01,
            "significant_at_05": "不显著" if is_greater_than_05 else "显著",
            "significant_at_01": "不显著" if is_greater_than_01 else "显著"
        },
        "interpretation": {
            "z_test_interpretation": f"独立样本Z检验{'不' if p_value_two_sided > 0.05 else ''}显著 (p={'>' if p_value_two_sided > 0.05 else '≤'}0.05)",
            "z_test_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided > 0.01 else ''}显著 (p={'>' if p_value_two_sided > 0.01 else '≤'}0.01)"
        }
    }


def perform_independent_sample_z_test_from_raw(raw_data_group1: List[float], raw_data_group2: List[float]) -> Dict:
    """
    根据给定参数执行基于原始数据的独立样本Z检验
    
    参数:
    - raw_data_group1: 第一组原始数据列表
    - raw_data_group2: 第二组原始数据列表
    
    返回:
    - 包含独立样本Z检验完整结果的字典
    """
    # 验证参数
    if not raw_data_group1:
        raise ValueError("第一组原始数据不能为空")
    
    if not raw_data_group2:
        raise ValueError("第二组原始数据不能为空")
    
    # 执行独立样本Z检验
    results = independent_sample_z_test_from_raw_data(raw_data_group1, raw_data_group2)
    
    return results


def cal_result_z_indep2(param: ZParamIndep2) -> Dict[str, Any]:
    """
    生成基于原始数据的独立样本Z检验分析的完整报告字典
    
    Args:
        param: ZParamIndep2对象，包含stats_data_list
    
    Returns:
        Dict[str, Any]: 包含Z检验分析指标的字典
    """
    raw_data_group1 = param.stats_data_list[0].data_list
    raw_data_group2 = param.stats_data_list[1].data_list

    # 执行独立样本Z检验
    results = perform_independent_sample_z_test_from_raw(raw_data_group1, raw_data_group2)
    
    # 构建结果字典
    result_dict = {
        "table_name": "独立样本Z检验分析（原始资料）",
        "input_parameters": results["input_parameters"],
        "sample_statistics": results["sample_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"组1数据点数: {len(raw_data_group1)}, 组2数据点数: {len(raw_data_group2)}"
    }
    
    return result_dict