"""
独立样本Z检验统计分析模块 (Independent Samples Z-Test Statistics Module)

本模块提供独立样本Z检验功能，用于检验两个独立样本的均数是否有显著性差异。
该检验适用于两组数据相互独立，总体标准差已知且样本量较大的情况。

【模块功能概述】:
1. Z统计量计算：计算独立样本Z检验的统计量
2. P值计算：基于标准正态分布计算双侧P值
3. 显著性检验：判断统计结果的显著性
4. 结果解释：提供统计结果的专业解释

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
from typing import Dict, Any
from scipy.stats import norm
from app.schemas.request_data.z_param import ZParamIndep1


def independent_sample_z_test(n1: int, mean1: float, std1: float, n2: int, mean2: float, std2: float) -> Dict[str, Any]:
    """
    独立样本Z检验
    
    参数:
    - n1: 第一组样本量
    - mean1: 第一组样本均数
    - std1: 第一组样本标准差
    - n2: 第二组样本量
    - mean2: 第二组样本均数
    - std2: 第二组样本标准差
    
    返回:
    - 包含独立样本Z检验统计量和结果的字典
    """
    # 参数验证
    if n1 <= 0 or n2 <= 0:
        raise ValueError("样本量必须大于0")
    
    if std1 <= 0 or std2 <= 0:
        raise ValueError("样本标准差必须大于0")
    
    # 计算Z统计量
    # Z = (mean1 - mean2) / sqrt((std1^2/n1) + (std2^2/n2))
    standard_error = math.sqrt((std1**2 / n1) + (std2**2 / n2))
    z_statistic = (mean1 - mean2) / standard_error
    
    # 计算双侧p值
    # P值 = 2 × P(Z ≥ |z|)
    p_value_two_sided = 2 * (1 - norm.cdf(abs(z_statistic)))
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "group1_n": int(n1),
            "group1_mean": float(mean1),
            "group1_std": float(std1),
            "group2_n": int(n2),
            "group2_mean": float(mean2),
            "group2_std": float(std2)
        },
        "test_statistics": {
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


def perform_independent_sample_z_test(n1: int, mean1: float, std1: float, n2: int, mean2: float, std2: float) -> Dict[str, Any]:
    """
    根据给定参数执行独立样本Z检验
    
    本函数封装了独立样本Z检验的核心计算逻辑，包括参数验证、统计量计算
    和结果整理。适用于需要快速获得Z检验结果的场景。
    
    【统计学原理】:
    独立样本Z检验基于以下公式计算Z统计量：
    Z = (mean1 - mean2) / sqrt((std1²/n1) + (std2²/n2))
    
    其中：
    - mean1, mean2: 两组的样本均数
    - std1, std2: 两组的样本标准差（作为总体标准差的估计）
    - n1, n2: 两组的样本量
    
    【临床应用场景】:
    1. 药物临床试验：比较新药组与安慰剂组的主要疗效指标
    2. 流行病学研究：比较暴露组与非暴露组的某项生理指标
    3. 诊断试验评价：比较两种诊断方法的检测结果
    4. 健康筛查：比较不同人群的某项健康指标
    
    【使用注意事项】:
    - 确保两组数据相互独立，不存在配对或重复测量关系
    - 样本量应足够大（建议每组n≥30），以保证Z检验的近似有效性
    - 当总体标准差未知时，应考虑使用独立样本t检验
    - 需检查数据是否满足正态性假设，或通过大样本保证中心极限定理适用
    
    Args:
        n1: 第一组样本量，必须为正整数
        mean1: 第一组样本均数，实数
        std1: 第一组样本标准差，必须为正实数
        n2: 第二组样本量，必须为正整数
        mean2: 第二组样本均数，实数
        std2: 第二组样本标准差，必须为正实数
    
    Returns:
        Dict[str, Any]: 包含Z检验结果的字典，结构如下：
            - input_parameters: 输入参数
            - test_statistics: 检验统计量（标准误、Z值、双侧P值）
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    
    Raises:
        ValueError: 当样本量≤0或标准差≤0时抛出
    """
    # 验证参数
    if n1 <= 0 or n2 <= 0:
        raise ValueError("样本量必须大于0")
    
    if std1 <= 0 or std2 <= 0:
        raise ValueError("样本标准差必须大于0")
    
    # 执行独立样本Z检验
    results = independent_sample_z_test(n1, mean1, std1, n2, mean2, std2)
    
    return results


def cal_result_z_indep1(param: ZParamIndep1) -> Dict[str, Any]:
    """
    生成独立样本Z检验分析的完整报告字典
    
    Args:
        param: ZParamIndep1对象，包含n1, mean1, std1, n2, mean2, std2
    
    Returns:
        Dict[str, Any]: 包含Z检验分析指标的字典
    """
    n1 = param.n1
    mean1 = param.mean1
    std1 = param.std1
    n2 = param.n2
    mean2 = param.mean2
    std2 = param.std2

    # 执行独立样本Z检验
    results = perform_independent_sample_z_test(n1, mean1, std1, n2, mean2, std2)
    
    # 构建结果字典
    result_dict = {
        "table_name": "独立样本Z检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"第一组(n={n1}, 均数={mean1:.4f}, 标准差={std1:.4f}) vs 第二组(n={n2}, 均数={mean2:.4f}, 标准差={std2:.4f})"
    }
    
    return result_dict