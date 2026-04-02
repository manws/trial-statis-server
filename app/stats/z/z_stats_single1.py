"""
单样本Z检验统计分析模块 (Single Sample Z-Test Statistics Module)

本模块提供单样本Z检验功能，用于检验样本均数与已知总体均数之间是否存在显著性差异。
该检验适用于总体标准差已知且样本来自正态分布总体的情况。

【模块功能概述】:
1. Z统计量计算：计算单样本Z检验的统计量
2. P值计算：基于标准正态分布计算双侧P值
3. 显著性检验：判断统计结果的显著性
4. 结果解释：提供统计结果的专业解释

【临床应用价值】:
- 检验样本均数与标准值的差异：如检验某地区人群身高是否符合国家标准
- 检验治疗前后变化：如检验某种药物对血压的影响
- 质量控制：检验生产过程是否符合标准
- 诊断试验：检验某项指标是否偏离正常值

【统计方法选择指南】:
1. 单样本Z检验适用条件：
   - 数据来自正态分布总体
   - 总体标准差已知
   - 样本量较大（通常n≥30）
   - 数据独立

2. 临床应用场景：
   - 检验某项生化指标是否偏离正常范围
   - 评估干预措施的效果
   - 比较样本与历史对照的差异
   - 验证实验室检测结果的准确性

【结果解读注意事项】:
1. P值解释：在零假设（样本均数与总体均数无差异）成立的情况下，观察到当前或更极端结果的概率
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
- "我的数据适合用单样本Z检验吗？"
- "如何解释Z检验的P值？"
- "Z值的含义是什么？"
- "单样本Z检验的前提条件是什么？"
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


def single_sample_z_test(pop_mean: float, pop_std: float, n: int, sample_mean: float) -> Dict:
    """
    单样本Z检验
    
    参数:
    - pop_mean: 总体均数
    - pop_std: 总体标准差
    - n: 样本量
    - sample_mean: 样本均数
    
    返回:
    - 包含Z检验统计量和结果的字典
    """
    # 参数验证
    if n <= 0:
        raise ValueError("样本量必须大于0")
    
    if pop_std <= 0:
        raise ValueError("总体标准差必须大于0")
    
    # 计算标准误
    standard_error = pop_std / math.sqrt(n)
    
    # 计算Z统计量
    z_statistic = (sample_mean - pop_mean) / standard_error
    
    # 计算双侧P值
    # P值 = 2 × P(Z ≥ |z|)
    p_value_two_sided = 2 * (1 - norm.cdf(abs(z_statistic)))
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "population_mean": float(pop_mean),
            "population_std": float(pop_std),
            "sample_size": int(n),
            "sample_mean": float(sample_mean)
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
            "z_test_interpretation": f"单样本Z检验{'不' if p_value_two_sided > 0.05 else ''}显著 (p={'>' if p_value_two_sided > 0.05 else '≤'}0.05)",
            "z_test_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided > 0.01 else ''}显著 (p={'>' if p_value_two_sided > 0.01 else '≤'}0.01)"
        }
    }


def perform_single_sample_z_test(pop_mean: float, pop_std: float, n: int, sample_mean: float) -> Dict:
    """
    根据给定参数执行单样本Z检验
    
    参数:
    - pop_mean: 总体均数
    - pop_std: 总体标准差
    - n: 样本量
    - sample_mean: 样本均数
    
    返回:
    - 包含单样本Z检验完整结果的字典
    """
    # 验证参数
    if n <= 0:
        raise ValueError("样本量必须大于0")
    
    if pop_std <= 0:
        raise ValueError("总体标准差必须大于0")
    
    # 执行单样本Z检验
    results = single_sample_z_test(pop_mean, pop_std, n, sample_mean)
    
    return results


def cal_result_z_single1(pop_mean: float, pop_std: float, n: int, sample_mean: float) -> Dict[str, Any]:
    """
    生成单样本Z检验分析的完整报告字典
    
    此函数整合了单样本Z检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的Z检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解Z检验的特征。
    
    Args:
        pop_mean: 总体均数，浮点数类型
        pop_std: 总体标准差，浮点数类型
        n: 样本量，整数类型
        sample_mean: 样本均数，浮点数类型
    
    Returns:
        Dict[str, Any]: 包含Z检验分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"单样本Z检验分析"
            - input_parameters: 输入参数信息
                - population_mean: 总体均数，浮点数类型
                - population_std: 总体标准差，浮点数类型
                - sample_size: 样本量，整数类型
                - sample_mean: 样本均数，浮点数类型
            - test_statistics: 检验统计量
                - standard_error: 标准误，浮点数类型
                - z_value: Z值，浮点数类型
                - p_value_two_sided: 双侧P值，浮点数类型
            - significance_tests: 显著性检验结果
                - p_greater_than_0_05: P值是否大于0.05，布尔类型
                - p_greater_than_0_01: P值是否大于0.01，布尔类型
                - significant_at_05: 0.05水平显著性，字符串类型
                - significant_at_01: 0.01水平显著性，字符串类型
            - interpretation: 统计解释
                - z_test_interpretation: Z检验解释，字符串类型
                - z_test_interpretation_01: 0.01水平解释，字符串类型
            - remark: 备注信息，字符串类型
    """
    # 执行单样本Z检验
    results = perform_single_sample_z_test(pop_mean, pop_std, n, sample_mean)
    
    # 构建结果字典
    result_dict = {
        "table_name": "单样本Z检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"总体均数: {pop_mean}, 总体标准差: {pop_std}, 样本量: {n}, 样本均数: {sample_mean}"
    }
    
    return result_dict