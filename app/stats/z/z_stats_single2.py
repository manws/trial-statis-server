"""
基于原始数据的单样本Z检验统计分析模块 (Single Sample Z-Test Statistics Module - Based on Raw Data)

本模块提供基于原始数据的单样本Z检验功能，用于检验样本均数与已知总体均数之间是否存在显著性差异。
该检验适用于总体标准差已知且仅有原始数据的情况，模块将自动计算样本统计量。

【模块功能概述】:
1. 样本统计量计算：计算样本均数、标准差等描述性统计量
2. Z统计量计算：计算单样本Z检验的统计量
3. P值计算：基于标准正态分布计算双侧P值
4. 显著性检验：判断统计结果的显著性
5. 结果解释：提供统计结果的专业解释

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
from typing import Dict, List, Any
from scipy.stats import norm
import numpy as np
from app.schemas.request_data.z_param import ZParamSingle2


def single_sample_z_test_from_raw_data(pop_mean: float, pop_std: float, raw_data: List[float]) -> Dict:
    """
    基于原始数据的单样本Z检验核心计算函数
    
    本函数执行单样本Z检验的核心统计计算，包括样本统计量计算、Z统计量计算、
    P值计算和显著性判断。适用于总体标准差已知的情况。
    
    【统计学原理】:
    单样本Z检验用于检验样本均数与已知总体均数之间是否存在显著性差异。
    当总体标准差已知时，使用Z统计量：
    Z = (样本均数 - 总体均数) / (总体标准差 / √样本量)
    
    【计算公式】:
    - 标准误 (SE) = σ / √n
    - Z统计量 = (X̄ - μ) / SE
    - 双侧P值 = 2 × P(Z ≥ |z|)
    
    【临床应用示例】:
    例：已知某地区成年男性身高总体均数为170cm，标准差为8cm。
    现随机抽取100名成年男性，测得平均身高为172cm。
    问：该地区成年男性身高是否与已知总体均数有显著差异？
    
    【结果解读】:
    - P > 0.05：差异无统计学意义，尚不能认为样本均数与总体均数不同
    - P ≤ 0.05：差异有统计学意义，可以认为样本均数与总体均数不同
    - P ≤ 0.01：差异有高度统计学意义
    
    参数:
        pop_mean (float): 总体均数，已知的总体平均值
        pop_std (float): 总体标准差，已知的总体标准差（必须>0）
        raw_data (List[float]): 原始数据列表，包含所有样本观测值
    
    返回:
        Dict: 包含Z检验统计量和结果的字典，结构如下:
            - input_parameters: 输入参数信息
            - sample_statistics: 样本统计量
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    
    异常:
        ValueError: 当原始数据为空或总体标准差≤0时抛出
    
    注意事项:
        1. 本函数假设数据来自正态分布总体
        2. 要求总体标准差已知
        3. 建议样本量≥30以获得更可靠的结果
        4. 数据应相互独立
    """
    # 参数验证
    if not raw_data or len(raw_data) == 0:
        raise ValueError("原始数据不能为空")
    
    if pop_std <= 0:
        raise ValueError("总体标准差必须大于0")
    
    # 计算样本统计量
    n = len(raw_data)
    sample_mean = np.mean(raw_data)
    sample_std = np.std(raw_data, ddof=1)  # 样本标准差
    
    # 计算标准误（使用总体标准差）
    standard_error = pop_std / math.sqrt(n)
    
    # 计算Z统计量
    z_statistic = (sample_mean - pop_mean) / standard_error
    
    # 计算双侧P值
    # P值 = 2 × P(Z ≥ |z|)
    p_value_two_sided = 2 * (1 - norm.cdf(abs(z_statistic)))
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    # 参数验证
    if not raw_data or len(raw_data) == 0:
        raise ValueError("原始数据不能为空")
    
    if pop_std <= 0:
        raise ValueError("总体标准差必须大于0")
    
    # 计算样本统计量
    n = len(raw_data)
    sample_mean = np.mean(raw_data)
    sample_std = np.std(raw_data, ddof=1)  # 样本标准差
    
    # 计算标准误（使用总体标准差）
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
            "sample_size": int(n)
        },
        "sample_statistics": {
            "sample_mean": float(sample_mean),
            "sample_std": float(sample_std),
            "sample_min": float(np.min(raw_data)),
            "sample_max": float(np.max(raw_data))
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


def perform_single_sample_z_test_from_raw(pop_mean: float, pop_std: float, raw_data: List[float]) -> Dict:
    """
    执行基于原始数据的单样本Z检验
    
    本函数是对核心计算函数的封装，提供参数验证和调用接口。
    适用于需要直接传入参数进行单样本Z检验的场景。
    
    【函数职责】:
    1. 验证输入参数的有效性
    2. 调用核心计算函数执行统计检验
    3. 返回标准化的结果字典
    
    【参数说明】:
    - pop_mean: 已知的总体均数，用于与样本均数比较
    - pop_std: 已知的总体标准差，必须为正数
    - raw_data: 原始观测数据列表，不能为空
    
    【返回值结构】:
    返回字典包含以下键:
    - input_parameters: 输入的参数信息
    - sample_statistics: 计算得到的样本统计量
    - test_statistics: Z检验的统计量结果
    - significance_tests: 显著性检验的判断结果
    - interpretation: 对结果的统计学解释
    
    【使用示例】:
    ```python
    # 示例：检验某班级学生平均身高是否与全国平均水平有差异
    pop_mean = 170.0  # 全国平均身高
    pop_std = 8.0     # 全国身高标准差
    sample_data = [172, 168, 175, 170, 169, ...]  # 班级学生身高数据
    
    result = perform_single_sample_z_test_from_raw(pop_mean, pop_std, sample_data)
    print(f"Z值：{result['test_statistics']['z_value']}")
    print(f"P值：{result['test_statistics']['p_value_two_sided']}")
    ```
    
    参数:
        pop_mean (float): 总体均数
        pop_std (float): 总体标准差
        raw_data (List[float]): 原始数据列表
    
    返回:
        Dict: 包含单样本Z检验完整结果的字典
    
    异常:
        ValueError: 当原始数据为空或总体标准差≤0时抛出
    """
    # 验证参数
    if not raw_data or len(raw_data) == 0:
        raise ValueError("原始数据不能为空")
    
    if pop_std <= 0:
        raise ValueError("总体标准差必须大于0")
    
    # 执行单样本Z检验
    results = single_sample_z_test_from_raw_data(pop_mean, pop_std, raw_data)
    
    return results


def cal_result_z_single2(param: ZParamSingle2) -> Dict[str, Any]:
    """
    生成基于原始数据的单样本Z检验分析的完整报告字典
    
    Args:
        param: ZParamSingle2对象，包含pop_mean, pop_std, stats_data_list
    
    Returns:
        Dict[str, Any]: 包含Z检验分析指标的字典
    """
    pop_mean = param.pop_mean
    pop_std = param.pop_std
    raw_data = param.stats_data_list[0].data_list

    # 执行单样本Z检验
    results = perform_single_sample_z_test_from_raw(pop_mean, pop_std, raw_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "单样本Z检验分析（原始资料）",
        "input_parameters": results["input_parameters"],
        "sample_statistics": results["sample_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"总体均数：{pop_mean}, 总体标准差：{pop_std}, 数据点总数：{len(raw_data)}"
    }
    
    return result_dict