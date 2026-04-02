"""
基于原始数据的单样本t检验统计分析模块 (Single Sample T-Test Statistics Module - Based on Raw Data)

本模块提供基于原始数据的单样本t检验功能，用于检验样本均数与已知总体均数之间是否存在显著性差异。
该检验适用于仅有原始数据的情况，模块将自动计算样本统计量。

【模块功能概述】:
1. 样本统计量计算：计算样本均数、标准差等描述性统计量
2. t统计量计算：计算单样本t检验的统计量
3. P值计算：基于t分布计算双侧P值
4. 显著性检验：判断统计结果的显著性
5. 结果解释：提供统计结果的专业解释

【临床应用价值】:
- 检验样本均数与标准值的差异：如检验某地区人群身高是否符合国家标准
- 检验治疗前后变化：如检验某种药物对血压的影响
- 质量控制：检验生产过程是否符合标准
- 诊断试验：检验某项指标是否偏离正常值

【统计方法选择指南】:
1. 单样本t检验适用条件：
   - 数据来自正态分布总体
   - 总体标准差未知
   - 样本量较小（通常n<30）
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
4. 需注意t检验对正态分布的假设

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用单样本t检验吗？"
- "如何解释t检验的P值？"
- "t值的含义是什么？"
- "单样本t检验的前提条件是什么？"
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
from app.schemas.request_data.t_param import TParamSingle2
import numpy as np


def single_sample_t_test_from_raw_data(pop_mean: float, raw_data: List[float]) -> Dict:
    """
    基于原始数据的单样本t检验
    
    参数:
    - pop_mean: 总体均数（假设值）
    - raw_data: 原始数据列表
    
    返回:
    - 包含单样本t检验统计量和结果的字典
    """
    # 参数验证
    if len(raw_data) <= 1:
        raise ValueError("原始数据量必须大于1")
    
    # 计算样本统计量
    sample_size = len(raw_data)
    sample_mean = np.mean(raw_data)
    sample_std = np.std(raw_data, ddof=1)  # 使用样本标准差（ddof=1）
    
    if sample_std <= 0:
        raise ValueError("样本标准差必须大于0")
    
    # 计算自由度
    degrees_of_freedom = sample_size - 1
    
    # 计算t统计量
    # t = (样本均数 - 总体均数) / (样本标准差 / sqrt(n))
    t_statistic = (sample_mean - pop_mean) / (sample_std / math.sqrt(sample_size))
    
    # 计算双侧p值
    # P值 = 2 × P(T ≥ |t|)
    p_value_two_sided = 2 * (1 - t.cdf(abs(t_statistic), df=degrees_of_freedom))
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "population_mean": float(pop_mean),
            "raw_data": raw_data,
            "sample_size": int(sample_size)
        },
        "sample_statistics": {
            "sample_mean": float(sample_mean),
            "sample_std": float(sample_std),
            "sample_min": float(np.min(raw_data)),
            "sample_max": float(np.max(raw_data))
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
            "t_test_interpretation": f"单样本t检验{'不' if p_value_two_sided > 0.05 else ''}显著 (p={'>' if p_value_two_sided > 0.05 else '≤'}0.05)",
            "t_test_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided > 0.01 else ''}显著 (p={'>' if p_value_two_sided > 0.01 else '≤'}0.01)"
        }
    }


def perform_single_sample_t_test_from_raw(pop_mean: float, raw_data: List[float]) -> Dict:
    """
    根据给定参数执行基于原始数据的单样本t检验
    
    此函数封装了单样本t检验的核心计算逻辑，对输入参数进行验证后调用底层计算函数。
    它作为中间层，提供了参数验证和数据预处理功能。
    
    【统计学原理】:
    单样本t检验用于比较一个样本的均数与已知的总体均数是否有显著差异。
    当总体标准差未知且样本量较小时，使用t分布进行推断。
    
    零假设 (H0): 样本来自的总体均数等于已知的总体均数 (μ = μ0)
    备择假设 (H1): 样本来自的总体均数不等于已知的总体均数 (μ ≠ μ0)
    
    【计算公式】:
    t = (x̄ - μ0) / (s / √n)
    
    其中:
    - x̄: 样本均数
    - μ0: 假设的总体均数
    - s: 样本标准差
    - n: 样本量
    
    【临床应用场景】:
    1. 新生儿体重检验：检验某医院新生儿平均体重是否与全国标准值 (3.5kg) 有差异
    2. 血压控制评估：检验高血压患者服药后的平均收缩压是否降至目标值 (140mmHg) 以下
    3. 实验室质控：检验某批次试剂的检测结果是否与标准品标称值一致
    4. 生长发育评估：检验某地区儿童的平均身高是否达到国家同龄标准
    
    【前提条件】:
    1. 独立性：观测值之间相互独立
    2. 正态性：数据来自正态分布总体（小样本时尤为重要）
    3. 连续性：数据为连续型变量
    
    【注意事项】:
    - 样本量过小时 (<10)，正态性假设难以验证，需谨慎解读结果
    - 存在异常值时，建议使用非参数检验（如符号秩检验）
    - 大样本时 (n>100)，即使轻微偏离正态分布，t检验仍具有稳健性
    
    Args:
        pop_mean: 总体均数（假设值），浮点数类型。代表已知的或理论上的总体均数。
        raw_data: 原始数据列表，浮点数列表类型。包含所有观测值的集合。
    
    Returns:
        Dict: 包含单样本t检验完整结果的字典，结构如下:
            - input_parameters: 输入参数
            - sample_statistics: 样本统计量
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    
    Raises:
        ValueError: 当原始数据量不大于1时抛出异常
    """
    # 验证参数
    if len(raw_data) <= 1:
        raise ValueError("原始数据量必须大于1")
    
    # 执行单样本t检验
    results = single_sample_t_test_from_raw_data(pop_mean, raw_data)
    
    return results


def cal_result_t_single2(param: TParamSingle2) -> Dict[str, Any]:
    """
    生成基于原始数据的单样本t检验分析的完整报告字典
    
    此函数整合了基于原始数据的单样本t检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的t检验结果。报告包括输入参数、
    样本统计量、检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解t检验的特征。
    
    【函数定位】:
    本函数是单样本t检验分析的主入口函数，面向最终用户和报告生成系统。
    它整合了所有统计计算结果，并以统一的字典格式返回，便于后续处理和展示。
    
    【输出数据结构说明】:
    返回的字典包含以下核心部分:
    1. table_name: 报告名称，标识分析类型
    2. input_parameters: 输入参数，包括总体均数和样本量
    3. sample_statistics: 描述性统计量，包括均数、标准差、极值
    4. test_statistics: 推断性统计量，包括自由度、t值、P值
    5. significance_tests: 显著性判断结果
    6. interpretation: 专业统计解释
    7. remark: 附加备注信息
    
    【临床报告应用】:
    该函数的输出可直接用于:
    - 生成临床研究论文的统计结果部分
    - 制作临床试验报告的统计附录
    - 构建交互式统计分析报告
    - 支持医疗决策系统的统计推断模块
    
    【结果解读示例】:
    假设检验某新药降压效果，已知标准值为140mmHg:
    - 若p<0.05，可认为该药降压效果与标准值有统计学差异
    - 若样本均数<140且p<0.05，提示该药可能有降压作用
    - 需结合临床专业知识判断差异的临床意义
    
    【与其他函数的关系】:
    - 调用 perform_single_sample_t_test_from_raw() 执行核心计算
    - 整合并标准化输出格式
    - 添加报告元数据（表名、备注等）
    
    Args:
        pop_mean: 总体均数（假设值），浮点数类型。代表已知的或理论上的总体均数。
        raw_data: 原始数据列表，浮点数列表类型。包含所有观测值的集合。
    
    Returns:
        Dict[str, Any]: 包含t检验分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"单样本t检验分析（原始数据）"
            - input_parameters: 输入参数信息
                - population_mean: 总体均数，浮点数类型
                - raw_data: 原始数据列表，列表类型
                - sample_size: 样本量，整数类型
            - sample_statistics: 样本统计量
                - sample_mean: 样本均数，浮点数类型
                - sample_std: 样本标准差，浮点数类型
                - sample_min: 样本最小值，浮点数类型
                - sample_max: 样本最大值，浮点数类型
            - test_statistics: 检验统计量
                - degrees_of_freedom: 自由度，整数类型
                - t_value: t值，浮点数类型
                - p_value_two_sided: 双侧P值，浮点数类型
            - significance_tests: 显著性检验结果
                - p_greater_than_0_05: P值是否大于0.05，布尔类型
                - p_greater_than_0_01: P值是否大于0.01，布尔类型
                - significant_at_05: 0.05水平显著性，字符串类型 ("显著" 或 "不显著")
                - significant_at_01: 0.01水平显著性，字符串类型 ("显著" 或 "不显著")
            - interpretation: 统计解释
                - t_test_interpretation: t检验解释（0.05水平），字符串类型
                - t_test_interpretation_01: t检验解释（0.01水平），字符串类型
            - remark: 备注信息，字符串类型，包含总体均数和样本量摘要
    
    Example:
        >>> data = [138.5, 142.3, 139.8, 141.2, 137.9, 140.5, 139.1, 141.8, 138.7, 140.2]
        >>> result = cal_result_t_single2(pop_mean=140.0, raw_data=data)
        >>> print(result['table_name'])
        '单样本t检验分析（原始数据）'
        >>> print(result['test_statistics']['t_value'])
        0.123456  # 示例值
        >>> print(result['significance_tests']['significant_at_05'])
        '不显著'  # 示例值
    """
    # 从参数对象解构
    pop_mean = param.pop_mean
    raw_data = param.stats_data_list[0].data_list

    # 执行单样本t检验
    results = perform_single_sample_t_test_from_raw(pop_mean, raw_data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "单样本t检验分析（原始数据）",
        "input_parameters": results["input_parameters"],
        "sample_statistics": results["sample_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"总体均数：{pop_mean}, 原始数据点数：{len(raw_data)}"
    }
    
    return result_dict