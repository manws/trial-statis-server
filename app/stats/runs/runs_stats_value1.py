"""
游程检验统计分析模块 - 按平均数分类的数值变量 (Runs Test Statistics Module - Numerical Variables Classified by Mean)

本模块提供全面的游程检验统计分析功能，用于检验数值变量序列的随机性，
通过将数值数据按平均数上下分为二分类数据后进行游程检验。
游程检验在医学研究和质量控制中广泛应用，如检验实验数据的随机性、
验证数据收集过程是否受到系统性因素影响等。

【模块功能概述】:
1. 数据转换：将数值数据按平均数上下分为二分类数据
2. 游程计数：计算二分类数据序列中的游程数
3. 游程检验：执行游程检验，计算统计量和 P 值
4. 随机性评估：评估数据序列的随机性
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 随机性检验：验证数据序列是否随机
- 质量控制：监测数据收集过程的稳定性
- 实验设计验证：检验随机分配的有效性
- 时间序列分析：评估时间趋势的存在性

【统计方法选择指南】:
1. 游程检验适用条件：
   - 数据序列独立
   - 二分类变量
   - 样本量不宜过小（一般要求各类别个数不少于 10）
   - 适用于检验序列的随机性

2. 临床应用场景：
   - 检验实验数据的随机性
   - 验证随机化分组的效果
   - 监控实验室数据的质量
   - 评估测量过程的稳定性

【结果解读注意事项】:
1. 游程数解释：游程数过多或过少都可能表示非随机性
2. P 值解释：在零假设（数据随机）成立的情况下，观察到当前或更极端结果的概率
3. 显著性水平：通常使用α=0.05 作为判断随机性的标准
4. 临床意义：随机性是许多统计方法和实验设计的重要前提

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用游程检验吗？"
- "如何解释游程检验的结果？"
- "P 值小于 0.05 意味着什么？"
- "游程检验的前提条件是什么？"
- "游程检验与随机性有什么关系？"

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
from typing import List, Tuple, Dict, Any
from scipy import stats
from app.schemas.request_data.runs_param import RunsParamValue1


def convert_to_binary_by_mean(data: List[float]) -> List[int]:
    """
    将数值数据按平均数上下分为二分类数据
    
    在临床研究中，这用于将连续变量转换为二分类变量，
    便于进行游程检验等非参数检验。
    
    Args:
        data: 数值数据列表
        
    Returns:
        二分类数据列表：大于等于平均数为 1，小于平均数为 0
        
    Raises:
        ValueError: 当数据无效时抛出异常
    """
    if not data:
        raise ValueError("数据不能为空")
    
    if len(data) < 2:
        raise ValueError("数据长度必须至少为 2")
    
    mean_value = np.mean(data)
    binary_data = [1 if x >= mean_value else 0 for x in data]
    
    return binary_data


def count_runs(data: List[int]) -> Tuple[int, int, int]:
    """
    计算二分类数据序列中的游程数
    
    在游程检验中，这是计算统计量的基础步骤，
    游程是指序列中连续相同值的段落。
    
    Args:
        data: 二分类数据列表，只包含 0 和 1
        
    Returns:
        tuple: (游程数，0 的个数，1 的个数)
    """
    if not data:
        return 0, 0, 0
    
    runs = 1  # 至少有一个游程
    count_0 = 0
    count_1 = 0
    
    # 统计 0 和 1 的个数
    for value in data:
        if value == 0:
            count_0 += 1
        elif value == 1:
            count_1 += 1
    
    # 计算游程数
    for i in range(1, len(data)):
        if data[i] != data[i-1]:
            runs += 1
    
    return runs, count_0, count_1


def calculate_runs_test_value_by_mean(data: List[float]) -> Dict:
    """
    对数值变量按平均数上下分类后执行游程检验
    
    在临床研究中，这用于检验数值序列的随机性，
    通过将数据按平均数分为两类后计算游程数并进行检验。
    
    Args:
        data: 数值数据列表
        
    Returns:
        包含检验统计量和结果的字典
        
    Raises:
        ValueError: 当数据无效时抛出异常
    """
    if not data:
        raise ValueError("数据不能为空")
    
    original_data = data
    
    if len(original_data) < 2:
        raise ValueError("数据长度必须至少为 2")
    
    # 转换为二分类数据
    binary_data = convert_to_binary_by_mean(original_data)
    
    # 计算基本统计量
    mean_value = np.mean(original_data)
    std_value = np.std(original_data, ddof=1)
    
    # 计算游程数和各类别计数
    runs, count_0, count_1 = count_runs(binary_data)
    
    # 计算期望游程数和方差
    n = len(binary_data)
    if count_0 == 0 or count_1 == 0:
        expected_runs = 1
        var_runs = 0
    else:
        expected_runs = 1 + (2 * count_0 * count_1) / n
        var_runs = (2 * count_0 * count_1 * (2 * count_0 * count_1 - n)) / (n * n * (n - 1))
    
    # 计算标准化统计量
    if var_runs > 0:
        z_statistic = (runs - expected_runs) / np.sqrt(var_runs)
    else:
        z_statistic = 0
    
    # 计算双侧 P 值
    p_value_two_sided = 2 * (1 - stats.norm.cdf(abs(z_statistic)))
    
    # 确定 P 值范围
    if p_value_two_sided > 0.05:
        p_value_range = "P > 0.05"
    elif 0.01 < p_value_two_sided <= 0.05:
        p_value_range = "0.01 < P ≤ 0.05"
    else:
        p_value_range = "P ≤ 0.01"
    
    return {
        "input_parameters": {
            "sample_size": n,
            "original_mean": float(mean_value),
            "original_std": float(std_value),
            "count_above_mean": count_1,
            "count_below_mean": count_0
        },
        "conversion_info": {
            "classification_method": "按平均数上下分类",
            "threshold_value": float(mean_value),
            "above_threshold_count": count_1,
            "below_threshold_count": count_0
        },
        "run_statistics": {
            "observed_runs": runs,
            "expected_runs": expected_runs,
            "var_runs": var_runs,
            "std_runs": np.sqrt(var_runs) if var_runs > 0 else 0
        },
        "test_statistics": {
            "z_statistic": z_statistic,
            "p_value_two_sided": p_value_two_sided,
            "p_value_range": p_value_range
        },
        "significance_tests": {
            "significant_at_05": p_value_two_sided <= 0.05,
            "significant_at_01": p_value_two_sided <= 0.01
        },
        "interpretation": {
            "data_conversion": f"原始数据按平均数{mean_value:.4f}分为两类：≥{mean_value:.4f}为 1，<{mean_value:.4f}为 0",
            "runs_interpretation": f"观察到的游程数为{runs}，期望游程数为{expected_runs:.2f}",
            "significance_interpretation": f"在 0.05 显著性水平下{'显著' if p_value_two_sided <= 0.05 else '不显著'} {p_value_range}",
            "randomness_assessment": "数据序列呈现随机性" if p_value_two_sided > 0.05 else "数据序列可能存在非随机模式"
        }
    }


def cal_result_runs_value1(param: RunsParamValue1) -> Dict[str, Any]:
    """
    生成按平均数分类的数值变量游程检验分析的完整报告字典
    
    此函数整合了游程检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的游程检验结果。报告包括输入参数、
    数据转换信息、游程统计、检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解游程检验的特征。
    
    Args:
        data: 数值数据列表
    
    Returns:
        Dict[str, Any]: 包含游程检验分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 表格名称，固定为"按平均数分类的数值变量游程检验分析"
            - input_parameters: 输入参数信息，包含以下子字段：
                - sample_size: 样本总数
                - original_mean: 原始数据的平均值
                - original_std: 原始数据的标准差
                - count_above_mean: 大于等于平均值的数据个数
                - count_below_mean: 小于平均值的数据个数
            - conversion_info: 数据转换信息，包含以下子字段：
                - classification_method: 分类方法，固定为"按平均数上下分类"
                - threshold_value: 分类阈值，即平均值
                - above_threshold_count: 高于阈值的数据个数
                - below_threshold_count: 低于阈值的数据个数
            - run_statistics: 游程统计信息，包含以下子字段：
                - observed_runs: 观察到的游程数
                - expected_runs: 期望的游程数
                - var_runs: 游程数的方差
                - std_runs: 游程数的标准差
            - test_statistics: 检验统计量，包含以下子字段：
                - z_statistic: Z统计量
                - p_value_two_sided: 双侧检验的P值
                - p_value_range: P值所属的范围
            - significance_tests: 显著性检验结果，包含以下子字段：
                - significant_at_05: 是否在0.05显著性水平下显著
                - significant_at_01: 是否在0.01显著性水平下显著
            - interpretation: 统计解释，包含以下子字段：
                - data_conversion: 数据转换说明
                - runs_interpretation: 游程解释
                - significance_interpretation: 显著性解释
                - randomness_assessment: 随机性评估
            - remark: 备注信息，包含样本量
    """
    # 执行按平均数分类的游程检验
    # 从参数对象解构
    data = param.stats_data_list[0].data_list

    results = calculate_runs_test_value_by_mean(data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "按平均数分类的数值变量游程检验分析",
        "input_parameters": results["input_parameters"],
        "conversion_info": results["conversion_info"],
        "run_statistics": results["run_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"样本量：{results['input_parameters']['sample_size']}"
    }
    
    return result_dict