"""
游程检验统计分析模块 - 按中位数分类的数值变量 (Runs Test Statistics Module - Numerical Variables Classified by Median)

本模块提供全面的游程检验统计分析功能，用于检验数值变量序列的随机性，
通过将数值数据按中位数上下分为二分类数据后进行游程检验。
游程检验在医学研究和质量控制中广泛应用，如检验实验数据的随机性、
验证数据收集过程是否受到系统性因素影响等。

【模块功能概述】:
1. 数据转换：将数值数据按中位数上下分为二分类数据
2. 游程计数：计算二分类数据序列中的游程数
3. 游程检验：执行游程检验，计算统计量和P值
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
   - 样本量不宜过小（一般要求各类别个数不少于10）
   - 适用于检验序列的随机性

2. 临床应用场景：
   - 检验实验数据的随机性
   - 验证随机化分组的效果
   - 监控实验室数据的质量
   - 评估测量过程的稳定性

【结果解读注意事项】:
1. 游程数解释：游程数过多或过少都可能表示非随机性
2. P值解释：在零假设（数据随机）成立的情况下，观察到当前或更极端结果的概率
3. 显著性水平：通常使用α=0.05作为判断随机性的标准
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
- "P值小于0.05意味着什么？"
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

import math
from typing import Any, Dict, List, Tuple
from scipy.stats import norm
import numpy as np


def convert_to_binary_by_median(data: List[float]) -> List[int]:
    """
    将数值数据按中位数分为二分类数据
    
    在临床研究中，这用于将连续变量转换为二分类变量，
    便于进行游程检验等非参数检验。
    
    统计学原理：
    - 中位数是将数据集分为两半的数值，50%的数据小于中位数，50%的数据大于等于中位数
    - 通过中位数转换可以消除原始数据的分布特征，专注于序列的随机性
    
    临床应用场景：
    - 将连续的血药浓度数据转换为高低两组，检验给药时间的随机性
    - 将患者的血压测量值转换为正常/异常两类，评估测量顺序是否存在系统性偏差
    - 在质量控制中，将测量值转换为合格/不合格，检验生产过程的稳定性
    
    Args:
        data: 数值数据列表，代表连续的观测值序列
        
    Returns:
        List[int]: 二分类数据列表，其中：
            - 1 表示大于等于中位数的值
            - 0 表示小于中位数的值
            
    Raises:
        ValueError: 当数据为空或观测值少于2个时抛出异常
        
    Example:
        >>> data = [10.5, 12.3, 11.8, 13.2, 9.7]
        >>> convert_to_binary_by_median(data)
        [0, 1, 1, 1, 0]  # 假设中位数为11.8
    """
    if not data:
        raise ValueError("数据不能为空")
    
    if len(data) < 2:
        raise ValueError("数据至少需要2个观测值才能计算中位数")
    
    # 计算中位数
    median_value = np.median(data)
    
    # 转换为二分类数据
    binary_data = [1 if x >= median_value else 0 for x in data]
    
    return binary_data


def count_runs(data: List[int]) -> Tuple[int, int, int]:
    """
    计算二分类数据的游程个数
    
    在游程检验中，这是计算统计量的基础步骤，
    游程是指序列中连续相同值的段落。
    
    统计学原理：
    - 游程（Run）定义为序列中连续相同符号的最大子序列
    - 例如：序列 [1,1,0,0,0,1,0] 包含4个游程：[1,1], [0,0,0], [1], [0]
    - 游程数是检验序列随机性的关键统计量
    
    临床应用场景：
    - 检验临床试验中患者分组的随机性（如治疗组/对照组）
    - 评估实验室检测结果的交替模式是否随机
    - 监测药物不良反应发生的聚集性
    
    方法学局限性：
    - 对样本量敏感，小样本时检验效能较低
    - 仅能检测一阶自相关，无法识别高阶依赖关系
    - 要求数据独立性，不适用于存在明显时间趋势的序列
    
    Args:
        data: 二分类数据列表，只包含0和1
        
    Returns:
        Tuple[int, int, int]: 包含三个整数的元组：
            - 第一个值：游程个数
            - 第二个值：类别0的个数
            - 第三个值：类别1的个数
            
    Raises:
        ValueError: 当数据为空或包含非二分类值时抛出异常
        
    Example:
        >>> data = [1, 1, 0, 0, 0, 1, 0]
        >>> count_runs(data)
        (4, 4, 3)  # 4个游程，4个0，3个1
    """
    if not data:
        raise ValueError("数据不能为空")
    
    # 统计两类数据的数量
    count_0 = data.count(0)
    count_1 = data.count(1)
    
    if count_0 + count_1 != len(data):
        raise ValueError("数据必须是二分类的(只包含0和1)")
    
    if count_0 == 0 or count_1 == 0:
        return (1, count_0, count_1)  # 如果只有一类，游程数为1
    
    # 计算游程数
    runs = 1
    for i in range(1, len(data)):
        if data[i] != data[i-1]:
            runs += 1
    
    return (runs, count_0, count_1)


def calculate_runs_test_value_by_median(data: List[float]) -> Dict[str, Any]:
    """
    按中位数上下分类的数值变量游程检验
    
    在临床研究中，这用于检验数值序列的随机性，
    通过将数据按中位数分为两类后计算游程数并进行检验。
    
    统计学原理：
    - 零假设（H0）：数据序列是随机的
    - 备择假设（H1）：数据序列存在非随机模式
    - 检验统计量：基于观察到的游程数与期望游程数的差异
    - 大样本时使用正态近似计算Z统计量和P值
    
    数学公式：
    - 期望游程数：E(R) = 1 + (2*n1*n2) / (n1+n2)
    - 游程方差：Var(R) = [2*n1*n2*(2*n1*n2-n1-n2)] / [(n1+n2)^2*(n1+n2-1)]
    - Z统计量：Z = (R - E(R)) / sqrt(Var(R))
    
    临床应用场景：
    - 验证随机对照试验中的随机分配是否成功
    - 检验时间序列数据是否存在周期性或趋势
    - 评估测量仪器的稳定性（如重复测量的随机性）
    - 检测数据收集过程中是否存在系统性偏差
    
    结果解读指导：
    - P > 0.05：无充分证据拒绝随机性假设，数据序列呈现随机性
    - P ≤ 0.05：拒绝随机性假设，数据序列可能存在聚集性或周期性模式
    - 游程数过少：提示数据存在聚集性或正向自相关
    - 游程数过多：提示数据存在交替模式或负向自相关
    
    使用注意事项：
    - 要求数据序列独立，不适用于存在自相关的时序数据
    - 样本量应足够大（建议每类不少于10个观测值）
    - 对异常值不敏感，因为使用中位数进行分类
    - 仅能检测一阶依赖关系，可能需要结合其他检验方法
    
    Args:
        data: 数值数据列表，代表待检验的观测值序列
        
    Returns:
        Dict[str, Any]: 包含游程检验统计量和结果的字典，包括：
            - input_parameters: 输入参数信息（样本量、中位数、均值、标准差等）
            - conversion_info: 数据转换信息（分类方法、阈值、各类数量）
            - run_statistics: 游程统计信息（观察游程数、期望游程数、方差、标准差）
            - test_statistics: 检验统计量（Z统计量、双侧P值、P值范围）
            - significance_tests: 显著性检验结果（0.05和0.01水平的显著性）
            - interpretation: 统计解释（数据转换说明、游程解释、显著性解释、随机性评估）
            
    Raises:
        ValueError: 当数据为空或观测值少于2个时抛出异常
        
    Example:
        >>> data = [10.5, 12.3, 11.8, 13.2, 9.7, 14.1, 10.9]
        >>> result = calculate_runs_test_value_by_median(data)
        >>> result['test_statistics']['p_value_two_sided']
        0.xxxx  # 具体的P值
    """
    # 参数验证
    if not data:
        raise ValueError("数据不能为空")
    
    if len(data) < 2:
        raise ValueError("数据至少需要2个观测值")
    
    # 提取原始数据
    original_data = data
    
    # 按中位数转换为二分类数据
    binary_data = convert_to_binary_by_median(original_data)
    
    # 计算基本统计量
    median_value = np.median(original_data)
    mean_value = np.mean(original_data)
    std_value = np.std(original_data, ddof=1)  # 样本标准差
    count_above_median = binary_data.count(1)
    count_below_median = binary_data.count(0)
    
    # 计算游程统计量
    observed_runs, count_0, count_1 = count_runs(binary_data)
    total_n = count_0 + count_1
    
    # 计算期望游程数和方差
    expected_runs = 1 + (2 * count_0 * count_1) / total_n
    variance_runs = (2 * count_0 * count_1 * (2 * count_0 * count_1 - total_n)) / (total_n**2 * (total_n - 1))
    
    # 计算标准差
    std_runs = math.sqrt(max(variance_runs, 0))  # 避免负数开方
    
    # 计算近似Z统计量
    if std_runs > 0:
        z_statistic = (observed_runs - expected_runs) / std_runs
        # 计算双侧P值
        p_value_two_sided = 2 * (1 - norm.cdf(abs(z_statistic)))
    else:
        z_statistic = 0
        p_value_two_sided = 1.0
    
    # P值范围判断
    p_range = ""
    if p_value_two_sided > 0.05:
        p_range = "P > 0.05"
    elif p_value_two_sided > 0.01:
        p_range = "0.01 < P ≤ 0.05"
    else:
        p_range = "P ≤ 0.01"
    
    # 显著性判断
    is_significant_05 = p_value_two_sided <= 0.05
    is_significant_01 = p_value_two_sided <= 0.01
    
    return {
        "input_parameters": {
            "sample_size": int(total_n),
            "original_median": float(median_value),
            "original_mean": float(mean_value),
            "original_std": float(std_value),
            "count_above_median": int(count_above_median),
            "count_below_median": int(count_below_median)
        },
        "conversion_info": {
            "classification_method": "按中位数上下分类",
            "threshold_value": float(median_value),
            "above_threshold_count": int(count_above_median),
            "below_threshold_count": int(count_below_median)
        },
        "run_statistics": {
            "observed_runs": int(observed_runs),
            "expected_runs": float(expected_runs),
            "variance_runs": float(variance_runs),
            "std_runs": float(std_runs)
        },
        "test_statistics": {
            "z_statistic": float(z_statistic),
            "p_value_two_sided": float(p_value_two_sided),
            "p_value_range": p_range
        },
        "significance_tests": {
            "significant_at_05": is_significant_05,
            "significant_at_01": is_significant_01,
            "p_value_le_0_05": not is_significant_05,
            "p_value_le_0_01": not is_significant_01
        },
        "interpretation": {
            "data_conversion": f"原始数据按中位数{median_value:.4f}分为两类：≥中位数为1，<中位数为0",
            "runs_interpretation": f"观察到的游程数为{observed_runs}，期望游程数为{expected_runs:.2f}",
            "significance_interpretation": f"在0.05显著性水平下{'不' if not is_significant_05 else ''}显著 {p_range}",
            "randomness_assessment": "数据序列呈现随机性" if not is_significant_05 else "数据序列可能存在非随机模式"
        }
    }


def cal_result_runs_value2(data: List[float]) -> Dict[str, Any]:
    """
    生成按中位数分类的数值变量游程检验分析的完整报告字典
    
    此函数整合了游程检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的游程检验结果。报告包括输入参数、
    数据转换信息、游程统计、检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解游程检验的特征。
    
    函数特点：
    - 输出格式标准化：采用统一的字典结构，便于后续处理和展示
    - 信息完整性：包含从原始数据到最终解释的全流程信息
    - 临床友好性：提供易于理解的统计解释和临床意义说明
    - 可扩展性：便于添加新的统计指标或解释维度
    
    与其他统计方法的关系：
    - 与卡方检验比较：游程检验关注序列顺序，卡方检验关注频数分布
    - 与自相关分析比较：游程检验检测一阶依赖，自相关可检测多阶依赖
    - 与趋势检验比较：游程检验检测随机性，趋势检验检测单调变化
    
    在AI问答系统中的应用：
    - 可作为知识库的结构化数据来源
    - 支持自动生成统计报告
    - 便于进行多组数据的比较分析
    - 可用于教学演示和案例学习
    
    Args:
        data: 数值数据列表，代表待检验的观测值序列
    
    Returns:
        Dict[str, Any]: 包含游程检验分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"按中位数分类的数值变量游程检验分析"
            - input_parameters: 输入参数信息
                - sample_size: 样本量，整数类型
                - original_median: 原始数据中位数，浮点数类型
                - original_mean: 原始数据均值，浮点数类型
                - original_std: 原始数据标准差，浮点数类型
                - count_above_median: 大于等于中位数的数量，整数类型
                - count_below_median: 小于中位数的数量，整数类型
            - conversion_info: 数据转换信息
                - classification_method: 分类方法描述，固定为"按中位数上下分类"
                - threshold_value: 分类阈值（中位数），浮点数类型
                - above_threshold_count: 高于阈值的数量，整数类型
                - below_threshold_count: 低于阈值的数量，整数类型
            - run_statistics: 游程统计信息
                - observed_runs: 观察到的游程数，整数类型
                - expected_runs: 期望游程数，浮点数类型
                - variance_runs: 游程方差，浮点数类型
                - std_runs: 游程标准差，浮点数类型
            - test_statistics: 检验统计量
                - z_statistic: Z统计量，浮点数类型
                - p_value_two_sided: 双侧P值，浮点数类型
                - p_value_range: P值范围描述，字符串类型
            - significance_tests: 显著性检验结果
                - significant_at_05: 0.05水平是否显著，布尔类型
                - significant_at_01: 0.01水平是否显著，布尔类型
                - p_value_le_0_05: P值是否大于0.05，布尔类型
                - p_value_le_0_01: P值是否大于0.01，布尔类型
            - interpretation: 统计解释
                - data_conversion: 数据转换说明，字符串类型
                - runs_interpretation: 游程解释，字符串类型
                - significance_interpretation: 显著性解释，字符串类型
                - randomness_assessment: 随机性评估，字符串类型
            - remark: 备注信息（样本量），字符串类型
            
    Example:
        >>> data = [10.5, 12.3, 11.8, 13.2, 9.7, 14.1, 10.9]
        >>> result = cal_result_runs_value2(data)
        >>> print(result['table_name'])
        '按中位数分类的数值变量游程检验分析'
        >>> print(result['test_statistics']['p_value_range'])
        'P > 0.05'  # 或其他P值范围
    """
    # 执行按中位数分类的游程检验
    results = calculate_runs_test_value_by_median(data)
    
    # 构建结果字典
    result_dict = {
        "table_name": "按中位数分类的数值变量游程检验分析",
        "input_parameters": results["input_parameters"],
        "conversion_info": results["conversion_info"],
        "run_statistics": results["run_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"样本量: {results['input_parameters']['sample_size']}"
    }
    
    return result_dict