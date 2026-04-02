"""
临床频率统计分析模块 (Clinical Frequency Statistics Module)

本模块提供全面的频率统计分析功能，用于临床数据的频数分布分析，是探索性数据分析（EDA）的重要组成部分。
频率统计分析能够直观展示数据的分布规律和特征，帮助研究者了解数据的集中趋势、离散程度和分布形态，
为后续的统计推断提供重要依据。

【模块功能概述】:
1. 频数分布表：将连续型数据分组并计算各组段的频数
2. 频率计算：计算各组段占总数的百分比
3. 累积频数与累积频率：显示逐步累积的频数和频率
4. 组中值计算：计算每个组段的中心值
5. 区间划分：自动将数据划分为指定数量的区间
6. 直方图数据生成：为可视化提供数据支持
7. 分位数区间识别：确定数据落在哪个分位数区间

【临床应用价值】:
- 数据分布探索：了解临床指标（如血压、血糖、年龄等）的分布形态
- 异常值识别：通过频数分布识别极端值或异常数据
- 统计方法选择：根据数据分布形态选择合适的统计方法
- 结果可视化基础：为直方图等可视化提供数据基础
- 分组分析：为临床决策提供分组界限参考

【统计方法选择指南】:
1. 数据分组策略：
   - 组数选择：一般为数据量的平方根或Sturges公式确定
   - 组距确定：(最大值-最小值)/组数，保证各组段连续且互斥
   
2. 适用场景：
   - 连续型变量的分布特征分析
   - 大样本数据的初步探索
   - 数据质量评估（如识别异常值、缺失模式）

【结果解读注意事项】:
1. 组数选择：组数过少会掩盖数据细节，过多则无法达到简化数据的目的
2. 分组合理性：确保各组段连续且互斥，避免数据遗漏或重复计算
3. 频数分布特征：关注分布的对称性、峰度、是否有异常峰值等
4. 临床意义：结合专业知识解读频数分布的临床含义
5. 与理论分布对比：可与正态分布等理论分布对比，判断数据特征

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据应该如何分组才能最好地展示分布特征？"
- "频数分布表中的累积频率有什么意义？"
- "如何根据频数分布判断数据是否符合正态分布？"
- "频数分布表在临床研究中的作用是什么？"
- "如何设置合适的分组数量？"

【相关标准和规范】:
- CONSORT 声明：随机临床试验报告规范
- STROBE 声明：观察性研究报告规范
- ICH E9 指导原则：临床试验的统计学原则
- SAMPL 指南：医学研究报告中的统计方法描述规范

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from typing import List, Dict, Any
from decimal import Decimal, ROUND_HALF_UP
import math
from app.schemas.request_data.base_param import BaseParamFreq
from app.schemas.request_data.stats_data import StatsData

def calculate_frequency_stats(data: List[float], num_intervals: int = 10, decimal_places: int = 2) -> List[Dict[str, Any]]:
    """
    计算频率统计信息，生成频数分布表
    
    频数分布表是统计分析的基础工具，通过将数据按一定间隔分组，统计各组内的观测值数量，
    从而展示数据的分布特征。在临床研究中，频数分布表常用于探索变量的分布形态、
    识别异常值以及为后续统计分析提供基础。

    Args:
        data: 输入的数值列表，通常为连续型临床指标（如血压、血糖、年龄等）
        num_intervals: 分组数量，默认为10组，可根据数据特点调整
        decimal_places: 结果保留的小数位数，默认为2位

    Returns:
        List[Dict[str, Any]]: 包含每组统计信息的字典列表，每个字典包含:
            - interval: 组段标识（如"10.00-15.00"）
            - lower_bound: 组下限
            - upper_bound: 组上限
            - midpoint: 组中值
            - frequency: 频数（该组内观测值数量）
            - cumulative_frequency: 累积频数
            - relative_frequency: 频率（百分比）
            - cumulative_relative_frequency: 累积频率（百分比）

    Raises:
        ValueError: 当输入数据为空或分组数无效时抛出异常

    Example:
        >>> data = [1.2, 2.5, 3.7, 4.1, 5.3, 6.8, 7.2, 8.5, 9.1, 10.0]
        >>> result = calculate_frequency_stats(data, num_intervals=5)
        >>> print(result[0]['interval'])  # 输出首个组段
        >>> print(result[0]['frequency'])  # 输出首个组段的频数
    """
    if not data:
        raise ValueError("数据列表不能为空")
    
    if num_intervals <= 0:
        raise ValueError("组数必须大于0")
    
    # 对数据进行排序以便计算
    sorted_data = sorted(data)
    min_val = min(sorted_data)
    max_val = max(sorted_data)
    
    # 计算组距
    interval_width = (max_val - min_val) / num_intervals
    
    # 初始化结果列表
    result = []
    
    # 计算每组的统计信息
    for i in range(num_intervals):
        # 计算当前组的范围
        lower_bound = min_val + i * interval_width
        upper_bound = min_val + (i + 1) * interval_width
        
        # 确保最后一组包含最大值
        if i == num_intervals - 1:
            upper_bound = max_val
            
        # 计算组中值
        midpoint = (lower_bound + upper_bound) / 2
        
        # 计算当前组的频数
        frequency = 0
        for value in data:
            # 判断值是否在当前组内
            if i == num_intervals - 1:  # 最后一组包含右边界
                if lower_bound <= value <= upper_bound:
                    frequency += 1
            else:  # 其他组不包含右边界
                if lower_bound <= value < upper_bound:
                    frequency += 1
        
        # 四舍五入到指定小数位数
        rounded_lower = round_to_decimal(lower_bound, decimal_places)
        rounded_upper = round_to_decimal(upper_bound, decimal_places)
        rounded_midpoint = round_to_decimal(midpoint, decimal_places)
        
        result.append({
            "interval": f"{rounded_lower}-{rounded_upper}",
            "lower_bound": rounded_lower,
            "upper_bound": rounded_upper,
            "midpoint": rounded_midpoint,
            "frequency": frequency
        })
    
    # 计算累计频数和频率
    total_count = len(data)
    cumulative_freq = 0
    
    for i in range(len(result)):
        # 累计频数
        cumulative_freq += result[i]["frequency"]
        result[i]["cumulative_frequency"] = cumulative_freq
        
        # 频率 (%)
        relative_freq = (result[i]["frequency"] / total_count) * 100
        result[i]["relative_frequency"] = round_to_decimal(relative_freq, decimal_places)
        
        # 累计频率 (%)
        cumulative_relative_freq = (result[i]["cumulative_frequency"] / total_count) * 100
        result[i]["cumulative_relative_frequency"] = round_to_decimal(cumulative_relative_freq, decimal_places)
    
    return result


def round_to_decimal(value: float, decimal_places: int) -> float:
    """
    将浮点数四舍五入到指定的小数位数
    
    在临床数据分析中，为了结果的可读性和报告的规范性，通常需要将计算结果
    四舍五入到指定位数。此函数使用精确的十进制运算，避免浮点数运算误差。

    Args:
        value: 待四舍五入的浮点数值
        decimal_places: 保留的小数位数

    Returns:
        float: 四舍五入后的浮点数值
    """
    quantize_str = '0.' + '0' * decimal_places if decimal_places > 0 else '1'
    rounded_value = Decimal(str(value)).quantize(Decimal(quantize_str), rounding=ROUND_HALF_UP)
    return float(rounded_value)


def get_histogram_data(data_list: List[float], num_intervals: int = 10) -> Dict[str, List]:
    """
    为直方图可视化生成所需的数据
    
    此函数生成适合用于绘制直方图的数据结构，包括区间边界、频数、频率等信息，
    便于前端可视化组件直接使用。在临床数据分析中，直方图是展示数据分布的
    重要图形工具。

    Args:
        data_list: 包含浮点数的列表，代表一组连续的观测数据
        num_intervals: 分组数量，默认为10组

    Returns:
        Dict[str, List]: 包含直方图绘制所需数据的字典
            - bins: 区间边界值列表
            - frequencies: 各区间的频数列表
            - relative_frequencies: 各区间的频率百分比列表
            - bin_labels: 区间标签列表（用于显示）
    """
    if not data_list:
        raise ValueError("数据列表不能为空")
    
    if num_intervals <= 0:
        raise ValueError("组数必须大于0")
    
    # 获取频率统计数据
    freq_stats = calculate_frequency_stats(data_list, num_intervals)
    
    # 提取数据用于直方图
    bins = []
    frequencies = []
    relative_frequencies = []
    bin_labels = []
    
    # 构建bins数组，需要n+1个值来定义n个区间
    if freq_stats:
        bins.append(freq_stats[0]["lower_bound"])
        for stat in freq_stats:
            bins.append(stat["upper_bound"])
            frequencies.append(stat["frequency"])
            relative_frequencies.append(stat["relative_frequency"])
            bin_labels.append(stat["interval"])
    
    return {
        "bins": bins,
        "frequencies": frequencies,
        "relative_frequencies": relative_frequencies,
        "bin_labels": bin_labels
    }


def find_interval_for_value(value: float, data_list: List[float], num_intervals: int = 10) -> Dict[str, Any]:
    """
    确定给定值属于哪个区间
    
    此函数用于确定一个特定的数值属于哪个频数分布区间，这在临床决策支持中很有用，
    例如确定某个患者的指标值属于正常、异常还是危急范围。

    Args:
        value: 需要查找区间的数值
        data_list: 用于定义区间的参考数据列表
        num_intervals: 分组数量，默认为10组

    Returns:
        Dict[str, Any]: 包含区间信息的字典
            - interval_index: 区间索引
            - interval_info: 区间详细信息（边界、中值、频数等）
            - is_within_range: 值是否在数据范围内
    """
    if not data_list:
        raise ValueError("数据列表不能为空")
    
    if num_intervals <= 0:
        raise ValueError("组数必须大于0")
    
    # 获取频率统计数据
    freq_stats = calculate_frequency_stats(data_list, num_intervals)
    
    # 查找值所在的区间
    interval_index = -1
    is_within_range = True
    
    # 检查值是否超出数据范围
    min_val = min(data_list)
    max_val = max(data_list)
    
    if value < min_val or value > max_val:
        is_within_range = False
    
    # 遍历区间查找值所属区间
    for idx, stat in enumerate(freq_stats):
        lower_bound = stat["lower_bound"]
        upper_bound = stat["upper_bound"]
        
        # 最后一个区间包含上边界
        if idx == len(freq_stats) - 1:
            if lower_bound <= value <= upper_bound:
                interval_index = idx
                break
        else:
            if lower_bound <= value < upper_bound:
                interval_index = idx
                break
    
    # 构建结果
    result = {
        "interval_index": interval_index,
        "is_within_range": is_within_range
    }
    
    if interval_index != -1:
        result["interval_info"] = freq_stats[interval_index]
    else:
        result["interval_info"] = None
    
    return result


def get_frequency_summary(data_list: List[float]) -> Dict[str, float]:
    """
    生成频率分布的关键摘要统计量
    
    此函数提供关于频率分布的一些关键统计量，如众数（最频繁的区间）、
    数据分布的偏度指示等。这些指标有助于快速了解数据分布特征。

    Args:
        data_list: 包含浮点数的列表，代表一组连续的观测数据

    Returns:
        Dict[str, float]: 包含频率分布摘要统计量的字典
            - mode_interval_frequency: 众数区间的频数
            - mode_interval_midpoint: 众数区间的中值
            - mode_interval_bounds: 众数区间的边界 (lower, upper)
            - max_frequency: 最大频数
            - min_frequency: 最小频数
            - frequency_range: 频数范围
    """
    if not data_list:
        raise ValueError("数据列表不能为空")
    
    # 使用默认10个区间计算频率统计
    freq_stats = calculate_frequency_stats(data_list)
    
    # 找到众数区间（频数最大的区间）
    max_freq_item = max(freq_stats, key=lambda x: x["frequency"])
    min_freq_item = min(freq_stats, key=lambda x: x["frequency"])
    
    return {
        "mode_interval_frequency": max_freq_item["frequency"],
        "mode_interval_midpoint": max_freq_item["midpoint"],
        "mode_interval_bounds": (max_freq_item["lower_bound"], max_freq_item["upper_bound"]),
        "max_frequency": max_freq_item["frequency"],
        "min_frequency": min_freq_item["frequency"],
        "frequency_range": max_freq_item["frequency"] - min_freq_item["frequency"]
    }


def cal_result_freq(data_list: List[float], num_intervals: int = 10) -> Dict[str, Any]:
    """
    生成频率统计分析的完整报告字典
    
    此函数整合了频率统计的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的数据分布特征概览。
    报告包括各组段的频数、频率、累积频数和累积频率等信息，
    便于临床医生和研究人员快速理解数据的分布特征。

    Args:
        data_list: 包含浮点数的列表，代表一组连续的观测数据
        num_intervals: 分组数量，默认为10组，可根据数据特点调整

    Returns:
        Dict[str, Any]: 包含频率统计分析指标的字典，键为指标名称，值为对应的统计量
            - sample_size: 样本量，表示观测值的总数
            - intervals: 分组信息列表，包含每个区间的详细统计
              - interval: 区间范围字符串
              - lower_bound: 区间下界
              - upper_bound: 区间上界
              - midpoint: 区间中点值
              - frequency: 频数
              - relative_frequency: 频率(%)
              - cumulative_frequency: 累积频数
              - cumulative_relative_frequency: 累积频率(%)
            - total_count: 总计数

    Raises:
        ValueError: 当输入列表为空或分组数无效时抛出异常

    Example:
        >>> data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        >>> result = cal_result_freq(data, num_intervals=5)
        >>> print(result['sample_size'])  # 输出：10
        >>> print(result['intervals'][0]['frequency'])  # 输出第一个区间的频数
    """
    if not data_list:
        raise ValueError("数据列表不能为空")
    
    if num_intervals <= 0:
        raise ValueError("组数必须大于0")
    
    # 计算频率统计详情
    freq_details = calculate_frequency_stats(data_list, num_intervals)
    
    # 构建并返回完整的频率统计报告字典
    report_dict = {
        'sample_size': len(data_list),
        'intervals': freq_details,
        'total_count': len(data_list)
    }
    
    return report_dict