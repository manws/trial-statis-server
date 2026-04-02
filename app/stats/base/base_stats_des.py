"""
临床描述性统计模块 (Clinical Descriptive Statistics Module)

本模块提供全面的描述性统计分析功能，是临床研究数据探索性分析（Exploratory Data Analysis, EDA）的核心工具集。
描述性统计是推断性统计的基础，用于系统性地组织、总结和呈现数据的基本特征，帮助研究者理解数据的分布形态、
集中趋势和离散程度，为后续的假设检验和统计建模奠定基础。

【模块功能概述】:
1. 集中趋势测量：算术平均数、几何平均数、调和平均数、截尾均值、中位数、众数
2. 离散程度测量：方差、标准差、极差、标准误、变异系数、四分位间距
3. 位置参数：最大值、最小值、第一四分位数 (Q1)、第三四分位数 (Q3)
4. 推断统计：总体均值的 95% 置信区间估计
5. 综合报告：一键生成包含所有描述性统计指标的完整报告

【临床应用价值】:
- 基线特征描述：在随机对照试验 (RCT) 中描述研究人群的人口学特征和临床指标
- 疗效评估：总结治疗组与对照组的主要结局指标的分布特征
- 安全性分析：描述不良事件发生率、实验室检查异常值等安全性指标
- 数据质量评估：识别异常值、缺失值模式和数据分布特征
- 样本量估算：基于预试验的描述性统计结果进行正式研究的样本量计算

【统计方法选择指南】:
1. 正态分布数据：
   - 集中趋势：优先使用算术平均数 ± 标准差
   - 适用场景：血压、血糖、胆固醇等生理指标
   
2. 偏态分布数据：
   - 集中趋势：使用中位数 (四分位间距) 或几何平均数
   - 适用场景：生存时间、费用数据、抗体滴度、炎症标志物 (如 CRP、IL-6)
   
3. 存在极端值的数据：
   - 集中趋势：使用截尾均值或中位数，减少异常值影响
   - 离散程度：使用四分位间距而非标准差
   
4. 比较不同量纲数据的变异程度：
   - 使用变异系数 (CV)，消除测量单位影响
   - 适用场景：比较身高和体重的变异程度

【结果解读注意事项】:
1. 正态性检验：在使用参数统计方法前，应通过 Shapiro-Wilk 检验或 Kolmogorov-Smirnov 检验验证数据正态性
2. 异常值处理：识别并审慎处理异常值，区分真实极端值和录入错误
3. 缺失数据：明确缺失值处理方式（删除、插补等），避免引入偏倚
4. 样本量影响：小样本 (n<30) 时，置信区间较宽，估计精度较低
5. 效应量报告：除 P 值外，应报告效应量（如 Cohen's d）以评估临床意义

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据是偏态分布，应该用哪个统计量描述集中趋势？"
- "如何解释标准差和标准误的区别？"
- "什么情况下需要使用几何平均数？"
- "95% 置信区间的临床意义是什么？"
- "如何判断我的数据是否存在异常值？"

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
import math
from collections import Counter
from scipy import stats

from app.schemas.request_data.base_param import BaseParamDes
from app.schemas.request_data.stats_data import StatsData

# ============================================================================
# 集中趋势测量函数 (Measures of Central Tendency)
# ============================================================================

def cal_mean(data_list: List[float]) -> float:
    """
    计算算术平均数 (Arithmetic Mean)
    
    在临床研究中，算术平均数是描述数据集中趋势的重要指标，表示所有观测值的平均水平。
    它是最常用的集中趋势指标，但容易受到极端值的影响，因此在分析前应检查数据分布形态。
    
    Args:
        data_list: 包含数值的列表，通常为连续型变量（如血压、血糖水平等）
        
    Returns:
        float: 算术平均数，表示数据的中心位置
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute mean of an empty list")
    return sum(data_list) / len(data_list)


def cal_geometric_mean(data_list: List[float]) -> float:
    """
    计算几何平均数 (Geometric Mean)
    
    几何平均数适用于呈倍数关系的数据或对数正态分布的数据，在临床药代动力学研究中尤为重要，
    如计算血药浓度的几何平均值。对于比例数据或增长率数据，几何平均数比算术平均数更准确。
    
    Args:
        data_list: 包含正数值的列表，如抗体滴度、细菌计数等
        
    Returns:
        float: 几何平均数，适用于对数正态分布数据
        
    Raises:
        ValueError: 当输入列表为空或包含非正值时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute geometric mean of an empty list")
    if any(x <= 0 for x in data_list):
        raise ValueError("Geometric mean is only defined for positive numbers")
    
    log_sum = sum(math.log(x) for x in data_list)
    return math.exp(log_sum / len(data_list))


def cal_harmonic_mean(data_list: List[float]) -> float:
    """
    计算调和平均数 (Harmonic Mean)
    
    调和平均数适用于计算比率或速度的平均值，当数据中存在较小值时特别有用，
    在临床研究中可用于计算平均速率或某些特定指标的平均值。
    
    Args:
        data_list: 包含非零数值的列表
        
    Returns:
        float: 调和平均数，适用于比率数据的平均计算
        
    Raises:
        ValueError: 当输入列表为空或包含零值时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute harmonic mean of an empty list")
    if any(x == 0 for x in data_list):
        raise ValueError("Harmonic mean is not defined when any value is zero")
    
    reciprocal_sum = sum(1/x for x in data_list)
    return len(data_list) / reciprocal_sum


def cal_trimmed_mean(data_list: List[float], trim_ratio: float = 0.1) -> float:
    """
    计算截尾均值 (Trimmed Mean)
    
    截尾均值通过去除一定比例的极值来减少极端值对平均值的影响，
    在临床试验中，当数据可能存在异常值或严重偏态时，截尾均值比普通均值更能反映数据的集中趋势。
    默认截去上下各 5% 的极值（共 10%）。
    
    Args:
        data_list: 包含数值的列表
        trim_ratio: 截尾比例，范围应在 0 到 0.5 之间，默认为 0.1（即去除上下各 5% 的值）
        
    Returns:
        float: 截尾均值，减少极端值影响的集中趋势度量
        
    Raises:
        ValueError: 当输入列表为空、截尾比例不合法或截尾后列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute trimmed mean of an empty list")
    if not 0 <= trim_ratio < 0.5:
        raise ValueError("Trim ratio must be between 0 and 0.5")
    
    sorted_list = sorted(data_list)
    n = len(sorted_list)
    
    # 计算需要截去的元素数量
    trim_count = int(n * trim_ratio)
    
    # 截去首尾的极端值
    trimmed_list = sorted_list[trim_count:n-trim_count]
    
    if not trimmed_list:
        raise ValueError("Trimmed list is empty, reduce the trim ratio")
        
    return sum(trimmed_list) / len(trimmed_list)


def cal_median(data_list: List[float]) -> float:
    """
    计算中位数 (Median)
    
    中位数是描述数据集中趋势的稳健指标，不受极端值影响，适用于偏态分布的数据，
    在临床研究中经常用于描述生存时间、费用等右偏数据的中心位置。
    当数据呈对称分布时，中位数与均值接近；当数据偏态时，中位数更有代表性。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 中位数，表示数据排序后的中间值
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute median of an empty list")
    
    sorted_list = sorted(data_list)
    n = len(sorted_list)
    mid = n // 2

    if n % 2 == 0:
        return (sorted_list[mid - 1] + sorted_list[mid]) / 2.0
    else:
        return float(sorted_list[mid])


def cal_mode(data_list: List[float]) -> List[float]:
    """
    计算众数 (Mode)
    
    众数是数据中出现频率最高的值，可用来描述分类数据的集中趋势，
    在临床研究中可用于识别最常见的诊断结果、最频繁出现的症状等。
    数据可能有一个众数、多个众数或没有明确的众数。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        List[float]: 众数列表（可能存在多个众数）
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute mode of an empty list")
    
    counter = Counter(data_list)
    max_freq = max(counter.values())
    
    # 找到所有出现频率最高的值
    modes = [val for val, freq in counter.items() if freq == max_freq]
    return modes

# ============================================================================
# 离散程度测量函数 (Measures of Dispersion)
# ============================================================================

def cal_range(data_list: List[float]) -> float:
    """
    计算极差 (Range)
    
    极差是最大值与最小值之差，是最简单的离散程度测度，但在临床数据分析中需谨慎解释，
    因为它只考虑了两个极端值，不能反映中间数据的分布情况，且易受抽样波动影响。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 极差，表示数据的最大变化幅度
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute range of an empty list")
    return cal_max(data_list) - cal_min(data_list)


def cal_variance(data_list: List[float]) -> float:
    """
    计算方差 (Variance)
    
    方差是衡量数据离散程度的重要指标，表示观测值与均值差值平方的平均数，
    在临床研究中常用于评估疗效的一致性、生物标志物的变异性等。
    方差越大，说明数据越分散；方差越小，说明数据越集中在均值附近。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 样本方差（除以 n-1），表示数据的离散程度
        
    Raises:
        ValueError: 当数据少于 2 个值时抛出异常
    """
    if len(data_list) < 2:
        raise ValueError("Variance requires at least two values")
    
    mean_val = cal_mean(data_list)
    squared_diffs = [(x - mean_val) ** 2 for x in data_list]
    return sum(squared_diffs) / (len(data_list) - 1)  # 样本方差


def cal_std(data_list: List[float]) -> float:
    """
    计算标准差 (Standard Deviation)
    
    标准差是方差的平方根，与原始数据具有相同的单位，是临床研究报告中最常用的离散程度指标，
    常与均值一起报告（如 "均值±标准差"）。标准差反映了数据相对于均值的平均偏离程度，
    可用于估计数据的分布范围（如约 68% 的数据落在均值±1 个标准差范围内）。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 标准差，与原数据单位相同，表示数据的离散程度
        
    Raises:
        ValueError: 当数据少于 2 个值时抛出异常
    """
    if len(data_list) < 2:
        raise ValueError("Standard deviation requires at least two values")
    return math.sqrt(cal_variance(data_list))


def cal_se(data_list: List[float]) -> float:
    """
    计算标准误 (Standard Error)
    
    标准误衡量样本均值的精确度，反映抽样误差的大小，是构建置信区间和进行假设检验的基础，
    在临床试验中用于评估治疗效果估计的精度。标准误越小，样本均值对总体均值的估计越精确。
    标准误 = 标准差 / √n，可通过增加样本量来减小标准误。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 标准误，表示样本均值的变异性
        
    Raises:
        ValueError: 当数据少于 2 个值时抛出异常
    """
    if len(data_list) < 2:
        raise ValueError("Standard error requires at least two values")
    
    std = cal_std(data_list)
    n = len(data_list)
    return std / math.sqrt(n)


def cal_cv(data_list: List[float]) -> float:
    """
    计算变异系数 (Coefficient of Variation)
    
    变异系数是标准差与均值的比值，是无量纲的离散程度指标，便于比较不同单位或不同量级数据的相对变异程度，
    在临床实验室质量控制中广泛应用。变异系数越大，说明相对变异程度越大，数据的稳定性越差。
    注意：当均值接近 0 时，变异系数可能变得不稳定或无意义。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 变异系数，表示相对离散程度，无单位
        
    Raises:
        ValueError: 当输入列表为空或均值为 0 时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute CV of an empty list")
    
    mean_val = cal_mean(data_list)
    if mean_val == 0:
        raise ValueError("Coefficient of variation is undefined when mean is zero")
    
    return cal_std(data_list) / abs(mean_val)

# ============================================================================
# 位置参数函数 (Position Parameters)
# ============================================================================

def cal_q1(data_list: List[float]) -> float:
    """
    计算下四分位数 (First Quartile, Q1)
    
    下四分位数是位于数据序列 25% 位置的值，将数据分为较小的 25% 部分和其他 75% 部分，
    在临床研究中与中位数和上四分位数结合使用，可以描述数据分布的形状和对称性，
    并用于识别异常值（如用箱线图方法）。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 第一四分位数（25% 分位数）
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute Q1 of an empty list")
    
    sorted_list = sorted(data_list)
    n = len(sorted_list)
    
    # 计算下四分位数的位置 (n+1)/4
    pos = (n + 1) * 0.25
    lower_idx = int(pos) - 1
    upper_idx = lower_idx + 1
    
    # 如果位置正好是一个整数，则直接取该位置的值
    if pos == int(pos):
        return float(sorted_list[int(pos) - 1])
    # 否则进行插值
    elif upper_idx < n:
        fraction = pos - int(pos)
        return sorted_list[lower_idx] + fraction * (sorted_list[upper_idx] - sorted_list[lower_idx])
    else:
        return float(sorted_list[lower_idx])


def cal_q3(data_list: List[float]) -> float:
    """
    计算上四分位数 (Third Quartile, Q3)
    
    上四分位数是位于数据序列 75% 位置的值，将数据分为较小的 75% 部分和其他 25% 部分，
    与下四分位数配合使用，可以确定数据的四分位间距，用于描述中间 50% 数据的散布范围，
    在临床研究中常用于绘制箱线图和识别异常值。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 第三四分位数（75% 分位数）
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute Q3 of an empty list")
    
    sorted_list = sorted(data_list)
    n = len(sorted_list)
    
    # 计算上四分位数的位置 3*(n+1)/4
    pos = (n + 1) * 0.75
    lower_idx = int(pos) - 1
    upper_idx = lower_idx + 1
    
    # 如果位置正好是一个整数，则直接取该位置的值
    if pos == int(pos):
        return float(sorted_list[int(pos) - 1])
    # 否则进行插值
    elif upper_idx < n:
        fraction = pos - int(pos)
        return sorted_list[lower_idx] + fraction * (sorted_list[upper_idx] - sorted_list[lower_idx])
    else:
        return float(sorted_list[lower_idx])


def cal_iqr(data_list: List[float]) -> float:
    """
    计算四分位间距 (Interquartile Range, IQR)
    
    四分位间距是上四分位数与下四分位数之差（Q3-Q1），表示中间 50% 数据的分布范围，
    是一个稳健的离散程度指标，不受极端值影响。在临床研究中，IQR 常用于识别异常值：
    小于 Q1-1.5×IQR 或大于 Q3+1.5×IQR 的值被视为潜在异常值，
    大于 Q1-3×IQR 或小于 Q3+3×IQR 的值被视为极端异常值。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 四分位间距，表示中间 50% 数据的散布范围
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute IQR of an empty list")
    
    return cal_q3(data_list) - cal_q1(data_list)

# ============================================================================
# 极值函数 (Extreme Values)
# ============================================================================

def cal_max(data_list: List[float]) -> float:
    """
    计算最大值 (Maximum Value)
    
    最大值是数据集中的最大观测值，虽然不是稳健的统计量，但在临床研究中有重要意义，
    如记录最高血压、最大耐受剂量、最长生存期等。需要注意的是，最大值可能是个异常值，
    应结合其他统计量综合判断其临床意义和对分析的影响。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 数据集中的最大值
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute max of an empty list")
    return max(data_list)


def cal_min(data_list: List[float]) -> float:
    """
    计算最小值 (Minimum Value)
    
    最小值是数据集中的最小观测值，与最大值一样，虽不是稳健统计量，但在临床实践中很重要，
    如记录最低血糖、最小有效剂量、最短起效时间等。同样需要注意最小值可能是异常值，
    需要结合临床背景和数据分布进行判断，避免因录入错误或设备故障导致的异常值影响结论。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        float: 数据集中的最小值
        
    Raises:
        ValueError: 当输入列表为空时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute min of an empty list")
    return min(data_list)

# ============================================================================
# 推断统计函数 (Inferential Statistics)
# ============================================================================

def cal_mean_ci95(data_list: List[float]) -> tuple:
    """
    计算平均数的 95% 置信区间 (95% Confidence Interval for Mean)
    
    95% 置信区间提供了总体均值的可能范围，是推断统计的重要组成部分，
    在临床试验中用于估计治疗效果的真实值范围。如果区间不包含零（对于差值）或 1（对于比值），
    通常表明结果在α=0.05 的显著性水平下具有统计学意义。置信区间的宽度反映了估计的精确度，
    宽度越窄，估计越精确。使用 t 分布计算，适用于小样本或总体标准差未知的情况。
    
    Args:
        data_list: 包含数值的列表
        
    Returns:
        tuple: (下限，上限)，表示 95% 置信区间的边界值
        
    Raises:
        ValueError: 当输入列表为空或数据少于 2 个值时抛出异常
    """
    if not data_list:
        raise ValueError("Cannot compute CI of an empty list")
    
    n = len(data_list)
    if n < 2:
        raise ValueError("Confidence interval requires at least two values")
    
    mean_val = cal_mean(data_list)
    std_dev = cal_std(data_list)
    se = std_dev / math.sqrt(n)  # 标准误
    
    # 使用 t 分布计算 95% 置信区间
    # 自由度 = n - 1
    df = n - 1
    # 获取双侧 95% 置信区间的 t 临界值
    t_critical = stats.t.ppf(0.975, df)  # 0.975 对应双侧 95% 置信区间
    margin_error = t_critical * se
    lower_bound = round(mean_val - margin_error, 4)
    upper_bound = round(mean_val + margin_error, 4)
    
    return (lower_bound, upper_bound)

# ============================================================================
# 综合报告生成函数 (Comprehensive Report Generation)
# ============================================================================

def cal_result_des(param: BaseParamDes) -> Dict[str, Any]:
    """
    生成描述性统计量统计分析的完整报告字典
    
    此函数整合了描述性统计量的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的描述性统计分析结果。
    报告包括集中趋势、离散趋势、分布形状等信息，
    便于临床医生和研究人员快速理解数据的基本特征。
    
    【输出数据结构说明】:
    返回的字典包含以下键值对：
    - table_name: 表格名称
    - input_parameters: 输入参数信息
        - variable_count: 变量数量
        - variables: 变量名称列表
    - descriptive_statistics: 每个变量的描述性统计结果字典
    - summary: 统计摘要信息
        - total_sample_size: 总样本量
        - overall_mean: 总体均数
    - interpretation: 结果解读
        - statistical_interpretation: 统计学解读
        - data_quality_note: 数据质量提示
    - remark: 备注信息
    
    【临床应用场景】:
    1. 基线特征表（Table 1）生成：在 RCT 或观察性研究中描述研究人群特征
    2. 主要疗效指标总结：报告治疗前后关键指标的变化
    3. 安全性数据分析：描述实验室检查值、生命体征等安全性指标的分布
    4. 数据质量报告：识别数据分布特征、异常值和缺失模式
    5. 中期分析报告：在独立数据监查委员会（IDMC）会议中呈现累积数据
    
    Args:
        param: BaseParamDes 对象，包含 stats_data_list
        
    Returns:
        Dict[str, Any]: 包含描述性统计量统计分析指标的字典
        
    Raises:
        ValueError: 当输入参数无效或数据列表为空时抛出异常
    """
    # 从参数对象中提取数据
    stats_data_list = param.stats_data_list
    
    if not stats_data_list:
        raise ValueError("统计数据列表不能为空")
    
    # 计算每个数据集的描述性统计量
    all_stats = []
    all_descriptive_stats = {}
    
    for i, stats_data in enumerate(stats_data_list):
        if not isinstance(stats_data, StatsData):
            raise ValueError(f"第{i+1}个元素不是 StatsData 类型")
        
        if not stats_data.data_list:
            raise ValueError(f"第{i+1}个数据集为空")
        
        data_list = stats_data.data_list
        
        # 计算详细的描述性统计量
        stats_result = {
            'sample_size': len(data_list),
            'sum': round(sum(data_list), 4),
            'min': round(cal_min(data_list), 4),
            'max': round(cal_max(data_list), 4),
            'mean': round(cal_mean(data_list), 4),
            'mean_ci95': cal_mean_ci95(data_list),
            'geometric_mean': round(cal_geometric_mean(data_list), 4),
            'harmonic_mean': round(cal_harmonic_mean(data_list), 4),
            'trimmed_mean': round(cal_trimmed_mean(data_list), 4),
            'variance': round(cal_variance(data_list), 4),
            'std': round(cal_std(data_list), 4),
            'range': round(cal_range(data_list), 4),
            'standard_error': round(cal_se(data_list), 4),
            'cv': round(cal_cv(data_list), 4),
            'median': round(cal_median(data_list), 4),
            'q1': round(cal_q1(data_list), 4),
            'q3': round(cal_q3(data_list), 4),
            'iqr': round(cal_iqr(data_list), 4),
            'mode': cal_mode(data_list)
        }
        
        all_stats.append({
            'variable_name': stats_data.field_name,
            'statistics': stats_result
        })
        all_descriptive_stats[stats_data.field_name] = stats_result
    
    # 计算总体统计信息
    total_samples = sum(stats['statistics']['sample_size'] for stats in all_stats)
    total_sum = sum(stats['statistics']['mean'] * stats['statistics']['sample_size'] for stats in all_stats)
    overall_mean = total_sum / total_samples if total_samples > 0 else 0
    
    # 构建结果字典
    result_dict = {
        "table_name": "描述性统计量",
        "input_parameters": {
            "variable_count": len(stats_data_list),
            "variables": [data.field_name for data in stats_data_list]
        },
        "descriptive_statistics": all_descriptive_stats,
        "summary": {
            "total_sample_size": total_samples,
            "overall_mean": round(overall_mean, 4)
        },
        "interpretation": {
            "statistical_interpretation": f"共分析了{len(stats_data_list)}个变量，总计{total_samples}个观测值",
            "data_quality_note": "请检查数据中是否存在异常值或缺失值"
        },
        "remark": f"变量数={len(stats_data_list)}, 总样本量={total_samples}, 总体平均数={overall_mean:.4f}"
    }
    
    return result_dict