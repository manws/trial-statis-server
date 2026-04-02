"""
临床Wilcoxon符号秩检验统计分析模块 (Clinical Wilcoxon Signed-Rank Test Statistics Module)

本模块提供全面的Wilcoxon符号秩检验统计分析功能，用于临床数据中配对样本的比较，
是一种重要的非参数统计分析方法。Wilcoxon符号秩检验通过比较配对数据的差值的秩次，
来判断配对样本是否存在显著差异，广泛应用于医学研究中的前后对照实验、
配对设计研究、重复测量分析等领域。

【模块功能概述】:
1. 差值计算：计算配对样本的差值
2. 符号秩统计：计算正负差值的秩和
3. 检验统计量：计算Wilcoxon检验统计量W
4. 显著性检验：执行正态近似检验判断统计量的显著性
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 前后对照实验：比较患者治疗前后的指标变化
- 配对设计研究：比较不同处理在相同受试者上的效果
- 重复测量分析：分析同一受试者多次测量结果的差异
- 交叉设计分析：比较交叉试验中不同阶段的治疗效果

【统计方法选择指南】:
1. Wilcoxon符号秩检验适用条件：
   - 数据为连续型变量
   - 数据分布不要求正态性
   - 配对观测值相互独立
   - 差值分布对称

2. 临床应用场景：
   - 比较患者治疗前后血压的差异
   - 评估某种药物治疗前后血糖水平的变化
   - 分析手术前后某项指标的改变
   - 研究不同治疗方案在同一患者的疗效差异

【结果解读注意事项】:
1. 检验统计量解释：W值越小，表明配对样本的差异越显著
2. P值解释：P值小于显著性水平（如0.05）时拒绝原假设
3. 临床意义：统计显著性不等于临床意义，需结合专业知识进行解释

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用Wilcoxon符号秩检验吗？"
- "如何解释Wilcoxon符号秩检验的结果？"
- "Wilcoxon符号秩检验的前提条件是什么？"
- "如何判断差异是否具有统计学意义？"

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
from typing import Dict, List
from scipy.stats import norm, rankdata
import numpy as np


def wilcoxon_signed_rank_test(data1: List[float], data2: List[float]) -> Dict:
    """
    Wilcoxon符号秩检验
    
    对两组配对数据进行Wilcoxon符号秩检验，用于判断配对样本是否存在显著差异。
    该检验是一种非参数检验方法，不要求数据服从正态分布，适用于小样本或
    分布未知的配对数据比较。
    
    统计学原理：
    1. 计算每对数据的差值
    2. 去除差值为0的配对
    3. 对差值的绝对值进行秩次排列
    4. 分别计算正差值和负差值的秩和
    5. 使用较小的秩和作为检验统计量W
    6. 通过正态近似计算Z值和P值
    
    临床研究应用场景：
    - 比较患者治疗前后的血压、血糖等指标变化
    - 评估手术前后某项生理指标的改变
    - 分析同一受试者接受不同治疗方案的效果差异
    - 重复测量设计中时间点间的比较
    
    结果解读指导：
    - W统计量：较小的秩和，值越小表示差异越显著
    - Z值：标准化检验统计量，用于计算P值
    - P值：双侧检验的显著性概率，P<0.05表示差异具有统计学意义
    - 注意：统计学显著性不等同于临床意义，需结合专业背景判断
    
    方法学局限性：
    - 要求配对观测值相互独立
    - 差值分布应大致对称
    - 对于大量结（相同绝对值）的情况，检验效能可能降低
    - 小样本时正态近似可能不够准确
    
    Args:
        data1: List[float]，第一组配对数据，如治疗前测量值
        data2: List[float]，第二组配对数据，如治疗后测量值
    
    Returns:
        Dict: 包含检验结果的字典，包括：
            - input_parameters: 输入参数信息（总配对数、有效配对数等）
            - sample_statistics: 样本统计量（均值、标准差、均数差等）
            - rank_statistics: 秩统计量（正负秩和、检验统计量W等）
            - test_statistics: 检验统计量（Z值、双侧P值）
            - significance_tests: 显著性检验结果（0.05和0.01水平）
            - interpretation: 统计解释（中文描述）
    
    Raises:
        ValueError: 当数据为空、长度不等或所有差值都为0时抛出
    """
    # 参数验证
    if not data1 or not data2:
        raise ValueError("两组配对数据都不能为空")
    
    if len(data1) != len(data2):
        raise ValueError("两组配对数据的长度必须相等")
    
    # 计算配对差值
    paired_differences = [data1[i] - data2[i] for i in range(len(data1))]
    
    # 移除差值为0的数据点
    non_zero_diffs = [diff for diff in paired_differences if diff != 0]
    
    if len(non_zero_diffs) == 0:
        raise ValueError("所有差值都为0，无法进行检验")
    
    n = len(non_zero_diffs)
    
    # 计算绝对值和符号
    abs_diffs = [abs(diff) for diff in non_zero_diffs]
    signs = [1 if diff > 0 else -1 for diff in non_zero_diffs]
    
    # 计算秩
    ranks = rankdata(abs_diffs, method='average')
    
    # 计算正秩和和负秩和
    positive_ranks = [ranks[i] for i in range(n) if signs[i] > 0]
    negative_ranks = [ranks[i] for i in range(n) if signs[i] < 0]
    
    W_plus = sum(positive_ranks)
    W_minus = sum(negative_ranks)
    
    # 使用较小的秩和作为检验统计量
    W_statistic = min(W_plus, W_minus)
    
    # 计算期望值和方差（无结时）
    E_W = n * (n + 1) / 4
    Var_W = n * (n + 1) * (2 * n + 1) / 24
    
    # 如果存在结（相同绝对值），需要调整方差
    # 计算结的修正因子
    unique_abs_diffs, counts = np.unique(abs_diffs, return_counts=True)
    tie_correction = sum([count**3 - count for count in counts if count > 1])
    if tie_correction > 0:
        Var_W -= tie_correction / 48
    
    # 计算标准化检验统计量Z
    if Var_W > 0:
        Z_statistic = (W_statistic - E_W) / math.sqrt(Var_W)
    else:
        Z_statistic = 0.0
    
    # 计算双侧P值
    # P值 = 2 × P(Z ≥ |z|)
    p_value_two_sided = 2 * (1 - norm.cdf(abs(Z_statistic)))
    
    # 显著性判断
    is_less_than_05 = p_value_two_sided < 0.05
    is_less_than_01 = p_value_two_sided < 0.01
    
    return {
        "input_parameters": {
            "total_pairs": len(data1),
            "non_zero_pairs": n,
            "zero_differences": len(data1) - n
        },
        "sample_statistics": {
            "group1_mean": float(np.mean(data1)),
            "group1_std": float(np.std(data1)),
            "group2_mean": float(np.mean(data2)),
            "group2_std": float(np.std(data2)),
            "mean_difference": float(np.mean(paired_differences))
        },
        "rank_statistics": {
            "positive_ranks_sum": float(W_plus),
            "negative_ranks_sum": float(W_minus),
            "test_statistic_W": float(W_statistic),
            "expected_W": float(E_W),
            "variance_W": float(Var_W)
        },
        "test_statistics": {
            "z_value": float(Z_statistic),
            "p_value_two_sided": float(p_value_two_sided)
        },
        "significance_tests": {
            "p_less_than_0_05": is_less_than_05,
            "p_less_than_0_01": is_less_than_01,
            "significant_at_05": "显著" if is_less_than_05 else "不显著",
            "significant_at_01": "显著" if is_less_than_01 else "不显著"
        },
        "interpretation": {
            "wilcoxon_interpretation": f"Wilcoxon符号秩检验{'不' if p_value_two_sided >= 0.05 else ''}显著 (p={'≥' if p_value_two_sided >= 0.05 else '<'}0.05)",
            "wilcoxon_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided >= 0.01 else ''}显著 (p={'≥' if p_value_two_sided >= 0.01 else '<'}0.01)"
        }
    }


def cal_result_rs_paired(data1: List[float], data2: List[float]) -> Dict:
    """
    生成Wilcoxon符号秩检验统计分析的完整报告字典
    
    此函数整合了Wilcoxon符号秩检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，便于临床医生和研究人员快速理解
    Wilcoxon符号秩检验的特征。
    
    函数特点：
    - 输入简洁：仅需两个配对的浮点数列表
    - 输出完整：包含从原始数据到统计解释的全方位信息
    - 格式规范：符合临床研究报告的标准格式要求
    - 易于集成：返回字典格式，便于与其他系统对接
    
    临床应用示例：
    1. 药物疗效评估：比较患者用药前后的症状评分
    2. 诊断方法比较：同一患者使用两种检测方法的结果对比
    3. 康复效果分析：康复训练前后功能评分的变化
    4. 生活方式干预：干预前后生化指标的改变
    
    输出指标说明：
    - table_name: 分析报告的名称标识
    - input_parameters: 数据基本情况，包括有效样本量
    - sample_statistics: 描述性统计量，了解数据分布特征
    - rank_statistics: 秩次相关统计量，反映检验的核心计算过程
    - test_statistics: 推断性统计量，用于假设检验决策
    - significance_tests: 不同显著性水平下的判断结果
    - interpretation: 中文统计结论，便于非统计专业人员理解
    - remark: 补充说明信息
    
    Args:
        data1: List[float]，第一组配对数据，如基线测量值
        data2: List[float]，第二组配对数据，如随访测量值
    
    Returns:
        Dict: 包含Wilcoxon符号秩检验统计分析指标的字典，结构如下：
            {
                "table_name": str,              # 表格名称
                "input_parameters": dict,       # 输入参数
                "sample_statistics": dict,      # 样本统计量
                "rank_statistics": dict,        # 秩统计量
                "test_statistics": dict,        # 检验统计量
                "significance_tests": dict,     # 显著性检验
                "interpretation": dict,         # 统计解释
                "remark": str                   # 备注信息
            }
    
    Raises:
        ValueError: 当输入数据不符合检验要求时抛出异常
    """
    # 执行Wilcoxon符号秩检验
    results = wilcoxon_signed_rank_test(data1, data2)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Wilcoxon符号秩检验分析",
        "input_parameters": results["input_parameters"],
        "sample_statistics": results["sample_statistics"],
        "rank_statistics": results["rank_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"配对数据点数：{len(data1)}, 有效配对数：{results['input_parameters']['non_zero_pairs']}"
    }
    
    return result_dict