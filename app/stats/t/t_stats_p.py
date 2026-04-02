"""
t检验P值计算模块 (T-Test P-Value Calculation Module)

本模块提供基于t统计量和自由度计算P值的功能。
该模块适用于已知t统计量和自由度，需要计算相应P值的情况。

【模块功能概述】:
1. 双侧P值计算：基于t分布计算双侧检验的P值
2. 单侧P值计算：计算左右单侧检验的P值
3. 显著性判断：判断统计结果在特定显著性水平下的显著性
4. 结果解释：提供统计结果的专业解释

【临床应用价值】:
- 从已知的t统计量计算P值：当已有t值时快速获得显著性
- 验证统计结果：复核已有统计软件的计算结果
- 教学演示：展示P值的计算过程和意义
- 文献数据处理：处理文献中报告的t值和自由度

【统计方法选择指南】:
1. t分布P值计算适用条件：
   - 已知t统计量值
   - 已知自由度
   - 需要精确的P值

2. 临床应用场景：
   - 计算发表文献中t检验的P值
   - 验证其他软件的计算结果
   - 生成统计报告
   - 教学和培训

【结果解读注意事项】:
1. P值解释：在零假设成立的情况下，观察到当前或更极端结果的概率
2. 显著性水平：通常使用α=0.05作为判断标准
3. 双侧vs单侧：根据研究假设选择适当的检验类型
4. 自由度的重要性：影响t分布的形状和P值的计算

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "如何根据t值计算P值？"
- "双侧和单侧P值有什么区别？"
- "自由度如何影响P值？"
- "t分布的基本性质是什么？"
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
from scipy.stats import t


def calculate_p_value_from_t(t_value: float, df: int) -> Dict:
    """
    根据t值和自由度计算P值
    
    参数:
    - t_value: t统计量值
    - df: 自由度
    
    返回:
    - 包含P值计算结果的字典
    """
    # 参数验证
    if df <= 0:
        raise ValueError("自由度必须大于0")
    
    if not isinstance(t_value, (int, float)):
        raise ValueError("t值必须是数字")
    
    # 计算双侧P值
    # P值 = 2 × P(T ≥ |t|)
    p_value_two_sided = 2 * (1 - t.cdf(abs(t_value), df=df))
    
    # 计算单侧P值（右侧）
    # P值 = P(T ≥ t)
    p_value_one_sided_right = 1 - t.cdf(t_value, df=df)
    
    # 计算单侧P值（左侧）
    # P值 = P(T ≤ t)
    p_value_one_sided_left = t.cdf(t_value, df=df)
    
    # 显著性判断
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    return {
        "input_parameters": {
            "t_value": float(t_value),
            "degrees_of_freedom": int(df)
        },
        "p_values": {
            "two_sided": float(p_value_two_sided),
            "one_sided_right": float(p_value_one_sided_right),
            "one_sided_left": float(p_value_one_sided_left)
        },
        "significance_tests": {
            "p_greater_than_0_05": is_greater_than_05,
            "p_greater_than_0_01": is_greater_than_01,
            "significant_at_05": "不显著" if is_greater_than_05 else "显著",
            "significant_at_01": "不显著" if is_greater_than_01 else "显著"
        },
        "interpretation": {
            "two_sided_interpretation": f"双侧检验{'不' if p_value_two_sided > 0.05 else ''}显著 (p={'>' if p_value_two_sided > 0.05 else '≤'}0.05)",
            "two_sided_interpretation_01": f"在0.01水平下{'不' if p_value_two_sided > 0.01 else ''}显著 (p={'>' if p_value_two_sided > 0.01 else '≤'}0.01)"
        }
    }


def perform_t_p_calculation(t_value: float, df: int) -> Dict:
    """
    根据给定参数执行P值计算
    
    本函数基于t统计量和自由度计算各种类型的P值，包括双侧和单侧检验。
    适用于需要快速从已知t值获取显著性水平的场景。
    
    【统计学原理】:
    - t分布：当样本量较小且总体标准差未知时使用的概率分布
    - P值：在零假设成立的前提下，观察到当前统计量或更极端情况的概率
    - 双侧检验：考虑两个方向的极端情况，P值为单侧的两倍
    - 单侧检验：只考虑一个方向的极端情况
    
    【参数说明】:
    - t_value: t统计量，表示样本均值与总体均值的标准化差异
    - df: 自由度，通常为样本量减1，影响t分布的形状
    
    【返回值结构】:
    返回字典包含以下关键信息：
    - input_parameters: 输入的t值和自由度
    - p_values: 三种类型的P值（双侧、右侧单侧、左侧单侧）
    - significance_tests: 在0.05和0.01水平下的显著性判断
    - interpretation: 对结果的统计学解释
    
    【临床应用场景】:
    1. 文献数据再分析：从已发表的t值和自由度计算精确P值
    2. 结果验证：核对其他统计软件的输出结果
    3. 教学示范：展示t检验的计算过程
    4. 快速评估：在已有t值时快速判断显著性
    
    【注意事项】:
    - 自由度必须为正整数
    - t值可以为正数或负数
    - 双侧检验适用于无方向性假设
    - 单侧检验适用于有明确方向性的研究假设
    
    Args:
        t_value (float): t统计量值
        df (int): 自由度，必须大于0
    
    Returns:
        Dict: 包含P值计算完整结果的字典，结构如下：
            {
                "input_parameters": {"t_value": float, "degrees_of_freedom": int},
                "p_values": {"two_sided": float, "one_sided_right": float, "one_sided_left": float},
                "significance_tests": {"p_greater_than_0_05": bool, "p_greater_than_0_01": bool, 
                                       "significant_at_05": str, "significant_at_01": str},
                "interpretation": {"two_sided_interpretation": str, "two_sided_interpretation_01": str}
            }
    
    Raises:
        ValueError: 当df<=0或t_value不是数字时抛出异常
    
    Example:
        >>> result = perform_t_p_calculation(2.5, 20)
        >>> print(result["p_values"]["two_sided"])
        0.021315
    """
    # 验证参数
    if df <= 0:
        raise ValueError("自由度必须大于0")
    
    if not isinstance(t_value, (int, float)):
        raise ValueError("t值必须是数字")
    
    # 执行P值计算
    results = calculate_p_value_from_t(t_value, df)
    
    return results


def cal_result_t_p(t_value: float, df: int) -> Dict[str, Any]:
    """
    生成t检验P值计算的完整报告字典
    
    此函数整合了t检验P值计算的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的P值计算结果。报告包括输入参数、
    P值计算结果和统计解释等信息，
    便于临床医生和研究人员快速理解P值的特征。
    
    Args:
        t_value: t统计量值，浮点数类型
        df: 自由度，整数类型
    
    Returns:
        Dict[str, Any]: 包含P值计算指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"t检验P值计算"
            - input_parameters: 输入参数信息
                - t_value: t值，浮点数类型
                - degrees_of_freedom: 自由度，整数类型
            - p_values: P值计算结果
                - two_sided: 双侧P值，浮点数类型
                - one_sided_right: 右侧单侧P值，浮点数类型
                - one_sided_left: 左侧单侧P值，浮点数类型
            - significance_tests: 显著性检验结果
                - p_greater_than_0_05: P值是否大于0.05，布尔类型
                - p_greater_than_0_01: P值是否大于0.01，布尔类型
                - significant_at_05: 0.05水平显著性，字符串类型
                - significant_at_01: 0.01水平显著性，字符串类型
            - interpretation: 统计解释
                - two_sided_interpretation: 双侧检验解释，字符串类型
                - two_sided_interpretation_01: 0.01水平解释，字符串类型
            - remark: 备注信息，字符串类型
    """
    # 从参数对象解构
    t_value = param.t_value
    df = param.df

    # 执行P值计算
    results = perform_t_p_calculation(t_value, df)
    
    # 构建结果字典
    result_dict = {
        "table_name": "t检验P值计算",
        "input_parameters": results["input_parameters"],
        "p_values": results["p_values"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"t值: {t_value}, 自由度: {df}"
    }
    
    return result_dict