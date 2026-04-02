"""
单样本t检验统计分析模块 (Single Sample T-Test Statistics Module)

本模块提供单样本t检验功能，用于检验样本均数与已知总体均数之间是否存在显著性差异。
该检验适用于样本来自正态分布总体，且总体标准差未知的情况。

【模块功能概述】:
1. t统计量计算：计算单样本t检验的统计量
2. P值计算：基于t分布计算双侧P值
3. 显著性检验：判断统计结果的显著性
4. 结果解释：提供统计结果的专业解释

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
from typing import Dict, Any
from scipy.stats import t
from app.schemas.request_data.t_param import TParamSingle1


def single_sample_t_test(pop_mean: float, n: int, sample_mean: float, sample_std: float) -> Dict[str, Any]:
    """
    单样本t检验核心计算函数
    
    本函数实现单样本t检验的核心统计算法，计算t统计量和双侧P值。
    单样本t检验用于比较样本均数与已知总体均数之间的差异是否具有统计学意义。
    
    【统计学原理】:
    单样本t检验的零假设(H0)为：样本所代表的总体均数等于已知的总体均数(μ = μ0)。
    备择假设(H1)为：样本所代表的总体均数不等于已知的总体均数(μ ≠ μ0)。
    
    t统计量计算公式：
    t = (x̄ - μ0) / (s / √n)
    
    其中：
    - x̄: 样本均数
    - μ0: 假设的总体均数
    - s: 样本标准差
    - n: 样本量
    
    自由度(df) = n - 1
    
    【适用场景】:
    1. 医学研究：检验某治疗方案的效果是否与标准值有差异
       例：检验新降压药治疗后患者收缩压是否与目标值120mmHg有显著差异
   
    2. 质量控制：检验生产线产品指标是否符合标准
       例：检验某批次药片有效成分含量是否符合标示量100mg
   
    3. 诊断标准验证：检验某项检测指标是否偏离正常参考范围
       例：检验糖尿病患者空腹血糖是否显著高于正常上限6.1mmol/L
   
    4. 流行病学调查：检验某地区人群特征是否与全国平均水平有差异
       例：检验某城市成年人平均BMI是否与全国平均值24有显著差异
    
    【前提条件】:
    1. 独立性：样本观测值相互独立
    2. 正态性：样本来自正态分布总体（小样本时尤为重要）
    3. 连续性：数据为连续型变量
    
    【结果解读】:
    - P < 0.05：拒绝零假设，认为样本均数与总体均数差异有统计学意义
    - P ≥ 0.05：不拒绝零假设，尚不能认为样本均数与总体均数有差异
    
    【注意事项】:
    1. 小样本(n<30)时需严格检验正态性假设
    2. 大样本时可适当放宽正态性要求（中心极限定理）
    3. 存在离群值时需谨慎使用，考虑非参数检验替代
    4. 统计显著性不等同于临床/实际意义
    
    Args:
        pop_mean (float): 总体均数（假设值），即零假设中的μ0
        n (int): 样本量，必须大于1
        sample_mean (float): 样本均数，即样本观测值的平均值
        sample_std (float): 样本标准差，必须大于0
    
    Returns:
        Dict[str, Any]: 包含单样本t检验统计量和结果的字典，结构如下：
            - input_parameters: 输入参数
                - population_mean: 总体均数
                - sample_size: 样本量
                - sample_mean: 样本均数
                - sample_std: 样本标准差
            - test_statistics: 检验统计量
                - degrees_of_freedom: 自由度
                - t_value: t统计量值
                - p_value_two_sided: 双侧P值
            - significance_tests: 显著性检验结果
                - p_greater_than_0_05: P值是否大于0.05
                - p_greater_than_0_01: P值是否大于0.01
                - significant_at_05: 0.05水平下的显著性结论
                - significant_at_01: 0.01水平下的显著性结论
            - interpretation: 统计解释
                - t_test_interpretation: 0.05水平的解释说明
                - t_test_interpretation_01: 0.01水平的解释说明
    
    Raises:
        ValueError: 当样本量≤1或样本标准差≤0时抛出异常
    
    Example:
        >>> result = single_sample_t_test(pop_mean=100, n=25, sample_mean=105.5, sample_std=8.2)
        >>> print(result['test_statistics']['t_value'])
        3.353...
        >>> print(result['test_statistics']['p_value_two_sided'])
        0.002...
    """
    # 参数验证：确保样本量大于1
    if n <= 1:
        raise ValueError("样本量必须大于1，否则无法计算样本标准差和进行统计推断")
    
    # 参数验证：确保样本标准差大于0
    if sample_std <= 0:
        raise ValueError("样本标准差必须大于0，标准差为0表示所有观测值相同，无法进行变异分析")
    
    # 计算自由度：单样本t检验的自由度为n-1
    # 自由度反映了可用于估计总体参数的独立信息数量
    degrees_of_freedom = n - 1
    
    # 计算t统计量
    # 公式：t = (样本均数 - 总体均数) / 标准误
    # 标准误 = 样本标准差 / sqrt(样本量)
    # t值表示样本均数与总体均数的差异相对于抽样误差的大小
    standard_error = sample_std / math.sqrt(n)
    t_statistic = (sample_mean - pop_mean) / standard_error
    
    # 计算双侧P值
    # P值表示在零假设成立的前提下，观察到当前或更极端结果的概率
    # 双侧检验考虑了两个方向的极端情况（大于或小于）
    # P = 2 × P(T ≥ |t|)，其中T服从自由度为df的t分布
    p_value_two_sided = 2 * (1 - t.cdf(abs(t_statistic), df=degrees_of_freedom))
    
    # 显著性判断：常用的两个显著性水平
    # α=0.05：通常作为统计学显著性的标准
    # α=0.01：更严格的显著性标准，用于需要更高置信度的情况
    is_greater_than_05 = p_value_two_sided > 0.05
    is_greater_than_01 = p_value_two_sided > 0.01
    
    # 构建并返回结果字典
    return {
        "input_parameters": {
            "population_mean": float(pop_mean),
            "sample_size": int(n),
            "sample_mean": float(sample_mean),
            "sample_std": float(sample_std)
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


def perform_single_sample_t_test(pop_mean: float, n: int, sample_mean: float, sample_std: float) -> Dict[str, Any]:
    """
    执行单样本t检验并返回完整结果
    
    本函数作为单样本t检验的主入口函数，负责参数验证并调用核心计算函数。
    它提供了更友好的接口，适合在业务逻辑层直接调用。
    
    【函数职责】:
    1. 参数有效性验证：检查样本量和标准差的合理性
    2. 调用核心算法：委托single_sample_t_test进行实际计算
    3. 异常处理：对无效输入提供清晰的错误提示
    4. 结果返回：返回标准化的检验结果字典
    
    【应用场景】:
    1. 临床研究数据分析：快速评估治疗效果
    2. 质量检测报告生成：自动化产品质量评估
    3. 教学演示：展示统计检验过程
    4. 科研数据处理：批量分析实验数据
    
    【与其他函数的关系】:
    - 依赖 single_sample_t_test 进行核心统计计算
    - 被 cal_result_t_single1 调用以生成完整报告
    
    Args:
        pop_mean (float): 总体均数（假设值），代表理论值或标准值
        n (int): 样本量，即参与分析的观测值数量
        sample_mean (float): 样本均数，由实际观测数据计算得到
        sample_std (float): 样本标准差，反映数据的离散程度
    
    Returns:
        Dict[str, Any]: 包含完整检验结果的字典，结构与single_sample_t_test返回值相同
            包含输入参数、检验统计量、显著性检验结果和统计解释四个部分
    
    Raises:
        ValueError: 当输入参数不满足检验要求时抛出
            - 样本量≤1：无法估计总体方差
            - 标准差≤0：数据无变异或输入错误
    
    See Also:
        single_sample_t_test: 核心统计算法实现
        cal_result_t_single1: 生成标准化报告字典
    
    Example:
        >>> result = perform_single_sample_t_test(
        ...     pop_mean=75.0,
        ...     n=30,
        ...     sample_mean=78.5,
        ...     sample_std=6.2
        ... )
        >>> if result['significance_tests']['significant_at_05'] == '显著':
        ...     print("样本均数与总体均数存在显著差异")
    """
    # 参数验证：样本量必须大于1
    # 原因：自由度=n-1，若n=1则自由度为0，无法进行统计推断
    # 此外，单个观测值无法估计总体变异
    if n <= 1:
        raise ValueError("样本量必须大于1，单个观测值无法进行统计推断和方差估计")
    
    # 参数验证：样本标准差必须大于0
    # 原因：标准差为0意味着所有观测值完全相同，不存在变异
    # 这种情况下无法计算标准误，也无法进行假设检验
    if sample_std <= 0:
        raise ValueError("样本标准差必须大于0，标准差为0表示数据无变异，无法进行统计检验")
    
    # 调用核心计算函数执行单样本t检验
    # 将经过验证的参数传递给single_sample_t_test函数
    results = single_sample_t_test(pop_mean, n, sample_mean, sample_std)
    
    # 返回完整的检验结果
    return results


def cal_result_t_single1(param: TParamSingle1) -> Dict[str, Any]:
    """
    生成单样本t检验分析的完整报告字典

    Args:
        param: TParamSingle1对象，包含pop_mean, n, sample_mean, sample_std
    
    Returns:
        Dict[str, Any]: 包含t检验分析指标的字典
    """
    pop_mean = param.pop_mean
    n = param.n
    sample_mean = param.sample_mean
    sample_std = param.sample_std

    results = perform_single_sample_t_test(pop_mean, n, sample_mean, sample_std)
    
    result_dict = {
        "table_name": "单样本t检验分析",
        "input_parameters": results["input_parameters"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"总体均数：{pop_mean}, 样本量：{n}, 样本均数：{sample_mean:.4f}, 样本标准差：{sample_std:.4f}"
    }
    
    # 返回完整的报告字典
    # 该字典可直接用于JSON序列化、数据库存储或前端展示
    return result_dict