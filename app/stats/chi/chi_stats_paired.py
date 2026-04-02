"""
临床配对资料卡方检验统计分析模块 (Clinical Paired Chi-Square Test Statistics Module)

本模块提供全面的配对资料卡方检验统计分析功能，用于临床数据中配对分类变量之间的关联性检验，
特别是McNemar检验。该检验在医学研究中广泛应用，如比较同一受试者在不同时间点或
不同处理下的分类结果，评估诊断试验的一致性等。

【模块功能概述】:
1. McNemar检验：计算配对资料的卡方统计量和P值
2. 校正McNemar检验：对小样本数据进行连续性校正
3. Kappa一致性系数：评估配对分类变量的一致性
4. 优势比计算：计算并估计优势比及其置信区间
5. 不一致分析：分析配对数据中的不一致情况
6. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 诊断试验：评估同一诊断方法在不同时间点的一致性
- 治疗效果：比较同一患者治疗前后的分类结果
- 问卷评估：评估同一患者不同时间的心理量表结果
- 影像诊断：比较不同影像方法对同一病变的判断一致性

【统计方法选择指南】:
1. 配对卡方检验适用条件：
   - 配对设计的分类数据
   - 同一受试者在不同条件下的分类结果
   - 需要评估配对变量间的一致性或差异性

2. 临床应用场景：
   - 比较患者治疗前后的分类结果（如症状有无、疾病分级）
   - 评估两种诊断方法在同一批患者中的一致性
   - 分析同一诊断方法在不同时间点的可靠性

【结果解读注意事项】:
1. McNemar检验解释：检验配对变量间是否存在统计学显著差异
2. Kappa系数解释：评估一致性强度，而非单纯的关联性
3. 校正判断：当b+c<40时使用校正检验
4. P值解释：在零假设（无差异）成立的情况下，观察到当前或更极端结果的概率
5. 临床意义：统计学显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的配对分类数据适合用McNemar检验吗？"
- "如何解释McNemar检验的结果？"
- "Kappa系数如何评估一致性？"
- "何时需要使用校正检验？"
- "如何判断结果的临床意义？"

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
from scipy.stats import chi2, norm
import numpy as np
from app.schemas.request_data.chi_param import ChiSquareParamPaired

def paired_chi_square_test(a: int, b: int, c: int, d: int) -> Dict[str, Any]:
    """
    配对资料卡方检验（McNemar检验）
    
    在临床研究中，这用于检验配对分类变量之间的关联性，
    如比较同一受试者在不同时间点或不同处理下的分类结果。
    
    Args:
        a: 样本1阳性/成功，样本2阳性/成功的配对数
        b: 样本1阳性/成功，样本2阴性/失败的配对数
        c: 样本1阴性/失败，样本2阳性/成功的配对数
        d: 样本1阴性/失败，样本2阴性/失败的配对数
        
    Returns:
        Dict[str, Any]: 包含配对卡方检验各种统计量和结果的字典
        
    Raises:
        ValueError: 当配对资料中的值为负数时抛出异常
    """
    # 参数验证
    if any(x < 0 for x in [a, b, c, d]):
        raise ValueError("配对资料中的所有值必须为非负数")
    
    # 计算边际总计
    n1_positive = a + b  # 样本1阳性总数
    n1_negative = c + d  # 样本1阴性总数
    n2_positive = a + c  # 样本2阳性总数
    n2_negative = b + d  # 样本2阴性总数
    total_pairs = a + b + c + d  # 总配对数
    
    # 判断是否需要校正
    # 通常当b+c<25时建议使用校正，这里采用b+c<40的标准
    discordant_sum = b + c
    needs_correction = discordant_sum < 40
    
    # 计算McNemar卡方统计量（不校正）
    if b + c == 0:
        # 处理分母为0的情况
        chi_square = float('inf')
        p_value = 0.0
    else:
        chi_square = (abs(b - c) - 0)**2 / (b + c)
        p_value = 1 - chi2.cdf(chi_square, df=1)
    
    # 计算校正的McNemar卡方统计量
    if b + c == 0:
        corrected_chi_square = float('inf')
        corrected_p_value = 0.0
    else:
        # McNemar校正公式: (|b-c|-1)²/(b+c)
        corrected_chi_square = (abs(b - c) - 1)**2 / (b + c)
        corrected_p_value = 1 - chi2.cdf(corrected_chi_square, df=1)
    
    # 计算Kappa一致性系数
    # Kappa = (Po - Pe) / (1 - Pe)
    # 其中 Po = (a+d)/n，Pe = [(a+b)(a+c) + (c+d)(b+d)]/n²
    
    if total_pairs == 0:
        kappa = float('nan')
        kappa_se = float('nan')
        kappa_lower_ci = float('nan')
        kappa_upper_ci = float('nan')
    else:
        po = (a + d) / total_pairs  # 观察一致率
        pe = ((n1_positive * n2_positive) + (n1_negative * n2_negative)) / (total_pairs ** 2)  # 期望一致率
        
        if pe == 1:
            kappa = float('nan')
            kappa_se = float('nan')
            kappa_lower_ci = float('nan')
            kappa_upper_ci = float('nan')
        else:
            kappa = (po - pe) / (1 - pe)
            
            # 计算Kappa的标准误
            # 使用 delta method 近似
            if (1 - pe) ** 2 == 0:
                kappa_se = float('inf')
            else:
                numerator = po * (1 - po)
                denominator = total_pairs * (1 - pe) ** 2
                kappa_se = math.sqrt(numerator / denominator) if denominator > 0 else float('inf')
            
            # 计算95%置信区间
            if math.isfinite(kappa_se):
                z_critical = 1.96
                kappa_lower_ci = kappa - z_critical * kappa_se
                kappa_upper_ci = kappa + z_critical * kappa_se
            else:
                kappa_lower_ci = float('nan')
                kappa_upper_ci = float('nan')
    
    # 计算优势比(OR)及其置信区间
    # 对于配对资料，OR = b/c
    if c == 0:
        if b == 0:
            odds_ratio = float('nan')  # 无法计算
        else:
            odds_ratio = float('inf')  # b/c → ∞
        or_lower_ci = float('nan')
        or_upper_ci = float('nan')
    elif b == 0:
        odds_ratio = 0.0  # 0/c = 0
        or_lower_ci = float('nan')
        or_upper_ci = float('nan')
    else:
        odds_ratio = b / c
        
        # 计算OR的95%置信区间（使用自然对数变换）
        ln_or = math.log(odds_ratio)
        # 配对资料OR的方差: Var(ln OR) = 1/b + 1/c
        se_ln_or = math.sqrt(1/b + 1/c)
        z_critical = 1.96  # 95%置信水平
        
        ln_or_lower = ln_or - z_critical * se_ln_or
        ln_or_upper = ln_or + z_critical * se_ln_or
        
        or_lower_ci = math.exp(ln_or_lower)
        or_upper_ci = math.exp(ln_or_upper)
    
    return {
        "input_parameters": {
            "a": a, "b": b, "c": c, "d": d,
            "n1_positive": n1_positive, "n1_negative": n1_negative,
            "n2_positive": n2_positive, "n2_negative": n2_negative,
            "total_pairs": total_pairs
        },
        "discordant_analysis": {
            "b_value": b,
            "c_value": c,
            "discordant_sum": discordant_sum,
            "needs_correction": needs_correction,
            "correction_reason": "建议使用校正" if needs_correction else "无需校正"
        },
        "mcnemar_test": {
            "degrees_of_freedom": 1,
            "chi_square_value": chi_square,
            "p_value_two_sided": p_value
        },
        "corrected_mcnemar_test": {
            "degrees_of_freedom": 1,
            "chi_square_value": corrected_chi_square,
            "p_value_two_sided": corrected_p_value,
            "used_for_analysis": needs_correction
        },
        "kappa_coefficient": {
            "kappa_value": kappa,
            "kappa_se": kappa_se,
            "kappa_95_ci_lower": kappa_lower_ci,
            "kappa_95_ci_upper": kappa_upper_ci,
            "interpretation": interpret_kappa(kappa)
        },
        "odds_ratio": {
            "or_value": odds_ratio,
            "or_95_ci_lower": or_lower_ci,
            "or_95_ci_upper": or_upper_ci
        },
        "interpretation": {
            "mcnemar_interpretation": f"McNemar检验{'不' if p_value > 0.05 else ''}显著 (p={'>' if p_value > 0.05 else '≤'}0.05)",
            "corrected_interpretation": f"校正检验{'不' if corrected_p_value > 0.05 else ''}显著 (p={'>' if corrected_p_value > 0.05 else '≤'}0.05)" if needs_correction else "未使用校正检验"
        }
    }


def interpret_kappa(kappa: float) -> str:
    """
    解释Kappa一致性系数的意义
    
    Args:
        kappa: Kappa系数值
    
    Returns:
        解释字符串
    """
    if math.isnan(kappa):
        return "无法解释（计算异常）"
    elif kappa < 0:
        return "一致性差于随机水平"
    elif kappa < 0.20:
        return "一致性极弱"
    elif kappa < 0.40:
        return "一致性较弱"
    elif kappa < 0.60:
        return "一致性中等"
    elif kappa < 0.80:
        return "一致性较强"
    elif kappa < 0.90:
        return "一致性很强"
    else:
        return "一致性极强"


def perform_paired_chi_square_test(a: int, b: int, c: int, d: int) -> Dict[str, Any]:
    """
    根据给定参数执行配对资料卡方检验
    
    Args:
        a: 样本1阳性/成功，样本2阳性/成功的配对数
        b: 样本1阳性/成功，样本2阴性/失败的配对数
        c: 样本1阴性/失败，样本2阳性/成功的配对数
        d: 样本1阴性/失败，样本2阴性/失败的配对数
    
    Returns:
        Dict[str, Any]: 包含配对资料卡方检验完整结果的字典
    """
    # 验证参数
    if any(x < 0 for x in [a, b, c, d]):
        raise ValueError("配对资料中的所有值必须为非负数")
    
    # 执行配对资料卡方检验
    results = paired_chi_square_test(a, b, c, d)
    
    return results


def cal_result_chi_paired(a: int, b: int, c: int, d: int) -> Dict[str, Any]:
    """
    生成配对资料卡方检验统计分析的完整报告字典
    
    此函数整合了配对资料卡方检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的配对卡方检验结果。报告包括输入参数、
    不一致分析、McNemar检验、校正检验、Kappa一致性系数、优势比分析
    和统计解释等信息，便于临床医生和研究人员快速理解配对卡方检验的特征。
    
    Args:
        a: 样本1阳性/成功，样本2阳性/成功的配对数，整数类型
        b: 样本1阳性/成功，样本2阴性/失败的配对数，整数类型
        c: 样本1阴性/失败，样本2阳性/成功的配对数，整数类型
        d: 样本1阴性/失败，样本2阴性/失败的配对数，整数类型
    
    Returns:
        Dict[str, Any]: 包含配对资料卡方检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"配对资料卡方检验分析"
            - input_parameters: 输入参数信息
            - discordant_analysis: 不一致分析
            - mcnemar_test: McNemar检验结果
            - corrected_mcnemar_test: 校正McNemar检验结果
            - kappa_coefficient: Kappa一致性系数
            - odds_ratio: 优势比分析
            - interpretation: 统计解释
    """
    # 从参数对象解构
    a = param.a
    b = param.b
    c = param.c
    d = param.d

    # 执行配对资料卡方检验
    results = perform_paired_chi_square_test(a, b, c, d)
    
    # 构建结果字典
    result_dict = {
        "table_name": "配对资料卡方检验分析",
        "input_parameters": results["input_parameters"],
        "discordant_analysis": results["discordant_analysis"],
        "mcnemar_test": results["mcnemar_test"],
        "corrected_mcnemar_test": results["corrected_mcnemar_test"],
        "kappa_coefficient": results["kappa_coefficient"],
        "odds_ratio": results["odds_ratio"],
        "interpretation": results["interpretation"],
        "remark": f"配对数据: a={a}, b={b}, c={c}, d={d}"
    }
    
    return result_dict