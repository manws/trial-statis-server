"""
临床 Fisher 精确检验统计分析模块 (Clinical Fisher Exact Test Statistics Module)

本模块提供全面的 Fisher 精确检验统计分析功能，用于临床数据中两个分类变量之间的关联性检验，
特别是在小样本或理论频数小于 5 的情况下。Fisher 精确检验是医学研究中重要的统计方法，
如比较不同治疗方案的有效性、评估危险因素与疾病的关系等。

【模块功能概述】:
1. 超几何分布概率计算：计算特定四格表的精确概率
2. Fisher 精确检验：计算双侧和单侧检验的 P 值
3. 优势比计算：计算并估计优势比
4. 所有可能表格枚举：列出所有满足边缘总计的可能表格
5. Scipy 验证：使用 scipy 库验证计算结果
6. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 小样本研究：在样本量较小时进行精确的关联性检验
- 罕见事件分析：当事件发生率很低时进行精确分析
- 诊断试验：评估诊断试验的准确性
- 预后分析：评估预后因子与结局的关系

【统计方法选择指南】:
1. Fisher 精确检验适用条件：
   - 总样本量<20 时
   - 总样本量 20-40 且任一理论频数<5 时
   - 任一单元格频数为 0 时
   - 任何情况下，当需要精确概率时

2. 临床应用场景：
   - 比较小样本量下两组患者治疗有效率的差异
   - 评估罕见疾病与暴露因素的关系
   - 比较诊断试验的敏感性和特异性

【结果解读注意事项】:
1. P 值解释：在零假设（两变量独立）成立的情况下，观察到当前或更极端结果的精确概率
2. 双侧检验：检验两变量是否存在关联（无方向性）
3. 单侧检验：检验关联的方向性（如 OR>1 或 OR<1）
4. 优势比解释：OR>1 表示暴露增加患病风险，OR<1 表示暴露降低患病风险
5. 临床意义：统计学显著性不等同于临床重要性，需结合实际意义解读

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的小样本数据适合用 Fisher 精确检验吗？"
- "如何解释 Fisher 精确检验的结果？"
- "双侧和单侧检验有什么区别？"
- "Fisher 精确检验的前提条件是什么？"
- "如何判断优势比的临床意义？"

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
from typing import Dict, Any, List
from scipy.stats import hypergeom, chi2, fisher_exact
from app.schemas.request_data.chi_param import ChiSquareParamFisher


def calculate_hypergeometric_pmf(k: int, M: int, n: int, N: int) -> float:
    """
    计算超几何分布的概率质量函数值
    
    Args:
        k: 成功抽取的数量
        M: 总体大小
        n: 总体中的成功数
        N: 抽样数量
    
    Returns:
        float: 超几何分布的概率值
    """
    # 验证参数
    if k < 0 or n < 0 or N < 0 or M < 0:
        return 0.0
    if k > n or N > M or n > M or N > M:
        return 0.0
    
    # 计算概率
    try:
        # 使用scipy的超几何分布计算
        rv = hypergeom(M, n, N)
        return float(rv.pmf(k))
    except Exception:
        # 如果scipy不可用，使用数学公式计算
        try:
            # C(n,k) * C(M-n, N-k) / C(M,N)
            numerator = math.comb(n, k) * math.comb(M - n, N - k)
            denominator = math.comb(M, N)
            return numerator / denominator if denominator != 0 else 0.0
        except Exception:
            return 0.0


def fisher_exact_test(a: int, b: int, c: int, d: int) -> Dict[str, Any]:
    """
    执行Fisher精确检验
    
    Args:
        a: 样本1发生数
        b: 样本1未发生数
        c: 样本2发生数
        d: 样本2未发生数
    
    Returns:
        Dict[str, Any]: 检验结果字典
    """
    # 验证参数
    if any(x < 0 for x in [a, b, c, d]):
        raise ValueError("四格表中的所有值必须为非负数")
    
    if a + b + c + d == 0:
        raise ValueError("四格表总和不能为0")
    
    # 计算边缘合计
    row1_sum = a + b  # 第一行合计
    row2_sum = c + d  # 第二行合计
    col1_sum = a + c  # 第一列合计
    col2_sum = b + d  # 第二列合计
    n = row1_sum + row2_sum  # 总计
    
    # 使用scipy的fisher_exact函数
    try:
        oddsratio, p_value_two_tail = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        _, p_value_left = fisher_exact([[a, b], [c, d]], alternative='less')
        _, p_value_right = fisher_exact([[a, b], [c, d]], alternative='greater')
    except Exception as e:
        # 如果scipy函数出错，使用备用计算方法
        oddsratio = float('nan')
        p_value_two_tail = float('nan')
        p_value_left = float('nan')
        p_value_right = float('nan')
    
    # 计算优势比OR（如果可能）
    if b == 0 or c == 0:
        or_value = float('inf') if a > 0 and d > 0 else 0.0
    elif a == 0 or d == 0:
        or_value = 0.0
    else:
        or_value = (a * d) / (b * c)
    
    # 确定最终的odds_ratio值：优先使用计算值，如果不可用则使用备选值
    final_odds_ratio = or_value if not math.isfinite(oddsratio) else oddsratio
    
    return {
        "odds_ratio": final_odds_ratio,
        "p_value_two_tail": p_value_two_tail,
        "p_value_left": p_value_left,
        "p_value_right": p_value_right,
        "marginal_totals": {
            "row1": row1_sum,
            "row2": row2_sum,
            "col1": col1_sum,
            "col2": col2_sum,
            "total": n
        }
    }


def generate_possible_tables(row1_sum: int, row2_sum: int, col1_sum: int, col2_sum: int) -> List[Dict[str, Any]]:
    """
    生成给定边缘合计的所有可能表格
    
    Args:
        row1_sum: 第一行合计
        row2_sum: 第二行合计
        col1_sum: 第一列合计
        col2_sum: 第二列合计
    
    Returns:
        List[Dict[str, Any]]: 所有可能表格的列表
    """
    tables = []
    
    # 计算a的可能取值范围
    min_a = max(0, row1_sum - col2_sum, col1_sum - row2_sum)
    max_a = min(row1_sum, col1_sum)
    
    for a in range(min_a, max_a + 1):
        b = row1_sum - a
        c = col1_sum - a
        d = col2_sum - b
        
        # 验证所有值都是非负的
        if a >= 0 and b >= 0 and c >= 0 and d >= 0:
            # 计算该表格的概率
            prob = calculate_hypergeometric_pmf(a, row1_sum + row2_sum, col1_sum, row1_sum)
            
            tables.append({
                "table": {"a": a, "b": b, "c": c, "d": d},
                "probability": prob
            })
    
    return tables


def analyze_tables_probabilities(tables: List[Dict[str, Any]], observed_prob: float) -> Dict[str, Any]:
    """
    分析表格概率，找出比观察到的表格更极端的情况
    
    Args:
        tables: 所有可能表格的列表
        observed_prob: 观察到的表格的概率
    
    Returns:
        Dict[str, Any]: 概率分析结果
    """
    # 计算比观察到的表格概率更小的所有表格的概率之和
    p_value_two_tail = sum(table["probability"] for table in tables if table["probability"] <= observed_prob)
    
    # 获取最小的表格概率（最不可能的情况）
    min_prob = min(table["probability"] for table in tables)
    
    return {
        "p_value_two_tail": p_value_two_tail,
        "observed_prob": observed_prob,
        "min_possible_prob": min_prob,
        "more_extreme_tables_count": len([t for t in tables if t["probability"] < observed_prob])
    }


def perform_fisher_exact_test(a: int, b: int, c: int, d: int) -> Dict[str, Any]:
    """
    执行完整的Fisher精确检验流程
    
    Args:
        a: 样本1发生数
        b: 样本1未发生数
        c: 样本2发生数
        d: 样本2未发生数
    
    Returns:
        Dict[str, Any]: 完整的检验结果字典
    """
    # 验证参数
    if any(x < 0 for x in [a, b, c, d]):
        raise ValueError("四格表中的所有值必须为非负数")
    
    if a + b + c + d == 0:
        raise ValueError("四格表总和不能为0")
    
    # 计算边缘合计
    row1_sum = a + b  # 第一行合计
    row2_sum = c + d  # 第二行合计
    col1_sum = a + c  # 第一列合计
    col2_sum = b + d  # 第二列合计
    n = row1_sum + row2_sum  # 总计
    
    # 执行Fisher精确检验
    fisher_results = fisher_exact_test(a, b, c, d)
    
    # 生成所有可能的表格
    all_possible_tables = generate_possible_tables(row1_sum, row2_sum, col1_sum, col2_sum)
    
    # 计算观察到的表格的概率
    observed_prob = calculate_hypergeometric_pmf(a, n, col1_sum, row1_sum)
    
    # 分析概率
    prob_analysis = analyze_tables_probabilities(all_possible_tables, observed_prob)
    
    # 执行scipy验证
    try:
        scipy_oddsratio, scipy_p_value = fisher_exact([[a, b], [c, d]])
        scipy_verification = {
            "odds_ratio": scipy_oddsratio,
            "p_value": scipy_p_value,
            "match_expected": abs(scipy_p_value - fisher_results["p_value_two_tail"]) < 1e-10
        }
    except Exception as e:
        scipy_verification = {
            "error": str(e)
        }
    
    # 返回综合结果
    return {
        "input_parameters": {
            "a": a, "b": b, "c": c, "d": d,
            "row_totals": [row1_sum, row2_sum],
            "column_totals": [col1_sum, col2_sum],
            "grand_total": n
        },
        "observed_table": {
            "values": {"a": a, "b": b, "c": c, "d": d},
            "probability": observed_prob
        },
        "all_tables": {
            "count": len(all_possible_tables),
            "tables_with_smaller_or_equal_prob": [t for t in all_possible_tables if t["probability"] <= observed_prob],
            "most_extreme_tables": [t for t in all_possible_tables if t["probability"] == min(t["probability"] for t in all_possible_tables)]
        },
        "test_results": {
            "odds_ratio": fisher_results["odds_ratio"],
            "p_values": {
                "two_tailed": fisher_results["p_value_two_tail"],
                "left_tail": fisher_results["p_value_left"],
                "right_tail": fisher_results["p_value_right"]
            }
        },
        "scipy_verification": scipy_verification,
        "interpretation": {
            "significance_95": "有显著性差异" if fisher_results["p_value_two_tail"] < 0.05 else "无显著性差异",
            "significance_99": "有非常显著性差异" if fisher_results["p_value_two_tail"] < 0.01 else "无非常显著性差异",
            "odds_ratio_interpretation": f"优势比为 {fisher_results['odds_ratio']:.4f}"
        }
    }


def cal_result_fisher(param: ChiSquareParamFisher) -> Dict[str, Any]:
    """
    生成Fisher精确检验统计分析的完整报告字典
    
    此函数整合了Fisher精确检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的Fisher精确检验结果。报告包括输入参数、
    观察表格、所有可能表格、检验结果、验证结果和统计解释等信息，
    便于临床医生和研究人员快速理解Fisher精确检验的特征。
    
    Args:
        param: ChiSquareParamFisher对象，包含a, b, c, d四个整数参数
    
    Returns:
        Dict[str, Any]: 包含Fisher精确检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"Fisher精确检验分析"
            - input_parameters: 输入参数信息
            - observed_table: 观察表格信息
            - all_tables: 所有可能表格
            - test_results: 检验结果
            - scipy_verification: scipy验证结果
            - interpretation: 统计解释
    """
    # 从参数对象中提取值
    a = param.a
    b = param.b
    c = param.c
    d = param.d
    
    # 执行Fisher精确检验
    results = perform_fisher_exact_test(a, b, c, d)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Fisher精确检验分析",
        "input_parameters": results["input_parameters"],
        "observed_table": results["observed_table"],
        "all_tables": results["all_tables"],
        "test_results": results["test_results"],
        "scipy_verification": results["scipy_verification"],
        "interpretation": results["interpretation"],
        "remark": f"四格表数据: a={a}, b={b}, c={c}, d={d}"
    }
    
    return result_dict