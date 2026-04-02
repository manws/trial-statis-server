"""
临床Cochran's Q检验统计分析模块 (Clinical Cochran's Q Test Statistics Module)

本模块提供全面的Cochran's Q检验统计分析功能，用于临床数据中多个相关样本的二分类数据检验，
是一种重要的非参数统计分析方法。Cochran's Q检验用于判断多个相关样本（如不同治疗方法）的
二分类结果（成功/失败）是否存在显著差异，广泛应用于医学研究中的治疗效果比较、诊断方法
评价、干预措施评估等领域。

【模块功能概述】:
1. Cochran's Q统计量计算：计算Cochran's Q检验统计量
2. 显著性检验：执行卡方检验判断Q统计量的显著性
3. P值计算：计算检验的P值
4. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 治疗效果比较：比较多种治疗方法的效果差异
- 诊断方法评价：评估不同诊断方法的准确性差异
- 干预措施评估：分析不同干预措施的有效性差异
- 研究假设验证：验证多个相关样本是否存在显著差异

【统计方法选择指南】:
1. Cochran's Q检验适用条件：
   - 数据为二分类（通常编码为0和1）
   - 样本为相关样本（如同一组受试者接受多种处理）
   - 各处理组样本量相同
   - 观测值独立

2. 临床应用场景：
   - 比较三种或以上不同药物治疗同一疾病的疗效
   - 评估多个诊断测试的准确性差异
   - 分析不同手术方式的治疗效果差异
   - 研究多种康复方案的疗效对比

【结果解读注意事项】:
1. Q统计量解释：Q值越大，表明各处理组间差异越显著
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
- "我的数据适合用Cochran's Q检验吗？"
- "如何解释Cochran's Q检验的结果？"
- "Cochran's Q检验的前提条件是什么？"
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
from scipy.stats import chi2
import numpy as np


def cochran_q_test(data_matrix: List[List[int]]) -> Dict:
    """
    Cochran's Q检验（适用于二分类数据的多个相关样本检验）
    
    参数:
    - data_matrix: 数据矩阵，每一行代表一个受试者，每一列代表一种处理/条件
                   数据应为二分类（0或1）
    
    返回:
    - 包含检验统计量和结果的字典
    """
    # 参数验证
    if not data_matrix or len(data_matrix) == 0:
        raise ValueError("数据矩阵不能为空")
    
    # 验证数据矩阵的一致性
    n_subjects = len(data_matrix)
    n_treatments = len(data_matrix[0]) if data_matrix else 0
    
    if n_treatments < 2:
        raise ValueError("至少需要两种处理进行Cochran's Q检验")
    
    # 验证每个受试者的数据完整性
    for i, subject_data in enumerate(data_matrix):
        if len(subject_data) != n_treatments:
            raise ValueError(f"第{i+1}个受试者的数据长度不一致")
        # 验证数据是否为二分类
        for j, value in enumerate(subject_data):
            if value not in [0, 1]:
                raise ValueError(f"数据必须为二分类（0或1），第{i+1}个受试者的第{j+1}个处理值为{value}")
    
    # 计算统计量
    # Lj = 第j个处理的成功次数（1的个数）
    # Ci = 第i个受试者的总成功次数
    L = [0] * n_treatments  # 每个处理的成功次数
    C = [0] * n_subjects    # 每个受试者的总成功次数
    
    for i in range(n_subjects):
        for j in range(n_treatments):
            if data_matrix[i][j] == 1:
                L[j] += 1
                C[i] += 1
    
    # 计算Cochran's Q统计量
    # Q = (n-1) * Σ(Lj - L_bar)² / (ΣCi - ΣCi²/n)
    # 其中 L_bar = ΣLj/k，k为处理数
    
    L_sum = sum(L)
    L_bar = L_sum / n_treatments
    
    numerator = sum((Lj - L_bar) ** 2 for Lj in L)
    
    C_sum = sum(C)
    C_squared_sum = sum(ci ** 2 for ci in C)
    
    denominator = C_sum - C_squared_sum / n_subjects
    
    # 避免除零错误
    if denominator == 0:
        raise ValueError("分母为零，无法计算Q统计量")
    
    Q_statistic = (n_subjects - 1) * numerator / denominator
    
    # 自由度
    df = n_treatments - 1
    
    # 计算P值（近似卡方检验）
    if Q_statistic >= 0:
        p_value = 1 - chi2.cdf(Q_statistic, df)
    else:
        p_value = 1.0  # Q统计量不应为负，如果出现则设为最大P值
    
    # 显著性判断 (> 0.05 和 > 0.01)
    is_greater_than_05 = p_value > 0.05
    is_greater_than_01 = p_value > 0.01
    
    # 计算基本统计量
    treatment_stats = []
    for j in range(n_treatments):
        success_rate = L[j] / n_subjects if n_subjects > 0 else 0
        treatment_stats.append({
            "treatment_index": j + 1,
            "success_count": L[j],
            "success_rate": float(success_rate),
            "failure_count": n_subjects - L[j]
        })
    
    subject_stats = []
    for i in range(n_subjects):
        subject_stats.append({
            "subject_index": i + 1,
            "total_successes": C[i],
            "success_rate": float(C[i] / n_treatments) if n_treatments > 0 else 0
        })
    
    return {
        "input_parameters": {
            "n_subjects": n_subjects,
            "n_treatments": n_treatments,
            "total_observations": n_subjects * n_treatments
        },
        "treatment_statistics": treatment_stats,
        "subject_statistics": subject_stats,
        "test_statistics": {
            "q_value": float(Q_statistic),
            "degrees_of_freedom": df,
            "p_value_two_sided": float(p_value)
        },
        "significance_tests": {
            "p_greater_than_0_05": is_greater_than_05,
            "p_greater_than_0_01": is_greater_than_01,
            "not_significant_at_05": "不显著" if is_greater_than_05 else "显著",
            "not_significant_at_01": "不显著" if is_greater_than_01 else "显著"
        },
        "interpretation": {
            "cochran_q_interpretation": f"Cochran's Q检验在0.05水平下{'不' if is_greater_than_05 else ''}显著 (p={'>' if is_greater_than_05 else '≤'}0.05)",
            "cochran_q_interpretation_01": f"Cochran's Q检验在0.01水平下{'不' if is_greater_than_01 else ''}显著 (p={'>' if is_greater_than_01 else '≤'}0.01)"
        }
    }


def cal_result_cq(data_matrix: List[List[int]]) -> Dict:
    """
    生成Cochran's Q检验统计分析的完整报告字典
    
    此函数整合了Cochran's Q检验的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的检验结果。报告包括输入参数、
    检验统计量和统计解释等信息，
    便于临床医生和研究人员快速理解Cochran's Q检验的特征。
    
    Args:
        data_matrix: List[List[int]]，数据矩阵，每一行代表一个受试者，每一列代表一种处理/条件
                   数据应为二分类（0或1）
    
    Returns:
        Dict: 包含Cochran's Q检验统计分析指标的字典，键为指标名称，值为对应的统计量
            - input_parameters: 输入参数信息
            - treatment_statistics: 处理统计信息
            - subject_statistics: 受试者统计信息
            - test_statistics: 检验统计量
            - significance_tests: 显著性检验结果
            - interpretation: 统计解释
    """
    # 从参数对象解构
    data_matrix = [[int(x) for x in item.data_list] for item in param.stats_data_list]

    # 执行Cochran's Q检验
    results = cochran_q_test(data_matrix)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Cochran's Q检验分析",
        "input_parameters": results["input_parameters"],
        "treatment_statistics": results["treatment_statistics"],
        "subject_statistics": results["subject_statistics"],
        "test_statistics": results["test_statistics"],
        "significance_tests": results["significance_tests"],
        "interpretation": results["interpretation"],
        "remark": f"处理数: {results['input_parameters']['n_treatments']}, 受试者数: {results['input_parameters']['n_subjects']}"
    }
    
    return result_dict
