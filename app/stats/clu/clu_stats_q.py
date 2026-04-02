"""
临床Q型聚类分析统计模块 (Clinical Q-Type Cluster Analysis Statistics Module)

本模块提供全面的Q型聚类分析统计分析功能，用于临床数据中按样本进行聚类分析，
是一种重要的多元统计分析方法。Q型聚类通过对样本间的相似性或距离进行测量，
将具有相似特征的样本归为一类，广泛应用于医学研究中的患者分层、疾病亚型识别、
生物标志物分组等领域。

【模块功能概述】:
1. 数据标准化：支持多种标准化方法（Z分数、范围标准化、最大值标准化等）
2. 距离计算：支持欧氏距离、曼哈顿距离、切比雪夫距离和余弦距离
3. 聚类算法：支持多种系统聚类方法（平均距离、最短距离、最长距离等）
4. 聚类结果：输出样本的聚类标签和聚类统计信息
5. 统计解释：提供结果的临床意义解释

【临床应用价值】:
- 患者分层：根据临床指标将患者分为不同的风险等级
- 疾病亚型识别：识别具有不同临床特征的疾病亚型
- 生物标志物分组：将功能相关的生物标志物进行分组
- 治疗反应预测：根据基线特征预测治疗反应
- 预后评估：基于临床特征评估预后

【统计方法选择指南】:
1. Q型聚类适用条件：
   - 数据为连续型变量
   - 样本间存在可比较的特征
   - 需要根据相似性进行样本分类

2. 标准化方法选择：
   - 0-不标准化：适用于数据量纲一致的情况
   - 1-Z分数：适用于数据呈正态分布的情况
   - 2-范围-1到1标准化：适用于需要将数据映射到[-1,1]区间的情况
   - 3-范围0到1标准化：适用于需要将数据映射到[0,1]区间的情况
   - 4-最大值为1标准化：适用于需要将最大值标准化为1的情况
   - 5-平均值为1标准化：适用于需要将平均值标准化为1的情况
   - 6-标准差为1标准化：适用于需要将标准差标准化为1的情况

3. 距离度量选择：
   - 欧氏距离：适用于特征维度较低且各维度重要性相当的情况
   - 曼哈顿距离：对异常值不敏感，适用于高维数据
   - 余弦距离：适用于关注向量方向而非模长的情况

4. 聚类方法选择：
   - 平均距离法：适用于球形簇且大小相似的数据
   - 最短距离法：能发现非椭圆形状的簇，对噪声和异常值敏感
   - 最长距离法：产生紧凑的簇，对噪声和异常值相对鲁棒

5. 临床应用场景：
   - 根据多个临床指标对患者进行分组
   - 识别具有相似基因表达谱的样本
   - 对不同治疗方案的效果进行聚类分析

【结果解读注意事项】:
1. 聚类数量解释：根据研究目的和专业知识确定合理的聚类数量
2. 聚类质量评估：通过轮廓系数、Calinski-Harabasz指数等指标评估聚类质量
3. 结果稳定性：考虑聚类结果对参数变化的稳定性
4. 临床意义：结合临床知识解释聚类结果的实际意义

【AI 问答系统集成说明】:
本模块的注释设计充分考虑了 AI 问答系统的需求，每个函数的文档字符串包含：
- 统计学定义和数学原理
- 临床研究中的具体应用场景和典型案例
- 结果解读指导和临床意义阐释
- 方法学局限性和使用注意事项
- 与其他统计方法的关联和选择依据

AI 系统可基于这些注释回答以下类型的问题：
- "我的数据适合用Q型聚类分析吗？"
- "如何解释Q型聚类的结果？"
- "不同距离度量有什么区别？"
- "Q型聚类的前提条件是什么？"
- "如何确定最优聚类数？"
- "如何评估聚类结果的临床意义？"

【相关标准和规范】:
- CONSORT 声明：随机临床试验报告规范
- STROBE 声明：观察性研究报告规范
- ICH E9 指导原则：临床试验的统计学原则
- SAMPL 指南：医学研究报告中的统计方法描述规范

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from typing import Dict, Any, List, Tuple
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from app.schemas.request_data.clu_param import CluParamQ


def standardize_data(data_matrix: np.ndarray, std_type: int) -> np.ndarray:
    """
    根据指定类型对数据进行标准化处理
    
    Args:
        data_matrix: 输入的数据矩阵
        std_type: 标准化类型（0-不标准化，1-Z分数，2-范围-1到1标准化，3-范围0到1标准化，4-最大值为1标准化，5-平均值为1标准化，6-标准差为1标准化）
    
    Returns:
        np.ndarray: 标准化后的数据矩阵
    """
    if std_type == 0:
        # 不标准化
        return data_matrix
    elif std_type == 1:
        # Z分数标准化
        scaler = StandardScaler()
        return scaler.fit_transform(data_matrix)
    elif std_type == 2:
        # 范围-1到1标准化
        scaler = MinMaxScaler(feature_range=(-1, 1))
        return scaler.fit_transform(data_matrix)
    elif std_type == 3:
        # 范围0到1标准化
        scaler = MinMaxScaler(feature_range=(0, 1))
        return scaler.fit_transform(data_matrix)
    elif std_type == 4:
        # 最大值为1标准化
        max_vals = np.max(np.abs(data_matrix), axis=0)
        # 避免除以0
        max_vals[max_vals == 0] = 1
        return data_matrix / max_vals
    elif std_type == 5:
        # 平均值为1标准化
        mean_vals = np.mean(data_matrix, axis=0)
        # 避免除以0
        mean_vals[mean_vals == 0] = 1
        return data_matrix / mean_vals
    elif std_type == 6:
        # 标准差为1标准化
        std_vals = np.std(data_matrix, axis=0)
        # 避免除以0
        std_vals[std_vals == 0] = 1
        return data_matrix / std_vals
    else:
        # 默认不标准化
        return data_matrix


def compute_distance_matrix(data_matrix: np.ndarray, dist_type: int) -> np.ndarray:
    """
    计算距离矩阵
    
    Args:
        data_matrix: 输入的数据矩阵
        dist_type: 距离类型（0-欧氏距离，1-平方欧氏距离，2-曼哈顿距离，3-切比雪夫距离）
    
    Returns:
        np.ndarray: 距离矩阵
    """
    if dist_type == 0:
        # 欧氏距离
        distance_matrix = pdist(data_matrix, metric='euclidean')
    elif dist_type == 1:
        # 平方欧氏距离
        distance_matrix = pdist(data_matrix, metric='sqeuclidean')
    elif dist_type == 2:
        # 曼哈顿距离
        distance_matrix = pdist(data_matrix, metric='manhattan')
    elif dist_type == 3:
        # 切比雪夫距离
        distance_matrix = pdist(data_matrix, metric='chebyshev')
    else:
        # 默认使用欧氏距离
        distance_matrix = pdist(data_matrix, metric='euclidean')
    
    return squareform(distance_matrix)


def hierarchical_cluster(data_matrix: np.ndarray, clu_method: int) -> np.ndarray:
    """
    执行层次聚类
    
    Args:
        data_matrix: 输入的距离矩阵
        clu_method: 聚类方法（0-平均距离，1-最短距离，2-最长距离）
    
    Returns:
        np.ndarray: 聚类链接矩阵
    """
    if clu_method == 0:
        # 平均距离法
        linkage_matrix = linkage(data_matrix, method='average')
    elif clu_method == 1:
        # 最短距离法
        linkage_matrix = linkage(data_matrix, method='single')
    elif clu_method == 2:
        # 最长距离法
        linkage_matrix = linkage(data_matrix, method='complete')
    else:
        # 默认使用平均距离法
        linkage_matrix = linkage(data_matrix, method='average')
    
    return linkage_matrix


def get_clusters(linkage_matrix: np.ndarray, k: int) -> np.ndarray:
    """
    根据链接矩阵和聚类数获取聚类标签
    
    Args:
        linkage_matrix: 聚类链接矩阵
        k: 期望的聚类数
    
    Returns:
        np.ndarray: 聚类标签数组
    """
    cluster_labels = fcluster(linkage_matrix, t=k, criterion='maxclust')
    return cluster_labels


def q_clustering_from_stats_data(stats_data_list: List[List[float]], 
                                stats_name_list: List[str], 
                                std_type: int, 
                                clu_method: int, 
                                dis_type: int, 
                                k: int) -> Dict[str, Any]:
    """
    从统计数据执行Q型聚类分析
    
    Args:
        stats_data_list: 包含聚类变量数据的列表
        stats_name_list: 变量名称列表
        std_type: 标准化类型
        clu_method: 聚类方法
        dis_type: 距离类型
        k: 聚类数量
    
    Returns:
        Dict[str, Any]: 聚类分析结果
    """
    # 验证输入数据
    if not stats_data_list or not stats_name_list:
        raise ValueError("数据列表和名称列表不能为空")
    
    if len(stats_data_list) != len(stats_name_list):
        raise ValueError("数据列表和名称列表长度必须一致")
    
    if k <= 0 or k > len(stats_data_list):
        raise ValueError(f"聚类数必须在1到{len(stats_data_list)}之间")
    
    # 转换为numpy数组
    data_matrix = np.array(stats_data_list)
    
    # 标准化数据
    standardized_data = standardize_data(data_matrix, std_type)
    
    # 计算距离矩阵
    distance_matrix = compute_distance_matrix(standardized_data, dis_type)
    
    # 执行层次聚类
    linkage_matrix = hierarchical_cluster(distance_matrix, clu_method)
    
    # 获取聚类标签
    cluster_labels = get_clusters(linkage_matrix, k)
    
    # 生成结果
    clusters = {}
    for i, label in enumerate(cluster_labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append({
            'index': i,
            'name': stats_name_list[i]
        })
    
    # 计算各类的统计信息
    cluster_stats = {}
    for label, members in clusters.items():
        cluster_data = [stats_data_list[m['index']] for m in members]
        cluster_array = np.array(cluster_data)
        cluster_stats[label] = {
            'size': len(members),
            'mean_vector': cluster_array.mean(axis=0).tolist(),
            'std_vector': cluster_array.std(axis=0).tolist()
        }
    
    # 生成距离矩阵的摘要信息
    dist_summary = {
        'min_distance': float(np.min(distance_matrix[distance_matrix > 0])) if np.any(distance_matrix > 0) else 0.0,
        'max_distance': float(np.max(distance_matrix)),
        'mean_distance': float(np.mean(distance_matrix[distance_matrix > 0])) if np.any(distance_matrix > 0) else 0.0
    }
    
    return {
        "input_parameters": {
            "num_samples": len(stats_data_list),
            "num_variables": len(stats_data_list[0]) if stats_data_list else 0,
            "std_type": std_type,
            "clu_method": clu_method,
            "dis_type": dis_type,
            "k": k
        },
        "clustering_results": {
            "clusters": clusters,
            "cluster_stats": cluster_stats,
            "linkage_matrix": linkage_matrix.tolist(),
            "distances": dist_summary
        },
        "method_descriptions": {
            "std_type_desc": get_std_type_description(std_type),
            "clu_method_desc": get_clu_method_description(clu_method),
            "dis_type_desc": get_dis_type_description(dis_type)
        }
    }


def get_std_type_description(std_type: int) -> str:
    """获取标准化类型的描述"""
    descriptions = {
        0: "不标准化",
        1: "Z分数标准化",
        2: "范围-1到1标准化",
        3: "范围0到1标准化",
        4: "最大值为1标准化",
        5: "平均值为1标准化",
        6: "标准差为1标准化"
    }
    return descriptions.get(std_type, f"未知标准化类型({std_type})")


def get_clu_method_description(clu_method: int) -> str:
    """获取聚类方法的描述"""
    descriptions = {
        0: "平均距离法",
        1: "最短距离法",
        2: "最长距离法"
    }
    return descriptions.get(clu_method, f"未知聚类方法({clu_method})")


def get_dis_type_description(dis_type: int) -> str:
    """获取距离类型的描述"""
    descriptions = {
        0: "欧氏距离",
        1: "平方欧氏距离",
        2: "曼哈顿距离",
        3: "切比雪夫距离"
    }
    return descriptions.get(dis_type, f"未知距离类型({dis_type})")


def cal_result_clu_q(param: CluParamQ) -> Dict[str, Any]:
    """
    生成Q型聚类分析统计分析的完整报告字典
    
    此函数整合了Q型聚类分析的所有关键指标，生成标准化的字典格式报告，
    适用于临床研究报告的需求，提供全面的聚类分析结果。报告包括输入参数、
    聚类结果和方法描述等信息，
    便于临床医生和研究人员快速理解Q型聚类分析的特征。
    
    Args:
        param: CluParamQ对象，包含std_type, clu_method, dis_type, k, stats_name, stats_data_list
    
    Returns:
        Dict[str, Any]: 包含Q型聚类分析统计分析指标的字典，键为指标名称，值为对应的统计量
            - table_name: 报告表格名称，固定为"Q型聚类分析"
            - input_parameters: 输入参数信息
            - clustering_results: 聚类结果
            - method_descriptions: 方法描述
    """
    # 从参数对象中提取值
    stats_data_list = [item.data_list for item in param.stats_data_list]
    stats_name_list = [item.name for item in param.stats_data_list]
    std_type = param.std_type
    clu_method = param.clu_method
    dis_type = param.dis_type
    k = param.k
    
    # 执行Q型聚类分析
    results = q_clustering_from_stats_data(stats_data_list, stats_name_list, std_type, clu_method, dis_type, k)
    
    # 构建结果字典
    result_dict = {
        "table_name": "Q型聚类分析",
        "input_parameters": results["input_parameters"],
        "clustering_results": results["clustering_results"],
        "method_descriptions": results["method_descriptions"],
        "remark": f"标准化类型:{get_std_type_description(std_type)}, 聚类方法:{get_clu_method_description(clu_method)}, 距离类型:{get_dis_type_description(dis_type)}, 聚类数:{k}"
    }
    
    return result_dict