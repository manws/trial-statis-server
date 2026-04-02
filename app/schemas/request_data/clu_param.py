from pydantic import BaseModel
from typing import List, Optional
from app.schemas.request_data.stats_data import StatsData, StatsName

# 1 聚类分析 Q型聚类，相当于按行聚类 - 参数
class CluParamQ(BaseModel):
    # std_type : 标准化方法 
    # 0-不标准化 
    # 1-Z分数
    # 2-范围-1到1标准化
    # 3-范围0到1标准化
    # 4-最大值为1标准化
    # 5-平均值为1标准化
    # 6-标准差为1标准化
    std_type: int

    # clu_method : 聚类方法
    # 0-平均距离
    # 1-最短距离
    # 2-最长距离
    clu_method: int

    # dis_type : 类间距离
    # 0-欧氏距离
    # 1-平方欧氏距离
    # 2-曼哈顿距离
    # 3-切比雪夫距离
    dis_type: int

    # k : 期望聚类数
    k: int

    # clu_name_list : 聚类变量名称列表，长度应与stats_data_list中每个StatsData的data_list长度一致
    stats_name: StatsName

    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list: List[StatsData]

# 2 聚类分析又称R型聚类，相当于按列聚类 - 参数
class CluParamR(BaseModel):
    # std_type : 标准化方法 
    # 0-不标准化 
    # 1-Z分数
    # 2-范围-1到1标准化
    # 3-范围0到1标准化
    # 4-最大值为1标准化
    # 5-平均值为1标准化
    # 6-标准差为1标准化
    std_type: int

    # clu_method : 聚类方法
    # 0-平均距离
    # 1-最短距离
    # 2-最长距离
    clu_method: int

    # dis_type : 类间距离
    # 0-pearson相关系数距离
    # 1-夹角余弦距离
    dis_type: int

    # k : 期望聚类数
    k: int

    # clu_name_list : 聚类变量名称列表，长度应与stats_data_list中每个StatsData的data_list长度一致
    stats_name: StatsName

    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list: List[StatsData]