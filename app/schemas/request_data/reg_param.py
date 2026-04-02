from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 回归分析 一元线性回归 - 参数
class RegParam1(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 回归分析 多元线性回归 - 参数
class RegParamN(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#3 回归分析 二元逻辑回归 - 参数
class RegParamLog2(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#4 回归分析 多元逻辑回归 - 参数
class RegParamLogN(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]
