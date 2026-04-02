from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 相关性分析 Pearson直线相关 - 参数
class CorParamPearson(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 相关性分析 Spearman秩相关 - 参数
class CorParamSpearman(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#3 相关性分析 Kendall秩相关 - 参数
class CorParamKendall(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

