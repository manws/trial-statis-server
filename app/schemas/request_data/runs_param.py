from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 游程检验 - 二分类变量 - 参数
class RunsParamBC(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 游程检验 - 数值变量【平均数上下分类】 - 参数
class RunsParamValue1(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#3 游程检验 - 数值变量【中位数上下分类】 - 参数
class RunsParamValue2(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#4 游程检验 - 数值变量【给定值上下分类】 - 参数
class RunsParamValue3(BaseModel):
    # value : 给定值
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    value : float
    stats_data_list : list[StatsData]