from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 方差齐性 F检验（方差比检验） - 参数
class HVParamF(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 方差齐性 Bartlett检验 - 参数
class HVParamBartlett(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#3 方差齐性 Levene检验 - 参数
class HVParamLevene(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#4 方差齐性 Brown-Forsythe检验 - 参数
class HVParamBrownForsythe(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#5 方差齐性 Hartley检验 - 参数
class HVParamHartley(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]