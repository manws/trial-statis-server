from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 Z检验 - 单样本 - 参数
class ZParamSingle1(BaseModel):
    # pop_mean : 总体均数
    # pop_std : 总体标准差
    # n : 样本量
    # sample_mean : 样本均数
    pop_mean : float
    pop_std : float
    n : int
    sample_mean : float

#2 Z检验单样本原始资料 - 参数
class ZParamSingle2(BaseModel):
    # pop_mean : 总体均数
    # pop_std : 总体标准差
    # list[StatsData] : 样本数据列表
    pop_mean : float
    pop_std : float
    stats_data_list : list[StatsData]

#3 独立样本Z检验 - 参数
class ZParamIndep1(BaseModel):
    # n1 : 样本量
    # mean1 : 样本均数
    # std1 : 样本标准差
    # n2 : 样本量
    # mean2 : 样本均数
    # std2 : 样本标准差
    n1 : int
    mean1 : float
    std1 : float
    n2 : int
    mean2 : float
    std2 : float

#4 独立样本Z检验(原始资料) - 参数
class ZParamIndep2(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]