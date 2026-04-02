from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 T检验 - 参数
class TParamSingle1(BaseModel):
    # pop_mean : 总体均数
    # n : 样本量
    # sample_mean : 样本均数
    # sample_std : 样本标准差
    pop_mean : float
    n : int
    sample_mean : float
    sample_std : float


#2 T检验原始资料 - 参数
class TParamSingle2(BaseModel):
    # pop_mean : 总体均数
    # list[StatsData] : 样本数据列表
    pop_mean : float
    stats_data_list : list[StatsData]

#3 配对资料T检验 - 参数
class TParamPaired(BaseModel):
    # list[StatsData] : 配对数据列表，每个StatsData包含一对观测值
    stats_data_list : list[StatsData]

#4 独立资料T检验 - 参数
class TParamIndep(BaseModel):
    # list[StatsData] : 独立样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#5 T检验P值计算 - 参数
class TParamP(BaseModel):
    # t : t值
    # df : 自由度
    t_value : float
    df : int