from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 秩和检验 - 配对样本 - 参数
class RSParamPaired(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 秩和检验单样本 - 参数
class RSParamSingle(BaseModel):
    # mean : 单样本Wilcoxon秩和检验的中位数
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    mean : float
    stats_data_list : list[StatsData]

#3 秩和检验两独立样本 - 参数
class RSParamIndep(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#4 秩和检验等级资料 两样本Wilcoxon（ordinal data） - 参数
class RSParamOD1(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#5 秩和检验等级资料 多样本Kruskal-Wallis H（ordinal data） - 参数
class RSParamOD2(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#6 秩和检验Kruskal-Wallis H检验（ordinal data） - 参数
class RSParamKWH(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#7 秩和检验Friedman M检验（ordinal data） - 参数
class RSParamFM(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]