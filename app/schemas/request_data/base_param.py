from pydantic import BaseModel
from app.schemas.request_data.stats_data import StatsData, StatsName

#1 描述性统计 - 参数
class BaseParamDes(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 频数统计 - 参数
class BaseParamFreq(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]

#2 正态分布 - 参数
class BaseParamNormal(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]


#1.1 二项分布概率计算 - 参数
class BaseParamBino1(BaseModel):
    # total_posi_p:总体阳性概率
    # sample_size: 样本量
    total_posi_p : float
    sample_size : int

#1.2 总体率的区间估计 - 参数
class BaseParamBino2(BaseModel):
    # confidence_level:可信度
    # sample_size: 样本量
    # sample_posi_num: 样本阳性数
    confidence_level : float
    sample_size : int
    sample_posi_num : int

#1.3 样本率与总体率比较 - 参数
class BaseParamBino3(BaseModel):
    # total_posi_p:总体阳性概率
    # sample_size: 样本量
    # sample_posi_num: 样本阳性数
    total_posi_p : float
    sample_size : int
    sample_posi_num : int

#1.4 两样本率比较 - 参数
class BaseParamBino4(BaseModel):
    # sample1_num:样本1数
    # sample1_posi_num: 样本1阳性数
    # sample2_num: 样本2数
    # sample2_posi_num: 样本2阳性数
    sample1_num : int
    sample1_posi_num : int
    sample2_num : int
    sample2_posi_num : int

#2.1 泊松分布概率计算 - 参数
class BaseParamPois1(BaseModel):
    # total_avg:总体均数
    # total_p:总体率
    # sample_num: 样本数
    total_avg : float
    total_p : float
    sample_num : int

#2.2 总体均数的区间估计 - 参数
class BaseParamPois2(BaseModel):
    # confidence_level:可信度
    # total_avg: 总体均数
    # sample_posi_num: 样本阳性数
    confidence_level : float
    total_avg : float

#2.3 样本均数与总体均数比较 - 参数
class BaseParamPois3(BaseModel):
    # total_posi_p:总体阳性概率
    # sample_size: 样本量
    # sample_posi_num: 样本阳性数
    total_posi_p : float
    sample_size : int
    sample_posi_num : int

#2.4 两样本均数比较 - 参数
class BaseParamPois4(BaseModel):
    # sample1_num:观察单位数
    # sample1_posi_num: 样本1计数
    # sample2_num: 样本2数
    # sample2_tick: 样本2阳性数
    sample1_num : int
    sample1_tick : int
    sample2_num : int
    sample2_tick : int
