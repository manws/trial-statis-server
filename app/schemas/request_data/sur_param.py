from pydantic import BaseModel
from typing import List, Optional
from app.schemas.request_data.stats_data import StatsData, StatsName

# 1.1 生存分析 Kaplan-Meier方法 原始资料 - 参数
class SurParamKM1(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list: List[StatsData]

# 1.2 生存分析 Kaplan-Meier方法 汇总资料 - 参数
class SurParamKM2(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list: List[StatsData]

# 2 生存分析 寿命表法 - 汇总数据 - 参数
class SurParamLT(BaseModel):
    stats_name: StatsName
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list: List[StatsData]

# 3 生存分析生存率比较 - 参数
class SurParamLV(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list: List[StatsData]