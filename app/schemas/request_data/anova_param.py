from pydantic import BaseModel
from typing import List, Optional
from app.schemas.request_data.stats_data import StatsData

# 完全随机设计的方差分析参数
class AnovaParamCRD(BaseModel):
    stats_data_list: List[StatsData]

# 随机区组设计的方差分析参数
class AnovaParamRBD(BaseModel):
    stats_data_list: List[StatsData]