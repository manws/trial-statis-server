from pydantic import BaseModel

from app.schemas.request_data.stats_data import StatsData

#1 Cochran'Q 检验 - 参数
class CQParam(BaseModel):
    # list[StatsData] : 样本数据列表，每个StatsData代表一个组的数据
    stats_data_list : list[StatsData]