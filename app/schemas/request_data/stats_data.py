from pydantic import BaseModel
from typing import List

class StatsData(BaseModel):
    field_name : str
    data_list : List[float]

class StatsName(BaseModel):
    field_name : str
    name_list : List[str]
