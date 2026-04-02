from pydantic import BaseModel
from typing import Dict

class ResultData(BaseModel):
    code : int
    message : str
    result : Dict
