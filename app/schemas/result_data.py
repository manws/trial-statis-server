from pydantic import BaseModel
from typing import Dict, Optional

class ResultData(BaseModel):
    code : int
    message : str
    result : Optional[Dict] = None
