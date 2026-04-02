from fastapi import APIRouter, HTTPException, status
from app.schemas.result_data import ResultData

from app.schemas.request_data.cq_param import CQParam
from app.stats.cq.cq_stats_cq import cal_report_cq
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/cq", tags=["Stats-CQ"])

# reuse shared _wrap helper


@router.post("/q", response_model=ResultData)
async def cq_q(param: CQParam):
    try:
        report = cal_report_cq(param)
        return _wrap(report, data=param.model_dump(), message="CQ检验成功")
    except Exception as e:
        LoggerHelper.error(f"cq_q failed: {e}")
        return ResultData(code=500, message=str(e), result=None)
