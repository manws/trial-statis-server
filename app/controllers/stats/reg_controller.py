from fastapi import APIRouter, HTTPException, status
from app.schemas.result_data import ResultData

from app.schemas.request_data.reg_param import RegParam1, RegParamN, RegParamLog2, RegParamLogN
from app.stats.reg.reg_stats_1 import cal_result_reg_1
from app.stats.reg.reg_stats_n import cal_result_reg_n
from app.stats.reg.reg_stats_log2 import cal_result_reg_log2
from app.stats.reg.reg_stats_logn import cal_result_reg_logn
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/reg", tags=["Stats-Reg"])


# shared wrapper used for response formatting


@router.post("/1", response_model=ResultData)
async def reg_1(param: RegParam1):
    try:
        report = cal_result_reg_1(param)
        return _wrap(report, data=param.model_dump(), message="回归1号成功")
    except Exception as e:
        LoggerHelper.error(f"reg_1 failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/n", response_model=ResultData)
async def reg_n(param: RegParamN):
    try:
        report = cal_result_reg_n(param)
        return _wrap(report, data=param.model_dump(), message="幂函数回归n成功")
    except Exception as e:
        LoggerHelper.error(f"reg_n failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/log2", response_model=ResultData)
async def reg_log2(param: RegParamLog2):
    try:
        report = cal_result_reg_log2(param)
        return _wrap(report, data=param.model_dump(), message="对数2回归成功")
    except Exception as e:
        LoggerHelper.error(f"reg_log2 failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/logn", response_model=ResultData)
async def reg_logn(param: RegParamLogN):
    try:
        report = cal_result_reg_logn(param)
        return _wrap(report, data=param.model_dump(), message="对数回归logn成功")
    except Exception as e:
        LoggerHelper.error(f"reg_logn failed: {e}")
        return ResultData(code=500, message=str(e), result=None)
