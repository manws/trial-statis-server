from fastapi import APIRouter, HTTPException, status
from app.schemas.result_data import ResultData

from app.schemas.request_data.hv_param import (
    HVParamBartlett,
    HVParamBrownForsythe,
    HVParamF,
    HVParamHartley,
    HVParamLevene,
)
from app.stats.hv.hv_stats_bartlett import cal_report_hv_bartlett
from app.stats.hv.hv_stats_bf import cal_report_hv_bf
from app.stats.hv.hv_stats_f import perform_hv_f_test
from app.stats.hv.hv_stats_hartley import cal_report_hv_hartley
from app.stats.hv.hv_stats_levene import cal_report_hv_levene
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/hv", tags=["Stats-HV"])


# use shared _wrap helper from stats.utils


@router.post("/bartlett", response_model=ResultData)
async def hv_bartlett(param: HVParamBartlett):
    try:
        report = cal_report_hv_bartlett(param)
        return _wrap(report, data=param.model_dump(), message="Bartlett检验成功")
    except Exception as e:
        LoggerHelper.error(f"hv_bartlett failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/bf", response_model=ResultData)
async def hv_bf(param: HVParamBrownForsythe):
    try:
        report = cal_report_hv_bf(param)
        return _wrap(report, data=param.model_dump(), message="Brown-Forsythe检验成功")
    except Exception as e:
        LoggerHelper.error(f"hv_bf failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/f", response_model=ResultData)
async def hv_f(param: HVParamF):
    try:
        report = perform_hv_f_test(param)
        return _wrap(report, data=param.model_dump(), message="F检验成功")
    except Exception as e:
        LoggerHelper.error(f"hv_f failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/hartley", response_model=ResultData)
async def hv_hartley(param: HVParamHartley):
    try:
        report = cal_report_hv_hartley(param)
        return _wrap(report, data=param.model_dump(), message="Hartley检验成功")
    except Exception as e:
        LoggerHelper.error(f"hv_hartley failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/levene", response_model=ResultData)
async def hv_levene(param: HVParamLevene):
    try:
        report = cal_report_hv_levene(param)
        return _wrap(report, data=param.model_dump(), message="Levene检验成功")
    except Exception as e:
        LoggerHelper.error(f"hv_levene failed: {e}")
        return ResultData(code=500, message=str(e), result=None)
