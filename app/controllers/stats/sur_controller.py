from fastapi import APIRouter, HTTPException, status

from app.schemas.result_data import ResultData
from app.schemas.request_data.sur_param import SurParamKM1, SurParamKM2, SurParamLT, SurParamLV
from app.stats.sur.sur_stats_km1 import cal_result_sur_km1
from app.stats.sur.sur_stats_km2 import cal_result_sur_km2
from app.stats.sur.sur_stats_lt import cal_result_sur_lt
from app.stats.sur.sur_stats_lv import cal_result_sur_lv
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/sur", tags=["Stats-Sur"])

# shared wrapper for consistent responses


@router.post("/km1", response_model=ResultData)
async def sur_km1(param: SurParamKM1):
    try:
        report = cal_result_sur_km1(param)
        return _wrap(report, data=param.model_dump(), message="生存分析KM1报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"sur_km failed: {e}")
        return ResultData(code=500, message=str(e), result=None)
    
@router.post("/km2", response_model=ResultData)
async def sur_km2(param: SurParamKM2):
    try:
        report = cal_result_sur_km2(param)
        return _wrap(report, data=param.model_dump(), message="生存分析KM2报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"sur_km failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/lt", response_model=ResultData)
async def sur_lt(param: SurParamLT):
    try:
        report = cal_result_sur_lt(param)
        return _wrap(report, data=param.model_dump(), message="生命表生存分析报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"sur_lt failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/lv", response_model=ResultData)
async def sur_lv(param: SurParamLV):
    try:
        report = cal_result_sur_lv(param)
        return _wrap(report, data=param.model_dump(), message="生存率比较报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"sur_lv failed: {e}")
        return ResultData(code=500, message=str(e), result=None)
