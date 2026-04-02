from fastapi import APIRouter, HTTPException, status

from app.schemas.result_data import ResultData
from app.schemas.request_data.runs_param import RunsParamBC, RunsParamValue1, RunsParamValue2, RunsParamValue3
from app.stats.runs.runs_stats_bc import cal_report_runs_bc
from app.stats.runs.runs_stats_value1 import cal_report_runs_value1
from app.stats.runs.runs_stats_value2 import cal_report_runs_value2
from app.stats.runs.runs_stats_value3 import cal_report_runs_value3
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/runs", tags=["Stats-Runs"])

# endpoints mirror the example controller but accept user-provided data

@router.post("/bc", response_model=ResultData)
async def runs_bc(param: RunsParamBC):
    try:
        report = cal_report_runs_bc(param)
        return _wrap(report, data=param.model_dump(), message="游程检验BC报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"runs_bc failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/value1", response_model=ResultData)
async def runs_value1(param: RunsParamValue1):
    try:
        report = cal_report_runs_value1(param)
        return _wrap(report, data=param.model_dump(), message="游程检验value1报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"runs_value1 failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/value2", response_model=ResultData)
async def runs_value2(param: RunsParamValue2):
    try:
        report = cal_report_runs_value2(param)
        return _wrap(report, data=param.model_dump(), message="游程检验value2报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"runs_value2 failed: {e}")
        return ResultData(code=500, message=str(e), result=None)


@router.post("/value3", response_model=ResultData)
async def runs_value3(param: RunsParamValue3):
    try:
        report = cal_report_runs_value3(param)
        return _wrap(report, data=param.model_dump(), message="游程检验value3报告生成成功")
    except Exception as e:
        LoggerHelper.error(f"runs_value3 failed: {e}")
        return ResultData(code=500, message=str(e), result=None)
