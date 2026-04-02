"""
临床秩和检验统计分析模块控制器 (Clinical Rank Sum Test Statistics Analysis Module Controller)

本控制器提供全面的临床秩和检验统计分析功能接口，涵盖配对样本、单样本、独立样本等多种秩和检验方法，
为临床研究数据分析提供标准化的 API 服务。

【模块功能概述】:
1. 配对样本秩和检验：比较配对数据的差异（Wilcoxon 符号秩检验）
2. 单样本秩和检验：比较单样本与已知总体的差异
3. 独立样本秩和检验：比较两独立样本的差异（Mann-Whitney U 检验）
4. 多组独立样本检验：Kruskal-Wallis H 检验
5. Friedman 检验：多个相关样本的秩和检验

【API 设计原则】:
- 参数校验：使用 Pydantic 模型进行请求参数校验
- 错误处理：统一的异常捕获和错误响应格式
- 日志记录：详细的操作日志记录便于调试和监控
- 结果封装：使用统一的响应格式包装计算结果

【技术架构】:
- 采用 FastAPI 框架提供高性能异步 API 服务
- 与底层统计算法模块解耦，通过函数调用实现功能
- 遵循 RESTful API 设计规范，提供清晰的端点命名

【临床应用价值】:
- 为临床研究提供标准化的非参数检验工具
- 支持多种秩和检验的应用场景
- 提供直观的 API 接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from fastapi import APIRouter, HTTPException, status, Body
from app.schemas.result_data import ResultData

from app.schemas.request_data.rs_param import (
    RSParamPaired,
    RSParamSingle,
    RSParamIndep,
    RSParamOD1,
    RSParamOD2,
    RSParamKWH,
    RSParamFM,
)
from app.stats.rs.rs_stats_paired import cal_result_rs_paired
from app.stats.rs.rs_stats_single import cal_result_rs_single
from app.stats.rs.rs_stats_indep import cal_result_rs_indep
from app.stats.rs.rs_stats_od1 import cal_result_rs_od1
from app.stats.rs.rs_stats_od2 import cal_result_rs_od2
from app.stats.rs.rs_stats_kwh import cal_result_rs_kwh
from app.stats.rs.rs_stats_fm import cal_result_rs_fm
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/rs", tags=["Stats-RS"])


@router.post("/paired", response_model=ResultData)
async def rs_paired(param: RSParamPaired):
    """
    配对样本秩和检验接口 - Wilcoxon 符号秩检验
    
    功能：比较配对样本的差异
    参数：RSParamPaired - 包含两组配对数据的参数对象
    返回：包含配对样本秩和检验结果的字典
    临床应用：用于比较同一受试者治疗前后的指标变化
    """
    try:
        # 调用底层配对样本秩和检验算法函数
        result = cal_result_rs_paired(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="配对样本秩和检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_paired failed: {e}")
        return ResultData(code=500, message=f"配对样本秩和检验失败：{e}", result=None)


@router.post("/single", response_model=ResultData)
async def rs_single(param: RSParamSingle):
    """
    单样本秩和检验接口
    
    功能：比较单样本与已知总体的差异
    参数：RSParamSingle - 包含单样本数据和总体中位数的参数对象
    返回：包含单样本秩和检验结果的字典
    临床应用：用于比较单组数据与已知参考值的差异
    """
    try:
        # 调用底层单样本秩和检验算法函数
        result = cal_result_rs_single(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="单样本秩和检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_single failed: {e}")
        return ResultData(code=500, message=f"单样本秩和检验失败：{e}", result=None)


@router.post("/indep", response_model=ResultData)
async def rs_indep(param: RSParamIndep):
    """
    独立样本秩和检验接口 - Mann-Whitney U 检验
    
    功能：比较两独立样本的差异
    参数：RSParamIndep - 包含两组独立数据的参数对象
    返回：包含独立样本秩和检验结果的字典
    临床应用：用于比较两组独立受试者的指标差异
    """
    try:
        # 调用底层独立样本秩和检验算法函数
        result = cal_result_rs_indep(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="独立样本秩和检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_indep failed: {e}")
        return ResultData(code=500, message=f"独立样本秩和检验失败：{e}", result=None)


@router.post("/od1", response_model=ResultData)
async def rs_od1(param: RSParamOD1):
    """
    有序分组秩和检验接口 - Wilcoxon 秩和检验
    
    功能：比较两个有序分组的差异
    参数：RSParamOD1 - 包含有序分组数据的参数对象
    返回：包含有序分组秩和检验结果的字典
    临床应用：用于比较两个等级资料的差异
    """
    try:
        # 调用底层有序分组秩和检验算法函数
        result = cal_result_rs_od1(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="有序分组秩和检验 OD1 成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_od1 failed: {e}")
        return ResultData(code=500, message=f"有序分组秩和检验 OD1 失败：{e}", result=None)


@router.post("/od2", response_model=ResultData)
async def rs_od2(param: RSParamOD2):
    """
    多组有序分组秩和检验接口 - Kruskal-Wallis H 检验
    
    功能：比较多个有序分组的差异
    参数：RSParamOD2 - 包含多组有序分组数据的参数对象
    返回：包含多组有序分组秩和检验结果的字典
    临床应用：用于比较多个等级资料的差异
    """
    try:
        # 调用底层多组有序分组秩和检验算法函数
        result = cal_result_rs_od2(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="多组有序分组秩和检验 OD2 成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_od2 failed: {e}")
        return ResultData(code=500, message=f"多组有序分组秩和检验 OD2 失败：{e}", result=None)


@router.post("/kwh", response_model=ResultData)
async def rs_kwh(param: RSParamKWH):
    """
    Kruskal-Wallis H 检验接口
    
    功能：比较多个独立样本的差异
    参数：RSParamKWH - 包含多组独立数据的参数对象
    返回：包含 Kruskal-Wallis H 检验结果的字典
    临床应用：用于比较多组独立受试者的指标差异
    """
    try:
        # 调用底层 Kruskal-Wallis H 检验算法函数
        result = cal_result_rs_kwh(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Kruskal-Wallis H 检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_kwh failed: {e}")
        return ResultData(code=500, message=f"Kruskal-Wallis H 检验失败：{e}", result=None)


@router.post("/fm", response_model=ResultData)
async def rs_fm(param: RSParamFM):
    """
    Friedman 检验接口
    
    功能：比较多个相关样本的差异
    参数：RSParamFM - 包含多个相关样本数据的参数对象
    返回：包含 Friedman 检验结果的字典
    临床应用：用于比较同一受试者在多个时间点或处理条件下的指标差异
    """
    try:
        # 调用底层 Friedman 检验算法函数
        result = cal_result_rs_fm(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Friedman 检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回 400 状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回 500 状态码
        LoggerHelper.error(f"rs_fm failed: {e}")
        return ResultData(code=500, message=f"Friedman 检验失败：{e}", result=None)