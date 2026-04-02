"""
临床t检验统计分析模块控制器 (Clinical t-Test Statistics Analysis Module Controller)

本控制器提供全面的临床t检验统计分析功能接口，涵盖单样本、配对样本、独立样本等多种t检验方法，
为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. 单样本t检验: 比较样本均数与已知总体均数的差异
2. 配对样本t检验: 比较配对数据的差异
3. 独立样本t检验: 比较两独立样本的均数差异
4. t检验参数检验: 基于给定参数进行t检验

【API设计原则】:
- 参数校验: 使用Pydantic模型进行请求参数校验
- 错误处理: 统一的异常捕获和错误响应格式
- 日志记录: 详细的操作日志记录便于调试和监控
- 结果封装: 使用统一的响应格式包装计算结果

【技术架构】:
- 采用FastAPI框架提供高性能异步API服务
- 与底层统计算法模块解耦，通过函数调用实现功能
- 遵循RESTful API设计规范，提供清晰的端点命名

【临床应用价值】:
- 为临床研究提供标准化的t检验工具
- 支持多种t检验的应用场景
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from fastapi import APIRouter, HTTPException, status, Body
from app.schemas.result_data import ResultData

from app.schemas.request_data.t_param import (
    TParamSingle1,
    TParamSingle2,
    TParamPaired,
    TParamIndep,
    TParamP,
)
from app.stats.t.t_stats_paired import cal_result_t_paired
from app.stats.t.t_stats_indep import cal_result_t_indep
from app.stats.t.t_stats_single1 import cal_result_t_single1
from app.stats.t.t_stats_single2 import cal_result_t_single2
from app.stats.t.t_stats_p import cal_result_t_p
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/t", tags=["Stats-T"])


@router.post("/single1", response_model=ResultData)
async def t_single1(param: TParamSingle1):
    """
    单样本t检验接口 - 已知样本统计量
    
    功能: 比较样本均数与已知总体均数的差异
    参数: TParamSingle1 - 包含总体均数、样本量、样本均数和样本标准差的参数对象
    返回: 包含单样本t检验结果的字典
    临床应用: 用于比较单组数据与已知参考值的差异
    """
    try:
        # 调用底层单样本t检验算法函数
        result = cal_result_t_single1(param.pop_mean, param.n, param.sample_mean, param.sample_std)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="t单样本检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"t_single1 failed: {e}")
        return ResultData(code=500, message=f"t单样本检验失败: {e}", result=None)


@router.post("/single2", response_model=ResultData)
async def t_single2(param: TParamSingle2):
    """
    单样本t检验接口 - 基于原始数据
    
    功能: 基于原始数据比较样本均数与已知总体均数的差异
    参数: TParamSingle2 - 包含总体均数和原始数据列表的参数对象
    返回: 包含基于原始数据的单样本t检验结果的字典
    临床应用: 用于基于原始数据比较单组数据与已知参考值的差异
    """
    try:
        # 提取总体均数和原始数据
        pop_mean = param.pop_mean
        raw_data = param.stats_data_list[0].data_list
        
        # 调用底层基于原始数据的单样本t检验算法函数
        result = cal_result_t_single2(pop_mean, raw_data)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="t单样本检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"t_single2 failed: {e}")
        return ResultData(code=500, message=f"t单样本检验失败: {e}", result=None)


@router.post("/paired", response_model=ResultData)
async def t_paired(param: TParamPaired):
    """
    配对样本t检验接口
    
    功能: 比较配对数据的差异
    参数: TParamPaired - 包含两组配对数据的参数对象
    返回: 包含配对样本t检验结果的字典
    临床应用: 用于比较同一受试者治疗前后的指标变化
    """
    try:
        # 提取两组配对数据
        data1 = param.stats_data_list[0].data_list
        data2 = param.stats_data_list[1].data_list
        
        # 调用底层配对样本t检验算法函数
        result = cal_result_t_paired(data1, data2)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="配对检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"t_paired failed: {e}")
        return ResultData(code=500, message=f"配对检验失败: {e}", result=None)


@router.post("/indep", response_model=ResultData)
async def t_indep(param: TParamIndep):
    """
    独立样本t检验接口
    
    功能: 比较两独立样本的均数差异
    参数: TParamIndep - 包含两组独立数据的参数对象
    返回: 包含独立样本t检验结果的字典
    临床应用: 用于比较两组独立受试者的指标差异
    """
    try:
        # 提取两组独立数据
        data1 = param.stats_data_list[0].data_list
        data2 = param.stats_data_list[1].data_list
        
        # 调用底层独立样本t检验算法函数
        result = cal_result_t_indep(data1, data2)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="独立样本t检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"t_indep failed: {e}")
        return ResultData(code=500, message=f"独立样本t检验失败: {e}", result=None)


@router.post("/param", response_model=ResultData)
async def t_param(param: TParamP):
    """
    t检验参数检验接口
    
    功能: 基于给定参数进行t检验
    参数: TParamP - 包含t值、自由度和检验方向的参数对象
    返回: 包含t检验参数检验结果的字典
    临床应用: 用于基于已知统计量进行显著性检验
    """
    try:
        # 调用底层t检验参数检验算法函数
        result = cal_result_t_p(param.t_value, param.df, param.test_type)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="参数检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"t_param failed: {e}")
        return ResultData(code=500, message=f"参数检验失败: {e}", result=None)