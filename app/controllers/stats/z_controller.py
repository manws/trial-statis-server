"""
临床Z检验统计分析模块控制器 (Clinical Z-Test Statistics Analysis Module Controller)

本控制器提供全面的临床Z检验统计分析功能接口，涵盖单样本、独立样本等多种Z检验方法，
为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. 单样本Z检验: 比较样本均数与已知总体均数的差异
2. 独立样本Z检验: 比较两独立样本的均数差异

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
- 为临床研究提供标准化的Z检验工具
- 支持多种Z检验的应用场景
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from fastapi import APIRouter, HTTPException, status
from app.schemas.result_data import ResultData

from app.schemas.request_data.z_param import (
    ZParamSingle1,
    ZParamSingle2,
    ZParamIndep1,
    ZParamIndep2,
)
from app.stats.z.z_stats_single1 import cal_result_z_single1
from app.stats.z.z_stats_single2 import cal_result_z_single2
from app.stats.z.z_stats_indep1 import cal_result_z_indep1
from app.stats.z.z_stats_indep2 import cal_result_z_indep2
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/z", tags=["Stats-Z"])


@router.post("/single1", response_model=ResultData)
async def z_single1(param: ZParamSingle1):
    """
    单样本Z检验接口 - 已知样本统计量
    
    功能: 比较样本均数与已知总体均数的差异
    参数: ZParamSingle1 - 包含总体均数、总体标准差、样本量和样本均数的参数对象
    返回: 包含单样本Z检验结果的字典
    临床应用: 用于比较单组数据与已知参考值的差异
    """
    try:
        # 调用底层单样本Z检验算法函数
        result = cal_result_z_single1(param.pop_mean, param.pop_std, param.n, param.sample_mean)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="z单样本1检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"z_single1 failed: {e}")
        return ResultData(code=500, message=f"z单样本1检验失败: {e}", result=None)


@router.post("/single2", response_model=ResultData)
async def z_single2(param: ZParamSingle2):
    """
    单样本Z检验接口 - 基于原始数据
    
    功能: 基于原始数据比较样本均数与已知总体均数的差异
    参数: ZParamSingle2 - 包含总体均数、总体标准差和原始数据列表的参数对象
    返回: 包含基于原始数据的单样本Z检验结果的字典
    临床应用: 用于基于原始数据比较单组数据与已知参考值的差异
    """
    try:
        # 提取总体均数、总体标准差和原始数据
        pop_mean = param.pop_mean
        pop_std = param.pop_std
        raw_data = param.stats_data_list[0].data_list
        
        # 调用底层基于原始数据的单样本Z检验算法函数
        result = cal_result_z_single2(pop_mean, pop_std, raw_data)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="z单样本2检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"z_single2 failed: {e}")
        return ResultData(code=500, message=f"z单样本2检验失败: {e}", result=None)


@router.post("/indep1", response_model=ResultData)
async def z_indep1(param: ZParamIndep1):
    """
    独立样本Z检验接口 - 已知样本统计量
    
    功能: 比较两独立样本的均数差异
    参数: ZParamIndep1 - 包含两组样本的样本量、均数和标准差的参数对象
    返回: 包含独立样本Z检验结果的字典
    临床应用: 用于比较两组独立受试者的指标差异
    """
    try:
        # 提取两组样本的参数
        result = cal_result_z_indep1(
            param.n1, param.mean1, param.std1,
            param.n2, param.mean2, param.std2
        )
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="z独立样本1检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"z_indep1 failed: {e}")
        return ResultData(code=500, message=f"z独立样本1检验失败: {e}", result=None)


@router.post("/indep2", response_model=ResultData)
async def z_indep2(param: ZParamIndep2):
    """
    独立样本Z检验接口 - 基于原始数据
    
    功能: 基于原始数据比较两独立样本的均数差异
    参数: ZParamIndep2 - 包含两组原始数据的参数对象
    返回: 包含基于原始数据的独立样本Z检验结果的字典
    临床应用: 用于基于原始数据比较两组独立受试者的指标差异
    """
    try:
        # 提取两组原始数据
        data1 = param.stats_data_list[0].data_list
        data2 = param.stats_data_list[1].data_list
        
        # 调用底层基于原始数据的独立样本Z检验算法函数
        result = cal_result_z_indep2(data1, data2)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="z独立样本2检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"z_indep2 failed: {e}")
        return ResultData(code=500, message=f"z独立样本2检验失败: {e}", result=None)
