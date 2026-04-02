"""
临床卡方检验统计分析模块控制器 (Clinical Chi-Square Test Statistics Analysis Module Controller)

本控制器提供全面的临床卡方检验统计分析功能接口，涵盖四格表、列联表、配对卡方检验等多种卡方检验方法，
为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. 四格表卡方检验: 比较两个二分类变量之间的关联性
2. 列联表卡方检验: 比较多个分类变量之间的关联性
3. 配对卡方检验: 比较配对样本的分类差异
4. Fisher精确检验: 适用于小样本的精确概率检验
5. 优势比分析: 计算病例对照研究中的优势比及置信区间

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
- 为临床研究提供标准化的卡方检验工具
- 支持多种卡方检验的应用场景
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from fastapi import APIRouter, HTTPException, status, Body
from app.schemas.result_data import ResultData

from app.schemas.request_data.chi_param import (
    ChiSquareParam1_1,
    ChiSquareParam1_2,
    ChiSquareParam1_3,
    ChiSquareParamRC,
    ChiSquareParamPaired,
    ChiSquareParamRR,
    ChiSquareParamFisher,
    ChiSquareParamP,
)
from app.stats.chi.chi_stats_4 import cal_result_chi4
from app.stats.chi.chi_stats_rc import cal_result_chi_rc
from app.stats.chi.chi_stats_paired import cal_result_chi_paired
from app.stats.chi.chi_stats_rr import cal_result_chi_rr
from app.stats.chi.chi_stats_fisher import cal_result_fisher
from app.stats.chi.chi_stats_p import cal_result_chi_p
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/chi", tags=["Stats-Chi"])

# 使用统一的 _wrap 工具函数处理响应格式，和 example 控制器保持一致。


@router.post("/4_1", response_model=ResultData)
async def chi_4_1(param: ChiSquareParam1_1):
    """
    四格表卡方检验接口 - 标准四格表格式
    
    功能: 比较两个二分类变量之间的关联性
    参数: ChiSquareParam1_1 - 包含四格表数据(a, b, c, d)的参数对象
    返回: 包含四格表卡方检验结果的字典
    临床应用: 用于比较两组间的阳性率或发生率差异
    """
    try:
        # 调用底层四格表卡方检验算法函数
        result = cal_result_chi4(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="四格表卡方检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_4_1 failed: {e}")
        return ResultData(code=500, message=f"四格表卡方检验失败: {e}", result=None)


@router.post("/4_2", response_model=ResultData)
async def chi_4_2(param: ChiSquareParam1_2):
    """
    四格表卡方检验接口 - 发生数格式
    
    功能: 比较两个二分类变量之间的关联性（输入发生数和总数）
    参数: ChiSquareParam1_2 - 包含发生数和总数的参数对象
    返回: 包含四格表卡方检验结果的字典
    临床应用: 用于比较两组间的阳性率或发生率差异（按发生数和总数输入）
    """
    try:
        # 将发生数格式转换为四格表格式
        a = param.a
        b = param.n1 - param.a
        c = param.c
        d = param.n2 - param.c
        
        # 构造四格表参数对象并调用底层算法函数
        chi4_param = ChiSquareParam1_1(a=a, b=b, c=c, d=d)
        result = cal_result_chi4(chi4_param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="四格表卡方检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_4_2 failed: {e}")
        return ResultData(code=500, message=f"四格表卡方检验失败: {e}", result=None)


@router.post("/4_3", response_model=ResultData)
async def chi_4_3(param: ChiSquareParam1_3):
    """
    四格表卡方检验接口 - 发生率格式
    
    功能: 比较两个二分类变量之间的关联性（输入发生率和总数）
    参数: ChiSquareParam1_3 - 包含发生率和总数的参数对象
    返回: 包含四格表卡方检验结果的字典
    临床应用: 用于比较两组间的阳性率或发生率差异（按百分率和总数输入）
    """
    try:
        # 将发生率格式转换为四格表格式
        a = int(round(param.per1 / 100.0 * param.n1))
        c = int(round(param.per2 / 100.0 * param.n2))
        b = param.n1 - a
        d = param.n2 - c
        
        # 构造四格表参数对象并调用底层算法函数
        chi4_param = ChiSquareParam1_1(a=a, b=b, c=c, d=d)
        result = cal_result_chi4(chi4_param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="四格表卡方检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_4_3 failed: {e}")
        return ResultData(code=500, message=f"四格表卡方检验失败: {e}", result=None)


@router.post("/rc", response_model=ResultData)
async def chi_rc(param: ChiSquareParamRC):
    """
    RC列联表卡方检验接口
    
    功能: 比较R行C列分类变量之间的关联性
    参数: ChiSquareParamRC - 包含列联表数据的参数对象
    返回: 包含列联表卡方检验结果的字典
    临床应用: 用于比较多个分组间的分类变量差异
    """
    try:
        # 提取参数并调用底层列联表卡方检验算法函数
        result = cal_result_chi_rc(param.data_table)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="RC列联表卡方检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_rc failed: {e}")
        return ResultData(code=500, message=f"RC列联表卡方检验失败: {e}", result=None)


@router.post("/paired", response_model=ResultData)
async def chi_paired(param: ChiSquareParamPaired):
    """
    配对卡方检验接口
    
    功能: 比较配对样本的分类差异
    参数: ChiSquareParamPaired - 包含配对四格表数据的参数对象
    返回: 包含配对卡方检验结果的字典
    临床应用: 用于比较同一组受试者在不同时间或条件下的分类变化
    """
    try:
        # 提取参数并调用底层配对卡方检验算法函数
        result = cal_result_chi_paired(param.a, param.b, param.c, param.d)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="配对资料卡方检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_paired failed: {e}")
        return ResultData(code=500, message=f"配对资料卡方检验失败: {e}", result=None)


@router.post("/rr", response_model=ResultData)
async def chi_rr(param: ChiSquareParamRR):
    """
    RR值计算接口
    
    功能: 计算相对危险度及其置信区间
    参数: ChiSquareParamRR - 包含四格表数据用于RR计算的参数对象
    返回: 包含RR值计算结果的字典
    临床应用: 用于队列研究中评估暴露因素与结局的关联强度
    """
    try:
        # 提取参数并调用底层RR值计算算法函数
        result = cal_result_chi_rr(param.a, param.b, param.c, param.d)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="RR列表卡方检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_rr failed: {e}")
        return ResultData(code=500, message=f"RR列表卡方检验失败: {e}", result=None)


@router.post("/fisher", response_model=ResultData)
async def chi_fisher(param: ChiSquareParamFisher):
    """
    Fisher精确检验接口
    
    功能: 对于小样本进行精确概率检验
    参数: ChiSquareParamFisher - 包含四格表数据用于Fisher检验的参数对象
    返回: 包含Fisher精确检验结果的字典
    临床应用: 用于小样本或理论频数小于5的情况
    """
    try:
        # 提取参数并调用底层Fisher精确检验算法函数
        result = cal_result_fisher(param.a, param.b, param.c, param.d)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Fisher精确检验成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_fisher failed: {e}")
        return ResultData(code=500, message=f"Fisher精确检验失败: {e}", result=None)


@router.post("/p", response_model=ResultData)
async def chi_p(param: ChiSquareParamP):
    """
    P值计算接口
    
    功能: 根据卡方值和自由度计算P值
    参数: ChiSquareParamP - 包含卡方值和自由度的参数对象
    返回: 包含P值计算结果的字典
    临床应用: 用于已知卡方统计量时计算对应的P值
    """
    try:
        # 提取参数并调用底层P值计算算法函数
        result = cal_result_chi_p(param.chi_square_value, param.df)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="P值计算成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"chi_p failed: {e}")
        return ResultData(code=500, message=f"P值计算失败: {e}", result=None)