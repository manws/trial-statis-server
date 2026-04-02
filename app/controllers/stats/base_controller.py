"""
临床统计分析基础模块控制器 (Clinical Statistical Analysis Base Module Controller)

本控制器提供全面的临床统计分析基础功能接口，涵盖描述性统计、频率分析、正态性检验、
二项分布和泊松分布等多种统计方法，为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. 描述性统计: 计算均值、中位数、标准差等基本统计量
2. 频率分析: 生成频数分布表和直方图数据
3. 正态性检验: 使用矩法检验数据是否符合正态分布
4. 二项分布分析: 包括概率计算、区间估计、单样本和双样本比较
5. 泊松分布分析: 包括概率计算、区间估计、单样本和双样本比较

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
- 为临床研究提供标准化的统计分析工具
- 支持多种常见的统计方法满足不同分析需求
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from typing import List, Dict, Any
from fastapi import APIRouter, HTTPException, status, Body

from app.schemas.result_data import ResultData
from app.schemas.request_data.base_param import (
    BaseParamBino1,
    BaseParamBino2,
    BaseParamBino3,
    BaseParamBino4,
    BaseParamPois1,
    BaseParamPois2,
    BaseParamPois3,
    BaseParamPois4,
)
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

from app.stats.base.base_stats_des import cal_result_des
from app.stats.base.base_stats_freq import cal_result_freq
from app.stats.base.base_stats_normal import cal_result_normal
from app.stats.base.base_stats_bino1 import cal_result_bino1
from app.stats.base.base_stats_bino2 import cal_result_bino2
from app.stats.base.base_stats_bino3 import cal_result_bino3
from app.stats.base.base_stats_bino4 import cal_result_bino4
from app.stats.base.base_stats_pois1 import cal_result_pois1
from app.stats.base.base_stats_pois2 import cal_result_pois2
from app.stats.base.base_stats_pois3 import cal_result_pois3
from app.stats.base.base_stats_pois4 import cal_result_pois4


# 创建基础统计分析API路由器，设置路由前缀和标签
router = APIRouter(prefix="/base", tags=["Stats-Base"])

# 统计接口统一使用 `app.controllers.stats.utils._wrap` 进行返回包装，
# 以保持与示例接口相同的结果结构，并在必要时回显输入数据。


@router.post("/des", response_model=ResultData)
async def des(data: List[float] = Body(..., embed=True)):
    """
    描述性统计分析接口
    
    功能: 计算数据的基本统计特征，包括集中趋势、离散程度、分布形态等指标
    参数: 
      - data: List[float] - 输入的数值列表，代表一组连续的观测数据
    返回: 包含描述性统计指标的字典，如样本量、均值、标准差、四分位数等
    临床应用: 用于探索性数据分析，了解数据的基本分布特征
    """
    try:
        # 调用底层描述性统计算法函数计算结果
        result = cal_result_des(data)
        # 使用统一的包装函数返回结果
        return _wrap(result, message="描述性统计报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"descriptive_report failed: {e}")
        return ResultData(code=500, message=f"描述性统计失败: {e}", result=None)


@router.post("/freq", response_model=ResultData)
async def frequency(data: List[float] = Body(..., embed=True), num_intervals: int = Body(10, embed=True)):
    """
    频率分析接口
    
    功能: 计算数据的频数分布表，包括各组段的频数、频率和累积频率
    参数: 
      - data: List[float] - 输入的数值列表，代表一组连续的观测数据
      - num_intervals: int - 分组数量，默认为10组，可根据数据特点调整
    返回: 包含频率统计分析指标的字典，如各组段的频数、频率等
    临床应用: 用于了解数据的分布规律，为后续统计分析提供基础
    """
    try:
        # 调用底层频率统计算法函数计算结果
        result = cal_result_freq(data, num_intervals)
        # 使用统一的包装函数返回结果
        return _wrap(result, message="频数分析报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"frequency_report failed: {e}")
        return ResultData(code=500, message=f"频数分析失败: {e}", result=None)


@router.post("/normal", response_model=ResultData)
async def normality(data: List[float] = Body(..., embed=True)):
    """
    正态性检验接口
    
    功能: 使用矩法检验数据是否符合正态分布
    参数: 
      - data: List[float] - 待检验的浮点数列表
    返回: 包含正态性检验统计分析指标的字典，如偏度、峰度及其显著性检验结果
    临床应用: 用于验证数据是否符合正态分布假设，为后续参数检验提供依据
    """
    try:
        # 调用底层正态性检验算法函数计算结果
        result = cal_result_normal(data)
        # 使用统一的包装函数返回结果
        return _wrap(result, message="正态检验报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"normality_report failed: {e}")
        return ResultData(code=500, message=f"正态检验失败: {e}", result=None)


# 二项分布分析接口组
@router.post("/bino1", response_model=ResultData)
async def bino1(param: BaseParamBino1):
    """
    二项分布概率计算接口
    
    功能: 计算二项分布的概率值
    参数: BaseParamBino1 - 包含二项分布计算所需参数的对象
    返回: 二项分布概率计算结果的字典
    临床应用: 用于计算特定事件发生的概率，如治疗成功的概率
    """
    try:
        # 调用底层二项分布概率计算算法函数
        result = cal_result_bino1(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="二项分布概率计算报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"bino1 failed: {e}")
        return ResultData(code=500, message=f"二项分布概率计算失败: {e}", result=None)

@router.post("/bino2", response_model=ResultData)
async def bino2(param: BaseParamBino2):
    """
    总体率区间估计接口
    
    功能: 计算总体率的置信区间
    参数: BaseParamBino2 - 包含总体率区间估计所需参数的对象
    返回: 总体率区间估计结果的字典
    临床应用: 用于估计总体成功率的置信区间，如治愈率的置信区间
    """
    try:
        # 调用底层总体率区间估计算法函数
        result = cal_result_bino2(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="总体率区间估计报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"bino2 failed: {e}")
        return ResultData(code=500, message=f"总体率区间估计失败: {e}", result=None)

@router.post("/bino3", response_model=ResultData)
async def bino3(param: BaseParamBino3):
    """
    样本率与总体率比较接口
    
    功能: 比较样本率与总体率是否存在显著差异
    参数: BaseParamBino3 - 包含样本率与总体率比较所需参数的对象
    返回: 样本率与总体率比较结果的字典
    临床应用: 用于比较样本成功率与已知总体成功率的差异
    """
    try:
        # 调用底层样本率与总体率比较算法函数
        result = cal_result_bino3(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="样本率与总体率比较报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"bino3 failed: {e}")
        return ResultData(code=500, message=f"样本率与总体率比较失败: {e}", result=None)

@router.post("/bino4", response_model=ResultData)
async def bino4(param: BaseParamBino4):
    """
    两样本率比较接口
    
    功能: 比较两个样本率是否存在显著差异
    参数: BaseParamBino4 - 包含两样本率比较所需参数的对象
    返回: 两样本率比较结果的字典
    临床应用: 用于比较两种治疗方法的成功率是否存在差异
    """
    try:
        # 调用底层两样本率比较算法函数
        result = cal_result_bino4(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="两样本率比较报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"bino4 failed: {e}")
        return ResultData(code=500, message=f"两样本率比较失败: {e}", result=None)

# 泊松分布分析接口组
@router.post("/pois1", response_model=ResultData)
async def pois1(param: BaseParamPois1, x: int = Body(..., embed=True)):
    """
    泊松分布概率计算接口
    
    功能: 计算泊松分布的概率值
    参数:
      - param: BaseParamPois1 - 包含泊松分布计算所需参数的对象
      - x: int - 要计算概率的事件数
    返回: 泊松分布概率计算结果的字典
    临床应用: 用于计算稀有事件发生的概率，如不良反应发生概率
    """
    try:
        # 调用底层泊松分布概率计算算法函数
        result = cal_result_pois1(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="泊松分布概率计算报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"pois1 failed: {e}")
        return ResultData(code=500, message=f"泊松分布概率计算失败: {e}", result=None)

@router.post("/pois2", response_model=ResultData)
async def pois2(param: BaseParamPois2, observed_events: int = Body(..., embed=True)):
    """
    总体均数区间估计接口（泊松分布）
    
    功能: 计算泊松分布总体均数的置信区间
    参数:
      - param: BaseParamPois2 - 包含泊松分布区间估计所需参数的对象
      - observed_events: int - 观察到的事件数
    返回: 总体均数区间估计结果的字典
    临床应用: 用于估计稀有事件发生率的置信区间
    """
    try:
        # 调用底层泊松分布区间估计算法函数
        result = cal_result_pois2(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="总体均数区间估计报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"pois2 failed: {e}")
        return ResultData(code=500, message=f"总体均数区间估计失败: {e}", result=None)

@router.post("/pois3", response_model=ResultData)
async def pois3(param: BaseParamPois3):
    """
    样本均数与总体均数比较接口（泊松分布）
    
    功能: 比较样本均数与总体均数是否存在显著差异
    参数: BaseParamPois3 - 包含泊松分布单样本比较所需参数的对象
    返回: 样本均数与总体均数比较结果的字典
    临床应用: 用于比较观察到的事件率与已知总体事件率的差异
    """
    try:
        # 调用底层泊松分布单样本比较算法函数
        result = cal_result_pois3(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="样本均数与总体均数比较报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"pois3 failed: {e}")
        return ResultData(code=500, message=f"样本均数与总体均数比较失败: {e}", result=None)

@router.post("/pois4", response_model=ResultData)
async def pois4(param: BaseParamPois4):
    """
    两样本均数比较接口（泊松分布）
    
    功能: 比较两个样本的泊松分布均数是否存在显著差异
    参数: BaseParamPois4 - 包含泊松分布两样本比较所需参数的对象
    返回: 两样本均数比较结果的字典
    临床应用: 用于比较两种条件下稀有事件发生率的差异
    """
    try:
        # 调用底层泊松分布两样本比较算法函数
        result = cal_result_pois4(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="两样本均数比较报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"pois4 failed: {e}")
        return ResultData(code=500, message=f"两样本均数比较失败: {e}", result=None)