"""
临床相关分析统计分析模块控制器 (Clinical Correlation Analysis Statistics Module Controller)

本控制器提供全面的临床相关分析统计分析功能接口，涵盖Pearson相关、Spearman秩相关、Kendall秩相关等多种相关分析方法，
为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. Pearson相关分析: 分析两个连续变量之间的线性关系
2. Spearman秩相关分析: 分析两个变量之间的单调关系（非参数方法）
3. Kendall秩相关分析: 分析两个变量之间的秩序一致性（非参数方法）

【API设计原则】:
- 参数校验: 使用Pydantic模型进行请求参数校验
- 错误处理: 统一的异常捕获和错误响应格式
- 日志记录: 详细的日志记录便于调试和监控
- 结果封装: 使用统一的响应格式包装计算结果

【技术架构】:
- 采用FastAPI框架提供高性能异步API服务
- 与底层统计算法模块解耦，通过函数调用实现功能
- 遵循RESTful API设计规范，提供清晰的端点命名

【临床应用价值】:
- 为临床研究提供标准化的相关分析工具
- 支持多种相关分析的应用场景
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from fastapi import APIRouter, HTTPException, status, Body
from app.schemas.result_data import ResultData
from app.schemas.request_data.cor_param import CorParamPearson, CorParamSpearman, CorParamKendall
from app.stats.cor.cor_stats_pearson import cal_result_cor_pearson
from app.stats.cor.cor_stats_spearman import cal_result_cor_spearman
from app.stats.cor.cor_stats_kendall import cal_result_cor_kendall
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/cor", tags=["Stats-Cor"])


@router.post("/pearson", response_model=ResultData)
async def cor_pearson(param: CorParamPearson):
    """
    Pearson相关分析接口
    
    功能: 分析两个连续变量之间的线性关系
    参数: CorParamPearson - 包含两个变量数据列表的参数对象
    返回: 包含Pearson相关分析结果的字典
    临床应用: 用于分析两个连续变量（如年龄与血压）之间的线性关系
    """
    try:
        # 提取数据列表并调用底层Pearson相关分析算法函数
        data_list = [data.data_list for data in param.stats_data_list]
        result = cal_result_cor_pearson(data_list)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Pearson相关分析成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"cor_pearson failed: {e}")
        return ResultData(code=500, message=f"Pearson相关分析失败: {e}", result=None)


@router.post("/spearman", response_model=ResultData)
async def cor_spearman(param: CorParamSpearman):
    """
    Spearman秩相关分析接口
    
    功能: 分析两个变量之间的单调关系（非参数方法）
    参数: CorParamSpearman - 包含两个变量数据列表的参数对象
    返回: 包含Spearman秩相关分析结果的字典
    临床应用: 用于分析非正态分布数据或等级数据之间的单调关系
    """
    try:
        # 提取数据列表并调用底层Spearman相关分析算法函数
        data_list = [data.data_list for data in param.stats_data_list]
        result = cal_result_cor_spearman(data_list)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Spearman相关分析成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"cor_spearman failed: {e}")
        return ResultData(code=500, message=f"Spearman相关分析失败: {e}", result=None)


@router.post("/kendall", response_model=ResultData)
async def cor_kendall(param: CorParamKendall):
    """
    Kendall秩相关分析接口
    
    功能: 分析两个变量之间的秩序一致性（非参数方法）
    参数: CorParamKendall - 包含两个变量数据列表的参数对象
    返回: 包含Kendall秩相关分析结果的字典
    临床应用: 用于分析小样本或存在异常值的等级数据之间的秩序关系
    """
    try:
        # 提取数据列表并调用底层Kendall相关分析算法函数
        data_list = [data.data_list for data in param.stats_data_list]
        result = cal_result_cor_kendall(data_list)
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Kendall相关分析成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"cor_kendall failed: {e}")
        return ResultData(code=500, message=f"Kendall相关分析失败: {e}", result=None)