"""
临床聚类分析统计分析模块控制器 (Clinical Cluster Analysis Statistics Module Controller)

本控制器提供全面的临床聚类分析统计分析功能接口，涵盖Q型聚类和R型聚类等多种聚类分析方法，
为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. Q型聚类分析: 根据样本间的相似性对样本进行聚类
2. R型聚类分析: 根据变量间的相似性对变量进行聚类

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
- 为临床研究提供标准化的聚类分析工具
- 支持多种聚类分析的应用场景
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from fastapi import APIRouter, HTTPException, status, Body
from app.schemas.result_data import ResultData

from app.schemas.request_data.clu_param import CluParamQ, CluParamR
from app.stats.clu.clu_stats_q import cal_result_clu_q
from app.stats.clu.clu_stats_r import cal_result_clu_r
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

router = APIRouter(prefix="/clu", tags=["Stats-Clu"])


@router.post("/q", response_model=ResultData)
async def clu_q(param: CluParamQ):
    """
    Q型聚类分析接口
    
    功能: 根据样本间的相似性对样本进行聚类
    参数: CluParamQ - 包含Q型聚类分析所需参数的对象
    返回: 包含Q型聚类分析结果的字典
    临床应用: 用于根据多个临床指标对患者进行分层或识别疾病亚型
    """
    try:
        # 提取参数并调用底层Q型聚类分析算法函数
        data_list = [data.data_list for data in param.stats_data_list]
        stats_name_list = [data.field_name for data in param.stats_data_list]
        
        result = cal_result_clu_q(
            data_list,
            stats_name_list,
            param.std_type,
            param.clu_method,
            param.dis_type,
            param.k
        )
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="Q型聚类分析成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"clu_q failed: {e}")
        return ResultData(code=500, message=f"Q型聚类分析失败: {e}", result=None)


@router.post("/r", response_model=ResultData)
async def clu_r(param: CluParamR):
    """
    R型聚类分析接口
    
    功能: 根据变量间的相似性对变量进行聚类
    参数: CluParamR - 包含R型聚类分析所需参数的对象
    返回: 包含R型聚类分析结果的字典
    临床应用: 用于将功能相关的生物标志物或临床指标进行分组
    """
    try:
        # 提取参数并调用底层R型聚类分析算法函数
        data_list = [data.data_list for data in param.stats_data_list]
        stats_name_list = [data.field_name for data in param.stats_data_list]
        
        result = cal_result_clu_r(
            data_list,
            stats_name_list,
            param.std_type,
            param.clu_method,
            param.dis_type,
            param.k
        )
        # 使用统一的包装函数返回结果
        return _wrap(result, data=param.model_dump(), message="R型聚类分析成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"clu_r failed: {e}")
        return ResultData(code=500, message=f"R型聚类分析失败: {e}", result=None)