"""
临床方差分析统计分析模块控制器 (Clinical ANOVA Statistics Analysis Module Controller)

本控制器提供全面的临床方差分析统计分析功能接口，涵盖完全随机设计(CRD)和随机区组设计(RBD)方差分析，
为临床研究数据分析提供标准化的API服务。

【模块功能概述】:
1. 完全随机设计方差分析(CRD): 比较多个处理组均值是否存在显著差异
2. 随机区组设计方差分析(RBD): 在控制区组效应后比较多个处理组均值差异

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
- 为临床研究提供标准化的方差分析工具
- 支持多种实验设计的方差分析需求
- 提供直观的API接口便于前端应用集成
- 确保统计结果的准确性和一致性

作者：Trial-Statis-Server Team
版本：1.0.0
最后更新：2026-04-01
"""

from typing import List
from fastapi import APIRouter, HTTPException, status, Body

from app.schemas.result_data import ResultData
from app.schemas.request_data.anova_param import AnovaParamCRD, AnovaParamRBD
from app.schemas.request_data.stats_data import StatsData
from app.stats.anova.anova_stats_crd import cal_result_anova_crd
from app.stats.anova.anova_stats_rbd import cal_result_anova_rbd
from app.utils.LoggerHelper import LoggerHelper
from app.controllers.stats.utils import _wrap

# 创建方差分析API路由器，设置路由前缀和标签
router = APIRouter(prefix="/anova", tags=["Stats-ANOVA"])


@router.post("/crd", response_model=ResultData)
async def crd_anova(data_list: List[List[float]] = Body(..., embed=True)):
    """
    完全随机设计方差分析接口
    
    功能: 比较多个处理组的均值是否存在显著差异
    参数: 
      - data_list: List[List[float]] - 每个子列表代表一个处理组的数据
    返回: 包含完全随机设计方差分析结果的字典
    临床应用: 用于比较不同治疗方案或干预措施的效果差异
    """
    try:
        # 构造参数对象并调用底层算法
        param = AnovaParamCRD(stats_data_list=[StatsData(field_name=f"组{i+1}", data_list=d) for i, d in enumerate(data_list)])
        result = cal_result_anova_crd(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, message="CRD方差分析报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"crd_anova failed: {e}")
        return ResultData(code=500, message=f"CRD方差分析失败: {e}", result=None)


@router.post("/rbd", response_model=ResultData)
async def rbd_anova(data_list: List[List[float]] = Body(..., embed=True)):
    """
    随机区组设计方差分析接口
    
    功能: 在控制区组效应后比较多个处理组的均值差异
    参数: 
      - data_list: List[List[float]] - 每个子列表代表一个区组的数据，每个区组内包含各处理的观测值
    返回: 包含随机区组设计方差分析结果的字典
    临床应用: 用于比较不同治疗方案效果，同时控制可能的混杂因素（区组效应）
    """
    try:
        # 构造参数对象并调用底层算法
        param = AnovaParamRBD(stats_data_list=[StatsData(field_name=f"区组{i+1}", data_list=d) for i, d in enumerate(data_list)])
        result = cal_result_anova_rbd(param)
        # 使用统一的包装函数返回结果
        return _wrap(result, message="RBD方差分析报告生成成功")
    except ValueError as e:
        # 捕获参数验证错误并返回400状态码
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        # 记录错误日志并返回500状态码
        LoggerHelper.error(f"rbd_anova failed: {e}")
        return ResultData(code=500, message=f"RBD方差分析失败: {e}", result=None)