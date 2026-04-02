from typing import Any, Dict

from app.schemas.result_data import ResultData
from app.utils.LoggerHelper import LoggerHelper


def _wrap(report: Any, data: Any = None, message: str = "操作成功") -> ResultData:
    """统一包装统计接口返回格式。

    :param report: 通过stats计算出的report对象，需实现`model_dump()`。
    :param data: 请求体中的原始数据，可选。如果提供则会放入result.data字段。
    :param message: 返回信息。
    :return: ResultData 实例
    """
    try:
        result: Dict[str, Any] = {"report": report.model_dump()}
        if data is not None:
            # 如果传入的是 Pydantic 对象或列表，将其序列化为 dict
            try:
                if hasattr(data, "model_dump"):
                    result["data"] = data.model_dump()
                else:
                    # 支持列表、字典等基本类型
                    result["data"] = data
            except Exception:
                result["data"] = data
        return ResultData(code=200, message=message, result=result)
    except Exception as e:
        LoggerHelper.error(f"_wrap failed: {e}")
        return ResultData(code=500, message=f"包装失败: {e}", result=None)
