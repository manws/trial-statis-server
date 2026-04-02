from pydantic import BaseModel
from typing import List, Optional

# 简单随机参数
class RandomParamSimple(BaseModel):
    # total_sample_size: 总样本量
    total_sample_size: int
    # group_num: 组数
    group_num: int
    # group_names: 组名称列表，可选
    group_names: Optional[List[str]] = None

# 区组随机参数
class RandomParamBlock(BaseModel):
    # total_sample_size: 总样本量
    total_sample_size: int
    # group_num: 组数
    group_num: int
    # block_size: 区组大小
    block_size: int
    # group_names: 组名称列表，可选
    group_names: Optional[List[str]] = None

# 分层区组随机参数
class RandomParamSBlock(BaseModel):
    # total_sample_size: 总样本量
    total_sample_size: int
    # group_num: 组数
    group_num: int
    # block_size: 区组大小
    block_size: int
    # strata_list: 分层变量列表，每个strata是变量名
    # 变量都是二分变量：0或者1
    strata_list: List[str]
    # group_names: 组名称列表，可选
    group_names: Optional[List[str]] = None

# 最小化随机参数
class RandomParamMin(BaseModel):
    # total_sample_size: 总样本量
    total_sample_size: int
    # group_num: 组数
    group_num: int
    # covariates: 协变量列表，用于最小化
    covariates: List[str]
    # weights: 权重列表，对应协变量，可选，默认等权重
    weights: Optional[List[float]] = None
    # group_names: 组名称列表，可选
    group_names: Optional[List[str]] = None
