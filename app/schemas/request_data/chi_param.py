from pydantic import BaseModel

#1.1 卡方检验 - 参数
class ChiSquareParam1_1(BaseModel):
    # 样本1 a: 发生数
    # 样本1 b: 未发生数
    # 样本2 c: 发生数
    # 样本2 d: 未发生数
    a : int
    b : int
    c : int
    d : int

#1.2 卡方检验（发生数） - 参数
class ChiSquareParam1_2(BaseModel):
    # 样本1 n1: 样本数
    # 样本1 a: 发生数
    # 样本2 n2: 样本数
    # 样本2 c: 发生数
    n1 : int
    a : int
    n2 : int
    c : int

#1.3 卡方检验（发生率） - 参数
class ChiSquareParam1_3(BaseModel):
    # 样本1 n1: 样本数
    # 样本1 per1: 发生率百分数
    # 样本2 n2: 样本数
    # 样本2 per2: 发生率百分数
    n1 : int
    per1 : float
    n2 : int
    per2 : float

#2 卡方检验 - RC表分析 - 参数
class ChiSquareParamRC(BaseModel):
    # 行数
    row_num : int
    # 列数
    col_num : int
    # 数据列表，按行优先顺序排列
    data_list : list[int]

#3 卡方检验 - 配对资料 - 参数
class ChiSquareParamPaired(BaseModel):
    # a: 样本1阳性/成功，样本2阳性/成功的配对数
    # b: 样本1阳性/成功，样本2阴性/失败的配对数
    # c: 样本1阴性/失败，样本2阳性/成功的配对数
    # d: 样本1阴性/失败，样本2阴性/失败的配对数
    a : int
    b : int
    c : int
    d : int

#4 卡方检验 - RR列表 - 参数
class ChiSquareParamRR(BaseModel):
    # n : 分类数
    n : int
    # data_list : 数据列表，长度应为n*n，按行优先顺序排列
    data_list : list[int]

#5 卡方检验 - Fisher检测 - 参数
class ChiSquareParamFisher(BaseModel):
    # 样本1 a: 发生数
    # 样本1 b: 未发生数
    # 样本2 c: 发生数
    # 样本2 d: 未发生数
    a : int
    b : int
    c : int
    d : int

#6 卡方检验 - P值计算 - 参数
class ChiSquareParamP(BaseModel):
    # chi : 卡方值
    # df : 自由度
    chi : float
    df : int

