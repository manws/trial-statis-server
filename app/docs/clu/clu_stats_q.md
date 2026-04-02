# 临床Q型聚类分析统计模块 (Clinical Q-Type Cluster Analysis Statistics Module)

## 模块概述

临床Q型聚类分析统计模块提供全面的Q型聚类分析统计功能，用于临床数据中样本（个案）的聚类分析，是临床试验数据分析的重要组成部分。Q型聚类在医学研究中广泛应用，如对患者进行亚型分类、识别疾病亚群等。

## 模块功能

### 1. 标准化处理
- 支持多种标准化方法
- 消除变量量纲影响

### 2. 距离计算
- 计算样本间的多种距离
- 适应不同类型数据

### 3. 聚类算法
- 实现不同的聚类方法
- 提供灵活的聚类选项

### 4. 结果解释
- 提供聚类结果的临床意义解释
- 生成详细的分析报告

## 临床应用价值

- **患者分层**：根据临床指标对患者进行分类
- **疾病亚型识别**：识别具有相似特征的疾病亚型
- **治疗反应预测**：根据基线特征预测治疗反应
- **预后评估**：基于临床特征评估预后

## 统计方法选择指南

### 1. Q型聚类适用条件
- 样本量适中（通常<200个样本）
- 变量间存在一定相关性
- 数据质量良好

### 2. 标准化方法选择
- 不标准化：变量单位相同且量级相近
- Z分数：变量单位不同，需消除量纲影响
- 范围标准化：变量范围差异较大

### 3. 临床应用场景
- 根据基因表达谱对肿瘤进行分型
- 根据症状表现对疾病进行分类
- 根据生物标志物对患者进行风险分层

## 主要函数说明

### 标准化数据函数 - standardize_data

根据指定类型对数据进行标准化处理。在临床研究中，这用于消除不同变量间的量纲差异，使聚类结果更具可比性。

**参数**：
- `data_matrix`: 输入的数据矩阵
- `std_type`: 标准化类型（0-不标准化，1-Z分数，2-范围-1到1标准化，3-范围0到1标准化，4-最大值为1标准化，5-平均值为1标准化，6-标准差为1标准化）

**返回值**：
- `np.ndarray`: 标准化后的数据矩阵

### 计算距离矩阵函数 - compute_distance_matrix

计算距离矩阵。在临床研究中，这用于量化样本间的相似性或差异性，是聚类分析的基础。

**参数**：
- `data_matrix`: 输入的数据矩阵
- `dist_type`: 距离类型（0-欧氏距离，1-平方欧氏距离，2-曼哈顿距离，3-切比雪夫距离）

**返回值**：
- `np.ndarray`: 距离矩阵

### 层次聚类函数 - hierarchical_cluster

执行层次聚类。在临床研究中，这用于构建样本间的层级关系，形成聚类树状图。

**参数**：
- `data_matrix`: 输入的距离矩阵
- `clu_method`: 聚类方法（0-平均距离，1-最短距离，2-最长距离）

**返回值**：
- `np.ndarray`: 聚类链接矩阵

### 获取聚类标签函数 - get_clusters

根据链接矩阵和聚类数获取聚类标签。在临床研究中，这用于确定每个样本所属的聚类。

**参数**：
- `linkage_matrix`: 聚类链接矩阵
- `k`: 期望的聚类数

**返回值**：
- `np.ndarray`: 聚类标签数组

### Q型聚类分析函数 - q_clustering_from_stats_data

从统计数据执行Q型聚类分析。在临床研究中，这用于对样本进行聚类，识别具有相似特征的样本组。

**参数**：
- `stats_data_list`: 包含聚类变量数据的列表
- `stats_name_list`: 变量名称列表
- `std_type`: 标准化类型
- `clu_method`: 聚类方法
- `dis_type`: 距离类型
- `k`: 聚类数量

**返回值**：
- `Dict[str, Any]`: 聚类分析结果

### 结果计算函数 - cal_result_clu_q

生成Q型聚类分析统计分析的完整报告字典。此函数整合了Q型聚类分析的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的聚类分析结果。报告包括输入参数、聚类结果和方法描述等信息，便于临床医生和研究人员快速理解Q型聚类分析的特征。

**参数**：
- `param`: CluParamQ对象，包含std_type, clu_method, dis_type, k, stats_name, stats_data_list

**返回值**：
- 包含Q型聚类分析统计分析指标的字典，包括：
  - `table_name`: 报告表格名称
  - `input_parameters`: 输入参数信息
  - `clustering_results`: 聚类结果
  - `method_descriptions`: 方法描述
  - `remark`: 备注信息

## 结果解读注意事项

1. **聚类数确定**：结合专业知识和统计准则确定最优聚类数
2. **聚类稳定性**：评估聚类结果的稳定性和可重现性
3. **临床意义**：结合临床知识解释聚类结果的意义
4. **模型验证**：通过外部验证评估聚类的合理性
5. **结果报告**：清晰描述聚类的特征和临床含义

## 使用示例

```python
from app.stats.clu.clu_stats_q import cal_result_clu_q
from app.schemas.request_data.clu_param import CluParamQ
from app.schemas.request_data.stats_data import StatsData

# 创建统计参数对象
data1 = StatsData(name="Patient1", data_list=[1.2, 2.3, 3.4, 4.5])
data2 = StatsData(name="Patient2", data_list=[1.1, 2.2, 3.3, 4.4])
data3 = StatsData(name="Patient3", data_list=[5.1, 6.2, 7.3, 8.4])

# 创建参数对象
param = CluParamQ(
    std_type=1,  # Z分数标准化
    clu_method=0,  # 平均距离法
    dis_type=0,  # 欧氏距离
    k=2,  # 2个聚类
    stats_data_list=[data1, data2, data3]
)

# 计算Q型聚类分析
result = cal_result_clu_q(param)

# 查看结果
print(f"输入参数: {result['input_parameters']}")
print(f"聚类结果: {result['clustering_results']}")
print(f"方法描述: {result['method_descriptions']}")
```

## 相关标准和规范

- CONSORT声明：随机临床试验报告规范
- STROBE声明：观察性研究报告规范
- ICH E9指导原则：临床试验的统计学原则
- SAMPL指南：医学研究报告中的统计方法描述规范