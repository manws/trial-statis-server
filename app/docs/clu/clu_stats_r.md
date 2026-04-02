# 临床R型聚类分析统计模块 (Clinical R-Type Cluster Analysis Statistics Module)

## 模块概述

临床R型聚类分析统计模块提供全面的R型聚类分析统计功能，用于临床数据中变量的聚类分析，是临床试验数据分析的重要组成部分。R型聚类在医学研究中广泛应用，如对变量进行归类、识别相关变量组等。

## 模块功能

### 1. 标准化处理
- 支持多种标准化方法
- 消除变量量纲影响

### 2. 相似性计算
- 计算变量间的多种相似性
- 适应不同类型数据

### 3. 聚类算法
- 实现不同的聚类方法
- 提供灵活的聚类选项

### 4. 结果解释
- 提供聚类结果的临床意义解释
- 生成详细的分析报告

## 临床应用价值

- **变量归类**：根据相关性对变量进行分类
- **指标筛选**：识别高度相关的指标组
- **维度约简**：减少冗余变量
- **生物标志物分组**：对生物标志物进行功能分组

## 统计方法选择指南

### 1. R型聚类适用条件
- 变量间存在一定相关性
- 样本量大于变量数
- 数据质量良好

### 2. 标准化方法选择
- 不标准化：变量单位相同且量级相近
- Z分数：变量单位不同，需消除量纲影响
- 范围标准化：变量范围差异较大

### 3. 临床应用场景
- 根据基因表达相关性对基因进行分组
- 根据症状相关性对症状进行分类
- 根据生物标志物相关性进行功能分组

## 主要函数说明

### 标准化数据函数 - standardize_data

根据指定类型对数据进行标准化处理。在临床研究中，这用于消除不同变量间的量纲差异，使聚类结果更具可比性。

**参数**：
- `data_matrix`: 输入的数据矩阵
- `std_type`: 标准化类型（0-不标准化，1-Z分数，2-范围-1到1标准化，3-范围0到1标准化，4-最大值为1标准化，5-平均值为1标准化，6-标准差为1标准化）

**返回值**：
- `np.ndarray`: 标准化后的数据矩阵

### 计算相似性矩阵函数 - compute_similarity_matrix

计算相似性矩阵。在临床研究中，这用于量化变量间的相似性或相关性，是聚类分析的基础。

**参数**：
- `data_matrix`: 输入的数据矩阵
- `sim_type`: 相似性类型（0-Pearson相关系数距离，1-夹角余弦距离）

**返回值**：
- `np.ndarray`: 相似性矩阵

### 层次聚类函数 - hierarchical_cluster

执行层次聚类。在临床研究中，这用于构建变量间的层级关系，形成聚类树状图。

**参数**：
- `data_matrix`: 输入的距离矩阵
- `clu_method`: 聚类方法（0-平均距离，1-最短距离，2-最长距离）

**返回值**：
- `np.ndarray`: 聚类链接矩阵

### 获取聚类标签函数 - get_clusters

根据链接矩阵和聚类数获取聚类标签。在临床研究中，这用于确定每个变量所属的聚类。

**参数**：
- `linkage_matrix`: 聚类链接矩阵
- `k`: 期望的聚类数

**返回值**：
- `np.ndarray`: 聚类标签数组

### R型聚类分析函数 - r_clustering_from_stats_data

从统计数据执行R型聚类分析。在临床研究中，这用于对变量进行聚类，识别具有相似特征的变量组。

**参数**：
- `stats_data_list`: 包含聚类样本数据的列表
- `stats_name_list`: 样本名称列表
- `std_type`: 标准化类型
- `clu_method`: 聚类方法
- `dis_type`: 相似性类型
- `k`: 聚类数量

**返回值**：
- `Dict[str, Any]`: 聚类分析结果

### 结果计算函数 - cal_result_clu_r

生成R型聚类分析统计分析的完整报告字典。此函数整合了R型聚类分析的所有关键指标，生成标准化的字典格式报告，适用于临床研究报告的需求，提供全面的聚类分析结果。报告包括输入参数、聚类结果和方法描述等信息，便于临床医生和研究人员快速理解R型聚类分析的特征。

**参数**：
- `param`: CluParamR对象，包含std_type, clu_method, dis_type, k, stats_name, stats_data_list

**返回值**：
- 包含R型聚类分析统计分析指标的字典，包括：
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
from app.stats.clu.clu_stats_r import cal_result_clu_r
from app.schemas.request_data.clu_param import CluParamR
from app.schemas.request_data.stats_data import StatsData

# 创建统计参数对象
data1 = StatsData(name="Var1", data_list=[1.2, 2.3, 3.4, 4.5, 5.6])
data2 = StatsData(name="Var2", data_list=[1.1, 2.2, 3.3, 4.4, 5.5])
data3 = StatsData(name="Var3", data_list=[5.1, 6.2, 7.3, 8.4, 9.5])

# 创建参数对象
param = CluParamR(
    std_type=1,  # Z分数标准化
    clu_method=0,  # 平均距离法
    dis_type=0,  # Pearson相关系数距离
    k=2,  # 2个聚类
    stats_data_list=[data1, data2, data3]
)

# 计算R型聚类分析
result = cal_result_clu_r(param)

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