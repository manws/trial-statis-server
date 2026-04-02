# Trial-Statis-Server 统计算法模块整理

## 📊 模块概览

本项目包含以下统计分析模块，共计约60个核心算法文件：

## 1️⃣ 方差分析 (ANVOA)
**目录**: `app/stats/anvoa/`
**文件数**: 3个
**核心算法**:
- `anvoa_stats_one.py` - 单因素方差分析
- `anvoa_stats_two.py` - 双因素方差分析  
- `anvoa_stats_mult.py` - 多因素方差分析

**特点**: 支持事后多重比较(LSD/SNK/Dunnett)，统一使用StatsData数据结构

## 2️⃣ 基础统计 (BASE)
**目录**: `app/stats/base/`
**文件数**: 11个
**核心算法**:
- `des_stats.py` - 描述性统计分析
- `freq_stats.py` - 频数统计分析
- `normal_stats.py` - 正态分布分析
- **二项分布系列**:
  - `bino_stats1.py` - 二项分布基础计算
  - `bino_stats2.py` - 二项分布置信区间
  - `bino_stats3.py` - 二项分布假设检验
  - `bino_stats4.py` - 二项分布样本量计算
- **泊松分布系列**:
  - `possion_stats1.py` - 泊松分布基础计算
  - `possion_stats2.py` - 泊松分布置信区间
  - `possion_stats3.py` - 泊松分布假设检验
  - `possion_stats4.py` - 泊松分布样本量计算

## 3️⃣ 卡方检验 (CHI)
**目录**: `app/stats/chi/`
**文件数**: 6个
**核心算法**:
- `chi_stats_p.py` - 拟合优度检验
- `chi_stats_rc.py` - 四格表卡方检验
- `chi_stats_rr.py` - R×C列联表卡方检验
- `chi_stats_paired.py` - 配对四格表McNemar检验
- `chi_stats_4.py` - 四格表精确概率法
- `chi_stats_fisher.py` - Fisher精确概率法

## 4️⃣ 聚类分析 (CLU)
**目录**: `app/stats/clu/`
**文件数**: 2个
**核心算法**:
- `clu_stats_q.py` - Q型聚类分析(样本聚类)
- `clu_stats_r.py` - R型聚类分析(变量聚类)

**特点**: 支持多种距离度量和聚类方法，已重构为统一的StatsData设计理念

## 5️⃣ 相关分析 (COR)
**目录**: `app/stats/cor/`
**文件数**: 3个
**核心算法**:
- `cor_stats_pearson.py` - Pearson相关系数
- `cor_stats_spearman.py` - Spearman等级相关系数
- `cor_stats_kendall.py` - Kendall相关系数

## 6️⃣ 判别分析 (CQ)
**目录**: `app/stats/cq/`
**文件数**: 1个
**核心算法**:
- `cochran_q.py` - Cochran Q检验(分类数据的方差分析)

## 7️⃣ 方差齐性检验 (HV)
**目录**: `app/stats/hv/`
**文件数**: 5个
**核心算法**:
- `hv_stats_f.py` - F检验(两组方差比较)
- `hv_stats_bartlett.py` - Bartlett检验(多组方差齐性)
- `hv_stats_levene.py` - Levene检验(稳健方差齐性)
- `hv_stats_bf.py` - Brown-Forsythe检验
- `hv_stats_hartley.py` - Hartley检验(最大方差比)

## 8️⃣ 回归分析 (REG)
**目录**: `app/stats/reg/`
**文件数**: 4个
**核心算法**:
- `reg_stats_1.py` - 一元线性回归
- `reg_stats_n.py` - 多元线性回归
- `reg_stats_log2.py` - 二元 logistic 回归
- `reg_stats_logn.py` - 多元 logistic 回归

## 9️⃣ 秩和检验 (RS)
**目录**: `app/stats/rs/`
**文件数**: 7个
**核心算法**:
- `rs_stats_single.py` - 单样本符号秩检验(Wilcoxon)
- `rs_stats_paired.py` - 配对样本符号秩检验
- `rs_stats_indep.py` - 两独立样本秩和检验(Mann-Whitney U)
- `rs_stats_kwh.py` - 多组独立样本Kruskal-Wallis H检验
- `rs_stats_fm.py` - Friedman M检验(随机区组设计)
- `rs_stats_od1.py` - 有序分类资料秩和检验1
- `rs_stats_od2.py` - 有序分类资料秩和检验2

## 🔟 游程检验 (RUNS)
**目录**: `app/stats/runs/`
**文件数**: 4个
**核心算法**:
- `runs_stats_bc.py` - 二分类变量游程检验
- `runs_stats_value1.py` - 按平均数分组的游程检验
- `runs_stats_value2.py` - 按中位数分组的游程检验
- `runs_stats_value3.py` - 按自定义值分组的游程检验

## 1️⃣1️⃣ 生存分析 (SUR)
**目录**: `app/stats/sur/`
**文件数**: 4个
**核心算法**:
- `sur_stats_km1.py` - Kaplan-Meier生存分析(原始数据)
- `sur_stats_km2.py` - Kaplan-Meier生存分析(汇总数据)
- `sur_stats_lt.py` - 寿命表法生存分析
- `sur_stats_lv.py` - 生存率比较分析(Log-rank等)

**特点**: 已重构为统一的StatsData设计理念，支持多种生存分析方法

## 1️⃣2️⃣ T检验 (T)
**目录**: `app/stats/t/`
**文件数**: 5个
**核心算法**:
- `t_stats_single1.py` - 单样本T检验(已知总体均数)
- `t_stats_single2.py` - 单样本T检验(未知总体均数)
- `t_stats_paired.py` - 配对样本T检验
- `t_stats_indep.py` - 两独立样本T检验
- `t_stats_p.py` - T检验功效分析

**特点**: 已重构为统一的StatsData设计理念

## 1️⃣3️⃣ Z检验 (Z)
**目录**: `app/stats/z/`
**文件数**: 4个
**核心算法**:
- `z_stats_single1.py` - 单样本Z检验(总体方差已知)
- `z_stats_single2.py` - 单样本Z检验(大样本近似)
- `z_stats_indep1.py` - 两独立样本Z检验(总体方差已知)
- `z_stats_indep2.py` - 两独立样本Z检验(大样本近似)

---

## ⚠️ 当前存在的不足与改进建议

### 1. 文档完整性问题
- [ ] 缺少各模块的详细技术文档
- [ ] 缺少API接口说明文档
- [ ] 缺少算法原理和使用指南
- [ ] 缺少性能基准测试报告

### 2. 测试覆盖率不足
- [ ] 部分模块缺少单元测试
- [ ] 边界条件测试不够充分
- [ ] 缺少集成测试和压力测试
- [ ] 缺少回归测试机制

### 3. 代码质量优化
- [ ] 部分代码重复度较高
- [ ] 异常处理机制需要完善
- [ ] 日志记录不够详细
- [ ] 缺少代码性能分析和优化

### 4. 功能完善性
- [ ] 缺少部分高级统计方法(如贝叶斯分析)
- [ ] 缺少机器学习相关算法
- [ ] 缺少时间序列分析功能
- [ ] 缺少多变量分析高级方法

### 5. 用户体验改进
- [ ] 缺少可视化图表生成功能
- [ ] 缺少交互式参数调整界面
- [ ] 缺少结果导出多样化格式
- [ ] 缺少在线帮助和示例库

### 6. 系统架构优化
- [ ] 缺少异步处理机制
- [ ] 缺少缓存策略优化
- [ ] 缺少负载均衡考虑
- [ ] 缺少监控和告警机制

### 7. 安全性增强
- [ ] 输入参数安全校验需要加强
- [ ] 缺少访问控制和权限管理
- [ ] 缺少数据加密传输
- [ ] 缺少审计日志功能

---

## 📈 未来发展规划

### 短期目标(1-3个月)
1. 完善现有模块的文档和技术说明
2. 提升测试覆盖率至80%以上
3. 优化核心算法性能
4. 增加常用统计图表生成功能

### 中期目标(3-6个月)
1. 补充缺失的高级统计方法
2. 实现异步处理和批处理功能
3. 增强系统的可扩展性和可维护性
4. 建立完整的CI/CD流程

### 长期目标(6-12个月)
1. 打造企业级统计分析平台
2. 支持大规模数据处理
3. 提供丰富的可视化和交互功能
4. 建立活跃的开发者社区

---
*最后更新: 2026年2月*
*整理人: Lingma*