# lib1 - 核相关基因sgRNA文库构建

## 项目简介
本项目针对核相关组分基因构建sgRNA文库（lib1），核心目标是合成总计6000条sgRNA（含目标基因sgRNA + 100条负对照sgRNA），最终需满足目标基因的unique gene数约1950~2000个（与对照合计凑够1500/2000个基因维度），支撑后续功能实验。

## 核心数据处理流程
### 1. GO术语筛选与基因初筛
#### 筛选原则
- 聚焦**核相关组分类GO术语**，排除process类术语（如maintenance of protein location in nucleus）、原生生物相关术语（macronucleus）、含亚基的复合物术语（nuclear protein-containing complex）；
- 初期限定evidence为`ECO:0000352`（manual assertion），因基因数量不足取消该限制；
- 补充GO:0006974术语扩充基因池，初筛去重后得到1680个unique gene。

#### 候选术语集（两组核心方案）
| 组别 | 包含术语 | 问题说明 |
|------|----------|----------|
| 组1 | nuclear chromosome、nuclear envelop、protein-containing complex、condensed nuclear chromosome、nuclear membrane、nuclear microtubule | complex术语可能影响后续分析 |
| 组2 | nuclear chromosome、nuclear envelop、nuclear lumen、condensed nuclear chromosome、nuclear membrane、nuclear microtubule | lumen导致基因数从2000暴涨至37000，数量过多 |

### 2. sgRNA筛选与序列设计
- **来源与规则**：从Sabatini参考文件获取sgRNA，每个基因筛选**最多4条**（或按总数调整为3条）；
- **序列规则**：前缀固定为 `GC + PacI(TTAATTAA) + U6 恒定片段`，剔除内部含PacI位点的序列；
- **筛选策略**：先映射10条sgRNA → 过滤低质量序列 → 截断为3条最终保留。

### 3. 对照设置
- **正对照**：Lig4、ERCC6L2等基因需包含在最终基因列表中；
- **负对照**：100条sgRNA（非基因维度），用于补齐6000条合成总数（如目标基因sgRNA为5910条时，随机选90条负对照补足）。

### 4. 基因ID映射（关键修正）
- 初期Symbol映射存在缺失（如PAXX、Shld1/2/3），改用**EntrezID/Ensembl ID**（物种内唯一标识）；
- 借助`org.Hs.eg.db`包实现ID转换，解决包兼容报错问题；
- 手动处理未映射Symbol：A8MVJ9、HERVK_113、HERV-K104、L1RE1无匹配结果，仅6个Symbol未完成映射（不影响最终结果）。

### 5. 最终数据校准
- 确保目标基因unique gene数接近2000个，按“基因数×3/4条sgRNA”计算基础条数，结合100条负对照补足至6000条；
- 验证所有基因在Sabatini文件中的sgRNA映射情况，未匹配基因通过ID手动校准。

## 代码使用说明
### 环境依赖
```R
# 核心依赖包
install.packages(c("org.Hs.eg.db", "dplyr", "tidyr"))
```

### 运行步骤
1. 修改代码开头的文件路径/文件名（指向目标GO表、Sabatini参考文件）；
2. 运行代码，查看控制台输出的过程量（如未映射Symbol、sgRNA筛选数）；
3. 手动校准未映射的Symbol（如需）；
4. 合并负对照sgRNA，输出最终6000条sgRNA列表。

### 代码特性
- 整合GO筛选、ID映射、sgRNA筛选、对照合并全流程；
- 内置注释与运行时过程量检验；
- 支持修改输入路径/文件，复用流程处理不同数据集。

## 关键注意事项
1. **ID唯一性**：优先用EntrezID/Ensembl ID映射，避免Symbol命名变更导致匹配失败；
2. **序列规则**：严格剔除内部含PacI位点的序列，确保前缀符合设计要求；
3. **数量校准**：最终需凑够6000条sgRNA，其中100条为负对照（非基因维度）；
4. **对照基因**：正对照（Lig4、ERCC6L2）必须包含在最终基因列表中。

## 已解决问题
1. 部分基因（Shld3等）Symbol无法映射，改用EntrezID后解决；
2. `org.Hs.eg.db`包兼容报错问题处理；
3. 手动映射无检索结果的Symbol，不影响最终结果；
4. 代码整合与注释补充，提升流程可复用性与健壮性。

## 后续优化方向
1. 分析两组GO术语集（含complex/含lumen）的差异，确定最优术语组合；
2. 进一步验证sgRNA序列的特异性，降低脱靶风险；
3. 自动化处理未映射Symbol，减少手动校准成本。

## 最终状态
文库构建订单已发出，全流程数据处理无核心问题，可支撑后续实验开展。
