# Library 1
# CRISPR sgRNA 文库构建流程

## 项目简介
本 R 脚本用于构建针对特定基因本体（GO）术语的 CRISPR sgRNA 文库。流程包括基因 Symbol 与 Entrez ID 的映射、GO 数据处理、sgRNA 序列筛选与质控、最终文库的组装与标准化。

## 依赖环境
### R 包
- `readr`：读取 TSV/CSV 文件
- `dplyr`：数据清洗与操作
- `stringr`：字符串处理
- `readxl`：读取 Excel 文件
- `org.Hs.eg.db`：人类基因 ID 映射数据库
- `AnnotationDbi`：基因注释工具

### 安装依赖
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))
install.packages(c("readr", "dplyr", "stringr", "readxl"))
```

## 输入文件
请确保以下文件路径在脚本中正确配置：
1. **GO10.tsv**：包含目标 GO 术语（如 `GO:0051457`、`GO:0006974` 等）的基因列表
2. **GO9.tsv**：包含另一组 GO 术语（如 `GO:0031981`、`GO:0005880` 等）的基因列表
3. **control_raw.xlsx**：对照组 sgRNA 序列文件
4. **human_gene-KO_sabatini.xlsx**：Sabatini 人类 sgRNA 文库（需包含 `Human sgRNA library` 工作表）

## 输出文件
脚本运行后将生成以下文件（路径可在脚本中修改）：
1. **constructs10_raw.csv**：GO10 组原始 sgRNA 构建体
2. **constructs9_raw.csv**：GO9 组原始 sgRNA 构建体
3. **control_constructs_raw.csv**：对照组原始 sgRNA 构建体
4. **constructs10_clean.csv**：GO10 组质控后 sgRNA 构建体
5. **constructs9_clean.csv**：GO9 组质控后 sgRNA 构建体
6. **control_constructs_clean.csv**：对照组质控后 sgRNA 构建体
7. **Final_final.csv**：最终组装的标准化 sgRNA 文库
8. **final_geneid_list.csv**：最终选择的基因 ID 及其 GO 注释信息

## 使用说明
1. **修改文件路径**：打开脚本，将 `path_go10`、`path_go9`、`path_control`、`path_sabatini` 等输入/输出路径修改为本地实际路径。
2. **运行脚本**：在 R 或 RStudio 中运行完整脚本。

## 主要流程
1. **Sabatini 文库处理**：读取并清洗 Sabatini sgRNA 文库，将基因 Symbol 映射为 Entrez ID。
2. **GO 数据处理**：读取 GO10/GO9 文件，映射基因 ID，去重并合并 GO 注释。
3. **构建体生成**：为每个基因选取最多 10 条 sgRNA，添加前缀、后缀和酶切位点，生成 `final-construct` 序列。
4. **质控过滤**：
   - 检查酶切位点（XbaI、BstXI、BlpI）是否唯一
   - 去除含 `TTTT` 终止信号的序列
   - 每个基因保留 3 条合格 sgRNA
5. **文库组装**：
   - 从 GO10/GO9 中随机抽取目标数量的基因（共 5910 条实验 sgRNA）
   - 从对照组中抽取 90 条 sgRNA
   - 合并并标准化基因 Symbol
6. **结果输出**：保存最终文库和基因注释列表。

## 注意事项
1. **随机种子**：脚本设置 `set.seed(25)` 以保证结果可复现，如需调整可修改此参数。
2. **手动映射表**：`get_manual_mappings()` 函数包含部分无法通过数据库映射的基因，可根据实际情况补充。
3. **质控标准**：酶切位点和序列过滤规则可在 `qc_filter_and_trim()` 函数中调整。
4. **文件格式**：输入文件需严格遵循脚本预期格式（如 GO 文件需包含 `SYMBOL`、`GO TERM`、`GO NAME` 列）。
