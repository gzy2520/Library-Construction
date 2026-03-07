library(readr)
library(dplyr)
library(stringr)
library(readxl)

rm(list = ls())
set.seed(25) 

# 输入文件路径
# GO:0140513    nuclear protein-containing complex
# GO:0000794    condensed nuclear chromosome
# GO:0000228    nuclear chromosome
# GO:0006974    DNA damage response
# GO:0051457    maintenance of protein location in nucleus
path_go8 <- "/Users/gzy2520/Desktop/library/GO8.tsv"
# GO:0031981    nuclear lumen
# GO:0051457    maintenance of protein location in nucleus
# GO:0006974    DNA damage response
# GO:0005880    nuclear microtubule
# GO:0031965    nuclear membrane
# GO:0000228    nuclear chromosome
# GO:0005635    nuclear envelope
# GO:0000794    condensed nuclear chromosome
# GO:0140513    nuclear protein-containing complex
path_go9 <- "/Users/gzy2520/Desktop/library/GO9.tsv"
path_control <- "/Users/gzy2520/Desktop/library/control_raw.xlsx" 
path_sabatini <- "/Users/gzy2520/Desktop/library/human_gene-KO_sabatini.xlsx"

# 输出文件路径
# 中间文件
path_go8_constructs_raw <- "/Users/gzy2520/Desktop/library/constructs8_raw.csv"
path_go9_constructs_raw <- "/Users/gzy2520/Desktop/library/constructs9_raw.csv"
path_control_constructs_raw <- "/Users/gzy2520/Desktop/library/control_constructs_raw.csv"
path_go8_constructs_clean <- "/Users/gzy2520/Desktop/library/constructs8_clean.csv"
path_go9_constructs_clean <- "/Users/gzy2520/Desktop/library/constructs9_clean.csv"
path_control_constructs_clean <- "/Users/gzy2520/Desktop/library/control_constructs_clean.csv"
# 最终文库文件
path_final_library <- "/Users/gzy2520/Desktop/library/Final_final.csv"


# 构建模块序列 
prefix <- "gtaccgggcccgcTCTAGA"
BstXI <- "CCACCTTGTTGG" 
suffix <- "GTTTAAGAGCTAAGCTGGA"

# 文库组合参数 
# 注意：这里的exp_count应该是3的倍数，因为我们按基因（每个基因3条sgRNA）进行抽样
exp_count <- 5910 
ctrl_count <- 90
total_count <- exp_count + ctrl_count


# 辅助函数
# GO文件处理函数 
process_go_file <- function(infile, outfile) {
  message("处理 GO 数据: ", infile)
  raw <- read_tsv(infile, show_col_types = FALSE, comment = "")
  flt <- raw %>% transmute(symbol = .data[['SYMBOL']], go_id = .data[['GO TERM']], go_name = .data[['GO NAME']])
  genes_uniq <- flt %>%
    filter(!is.na(symbol)) %>% 
    group_by(symbol) %>%
    summarise(go_ids = paste(sort(unique(go_id)), collapse=";"), go_names = paste(sort(unique(go_name)), collapse=";"), .groups="drop") %>% 
    arrange(symbol)
  write_csv(genes_uniq, outfile)
  message("完成: ", outfile)
  return(genes_uniq)
}

# 从Sabatini库为目标基因创建constructs的函数
# 现在初始会为每个基因选择10条sgRNA
create_go_constructs <- function(input_df, sabatini_df, n_per_gene = 10) {
  sg_matched <- sabatini_df %>%
    semi_join(input_df, by = c("Symbol" = "symbol")) %>%
    group_by(Symbol) %>%
    slice_head(n = n_per_gene) %>%
    ungroup()
  
  constructs_df <- sg_matched %>%
    mutate(
      barcode = str_sub(`sgRNA sequence`, 1, -3),
      `final-construct` = paste0(prefix, barcode, BstXI, `sgRNA sequence`, suffix)
    ) %>%
    # 初始过滤，去除TTTT终止信号
    filter(!str_detect(`sgRNA sequence`, "TTTT")) %>%
    select(Symbol, `sgRNA ID`, `sgRNA sequence`, barcode, `final-construct`)
  
  return(constructs_df)
}

# 质控、过滤并截断至每个基因3条sgRNA的函数
qc_filter_and_trim <- function(df, group_name) {
  message(sprintf("开始对【%s】组进行质控过滤", group_name))
  
  # 酶切位点质控 
  XBA_PATTERN <- "TCTAGA"
  BSTXI_PATTERN <- "CCA.{6}TGG"
  BLPI_PATTERN <- "GCT.AGC"
  
  validation_results <- df %>%
    mutate(
      construct_upper = toupper(`final-construct`),
      xba_count = str_count(construct_upper, XBA_PATTERN),
      bstxi_count = str_count(construct_upper, BSTXI_PATTERN),
      blpi_count = str_count(construct_upper, BLPI_PATTERN)
    )
  
  good_constructs <- validation_results %>%
    filter(xba_count == 1 & bstxi_count == 1 & blpi_count == 1)
  
  bad_constructs <- validation_results %>%
    filter(xba_count != 1 | bstxi_count != 1 | blpi_count != 1)
  
  message(sprintf("初始数量: %d 条。质控后剩余: %d 条合格。", nrow(df), nrow(good_constructs)))
  if (nrow(bad_constructs) > 0) {
    warning(sprintf("在【%s】组中发现 %d 条不合格序列，已被剔除。", group_name, nrow(bad_constructs)))
  }
  
  # 截断并筛选保证每个基因有3条sgRNA 
  message(sprintf("对【%s】组进行截断，保证每个基因有3条sgRNA", group_name))
  
  final_sgrnas <- good_constructs %>%
    group_by(Symbol) %>%
    slice_head(n = 3) %>% # 每个基因最多取3条
    add_count(Symbol, name = "sgrna_count") %>% # 计算每个基因剩余的sgRNA数量
    filter(sgrna_count == 3) %>% # 只保留那些正好有3条合格sgRNA的基因
    ungroup() %>%
    select(-construct_upper, -xba_count, -bstxi_count, -blpi_count, -sgrna_count) # 移除临时列
  
  n_symbols_final = n_distinct(final_sgrnas$Symbol)
  message(sprintf("【%s】组最终保留 %d 个基因 (共 %d 条sgRNA)，这些基因均有3条合格的sgRNA。", 
                  group_name, n_symbols_final, nrow(final_sgrnas)))
  
  return(final_sgrnas)
}

# Control组的质控函数
validate_control_constructs <- function(df) {
  message("开始对【Control】组进行质控过滤")
  XBA_PATTERN <- "TCTAGA"
  BSTXI_PATTERN <- "CCA.{6}TGG"
  BLPI_PATTERN <- "GCT.AGC"
  
  validation_results <- df %>%
    mutate(
      construct_upper = toupper(`final-construct`),
      xba_count = str_count(construct_upper, XBA_PATTERN),
      bstxi_count = str_count(construct_upper, BSTXI_PATTERN),
      blpi_count = str_count(construct_upper, BLPI_PATTERN)
    )
  
  good_constructs <- validation_results %>%
    filter(xba_count == 1 & bstxi_count == 1 & blpi_count == 1) %>%
    select(-construct_upper, -xba_count, -bstxi_count, -blpi_count)
  
  message(sprintf("Control组初始数量: %d 条。质控后剩余: %d 条合格。", nrow(df), nrow(good_constructs)))
  return(good_constructs)
}


# ============================================================================
# 执行流程
# ============================================================================

# 1. 生成所有原始 constructs 
message("生成所有原始 constructs (每个基因10条)")
go8_re <- process_go_file(path_go8, "/Users/gzy2520/Desktop/library/GO8_re.csv")
go9_re <- process_go_file(path_go9, "/Users/gzy2520/Desktop/library/GO9_re.csv")

sabatini_lib <- read_excel(path_sabatini, sheet = "Human sgRNA library", skip = 1) %>%
  mutate(`sgRNA sequence` = toupper(str_squish(`sgRNA sequence`)), Symbol = str_squish(Symbol)) %>%
  filter(nchar(`sgRNA sequence`) == 20, !is.na(Symbol), !grepl("[^ACGT]", `sgRNA sequence`)) %>%
  distinct(`sgRNA ID`, .keep_all = TRUE)

go8_constructs_raw <- create_go_constructs(go8_re, sabatini_lib, n_per_gene = 10)
go9_constructs_raw <- create_go_constructs(go9_re, sabatini_lib, n_per_gene = 10)
write_csv(go8_constructs_raw, path_go8_constructs_raw)
write_csv(go9_constructs_raw, path_go9_constructs_raw)

control_list <- read_excel(path_control, col_names = FALSE) %>%
  rename(`sgRNA sequence` = last_col()) %>%
  mutate(`sgRNA sequence` = toupper(str_squish(`sgRNA sequence`))) %>%
  filter(!is.na(`sgRNA sequence`), nchar(`sgRNA sequence`) == 20, !grepl("[^ACGT]", `sgRNA sequence`)) %>%
  distinct() 
control_constructs_raw <- control_list %>%
  mutate(Symbol = "Control", `sgRNA ID` = paste0("CTRL_", row_number()), barcode = str_sub(`sgRNA sequence`, 1, -3), `final-construct` = paste0(prefix, barcode, BstXI, `sgRNA sequence`, suffix)) %>%
  filter(!str_detect(`sgRNA sequence`, "TTTT")) %>%
  select(Symbol, `sgRNA ID`, `sgRNA sequence`, barcode, `final-construct`)
write_csv(control_constructs_raw, path_control_constructs_raw)
message(sprintf("原始构建完成：GO8 %d 条, GO9 %d 条, Control %d 条。", nrow(go8_constructs_raw), nrow(go9_constructs_raw), nrow(control_constructs_raw)))

# 2. 对所有 constructs 进行质控、过滤和截断
message("对所有 constructs 进行质控、过滤和截断")
go8_final_pool <- qc_filter_and_trim(go8_constructs_raw, "GO8")
go9_final_pool <- qc_filter_and_trim(go9_constructs_raw, "GO9")
control_constructs_clean <- validate_control_constructs(control_constructs_raw) # Control组逻辑不变
write_csv(go8_final_pool, path_go8_constructs_clean)
write_csv(go9_final_pool, path_go9_constructs_clean)
write_csv(control_constructs_clean, path_control_constructs_clean)

# 3. 用过滤后的 constructs 按基因Symbol组合最终文库
message("按基因Symbol组合最终文库")

# 确定GO8和GO9中可用的基因列表
go8_symbols_available <- unique(go8_final_pool$Symbol)
go9_symbols_available <- unique(go9_final_pool$Symbol)
n_go8_symbols <- length(go8_symbols_available)
message(sprintf("GO8中可用于抽样的基因有 %d 个。", n_go8_symbols))

# 计算需要的目标基因数量
n_exp_symbols_needed <- exp_count / 3
message(sprintf("目标实验sgRNA数量 %d, 对应需要 %d 个基因。", exp_count, n_exp_symbols_needed))

# 选择GO8的全部基因
final_go8_symbols <- go8_symbols_available

if (length(final_go8_symbols) >= n_exp_symbols_needed) {
  warning(sprintf("GO8的合格基因数 (%d) 已达到或超过目标基因数 (%d)，将从GO8中随机抽样。", length(final_go8_symbols), n_exp_symbols_needed))
  final_go8_symbols <- sample(final_go8_symbols, n_exp_symbols_needed)
  final_go9_symbols <- character(0) # 不需要GO9
} else {
  # GO8基因数不足，需要从GO9补充
  n_go9_symbols_needed <- n_exp_symbols_needed - length(final_go8_symbols)
  message(sprintf("GO8基因数不足，需要从GO9中补充 %d 个基因。", n_go9_symbols_needed))
  
  # 从GO9中筛选出不与GO8重复的基因
  go9_symbols_unique <- setdiff(go9_symbols_available, final_go8_symbols)
  
  if (length(go9_symbols_unique) < n_go9_symbols_needed) {
    warning(sprintf("GO9中不重复的合格基因只有 %d 个，不足所需的 %d 个。将使用所有可用的。", length(go9_symbols_unique), n_go9_symbols_needed))
    final_go9_symbols <- go9_symbols_unique
  } else {
    # 从不重复的GO9基因中随机抽样
    final_go9_symbols <- sample(go9_symbols_unique, n_go9_symbols_needed)
  }
}

# 根据选定的基因列表，从sgRNA池中提取所有对应的sgRNA
final_go8_sgrnas <- go8_final_pool %>% filter(Symbol %in% final_go8_symbols)
final_go9_sgrnas <- go9_final_pool %>% filter(Symbol %in% final_go9_symbols)
message(sprintf("最终从GO8选择了 %d 个基因 (%d 条sgRNA)。", n_distinct(final_go8_sgrnas$Symbol), nrow(final_go8_sgrnas)))
message(sprintf("最终从GO9选择了 %d 个基因 (%d 条sgRNA)。", n_distinct(final_go9_sgrnas$Symbol), nrow(final_go9_sgrnas)))

# 从合格的Control中抽样
n_control_available <- nrow(control_constructs_clean)
if (n_control_available < ctrl_count){
  warning(sprintf("合格的Control组只有 %d 条，不足 %d 条，将使用所有可用条目。", n_control_available, ctrl_count))
  final_control <- control_constructs_clean
} else {
  final_control <- control_constructs_clean %>% slice_sample(n = ctrl_count)
}
message(sprintf("从合格的Control组中选择了 %d 条。", nrow(final_control)))

# 合并所有部分
final_library <- bind_rows(
  mutate(final_go8_sgrnas, Group = "GO8"),
  mutate(final_go9_sgrnas, Group = "GO9_unique"),
  mutate(final_control, Group = "Control")
)

# 4. 保存最终文库
message("保存最终文库")
final_total <- nrow(final_library)
message(sprintf("最终文库组装完成，总计 %d 条。", final_total))
if (final_total < (n_exp_symbols_needed*3 + ctrl_count)) {
  warning(sprintf("最终文库大小为 %d, 未达到目标 %d。", final_total, (n_exp_symbols_needed*3 + ctrl_count)))
}

write_csv(final_library, path_final_library)
message("最终文件已保存到: ", path_final_library)
print(table(final_library$Group))



message("生成最终选择的 Symbol 列表")

# 1. 创建一个包含所有最终选择的 symbol 及其来源（GO8或GO9）的数据框
df_go8_selected <- tibble(symbol = final_go8_symbols, source_group = "GO8")
df_go9_selected <- tibble(symbol = final_go9_symbols, source_group = "GO9_unique")
all_selected_symbols_df <- bind_rows(df_go8_selected, df_go9_selected)

# 2. 合并 go8_re 和 go9_re 以创建一个全面的GO注释查找表
all_go_annotations <- bind_rows(go8_re, go9_re) %>%
  distinct(symbol, .keep_all = TRUE)

# 3. 将选择的 symbol 列表与GO注释表进行左连接，以获取每个symbol的GO信息
final_symbol_list_with_info <- all_selected_symbols_df %>%
  left_join(all_go_annotations, by = "symbol")

# 4. 保存这个新的数据框到CSV文件
write_csv(final_symbol_list_with_info, "/Users/gzy2520/Desktop/library/final_symbol_list.csv")
message(sprintf("最终选择的 %d 个基因及其GO信息已保存到: /Users/gzy2520/Desktop/library/final_symbol_list.csv", 
                nrow(final_symbol_list_with_info)))
