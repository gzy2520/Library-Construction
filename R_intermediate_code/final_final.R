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
path_go8_constructs_raw <- "/Users/gzy2520/Desktop/library/constructs8.csv"
path_go9_constructs_raw <- "/Users/gzy2520/Desktop/library/constructs9.csv"
path_control_constructs_raw <- "/Users/gzy2520/Desktop/library/control_constructs.csv"
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
exp_count <- 5910
ctrl_count <- 90
total_count <- exp_count + ctrl_count


# 辅助函数

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

# 用于质控和过滤的函数
validate_and_filter_constructs <- function(df, group_name) {
  message(sprintf("开始对【%s】组进行质控过滤", group_name))
  
  # 定义检查模式
  XBA_PATTERN <- "TCTAGA"
  BSTXI_PATTERN <- "CCA.{6}TGG"
  BLPI_PATTERN <- "GCT.AGC"
  
  # 检查并添加计数列
  validation_results <- df %>%
    mutate(
      construct_upper = toupper(`final-construct`),
      xba_count = str_count(construct_upper, XBA_PATTERN),
      bstxi_count = str_count(construct_upper, BSTXI_PATTERN),
      blpi_count = str_count(construct_upper, BLPI_PATTERN)
    )
  
  # 分离合格与不合格的条目
  good_constructs <- validation_results %>%
    filter(xba_count == 1 & bstxi_count == 1 & blpi_count == 1) %>%
    select(-construct_upper, -xba_count, -bstxi_count, -blpi_count) # 移除检查用的临时列
  
  bad_constructs <- validation_results %>%
    filter(xba_count != 1 | bstxi_count != 1 | blpi_count != 1)
  
  message(sprintf("初始数量: %d 条。质控后剩余: %d 条合格。", nrow(df), nrow(good_constructs)))
  
  if (nrow(bad_constructs) > 0) {
    warning(sprintf("在【%s】组中发现 %d 条不合格序列，已被剔除。", group_name, nrow(bad_constructs)))
    message("以下是被剔除的不合格序列详情：")
    print(bad_constructs %>% select(Symbol, `sgRNA sequence`, xba_count, bstxi_count, blpi_count))
  } else {
    message("所有序列均通过质控。")
  }
  
  return(good_constructs)
}

# ============================================================================
# 执行流程
# ============================================================================

# 1 生成所有原始 constructs 
message("生成所有原始 constructs ")
go8_re <- process_go_file(path_go8, "/Users/gzy2520/Desktop/library/GO8_re.csv")
go9_re <- process_go_file(path_go9, "/Users/gzy2520/Desktop/library/GO9_re.csv")

sabatini_lib <- read_excel(path_sabatini, sheet = "Human sgRNA library", skip = 1) %>%
  mutate(`sgRNA sequence` = toupper(str_squish(`sgRNA sequence`)), Symbol = str_squish(Symbol)) %>%
  filter(nchar(`sgRNA sequence`) == 20, !is.na(Symbol), !grepl("[^ACGT]", `sgRNA sequence`)) %>%
  distinct(`sgRNA ID`, .keep_all = TRUE)

go8_constructs_raw <- create_go_constructs(go8_re, sabatini_lib)
go9_constructs_raw <- create_go_constructs(go9_re, sabatini_lib)
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

# 2 对所有原始 constructs 进行质控和过滤
message("对所有 constructs 进行质控过滤 ")
go8_constructs_clean <- validate_and_filter_constructs(go8_constructs_raw, "GO8")
go9_constructs_clean <- validate_and_filter_constructs(go9_constructs_raw, "GO9")
control_constructs_clean <- validate_and_filter_constructs(control_constructs_raw, "Control")
write_csv(go8_constructs_clean, path_go8_constructs_clean)
write_csv(go9_constructs_clean, path_go9_constructs_clean)
write_csv(control_constructs_clean, path_control_constructs_clean)

# 3 用过滤后的 constructs 组合最终文库 
message("用过滤后的 constructs 组合最终文库")
# 取 GO8 的全部合格结果
final_go8 <- go8_constructs_clean
n_go8 <- nrow(final_go8)
message(sprintf("合格的GO8共 %d 条，全部用于最终文库。", n_go8))

if (n_go8 >= exp_count) {
  warning(sprintf("合格的GO8有 %d 条,  >= 目标 %d. 将从GO8中随机抽样。", n_go8, exp_count))
  final_go8 <- final_go8 %>% slice_sample(n = exp_count)
  n_go9_needed <- 0
} else {
  n_go9_needed <- exp_count - n_go8
  message(sprintf("需要从合格的GO9中补充 %d 条。", n_go9_needed))
}

# 从合格的GO9中筛选和抽样
if (n_go9_needed > 0) {
  go8_symbols <- unique(final_go8$Symbol)
  go9_candidates <- go9_constructs_clean %>% filter(!Symbol %in% go8_symbols)
  n_go9_candidates <- nrow(go9_candidates)
  
  if (n_go9_candidates < n_go9_needed) {
    warning(sprintf("合格且唯一的GO9基因只有 %d 条，不足 %d 条。将使用所有可用条目。", n_go9_candidates, n_go9_needed))
    final_go9 <- go9_candidates
  } else {
    final_go9 <- go9_candidates %>% slice_sample(n = n_go9_needed)
  }
  message(sprintf("从合格的GO9中选择了 %d 条。", nrow(final_go9)))
} else {
  final_go9 <- tibble() 
}

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
  mutate(final_go8, Group = "GO8"),
  mutate(final_go9, Group = "GO9_unique"),
  mutate(final_control, Group = "Control")
)


# 4 保存最终文库
message("保存最终文库 ")
final_total <- nrow(final_library)
message(sprintf("最终文库组装完成，总计 %d 条。", final_total))
if (final_total != total_count && final_total < total_count) {
  warning(sprintf("最终文库大小为 %d, 未达到目标 %d。", final_total, total_count))
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
#    使用 distinct(symbol, .keep_all = TRUE) 来确保每个symbol只出现一次，防止重复
all_go_annotations <- bind_rows(go8_re, go9_re) %>%
  distinct(symbol, .keep_all = TRUE)

# 3. 将选择的 symbol 列表与GO注释表进行左连接，以获取每个symbol的GO信息
final_symbol_list_with_info <- all_selected_symbols_df %>%
  left_join(all_go_annotations, by = "symbol")

# 4. 保存这个新的数据框到CSV文件
write_csv(final_symbol_list_with_info, "/Users/gzy2520/Desktop/library/final_symbol_list.csv")
message(sprintf("最终选择的 %d 个基因及其GO信息已保存到: /Users/gzy2520/Desktop/library/final_symbol_list.csv", 
                nrow(final_symbol_list_with_info)))

