library(readr)
library(dplyr)
library(stringr)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)

rm(list = ls())
set.seed(25) 

# 输入文件路径
# GO:0051457    maintenance of protein location in nucleus
# GO:0006974    DNA damage response
# GO:0005635    nuclear envelope
# GO:0140513    nuclear protein-containing complex
path_go10 <- "/Users/gzy2520/Desktop/library/GO10.tsv"
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
path_go10_constructs_raw <- "/Users/gzy2520/Desktop/library/constructs10_raw.csv"
path_go9_constructs_raw <- "/Users/gzy2520/Desktop/library/constructs9_raw.csv"
path_control_constructs_raw <- "/Users/gzy2520/Desktop/library/control_constructs_raw.csv"
path_go10_constructs_clean <- "/Users/gzy2520/Desktop/library/constructs10_clean.csv"
path_go9_constructs_clean <- "/Users/gzy2520/Desktop/library/constructs9_clean.csv"
path_control_constructs_clean <- "/Users/gzy2520/Desktop/library/control_constructs_clean.csv"
path_final_library <- "/Users/gzy2520/Desktop/library/Final_final.csv"

# 模块序列 
prefix <- "gtaccgggcccgcTCTAGA"
BstXI <- "CCACCTTGTTGG" 
suffix <- "GTTTAAGAGCTAAGCTGGA"

# 文库组合参数 
exp_count <- 5910 
ctrl_count <- 90
total_count <- exp_count + ctrl_count


# ============================================================================
# 辅助函数
# ============================================================================

# 手动映射表（基于NCBI Gene数据库查询）
get_manual_mappings <- function() {
  manual_map <- tribble(
    ~symbol,          ~gene_id,    
    # 手动映射
    "CRIPAK",         "285464",    
    "FAM231A",        "653203",     
    "KIAA1045",       "79789",     
    "KIAA1804",       "168850",    
    "PRAMEF16",       "653277",
    "PRAMEF3",        "441871",    
    "SMCR9",          "152078",    
    "LOC100127983",   "100127983",
    "LOC100130480",   "100130480",
    "LOC100133267",   "100133267",
    "LOC100287177",   "100287177",
    "LOC100505679",   "100505679", 
    "LOC100507003",   "100507003",
    "LOC100653515",   "100653515",
    "LOC100996485",   "100996485", 
    "LOC147646",      "147646", 
    "LOC643669",      "643669", 
    "LOC81691",       "81691", 
    "SMCR9",          "101048113",
    "SIX6OS1",        "100861532", 
    "CUSTOS",         "100507165",
    "ZHX1-C8ORF76",  "100533106"
    # A8MVJ9, HERVK_113, HERV-K104, L1RE1 人工搜索也未找到
  )
  
  return(manual_map)
}

# Symbol到Entrez ID映射函数（包含手动映射）
map_symbol_to_id <- function(symbols) {
  message("映射 symbol 到 Entrez ID...")
  
  # 使用org.Hs.eg.db进行映射
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # 处理可能存在的别名
  na_symbols <- names(entrez_ids)[is.na(entrez_ids)]
  if (length(na_symbols) > 0) {
    message(sprintf("  %d 个symbol未找到ID，尝试ALIAS映射...", length(na_symbols)))
    
    # 使用tryCatch处理可能的错误
    alias_ids <- tryCatch({
      mapIds(org.Hs.eg.db,
             keys = na_symbols,
             column = "ENTREZID",
             keytype = "ALIAS",
             multiVals = "first")
    }, error = function(e) {
      message("逐个尝试...")
      result <- setNames(rep(NA_character_, length(na_symbols)), na_symbols)
      for (sym in na_symbols) {
        tryCatch({
          id <- mapIds(org.Hs.eg.db,
                       keys = sym,
                       column = "ENTREZID",
                       keytype = "ALIAS",
                       multiVals = "first")
          if (!is.na(id)) {
            result[sym] <- id
          }
        }, error = function(e2) {
          # 单个symbol映射失败，保持NA
        })
      }
      return(result)
    })
    
    # 更新entrez_ids
    entrez_ids[na_symbols] <- alias_ids
  }
  
  # 创建映射表
  symbol_to_id_df <- tibble(
    symbol = names(entrez_ids),
    gene_id = as.character(entrez_ids)
  )
  
  # 应用手动映射 
  manual_map <- get_manual_mappings()
  
  # 找出仍未映射且在手动映射表中的symbol
  still_unmapped <- symbol_to_id_df %>%
    filter(is.na(gene_id), symbol %in% manual_map$symbol)
  
  if (nrow(still_unmapped) > 0) {
    message(sprintf("  应用手动映射补充 %d 个symbol...", nrow(still_unmapped)))
    
    # 更新映射
    for (i in 1:nrow(still_unmapped)) {
      sym <- still_unmapped$symbol[i]
      manual_id <- manual_map$gene_id[manual_map$symbol == sym][1]
      symbol_to_id_df$gene_id[symbol_to_id_df$symbol == sym] <- manual_id
      message(sprintf("    %s -> %s ", 
                      sym, manual_id))
    }
  }
  
  n_mapped <- sum(!is.na(symbol_to_id_df$gene_id))
  n_unmapped <- sum(is.na(symbol_to_id_df$gene_id))
  message(sprintf("  最终映射结果: 成功 %d, 未映射 %d", n_mapped, n_unmapped))
  
  # 显示仍未映射的symbol
  unmapped_symbols <- symbol_to_id_df$symbol[is.na(symbol_to_id_df$gene_id)]
  message(sprintf("  仍未映射的symbols: %s", paste(unmapped_symbols, collapse=", ")))
  
  return(symbol_to_id_df)
}

# 根据gene_id获取官方标准symbol
get_official_symbol <- function(gene_ids) {
  message("根据gene_id获取官方标准symbol...")
  
  # 获取唯一的gene_id
  unique_gene_ids <- unique(gene_ids)
  
  # 过滤掉特殊值（如"Control"）和NA
  real_gene_ids <- unique_gene_ids[unique_gene_ids != "Control" & !is.na(unique_gene_ids)]
  
  if (length(real_gene_ids) == 0) {
    return(tibble(gene_id = unique_gene_ids, official_symbol = unique_gene_ids))
  }
  
  # 使用mapIds获取官方symbol
  official_symbols <- mapIds(org.Hs.eg.db,
                             keys = real_gene_ids,
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")
  
  # 创建完整的映射表（包括Control等特殊值）
  id_to_symbol_df <- tibble(gene_id = unique_gene_ids) %>%
    mutate(
      official_symbol = case_when(
        gene_id == "Control" ~ "Control",
        is.na(gene_id) ~ NA_character_,
        gene_id %in% names(official_symbols) ~ as.character(official_symbols[gene_id]),
        TRUE ~ gene_id  # 如果找不到，使用gene_id本身
      )
    )
  
  n_real_ids <- length(real_gene_ids)
  n_mapped <- sum(!is.na(official_symbols))
  message(sprintf("  成功获取官方symbol: %d/%d", n_mapped, n_real_ids))
  
  return(id_to_symbol_df)
}

# GO文件处理函数（改为基于ID）
process_go_file <- function(infile, outfile) {
  message("处理 GO 数据: ", infile)
  raw <- read_tsv(infile, show_col_types = FALSE, comment = "")
  
  # 提取symbol, go_id, go_name
  flt <- raw %>% 
    transmute(
      symbol = .data[['SYMBOL']], 
      go_id = .data[['GO TERM']], 
      go_name = .data[['GO NAME']]
    ) %>%
    filter(!is.na(symbol))
  
  # 映射symbol到gene_id
  symbols <- unique(flt$symbol)
  symbol_to_id <- map_symbol_to_id(symbols)
  
  # 合并ID信息
  flt_with_id <- flt %>%
    left_join(symbol_to_id, by = "symbol") %>%
    filter(!is.na(gene_id))  # 只保留能映射到ID的
  
  message(sprintf("  原始symbol数: %d, 映射后保留: %d", 
                  length(symbols), n_distinct(flt_with_id$gene_id)))
  
  # 按gene_id进行去重和分组
  genes_uniq <- flt_with_id %>%
    group_by(gene_id) %>%
    summarise(
      symbol = first(symbol),  # 保留第一个symbol作为代表
      all_symbols = paste(unique(symbol), collapse=";"),  # 保留所有symbol（处理别名）
      go_ids = paste(sort(unique(go_id)), collapse=";"), 
      go_names = paste(sort(unique(go_name)), collapse=";"),
      .groups="drop"
    ) %>% 
    arrange(gene_id)
  
  write_csv(genes_uniq, outfile)
  message("完成: ", outfile, " (共 ", nrow(genes_uniq), " 个唯一基因ID)")
  return(genes_uniq)
}

# 从Sabatini库为目标基因创建constructs的函数（改为基于ID）
create_go_constructs <- function(input_df, sabatini_df, n_per_gene = 10) {
  # 通过gene_id进行匹配
  sg_matched <- sabatini_df %>%
    semi_join(input_df, by = "gene_id") %>%
    group_by(gene_id) %>%
    slice_head(n = n_per_gene) %>%
    ungroup()
  
  constructs_df <- sg_matched %>%
    mutate(
      barcode = str_sub(`sgRNA sequence`, 1, -3),
      `final-construct` = paste0(prefix, barcode, BstXI, `sgRNA sequence`, suffix)
    ) %>%
    # 初始过滤，去除TTTT终止信号
    filter(!str_detect(`sgRNA sequence`, "TTTT")) %>%
    select(gene_id, Symbol, `sgRNA ID`, `sgRNA sequence`, barcode, `final-construct`)
  
  return(constructs_df)
}

# 质控、过滤并截断至每个基因3条sgRNA的函数（改为基于ID）
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
  
  # 按gene_id截断并筛选
  message(sprintf("对【%s】组进行截断，保证每个基因有3条sgRNA", group_name))
  
  final_sgrnas <- good_constructs %>%
    group_by(gene_id) %>%
    slice_head(n = 3) %>%
    add_count(gene_id, name = "sgrna_count") %>%
    filter(sgrna_count == 3) %>%
    ungroup() %>%
    select(-construct_upper, -xba_count, -bstxi_count, -blpi_count, -sgrna_count)
  
  n_genes_final = n_distinct(final_sgrnas$gene_id)
  message(sprintf("【%s】组最终保留 %d 个基因 (共 %d 条sgRNA)，这些基因均有3条合格的sgRNA。", 
                  group_name, n_genes_final, nrow(final_sgrnas)))
  
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

# 1. 处理Sabatini库，添加gene_id映射
message("==== 处理Sabatini库 ====")
sabatini_lib_raw <- read_excel(path_sabatini, sheet = "Human sgRNA library", skip = 1) %>%
  mutate(`sgRNA sequence` = toupper(str_squish(`sgRNA sequence`)), 
         Symbol = str_squish(Symbol)) %>%
  filter(nchar(`sgRNA sequence`) == 20, 
         !is.na(Symbol), 
         !grepl("[^ACGT]", `sgRNA sequence`)) %>%
  distinct(`sgRNA ID`, .keep_all = TRUE)

# 映射Sabatini库的symbol到gene_id
sabatini_symbols <- unique(sabatini_lib_raw$Symbol)
sabatini_symbol_to_id <- map_symbol_to_id(sabatini_symbols)

sabatini_lib <- sabatini_lib_raw %>%
  left_join(sabatini_symbol_to_id, by = c("Symbol" = "symbol")) %>%
  filter(!is.na(gene_id))  # 只保留能映射到ID的

message(sprintf("Sabatini库处理结果: 原始 %d 条, 映射后保留 %d 条 (涉及 %d 个基因ID)", 
                nrow(sabatini_lib_raw), nrow(sabatini_lib), n_distinct(sabatini_lib$gene_id)))

# 2. 生成所有原始 constructs 
message("\n==== 生成所有原始 constructs (每个基因10条) ====")
go10_re <- process_go_file(path_go10, "/Users/gzy2520/Desktop/library/GO10_re.csv")
go9_re <- process_go_file(path_go9, "/Users/gzy2520/Desktop/library/GO9_re.csv")

go10_constructs_raw <- create_go_constructs(go10_re, sabatini_lib, n_per_gene = 10)
go9_constructs_raw <- create_go_constructs(go9_re, sabatini_lib, n_per_gene = 10)
write_csv(go10_constructs_raw, path_go10_constructs_raw)
write_csv(go9_constructs_raw, path_go9_constructs_raw)

# Control组处理
control_list <- read_excel(path_control, col_names = FALSE) %>%
  rename(`sgRNA sequence` = last_col()) %>%
  mutate(`sgRNA sequence` = toupper(str_squish(`sgRNA sequence`))) %>%
  filter(!is.na(`sgRNA sequence`), 
         nchar(`sgRNA sequence`) == 20, 
         !grepl("[^ACGT]", `sgRNA sequence`)) %>%
  distinct() 

control_constructs_raw <- control_list %>%
  mutate(
    gene_id = "Control",
    Symbol = "Control", 
    `sgRNA ID` = paste0("CTRL_", row_number()), 
    barcode = str_sub(`sgRNA sequence`, 1, -3), 
    `final-construct` = paste0(prefix, barcode, BstXI, `sgRNA sequence`, suffix)
  ) %>%
  filter(!str_detect(`sgRNA sequence`, "TTTT")) %>%
  select(gene_id, Symbol, `sgRNA ID`, `sgRNA sequence`, barcode, `final-construct`)

write_csv(control_constructs_raw, path_control_constructs_raw)
message(sprintf("原始构建完成：GO10 %d 条 (%d 个基因), GO9 %d 条 (%d 个基因), Control %d 条。", 
                nrow(go10_constructs_raw), n_distinct(go10_constructs_raw$gene_id),
                nrow(go9_constructs_raw), n_distinct(go9_constructs_raw$gene_id),
                nrow(control_constructs_raw)))

# 3. 对所有 constructs 进行质控、过滤和截断
message("\n==== 对所有 constructs 进行质控、过滤和截断 ====")
go10_final_pool <- qc_filter_and_trim(go10_constructs_raw, "GO10")
go9_final_pool <- qc_filter_and_trim(go9_constructs_raw, "GO9")
control_constructs_clean <- validate_control_constructs(control_constructs_raw)
write_csv(go10_final_pool, path_go10_constructs_clean)
write_csv(go9_final_pool, path_go9_constructs_clean)
write_csv(control_constructs_clean, path_control_constructs_clean)

# 4. 用过滤后的 constructs 按基因ID组合最终文库
message("\n==== 按基因ID组合最终文库 ====")

# 确定GO10和GO9中可用的基因ID列表
go10_ids_available <- unique(go10_final_pool$gene_id)
go9_ids_available <- unique(go9_final_pool$gene_id)
n_go10_ids <- length(go10_ids_available)
message(sprintf("GO10中可用于抽样的基因ID有 %d 个。", n_go10_ids))

# 计算需要的目标基因数量
n_exp_ids_needed <- exp_count / 3
message(sprintf("目标实验sgRNA数量 %d, 对应需要 %d 个基因。", exp_count, n_exp_ids_needed))

# 选择GO10的全部基因ID
final_go10_ids <- go10_ids_available

if (length(final_go10_ids) >= n_exp_ids_needed) {
  warning(sprintf("GO10的合格基因数 (%d) 已达到或超过目标基因数 (%d)，将从GO10中随机抽样。", 
                  length(final_go10_ids), n_exp_ids_needed))
  final_go10_ids <- sample(final_go10_ids, n_exp_ids_needed)
  final_go9_ids <- character(0)
} else {
  # GO10基因数不足，需要从GO9补充
  n_go9_ids_needed <- n_exp_ids_needed - length(final_go10_ids)
  message(sprintf("GO10基因数不足，需要从GO9中补充 %d 个基因。", n_go9_ids_needed))
  
  # 从GO9中筛选出不与GO10重复的基因ID
  go9_ids_unique <- setdiff(go9_ids_available, final_go10_ids)
  
  if (length(go9_ids_unique) < n_go9_ids_needed) {
    warning(sprintf("GO9中不重复的合格基因只有 %d 个，不足所需的 %d 个。将使用所有可用的。", 
                    length(go9_ids_unique), n_go9_ids_needed))
    final_go9_ids <- go9_ids_unique
  } else {
    # 从不重复的GO9基因ID中随机抽样
    final_go9_ids <- sample(go9_ids_unique, n_go9_ids_needed)
  }
}

# 根据选定的基因ID列表，从sgRNA池中提取所有对应的sgRNA
final_go10_sgrnas <- go10_final_pool %>% filter(gene_id %in% final_go10_ids)
final_go9_sgrnas <- go9_final_pool %>% filter(gene_id %in% final_go9_ids)
message(sprintf("最终从GO10选择了 %d 个基因ID (%d 条sgRNA)。", 
                n_distinct(final_go10_sgrnas$gene_id), nrow(final_go10_sgrnas)))
message(sprintf("最终从GO9选择了 %d 个基因ID (%d 条sgRNA)。", 
                n_distinct(final_go9_sgrnas$gene_id), nrow(final_go9_sgrnas)))

# 从合格的Control中抽样
n_control_available <- nrow(control_constructs_clean)
if (n_control_available < ctrl_count){
  warning(sprintf("合格的Control组只有 %d 条，不足 %d 条，将使用所有可用条目。", 
                  n_control_available, ctrl_count))
  final_control <- control_constructs_clean
} else {
  final_control <- control_constructs_clean %>% slice_sample(n = ctrl_count)
}
message(sprintf("从合格的Control组中选择了 %d 条。", nrow(final_control)))

# 合并所有部分
final_library <- bind_rows(
  mutate(final_go10_sgrnas, Group = "GO10"),
  mutate(final_go9_sgrnas, Group = "GO9_unique"),
  mutate(final_control, Group = "Control")
)

# 5. 将Symbol替换为官方标准symbol
message("\n==== 将Symbol替换为官方标准symbol ====")

# 获取所有gene_id的官方symbol映射
all_gene_ids <- unique(final_library$gene_id)
id_to_official_symbol <- get_official_symbol(all_gene_ids)

# 更新final_library中的Symbol列
final_library_standardized <- final_library %>%
  select(-Symbol) %>%  # 移除原来的Symbol列
  left_join(id_to_official_symbol, by = "gene_id") %>%
  rename(Symbol = official_symbol) %>%  # 将official_symbol重命名为Symbol
  select(gene_id, Symbol, `sgRNA ID`, `sgRNA sequence`, barcode, `final-construct`, Group)  # 重新排列列顺序

message(sprintf("已将 %d 条记录的symbol替换为官方标准symbol", nrow(final_library_standardized)))

# 6. 保存最终文库（标准化后）
message("\n==== 保存最终文库 ====")
final_total <- nrow(final_library_standardized)
message(sprintf("最终文库组装完成，总计 %d 条。", final_total))
if (final_total < (n_exp_ids_needed*3 + ctrl_count)) {
  warning(sprintf("最终文库大小为 %d, 未达到目标 %d。", 
                  final_total, (n_exp_ids_needed*3 + ctrl_count)))
}

write_csv(final_library_standardized, path_final_library)
message("最终文件已保存到: ", path_final_library)
print(table(final_library_standardized$Group))

# 7. 生成最终选择的基因ID列表（标准化symbol）
message("\n==== 生成最终选择的基因ID列表 ====")

# 创建一个包含所有最终选择的gene_id及其来源的数据框
df_go10_selected <- tibble(gene_id = final_go10_ids, source_group = "GO10")
df_go9_selected <- tibble(gene_id = final_go9_ids, source_group = "GO9_unique")
all_selected_ids_df <- bind_rows(df_go10_selected, df_go9_selected)

# 合并go10_re和go9_re以创建一个全面的GO注释查找表（基于gene_id）
all_go_annotations <- bind_rows(go10_re, go9_re) %>%
  distinct(gene_id, .keep_all = TRUE)

# 将选择的gene_id列表与GO注释表进行左连接
final_geneid_list_with_info <- all_selected_ids_df %>%
  left_join(all_go_annotations, by = "gene_id")

# 标准化symbol列
final_geneid_list_with_info_standardized <- final_geneid_list_with_info %>%
  select(-symbol) %>%  # 移除原来的symbol列
  left_join(id_to_official_symbol, by = "gene_id") %>%
  rename(symbol = official_symbol) %>%  # 将official_symbol重命名为symbol
  select(gene_id, symbol, all_symbols, source_group, go_ids, go_names)  # 重新排列列顺序

# 保存
write_csv(final_geneid_list_with_info_standardized, "/Users/gzy2520/Desktop/library/final_geneid_list.csv")
message(sprintf("最终选择的 %d 个基因ID及其GO信息已保存到: /Users/gzy2520/Desktop/library/final_geneid_list.csv", 
                nrow(final_geneid_list_with_info_standardized)))

message("\n==== 全部完成！ ====")