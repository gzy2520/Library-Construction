# 依赖包
library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(writexl)

# 常量（与你脚本一致）
BstXI   <- "CCACCTTGTTGG"          # BstXI
BUF    <- "GC"                    # 缓冲
CONST5 <- "gatccgacgcgccatctctag" # U6 前引导片段（去掉 XbaI 后）

# 路径
p_control <- "/Users/gzy2520/Desktop/library/control.xlsx"                 # 表1
p_annot   <- "/Users/gzy2520/Desktop/library/human_gene-KO_sabatini.xlsx"  # 表2
p_out     <- "/Users/gzy2520/Desktop/library/control_BstXI.xlsx"           # 输出（改为control专用）

# 工具函数（与之前一致）
add5G    <- function(s) ifelse(substr(s, 1, 1) == "G", "", "G")
has_BstXI <- function(s) grepl(BstXI, toupper(s), fixed = TRUE)

# 2) 读 Sabatini 
df2 <- readxl::read_excel(p_annot, skip = 1) %>%
  dplyr::rename_with(~stringr::str_squish(.x)) %>%
  dplyr::mutate(`sgRNA sequence` = stringr::str_to_upper(stringr::str_squish(`sgRNA sequence`))) %>%
  dplyr::distinct(`sgRNA sequence`, .keep_all = TRUE)

# 1) 读 control
df1_raw <- readxl::read_excel(p_control, col_names = T)

control_guides <- df1_raw %>%
  dplyr::transmute(`sgRNA sequence` = stringr::str_to_upper(stringr::str_squish(`sgRNA sequence`))) %>%
  dplyr::filter(!is.na(`sgRNA sequence`) & `sgRNA sequence` != "")

# 用序列去注释表拿信息
ctrl_anno <- control_guides %>%
  dplyr::left_join(df2, by = "sgRNA sequence")

# 设计 final-construct
ctrl_design <- ctrl_anno %>%
  dplyr::mutate(spacer = toupper(`sgRNA sequence`)) %>%
  dplyr::filter(nchar(spacer) == 20, !grepl("[^ACGT]", spacer)) %>%
  dplyr::mutate(
    u6_lead         = add5G(spacer),
    insert          = paste0(u6_lead, spacer),
    bad_BstXI        = has_BstXI(paste0(CONST5, insert)),
    bad_polyT       = grepl("TTTT", insert, fixed = TRUE)
  ) %>%
  dplyr::filter(!bad_BstXI, !bad_polyT) %>%
  dplyr::mutate(
    `final-construct` = paste0(BUF, BstXI, CONST5, insert),
    add_5prime_G      = u6_lead != "",
    add_5prime        = add_5prime_G,                               # 同步一列方便兼容
    oligo_len         = nchar(`final-construct`),
    gc_percent        = round(100 * stringr::str_count(toupper(`final-construct`), "[GC]") / oligo_len, 2)
  ) %>%
  # 保留关心的列
  dplyr::select(any_of(c(
    "sgRNA sequence","final-construct","add_5prime_G","add_5prime","oligo_len","gc_percent"
  )))

# 导出control表
writexl::write_xlsx(ctrl_design, p_out)

