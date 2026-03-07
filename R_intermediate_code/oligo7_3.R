# 依赖包
library(readr)
library(dplyr)
library(stringr)
library(readxl)

BstXI  <- "CCACCTTGTTGG"  # 酶切位点
BUF   <- "GC"  # 缓冲区
CONST5 <- "gatccgacgcgccatctctag"   # U6P 片段

# 读基因名单
genes <- read_csv("/Users/gzy2520/Desktop/library/GO7_re.csv", show_col_types = FALSE)

colnames(genes)[colnames(genes) == "symbol"] <- "Symbol"


# 2) 读 Sabatini
sg <- readxl::read_excel(
  "/Users/gzy2520/Desktop/library/human_gene-KO_sabatini.xlsx",
  sheet = "Human sgRNA library",
  skip = 1,
  col_names = TRUE
)

# 仅保留在genes里的基因；清洗 spacer；每基因同样取3条
sg2 <- sg %>%
  semi_join(genes, by = "Symbol") %>%  # 保留 genes 中存在的基因
  mutate(spacer = toupper(`sgRNA sequence`)) %>%  # 将 sgRNA 序列转换为大写
  filter(nchar(spacer) == 20, !grepl("[^ACGT]", spacer)) %>%  # 保证 spacer 为 20nt，且只包含 ACGT
  group_by(Symbol) %>%  # 按基因分组
  slice_head(n = 3) %>%  # 取3条sgRNA
  ungroup()  # 解除分组

# ==================== 修改开始 ====================

# 辅助函数，检查BstXI位点
has_BstXI <- function(s) grepl(BstXI, toupper(s), fixed = TRUE)

# 基因序号表
design <- sg2 %>%
  # 这里直接使用原始的spacer作为insert，不再判断和添加'G'
  dplyr::mutate(
    insert    = spacer, 
    bad_BstXI  = has_BstXI(paste0(CONST5, insert)),
    bad_polyT = grepl("TTTT", insert, fixed = TRUE)
  ) %>%
  dplyr::filter(!bad_BstXI, !bad_polyT) %>%
  dplyr::mutate(
    `sgRNA sequence`  = spacer,
    # final-construct现在只包含原始spacer，为后续插入barcode做准备
    `final-construct` = paste0(BUF, BstXI, CONST5, insert),
    # 因为我们不再为sgRNA加G，所以此列固定为FALSE
    add_5prime_G       = FALSE, 
    oligo_len          = nchar(`final-construct`),
    gc_percent         = round(100 * stringr::str_count(toupper(`final-construct`), "[GC]") / oligo_len, 2)
  )

# ==================== 修改结束 ====================

out <- design %>%
  mutate(index = match(Symbol, unique(Symbol))) %>%
  ungroup() %>%
  dplyr::select(
    index,
    Symbol,
    `sgRNA ID`,
    Chromosome,
    `sgRNA location`,
    `Genomic strand targeted`,
    `sgRNA sequence`,
    `final-construct`,
    add_5prime_G,
    oligo_len,
    gc_percent,
  )

# 导出
outfile <- "/Users/gzy2520/Desktop/library/oligo7_3.csv"
readr::write_csv(out, outfile)

message("脚本1完成：已生成不含额外'G'的sgRNA constructs。")