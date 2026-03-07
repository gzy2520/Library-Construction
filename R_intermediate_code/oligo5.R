# 依赖包
library(readr)  
library(dplyr)  
library(stringr) 
library(readxl)  

# PacI；U6 前引导片段
PACI  <- "TTAATTAA"  # 酶切位点
BUF   <- "GC"  # 缓冲区
CONST5 <- "gatccgacgcgccatctctag"   # XbaI 去掉后留下的 U6P 片段

# 读基因名单
genes <- read_csv("/Users/gzy2520/Desktop/library/GO5_re.csv", show_col_types = FALSE) 

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

# 设计寡核苷酸：PacI + 恒定区 + U6 兼容的 spacer（首位不是 G 则加 5'G）
add5G    <- function(s) ifelse(substr(s,1,1) == "G", "", "G")
has_pacI <- function(s) grepl(PACI, toupper(s), fixed = TRUE)

# 基因序号表
design <- sg2 %>%
  dplyr::mutate(spacer = toupper(`sgRNA sequence`)) %>%
  dplyr::mutate(
    u6_lead   = add5G(spacer),
    insert    = paste0(u6_lead, spacer),
    bad_pacI  = has_pacI(paste0(CONST5, insert)),
    bad_polyT = grepl("TTTT", insert, fixed = TRUE)
  ) %>%
  dplyr::filter(!bad_pacI, !bad_polyT) %>%
  dplyr::mutate(
    `sgRNA sequence`  = spacer,
    `final-construct` = paste0(BUF, PACI, CONST5, insert),
    add_5prime_G       = u6_lead != "",
    oligo_len          = nchar(`final-construct`),
    gc_percent         = round(100 * stringr::str_count(toupper(`final-construct`), "[GC]") / oligo_len, 2)
  )

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
outfile <- "/Users/gzy2520/Desktop/library/oligo5.csv"
readr::write_csv(out, outfile)
