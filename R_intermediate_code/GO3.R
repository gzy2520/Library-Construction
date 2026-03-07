#七个
library(readr); library(dplyr); library(stringr)

infile <- "/Users/gzy2520/Desktop/library/GO3.tsv"

raw <- read_tsv(infile, show_col_types = FALSE, comment = "")

# 过滤术语、去 NOT、只留人；同时保留原表里有用字段
flt <- raw %>%
  transmute(
    symbol      = .data[['SYMBOL']],
    go_id       = .data[['GO TERM']],
    go_name     = .data[['GO NAME']],
    eco_id      = .data[['ECO ID']],
    evidence    = .data[['GO EVIDENCE CODE']],
    reference   = .data[['REFERENCE']],
    with_from   = .data[['WITH/FROM']],
    taxon_id    = .data[['TAXON ID']],
    gp_db       = .data[['GENE PRODUCT DB']],
    gp_id       = .data[['GENE PRODUCT ID']],
    assigned_by = .data[['ASSIGNED BY']]
  )

# 1) 最终“去重名单”（库用的唯一基因列表）
uniq_symbols <- flt %>% distinct(symbol) %>% arrange(symbol)
write_csv(uniq_symbols, "library/genes_union_symbols.csv")

# 2) “去重＋补回信息”（一条/基因，汇总原始字段）
genes_master <- flt %>%
  group_by(symbol) %>%
  summarise(
    go_ids        = paste(sort(unique(go_id)), collapse=";"),
    go_names      = paste(sort(unique(go_name)), collapse=";"),
    evidence_set  = paste(sort(unique(evidence)), collapse=";"),
    eco_set       = paste(sort(unique(eco_id)), collapse=";"),
    refs          = paste(sort(unique(reference)), collapse=";"),
    gene_products = paste(sort(unique(gp_id)), collapse=";"),
    assigned_by   = paste(sort(unique(assigned_by)), collapse=";"),
    n_annotations = dplyr::n(),
    .groups="drop"
  ) %>% arrange(symbol)

write_csv(genes_master, "/Users/gzy2520/Desktop/library/GO3_re.csv")
