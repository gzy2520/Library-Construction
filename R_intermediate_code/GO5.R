#五术语加DDR
library(readr); library(dplyr); library(stringr)

infile <- "/Users/gzy2520/Desktop/library/GO2.tsv"
dir.create("library", showWarnings = FALSE)

raw <- read_tsv(infile, show_col_types = FALSE, comment = "")

# 同时保留原表里有用字段
flt <- raw %>%
  transmute(
    symbol      = .data[['SYMBOL']],
    go_id       = .data[['GO TERM']],
    go_name     = .data[['GO NAME']],
    qualifier   = .data[['QUALIFIER']],
    eco_id      = .data[['ECO ID']],
    evidence    = .data[['GO EVIDENCE CODE']],
    reference   = .data[['REFERENCE']],
    with_from   = .data[['WITH/FROM']],
    taxon_id    = .data[['TAXON ID']],
    gp_db       = .data[['GENE PRODUCT DB']],
    gp_id       = .data[['GENE PRODUCT ID']],
    assigned_by = .data[['ASSIGNED BY']]
  ) 

#去重＋补回信息
genes_uniq <- flt %>%
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

write_csv(genes_uniq, "/Users/gzy2520/Desktop/library/GO2_re.csv")
