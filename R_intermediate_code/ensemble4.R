suppressPackageStartupMessages({
  library(readr); library(AnnotationDbi); library(org.Hs.eg.db)
})

# 读结果表
dat <- readr::read_csv("/Users/gzy2520/Desktop/library/oligo4.csv", show_col_types = FALSE)

# 映射
map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(na.omit(dat$Symbol)),
  columns = c("SYMBOL","ENSEMBL"),
  keytype = "SYMBOL"
)
# 只保留有效 ENSG，并按 SYMBOL 去重取第一条
map1 <- map[!is.na(map$ENSEMBL) & grepl("^ENSG", map$ENSEMBL), c("SYMBOL","ENSEMBL")]
map1 <- map1[!duplicated(map1$SYMBOL), ]

# 改列名
names(map1)[names(map1)=="SYMBOL"]  <- "Symbol"
names(map1)[names(map1)=="ENSEMBL"] <- "ID"

dat1 <- merge(dat, map1, by = "Symbol", all.x = TRUE, sort = FALSE)

# 挪位置
ord <- c("index","Symbol","ID",
         setdiff(names(dat1), c("index","Symbol","ID")))
dat1 <- dat1[ , ord]

# 导出
readr::write_csv(dat1, "/Users/gzy2520/Desktop/library/oligo4_ID1.csv")
