library(readr); library(readxl); library(writexl); library(dplyr); library(stringr)

BstXI   <- "CCACCTTGTTGG"
BUF     <- "GC"
CONST5  <- "gatccgacgcgccatctctag"
bc_len  <- 12
set.seed(25)

exp_path    <- "/Users/gzy2520/Desktop/library/oligo7_3.csv"
control_xlsx<- "/Users/gzy2520/Desktop/library/control_BstXI.xlsx"

exp_out  <- sub("\\.csv$", "_with_barcode.csv", exp_path)
ctrl_out <- if(file.exists(control_xlsx)) sub("\\.xlsx$", "_with_barcode.csv", control_xlsx) else sub("\\.csv$", "_with_barcode.csv", control_csv)
map_out  <- sub("_with_barcode.csv","_barcode_map.csv", exp_out)

rand_base <- function(n) paste0(sample(c("A","C","G","T"), n, replace=TRUE), collapse="")
ok_barcode <- function(bc, forbidden_pattern){
  if(nchar(bc) != bc_len) return(FALSE)
  if(grepl("AAAAA|CCCCC|GGGGG|TTTTT", bc)) return(FALSE)
  gcpct <- 100 * (stringr::str_count(bc, "[GC]") / nchar(bc))
  if(gcpct < 40 | gcpct > 60) return(FALSE)
  if(grepl(forbidden_pattern, bc, fixed=TRUE)) return(FALSE)
  rc <- chartr("ACGT","TGCA", paste0(rev(strsplit(bc,"")[[1]]), collapse=""))
  if(grepl(forbidden_pattern, rc, fixed=TRUE)) return(FALSE)
  TRUE
}

# ==================== 修改开始 ====================
generate_barcodes <- function(n, forbidden_pattern, existing = character(0)){
  bcs <- character(0); tries <- 0; maxtries <- max(5000, n*200)
  while(length(bcs) < n && tries < maxtries){
    tries <- tries + 1
    # 核心修改：强制barcode第一个碱基为'G'，以满足U6启动子要求
    cand <- paste0("G", rand_base(bc_len - 1))
    
    if(!ok_barcode(cand, forbidden_pattern)) next
    if(cand %in% bcs) next
    if(cand %in% existing) next
    bcs <- c(bcs, cand)
  }
  if(length(bcs) < n) stop("无法生成足够数量的合格barcode，请检查参数或增加尝试次数。")
  bcs
}
# ==================== 修改结束 ====================

insert_barcode_into_final <- function(df){
  prefix <- paste0(BUF, BstXI, CONST5)
  df <- df %>% rowwise() %>%
    mutate(fc = `final-construct`,
           # 插入位置是在CONST5之后，这正是我们想要的
           pos = regexpr(CONST5, fc, fixed=TRUE)[1],
           pos = ifelse(pos > 0, pos + nchar(CONST5) - 1, nchar(prefix)),
           pre = substring(fc, 1, pos),
           post = substring(fc, pos + 1),
           new_final = ifelse(is.na(barcode) | barcode == "", fc, paste0(pre, barcode, post))) %>%
    ungroup()
  if(any(stringr::str_count(df$new_final, fixed(BstXI)) > 1, na.rm=TRUE)){
    idxs <- which(stringr::str_count(df$new_final, fixed(BstXI)) > 1)
    stop("barcode 导致额外 BstXI 位点，行: ", paste(idxs, collapse=","), ". 请检查。")
  }
  df$`final-construct` <- df$new_final
  df <- df %>% select(-fc, -pos, -pre, -post, -new_final)
  df
}

# exp
if(!file.exists(exp_path)) stop("找不到 experimental 文件: ", exp_path)
exp_df <- read_csv(exp_path, show_col_types = FALSE)


# ctrl
ctrl_df <- read_excel(control_xlsx)

#  existing barcodes
existing_bcs <- character(0)
if("barcode" %in% names(exp_df)) existing_bcs <- unique(na.omit(exp_df$barcode))
if("barcode" %in% names(ctrl_df)) existing_bcs <- unique(c(existing_bcs, na.omit(ctrl_df$barcode)))

assign_indices <- function(df){
  if(!("barcode" %in% names(df))) {
    df$barcode <- NA_character_
  }
  which(is.na(df$barcode) | df$barcode == "")
}

idx_exp <- assign_indices(exp_df)
idx_ctrl <- assign_indices(ctrl_df)

Nexp <- length(idx_exp)
Nctrl <- length(idx_ctrl)
total_needed <- Nexp + Nctrl
if(total_needed > 0){
  new_bcs <- generate_barcodes(total_needed, BstXI, existing_bcs)
  if(length(new_bcs) < total_needed) stop("生成的合格barcode数量不足，请重试。")
  exp_df$barcode[idx_exp] <- new_bcs[1:Nexp]
  if(Nctrl>0) ctrl_df$barcode[idx_ctrl] <- new_bcs[(Nexp+1):total_needed]
}

# 插入
exp_df_final  <- insert_barcode_into_final(exp_df)
ctrl_df_final <- insert_barcode_into_final(ctrl_df)


write_csv(exp_df_final, exp_out)
write_csv(ctrl_df_final, ctrl_out)

# barcode map
map_exp <- exp_df_final %>% transmute(source="experimental", barcode=barcode, final=`final-construct`)
map_ctrl <- ctrl_df_final %>% transmute(source="control", barcode=barcode, final=`final-construct`)
write_csv(bind_rows(map_exp, map_ctrl), map_out)
message("脚本2完成：已成功添加首位为'G'的barcode。最终文件已生成：\n", exp_out, "\n", ctrl_out, "\nBarcode映射文件：", map_out)