library(readr)
library(readxl)
library(dplyr)

experimental_file <- "/Users/gzy2520/Desktop/library/oligo7_2.csv"
control_file      <- "/Users/gzy2520/Desktop/library/control_BstXI.xlsx"

outfile <- "/Users/gzy2520/Desktop/library/final7_2.csv"

# 目标抽样数量
num_experimental <- 5910
num_control      <- 90

set.seed(25) 

experimental_oligos <- read_csv(experimental_file, show_col_types = FALSE)

# 随机抽取5910条
sampled_experimental <- experimental_oligos %>%
  slice_sample(n = num_experimental) %>%
  dplyr::select(`final-construct`, add_5prime_G) %>%
  # 添加分组标签
  mutate(Group = "experimental")

control_oligos <- read_excel(control_file) 

# 随机抽取90条
sampled_control <- control_oligos %>%
  slice_sample(n = num_control) %>%
  dplyr::select(`final-construct`, add_5prime_G) %>%
  # 添加分组标签
  mutate(Group = "control")

# 将两个数据框垂直合并
final_oligo_list <- bind_rows(sampled_experimental, sampled_control)

# 为了方便查看，可以将对照组排在前面
final_oligo_list_sorted <- final_oligo_list %>%
  arrange(Group) 

# 写入最终的CSV文件
write_csv(final_oligo_list_sorted, outfile)
