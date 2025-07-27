
args <- commandArgs(trailingOnly = TRUE)
ID <- args[1]
group <- args[2]
file <- args[3]
outdir <- args[4]

library(tidyr)

CTCF_dat <- read.table(file, head = FALSE)

CTCF_dat$V6 <- paste(CTCF_dat$V1, CTCF_dat$V2, sep = "_")
CTCF_dat <- CTCF_dat[, c(4:6)]

colnames(CTCF_dat) <- c("POS", "CTCF_SCORE", "GENE")

CTCF_dat$POS <- CTCF_dat$POS - 1501

CTCF_dat <- spread(CTCF_dat, key = "GENE", value = "CTCF_SCORE")


#normalize_column <- function(x) {
  # 计算特定元素的平均值M
#  M <- mean(x[c(1:500, 2502:3001)])
  
  # 如果M为0，则返回NULL，否则返回标准化后的列
#  if (M != 0) {
#    return(x / M)
#  } else {
#    return(NULL)
#  }
#}

# 应用这个函数到数据框的每一列
#CTCF_dat[] <- lapply(CTCF_dat, normalize_column)

# 移除值为NULL的列（即原来平均值为0的列）
#CTCF_dat  <- CTCF_dat[, colSums(is.na(CTCF_dat)) != nrow(CTCF_dat)]
CTCF_dat$Coverage=rowMeans(CTCF_dat[,2:ncol(CTCF_dat)])
CTCF_dat=CTCF_dat[,c("POS","Coverage")]
M=mean(CTCF_dat$Coverage[c(1:500, 2502:3001)])
CTCF_dat$Coverage=CTCF_dat$Coverage/M






write.table(CTCF_dat$Coverage, paste0(outdir, "/", ID, "_ctcf_mean_normalized_new.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

