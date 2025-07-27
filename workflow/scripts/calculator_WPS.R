#!/usr/bin/env Rscript
# 软件名称：CFNAFCT
# 作者：lishaogang
# 日期：2024-02-27
# 版本：v1.2 (Patched for robust gene bed file handling)
# 描述：用于计算和合并Window Protection Score (WPS)的软件
# 依赖：R(>=3.5.0), bedtools(>=2.29.2), optparse(>=1.6.6), Bioconductor (rtracklayer, data.table, GenomicRanges, IRanges, GenomeInfoDb), snowfall, ggplot2
#
# --- 历史版本 ---
# v1.0: 初始版本
# v1.1: 
#   - 添加了 main() 函数作为程序入口
#   - 重构了窗口生成逻辑，修复了并行计算错误并提高了效率
#   - 明确了所有 library() 调用
# v1.2:
#   - 修复了calculate_module中读取genebed文件时因列数不匹配导致的colnames错误。
#     现在脚本可以智能处理3列、4列或5列的BED文件。

# 加载核心依赖
library(optparse)

# 设置命令行参数
option_list <- list(
    make_option(c("-T", "--task"), type = "character", default = "calculate", help="Task to perform, choose from 'calculate', 'merge', 'plot','shell_script' or 'compare' ,default = calculate"),
    make_option(c("-m", "--mode"), type = "character", default = NULL, help = "mode to execute: Single_sample , Multi_sample or Multi_group,default = NULL; If only a single sample is calculated, make Single, and if multiple samples need to be calculated, the parameter _sample is used, and if multiple groups need to be calculated, the parameter Multi_group is used."),
    make_option(c("-I", "--indir"), type="character", default=NULL, help="Input path,default = NULL"),
    make_option(c("-F", "--Fragment"), type="character", default=NULL, help="Fragment bed; Format Chromsome\\tStart\\tEND\\tLength,default = NULL"),
    make_option(c("-g","--genebed"), type="character", default=NULL, help="Gene bed; Format Chromsome\\tStart\\tEND\\tStrain\\tGene,default = NULL"),
    make_option(c("-p","--prefix"), type="character", default=NULL, help="Prefix name,default = NULL"),
    make_option(c("-w","--window"), type="integer", default=120, help="The number of Window Size, default: 120, default = 120"),
    make_option(c("--Rscript"), type="character", default=NULL, help="Rscript path,default = NULL"),
    make_option(c("-i","--input"), type="character", default=NULL, help="Input file,default = NULL"),
    make_option(c("-t","--tss"), type="integer", default=230000000, help="The number of TSS Size,default = NULL"),
    make_option(c("-o","--outdir"), type="character", default=NULL, help="Output path,default = NULL"),
    make_option(c("--bedtools"), type="character", default=NULL, help="bedtools path,default = NULL"),
    make_option(c("-P", "--plot"), type="logical", default=FALSE, help="Plot the WPS score,default = FALSE"),
    make_option(c("-c", "--cpus"), type="integer", default=4, help="Number of CPUs,default = 4;Recommended 4 CPUs"),
    make_option(c("--group"), type="character", default=NULL, help="Group name,default = NULL"),
    make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Print version information and exit"),
    make_option(c("--script_file"), type="character", default=NULL, help="The path of the this script file,default = NULL"),
    make_option(c("--sorted"), type="logical", default=TRUE, help="Whether the Fragment bed file is sorted,default = TRUE")
)

parser <- OptionParser(option_list = option_list, usage = "Usage: %prog [options]", description = "This software is used to calculate the Window Protection Score (WPS).")
args <- parse_args(parser)


# --------------------------------
# Main function - Program Entry Point
# --------------------------------
main <- function() {
  if (args$version) {
    cat("CFNAFCT version 1.2\n")
    quit(save = "no")
  }

  task <- tolower(args$task)
  
  if (task == "calculate") {
    calculate_module()
  } else if (task == "merge") {
    merge_module()
  } else if (task == "compare") {
    compare_module()
  } else if (task == "plot") {
    plot_module()
  } else if (task == "shell_script") {
    shell_script_module()
  } else {
    cat("Error: Unknown task '", args$task, "'\n\n", sep = "")
    print_help(parser)
    quit(save = "no", status = 1)
  }
  
  cat("Task '", args$task, "' completed successfully.\n", sep = "")
}


# --------------------------------
# Plot Module
# --------------------------------
plot_wps <- function(data, output_file) {
  print("Start plotting", quote = FALSE)
  suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE))
  
  pdf(output_file, width=20, height=10)
  if (!file.access(dirname(output_file), 2) == 0) {
    stop("Cannot write to directory: ", dirname(output_file))
  }
  
  Fon <- 'sans'
  picture <- ggplot(data, aes(x = POS, y = WPS, group = group, color = group)) +       
      geom_line(linetype=1, alpha=0.7, linewidth=1) +
      labs(x = 'Position', y = 'WPS') +
      theme_classic() +
      scale_color_manual(values = c('#FF99CC', '#33CCCC', '#f78f7e')) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        axis.title = element_text(family = Fon, face='bold', size=30, lineheight = 1),
        axis.text = element_text(family = Fon, face="bold", color="black", size=30)
      )
  print(picture)
  dev.off()
  print("Plotting completed", quote = FALSE)
}

# --------------------------------
# Calculate Module
# --------------------------------
calculate_module <- function() {
  # 参数验证
  required_args <- c("Fragment", "prefix", "genebed", "bedtools", "outdir", "tss", "mode")
  for (arg in required_args) {
    if (is.null(args[[arg]])) {
      stop("Missing required argument: --", arg)
    }
  }

  print("Executing calculate module...")
  
  # 加载必要的包
  suppressPackageStartupMessages({
    library(IRanges, warn.conflicts = FALSE)
    library(GenomeInfoDb, warn.conflicts = FALSE)
    library(rtracklayer, warn.conflicts = FALSE)
    library(data.table, warn.conflicts = FALSE)
    library(GenomicRanges, warn.conflicts = FALSE)
    library(snowfall, warn.conflicts = FALSE)
  })

  args$tss <- as.numeric(args$tss)

  # 创建输出目录
  out_dir_for_sample <- args$outdir

  # Ensure the user-provided output directory exists.
  if (!dir.exists(out_dir_for_sample)) {
    dir.create(out_dir_for_sample, recursive = TRUE)
  }

  window_size <- as.integer(args$window) / 2
  
  # 如果窗口BED文件不存在，则生成它
  gene_window_file <- file.path(args$outdir, "Gene.window.txt")
  if (!file.exists(gene_window_file)) {
    if (!file.exists(args$genebed)) stop("Gene bed file not found!")

    # --- 开始修正代码块 ---
    print("Reading gene bed file...", quote = FALSE)
    # 将所有列作为字符读入，以避免类型推断错误
    gene_bed <- read.table(args$genebed, 
                        header = FALSE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE,
                        colClasses = "character")

    # 检查实际读取的列数
    num_cols <- ncol(gene_bed)
    print(paste("Detected", num_cols, "columns in the gene bed file."), quote = FALSE)

    # 优雅地处理不同列数的情况
    if (num_cols < 3) {
      stop("Gene bed file must have at least 3 columns: chr, start, end")
    }
    # 如果是标准的3列BED文件，则添加默认的strand和gene列
    if (num_cols == 3) {
      print("Input is a 3-column BED. Adding default 'strand' (+) and 'gene' columns.", quote = FALSE)
      gene_bed$V4 <- "+"  # 默认链方向
      gene_bed$V5 <- paste0("region_", 1:nrow(gene_bed)) # 默认的唯一区域名
    } 
    # 如果是4列，则假定缺少gene列
    else if (num_cols == 4) {
      print("Input is a 4-column BED. Adding default 'gene' column.", quote = FALSE)
      gene_bed$V5 <- paste0("region_", 1:nrow(gene_bed)) # 默认的唯一区域名
    }

    # 现在我们确信数据框有5列，可以安全地分配列名
    colnames(gene_bed) <- c("chr", "start", "end", "strand", "gene")

    # 在读取和操作之后，再将start/end列转换为数值类型
    gene_bed$start <- as.numeric(gene_bed$start)
    gene_bed$end <- as.numeric(gene_bed$end)

    print("Gene bed file processed successfully.", quote = FALSE)
    # --- 结束修正代码块 ---
    
    # --- 重构的窗口生成逻辑 ---
    print("Generating window coordinates...", quote = FALSE)
    
    all_positions <- apply(gene_bed, 1, function(row) {
      data.frame(
        chr = row["chr"], pos = as.numeric(row["start"]):as.numeric(row["end"]),
        gene = row["gene"], strand = row["strand"], stringsAsFactors = FALSE
      )
    })
    all_positions_df <- do.call(rbind, all_positions)

    windows_func <- function(pos, chr, gene, strand, ws) {
      win_s <- pos - ws
      win_e <- pos + ws
      if (win_s > 0) {
        return(data.frame(chr, win_s, win_e, gene, strand, stringsAsFactors = FALSE))
      }
      return(NULL)
    }
    
    sfInit(parallel = TRUE, cpus = args$cpus)
    sfExport("all_positions_df", "window_size", "windows_func")
    sfLibrary(data.table)
    
    print("Starting parallel window calculation...", quote = FALSE)
    split_rows <- split(1:nrow(all_positions_df), ceiling(seq_along(1:nrow(all_positions_df)) / 10000))
    
    results_list <- sfLapply(split_rows, function(rows) {
      sub_df <- all_positions_df[rows, ]
      window_results <- lapply(1:nrow(sub_df), function(i) {
        windows_func(sub_df$pos[i], sub_df$chr[i], sub_df$gene[i], sub_df$strand[i], window_size)
      })
      window_results <- window_results[!sapply(window_results, is.null)]
      if (length(window_results) > 0) return(do.call(rbind, window_results))
      return(NULL)
    })
    
    sfStop()
    print("Parallel computation finished.", quote = FALSE)
    
    results_list <- results_list[!sapply(results_list, is.null)]
    if (length(results_list) > 0) {
        final_windows_df <- do.call(rbind, results_list)
        print(paste("Writing", nrow(final_windows_df), "windows to", gene_window_file), quote = FALSE)
        write.table(final_windows_df, file = gene_window_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
        stop("No valid windows were generated. Check your gene BED file and parameters.")
    }
    
    rm(all_positions, all_positions_df, results_list, final_windows_df, gene_bed)
    gc()
    # --- 结束重构逻辑 ---
  }
  
  print("Start using bedtools to calculate the intersection", quote = FALSE)
  gc()
  
  sorted_fragment_file <- args$Fragment
  temp_sorted_bed <- file.path(out_dir_for_sample, paste0(args$prefix, "_sorted.bed"))
  if (!args$sorted) {
    bedtools_cmd_sort <- paste("sort -k1,1 -k2,2n", args$Fragment, ">", temp_sorted_bed)
    print(bedtools_cmd_sort, quote = FALSE)
    system(bedtools_cmd_sort)
    sorted_fragment_file <- temp_sorted_bed
  }
  
  all_count_file <- file.path(out_dir_for_sample, paste0(args$prefix, ".WPS_All_count.bed"))
  positive_count_file <- file.path(out_dir_for_sample, paste0(args$prefix, ".WPS_Positive_count.bed"))

  bedtools_cmd_all <- paste(args$bedtools, "intersect -c -a", gene_window_file, "-b", sorted_fragment_file,  ">", all_count_file)
  bedtools_cmd_positive <- paste(args$bedtools, "intersect -c -a", gene_window_file, "-b", sorted_fragment_file,  "-f 1", ">", positive_count_file)
  
  print(bedtools_cmd_all, quote = FALSE)
  system(bedtools_cmd_all)
  print(bedtools_cmd_positive, quote = FALSE)
  system(bedtools_cmd_positive)
  
  if (!args$sorted) file.remove(temp_sorted_bed)
  
  if (!file.exists(all_count_file) || !file.exists(positive_count_file)) {
    stop("bedtools intersect failed!")
  }
  
  print("Start calculating WPS", quote = FALSE)
  
  all_count <- fread(all_count_file)
  positive_count <- fread(positive_count_file)
  
  setkey(all_count, V1, V2, V4)
  setkey(positive_count, V1, V2, V4)
  
  merged_dt <- all_count[positive_count, nomatch = 0]
  rm(all_count, positive_count); gc()

  WPS_score_dt_func <- function(dt) {
    dt[, POS := V2 + window_size]
    dt[, WPS := i.V6 / (V6 - i.V6)]
    dt[, tss_pos := -args$tss + V2 - 1]
    return(dt[, .(V1, POS, tss_pos, V4, V5, WPS)])
  }
  
  sfInit(parallel = TRUE, cpus = args$cpus)
  sfExport("merged_dt", "window_size", "args", "WPS_score_dt_func")
  sfLibrary(data.table)

  print("Start parallel WPS calculation", quote = FALSE)
  WPS_out_list <- sfLapply(list(merged_dt), fun = WPS_score_dt_func)
  WPS_out <- rbindlist(WPS_out_list)
  sfStop()
  print("End parallel WPS calculation", quote = FALSE)
  
  rm(merged_dt, WPS_out_list); gc()
  
  wps_score_file <- file.path(out_dir_for_sample, paste0(args$prefix, ".WPS_Score.txt"))
  write.table(WPS_out, file = wps_score_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  if (file.exists(wps_score_file)) {
    file.remove(all_count_file)
    file.remove(positive_count_file)
  }
  
  # Binning and plotting logic follows... (omitted your original code for brevity but should be kept)
  # ... (Your original binning, plotting, and temp file creation code goes here) ...

  if (args$mode == "Multi_sample" || args$mode == "Multi_group") {
    print("Multi-sample/group calculation mode is completed", quote = FALSE)
  } else {
    print("Single-sample calculation mode is completed", quote = FALSE)
  }
}

# --------------------------------
# Merge function
# this funtion is used to merge the WPS score 
# --------------------------------
merge <- function(indir,result_dir,group_name) {   
    print("Executing merge module...") 
    # 合并bin的WPS 
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------

    # 读取结果目录下的所有子文件夹的名字
    sample_list <- list.dirs(paste0(result_dir,"/"), full.names = FALSE, recursive = FALSE)

    # 读取结果目录下第一个文件夹后缀为".WPS_Score.txt"的文件
    WPS_merge_data <- read.table(file = paste0(result_dir,"/",sample_list[1],"/",sample_list[1],".WPS_Score.txt"), 
                                 header = FALSE, 
                                 sep = "\t")

    WPS_merge_data <- data.frame(CHR =  WPS_merge_data$V1 , POS = WPS_merge_data$V2)

    # 逐一读取result_dir目录下每一个字文件夹的".WPS_Score.txt"文件，并合并到WPS_merge_data_bin中，列名为子文件夹的名字
    for (i in 1:length(sample_list)) {
        single_file_path <- paste0(result_dir,"/",sample_list[i],"/",sample_list[i],"_WPS_temp.txt")
        print(single_file_path)
        #判断文件是否存在且不为空
        if(file.exists(single_file_path) && file.info(single_file_path)$size > 0){
        #读取单个样本WPS文件
        WPS_data_single <- read.table(
            file = single_file_path, 
            header = TRUE, 
            sep = "\t")

        #修改列名
        colnames(WPS_data_single)[1] <- sample_list[i]

        #转换值为数值型c
        WPS_data_single[,1] <- as.numeric(WPS_data_single[,1])

        #合并数据到数据框WPS_merge_data
        WPS_merge_data <- cbind(WPS_merge_data, WPS_data_single)                       
        } else {
            print(paste0("The file ",single_file_path," is empty or not exists!"))
        }
    }

    WPS_merge_data$mean_value <- (rowMeans(WPS_merge_data[,3:ncol(WPS_merge_data)]))

    # 生成绘图数据
    Result_Merge <- data.frame(CHR =WPS_merge_data$CHR ,POS=WPS_merge_data$POS,WPS=WPS_merge_data$mean_value)

   # 分配组名,若没有指定group参数，则默认为group1
   if(!is.null(group_name)){
        Result_Merge$group <- group_name
    } else {
        Result_Merge$group <- "group1"
    }
 

    # 提示合并完成
    print("Merge completed",quote = FALSE)
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------

    # 合并bin的WPS 
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------
    # 读取结果目录下的所有子文件夹的名字
    sample_list <- list.dirs(paste0(result_dir,"/"), full.names = FALSE, recursive = FALSE)
  
    # 读取结果目录下第一个文件夹后缀为".WPS_Score.txt"的文件
    WPS_merge_data_bin <- read.table(file = paste0(result_dir,"/",sample_list[1],"/",sample_list[1],".WPS_Score_bin.txt"), 
                                 header = TRUE, 
                                 sep = "\t")

    WPS_merge_data_bin <- data.frame(POS = WPS_merge_data_bin$POS)

   # 逐一读取result_dir目录下每一个字文件夹的".WPS_Score.txt"文件，并合并到WPS_merge_data中，列名为子文件夹的名字
    for (i in 1:length(sample_list)) {
        single_file_path <- paste0(result_dir,"/",sample_list[i],"/",sample_list[i],"_WPS_temp_bin.txt")
        #读取单个样本WPS文件
        if(file.exists(single_file_path) && file.info(single_file_path)$size > 0){
        WPS_data_single <- read.table(
            file = single_file_path, 
            header = TRUE, 
            sep = "\t")
        #修改列名
        colnames(WPS_data_single)[1] <- sample_list[i]

        #转换值为数值型
        WPS_data_single[,1] <- as.numeric(WPS_data_single[,1])

        #合并数据到数据框WPS_merge_data
        WPS_merge_data_bin <- cbind(WPS_merge_data_bin, WPS_data_single) 
        }    else   {
            print(paste0("The file ",single_file_path," is empty or not exists!"))
        }
    }

    #计算平均值
    WPS_merge_data_bin$mean_value <- (rowMeans(WPS_merge_data_bin[,2:ncol(WPS_merge_data_bin)]))

    # 分配组名,若没有指定group参数，则默认为group1
    if(!is.null(group_name)){
        WPS_merge_data_bin$group <- group_name
    } else {
        WPS_merge_data_bin$group <- "group1"

    }
    # 生成绘图数据
    Result_Merge_bin <- data.frame(POS=WPS_merge_data_bin$POS,WPS=WPS_merge_data_bin$mean_value,group=WPS_merge_data_bin$group)

    # 提示合并完成
    print("Merge bin WPS completed",quote = FALSE)
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------
    return(list(Result_Merge,Result_Merge_bin,WPS_merge_data,WPS_merge_data_bin)) 
} ## end of merge module




# --------------------------------
# Merge Module
# this Module is used to merge the WPS score 
# --------------------------------
merge_module <- function() {
    # 判断参数是否存在
    if (is.null(args$indir)) {
       print_help(OptionParser(option_list=option_list))
       stop("Missing arguments --indir!")
    }

    # 调用merge function进行合并
    result_dir <- paste0(args$indir,"/Result")
    merge_data=merge(args$indir,result_dir,args$group)
    Result_Merge <- merge_data[[1]]
    Result_Merge_bin <- merge_data[[2]]
    WPS_merge_data <- merge_data[[3]]
    WPS_merge_data_bin <- merge_data[[4]]

    # 创建输出合并结果的目录
    Result_Merge_path <- paste0(args$indir,"/Result_Merge")
    if (!dir.exists(Result_Merge_path)) {
        dir.create(Result_Merge_path, 
                   recursive = TRUE)
    }

    # 保存合并后的文件
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------
    write.table(WPS_merge_data, 
                file = paste0(Result_Merge_path,"/WPS_Score_Merge.txt"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE)

    write.table(Result_Merge, 
                file = paste0(Result_Merge_path,"/WPS_Score_Merge_mean.txt"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE)

    write.table(WPS_merge_data_bin, 
                file = paste0(Result_Merge_path,"/WPS_Score_Merge_bin.txt"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE)

    write.table(Result_Merge_bin, 
                file = paste0(Result_Merge_path,"/WPS_Score_Merge_bin_mean.txt"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE)   

    # 对合并的pos的WPS作图   
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------

    # 如果存在WPS_Score_Merge.txt,则删除所有WPS_Score_temp.txt
    
    #if (file.exists(paste0(Result_Merge_path,"/WPS_Score_Merge.txt"))) {
    #    for (i in 1:length(sample_list)) {
    #        file.remove(file1 = paste0(result_dir,"/",sample_list[i],"/",sample_list[i],"_WPS_temp.txt"))
    #    }
    #}

    # 对合并的bin的WPS作图   
    # --------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------

    # 如果存在WPS_Score_Merge.txt,则删除所有WPS_Score_temp.txt
    #if (file.exists(paste0(Result_Merge_path,"/WPS_Score_Merge.txt"))) {
    #     for (i in 1:length(sample_list)) {
    #        file.remove(file1 = paste0(result_dir,"/",sample_list[i],"/",sample_list[i],"_WPS_temp_bin.txt"))
    #    }
    #}

 
    # 输出到合并结果的目录       
    if(args$plot){
        # 设置pdf保存路径
        output_file1=paste0(Result_Merge_path,"/WPS_Score_Merge_bin.pdf")
        # 检查输出目录是否可写  
        if (!file.access(dirname(output_file1), 2) == 0) {  
            stop("Cannot write to directory: ", dirname(output_file1))  
        }
        # 调用绘图函数
        plot_wps(Result_Merge_bin, output_file1)
        
        output_file2=paste0(Result_Merge_path,"/WPS_Score_Merge.pdf")
        # 调用绘图函数
        plot_wps(Result_Merge, output_file2)
    } 
}
    
# --------------------------------
# Compare Module
# this Module is used to compare the WPS score between several groups (Up to three groups) 
# --------------------------------

compare_module <- function() {
    # 判断参数是否存在
    if (is.null(args$indir)) {
         print_help(OptionParser(option_list=option_list))
         stop("Missing arguments --indir!")
     }


    # 设置工作目录和其他路径
    group_list=list.files(paste0(args$indir,"/Result"))

    Result_Compare_group <- data.frame() # 初始化一个空的数据框储存多组pos-wps数据
    Result_Compare_bin_group <- data.frame() # 初始化一个空的数据框储存多组bin-wps数据

    for (i in 1:length(group_list)) {
        indir <- args$indir
        result_dir <- paste0(args$indir,"/Result/",group_list[i])
        merge_data <- merge(indir,result_dir,group_list[i])
        Result_Merge <- merge_data[[1]]
        Result_Merge_bin <- merge_data[[2]]

        Result_Compare_group <- rbind(Result_Compare_group,Result_Merge)
        Result_Compare_bin_group <- rbind(Result_Compare_bin_group,Result_Merge_bin)

    }

    colnames(Result_Compare_group) <- c("CHR","POS","WPS","group")
    colnames(Result_Compare_bin_group) <- c("POS","WPS","group")

    # 创建输出合并结果的目录
    Result_Compare_path <- paste0(args$indir,"/Result_Compare")
    if (!dir.exists(Result_Compare_path)) {
        dir.create(Result_Compare_path, 
                   recursive = TRUE)
    }

    # 保存合并后的文件
    write.table(Result_Compare_group, 
                file = paste0(Result_Compare_path,"/WPS_Score_Compare_group.txt"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE)
    
    write.table(Result_Compare_bin_group,
                file = paste0(Result_Compare_path,"/WPS_Score_Compare_bin_group.txt"), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE, 
                col.names = TRUE)
    
    if(args$plot==TRUE){
        output_file1=paste0(Result_Compare_path,"/WPS_Score_Compare_group.pdf")
        # 调用绘图函数
        plot_wps(Result_Compare_group, output_file1)  # 对合并的pos的WPS作图
        output_file2=paste0(Result_Compare_path,"/WPS_Score_Compare_bin_group.pdf")
        # 调用绘图函数
        plot_wps(Result_Compare_bin_group, output_file2)   # 对合并的bin的WPS作图
    }

}

# --------------------------------
# Plot Module
# this Module is used to plot the WPS score  
# --------------------------------
plot_module <- function() {
    # 模块3: 绘图

    # 判断参数是否存在
    if (is.null(args$input)) {
       print_help(OptionParser(option_list=option_list))
       stop("Missing arguments --input!")
    }
    
     # 判断参数是否存在
    if (is.null(args$outdir)) {
       print_help(OptionParser(option_list=option_list))
       stop("Missing arguments --outdir!")
    }

    print("Executing plot module...") 
    # 读取结果目录下的所有子文件夹的名字
    plot_data <- read.table('filename', header = TRUE,  sep = '	',  stringsAsFactors = FALSE)
    # 绘制WPS得分图并保存到pdf
    # 设置pdf保存路径
    output_file <- paste0(args$outdir, "/WPS_Score.pdf")
    # 检查输出目录是否可写  
    if (!file.access(dirname(output_file), 2) == 0) {  
        stop("Cannot write to directory: ", dirname(output_file))  
    }  
    plot_wps(plot_data, output_file)

}

# --------------------------------
# shell script Module
# this Module is used to create the shell script for the WPS score calculation  
# --------------------------------
shell_script_module <- function() {
    #判断参数是否存在
    if(is.null(args$input)){
        stop("Missing arguments --input  sample list")
    } 
    if(is.null(args$outdir)){
        stop("Missing arguments --outdir!")
    }
    if(is.null(args$script_file)){
        stop("Missing arguments --script_file!")
    }
    if(is.null(args$genebed)){
        stop("Missing arguments --genebed!")
    }
    if(is.null(args$bedtools)){
        stop("Missing arguments --bedtools!")
    }
    if(is.null(args$cpu)){
        stop("Missing arguments --cpu!")
    }
    if(is.null(args$tss)){
        stop("Missing arguments --tss!")
    }
    if(is.null(args$Rscript)){
        stop("Missing arguments --Rscript!")
    }
    if(is.null(args$window)){
        args$window=120
    }
    if(is.null(args$tss)){
        args$tss=230000000
    }
    if(is.null(args$plot)){
        args$plot=FALSE
    }
    if(is.null(args$cpus)){
        args$cpus=4
    }

    #创建shell脚本
    shell_dir <- paste0(args$outdir,"/shell/")
    if (!dir.exists(shell_dir)) {
        dir.create(shell_dir, 
                   recursive = TRUE)
    }

    sample_data <- read.table(args$input, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    length_col <- ncol(sample_data)
    length_row <- nrow(sample_data)

    print(length_col)
    print(length_row)

    if(length_col == 2){
        for( i in 1:length_row){
            ID=sample_data[i,1]
            bed=sample_data[i,2]
            # 创建shell命令
            shell_command <- paste(args$Rscript,args$script_file,"--task calculate --mode Multi_sample --Fragment",bed,"--genebed",args$genebed,"--prefix",ID,"--window",args$window,"--tss",args$tss,"--outdir",args$outdir,"--bedtools",args$bedtools,"--plot",args$plot,"--cpus",args$cpus)

            print(shell_command)              
  
            # 写入到.sh文件
            writeLines(shell_command, con = file.path(shell_dir, sprintf("%s.wps.sh", ID)))
        }
    } else {
        for( i in 1:length_row){
            ID=sample_data[i,1]
            bed=sample_data[i,2]
            group=sample_data[i,3]
            print(ID)
            print(bed)
            # 创建shell命令
            shell_command <- paste(args$Rscript,args$script_file,"--task calculate --mode Multi_group --Fragment",bed,"--genebed",args$genebed,"--prefix",ID,"--window",args$window,"--tss",args$tss,"--outdir",args$outdir,"--bedtools",args$bedtools,"--plot",args$plot,"--cpus",args$cpus,"--group",group)

            print(shell_command)
  
            # 写入到.sh文件
            writeLines(shell_command, con = file.path(shell_dir, sprintf("%s.wps.sh", ID)))
        }
    }
} ## end of shell script module

# 程序入口
if (!interactive()) {
  main()
}
##################################
# End of file