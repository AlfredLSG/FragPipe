#!/usr/bin/env Rscript

# --- Load Libraries ---
suppressPackageStartupMessages(library(data.table))
# We no longer need tidyr for the memory-intensive spread()
# suppressPackageStartupMessages(library(tidyr)) 

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript TSS_normaolize.R <ID> <group> <input_file> <outdir>", call. = FALSE)
}
ID <- args[1]
group <- args[2]
file <- args[3]
outdir <- args[4]

# --- Main Logic ---
cat(paste("Processing sample:", ID, "\n"))
cat(paste("Input file:", file, "\n"))

tryCatch({
  TSS_dat <- fread(file, header = FALSE, col.names = c("chr", "start", "end", "GENE", "POS", "Coverage"))
}, error = function(e) {
  stop(paste("Failed to read the coverage file:", file, ". Error:", e$message))
})

if (nrow(TSS_dat) == 0) {
  cat("Warning: The coverage file is empty. No output will be generated.\n")
  file.create(file.path(outdir, paste0(ID, "_TSS_mean.txt")))
  quit(save = "no", status = 0)
}

# --- MEMORY-EFFICIENT DATA TRANSFORMATION ---
cat("Calculating mean coverage profile using a memory-efficient method...\n")

# Adjust POS to be centered around 0
TSS_dat[, POS := POS - 2501]

# The key to memory efficiency:
# Directly calculate the mean coverage for each position group,
# without creating a huge intermediate wide-format table.
# 'by = POS' groups the data by the POS column.
# '.(Mean_Coverage = mean(Coverage, na.rm = TRUE))' calculates the mean for each group.
TSS_profile <- TSS_dat[, .(Mean_Coverage = mean(Coverage, na.rm = TRUE)), by = POS]

# Sort the results by position for a clean profile
setorder(TSS_profile, POS)

# --- Write Output ---
output_file <- file.path(outdir, paste0(ID, "_TSS_mean.txt"))
cat(paste("Writing mean coverage profile to:", output_file, "\n"))

# The output is a single column of mean coverage values
# We extract just the Mean_Coverage column to match the original script's output format.
write.table(TSS_profile$Mean_Coverage, output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("Processing finished successfully.\n")