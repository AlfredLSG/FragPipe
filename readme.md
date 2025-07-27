# FragPipe: A Snakemake Pipeline for cfDNA Fragmentomics Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/Made%20with-Snakemake-blue.svg)](https://snakemake.readthedocs.io)

**FragPipe** is a reproducible, end-to-end Snakemake pipeline designed for **whole-genome sequencing (WGS)** and **whole-genome bisulfite sequencing (WGBS)** data from **cell-free DNA (cfDNA)**. It covers the entire workflow from raw FASTQ files to comprehensive **fragmentomics** analysis, including Window Protection Score (WPS), Transcription Start Site (TSS) enrichment, CTCF binding site enrichment, and fragment end motif analysis.

This pipeline aims to provide an automated and standardized solution for exploring the biological significance of cfDNA fragmentation patterns in liquid biopsy research.



---

## Key Features

*   **Automated Workflow**: Fully automated process from raw FASTQ files to final analysis results.
*   **Reproducibility**: Utilizes Conda for strict environment management, ensuring consistent results on any machine.
*   **Flexibility**: Supports both WGS and WGBS data types and allows for custom sequencing adapters per sample.
*   **Comprehensive Fragmentomics Analysis**:
    *   **Quality Control**: Quality control and adapter trimming using `fastp`.
    *   **Alignment**: Alignment using `BWA-MEM` (for WGS) or `bwameth` (for WGBS).
    *   **BAM Processing**: Sorting and marking PCR duplicates with `Picard MarkDuplicates`.
    *   **Window Protection Score (WPS)**: Calculation of fragment protection scores at specific genomic loci.
    *   **TSS Enrichment**: Assessment of fragment enrichment patterns around Transcription Start Sites.
    *   **CTCF Enrichment**: Evaluation of fragment enrichment patterns around CTCF binding sites.
    *   **End Motif Analysis**: Analysis of nucleotide preferences at cfDNA fragment ends.
*   **Tidy Output Structure**: All output files are neatly organized into a dedicated directory for each sample.

## Workflow Diagram
![alt text](fragpipe_dag-1.svg)


---

## Deployment and Installation

Follow these steps to deploy and prepare the FragPipe environment.

### Step 1: Clone the Repository

First, clone this repository to your local machine.
```bash
git clone https://github.com/your_username/FragPipe.git
cd FragPipe
```

### Step 2: Install Conda

This pipeline relies on Conda for environment management. If you don't have Conda installed, we recommend installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html), a lightweight distribution.

For faster environment creation, it is highly recommended to install Mamba after installing Conda:
```bash
conda install -n base -c conda-forge mamba
```

### Step 3: Create the Conda Environment

This repository provides an `environment.yaml` file that defines all necessary software dependencies.
```bash
# Recommended: Use Mamba (much faster)
mamba env create -f workflow/environment.yaml

# Alternative: Use Conda (slower)
# conda env create -f workflow/environment.yaml
```
This will create a Conda environment named `FragPipe-env` (or as defined in the yaml file). Always activate this environment before running the pipeline:
```bash
conda activate FragPipe-env
```

### Step 4: Download and Prepare Reference Genome

A convenience script is provided to automatically download the human reference genome (hg38) and build the required BWA and Samtools indexes.
```bash
# 1. Grant execution permissions to the script
chmod +x download_hg38_build_index.sh

# 2. Run the script (ensure the Conda environment is activated)
# This will download a ~1GB file and may take 30-60 minutes to build indexes.
./download_hg38_build_index.sh
```

---

## Usage

### Step 1: Prepare Input Files

1.  Place your raw **FASTQ files** in a directory of your choice.
2.  Edit the **`Projects/Myproject1/samples.tsv`** file. This is a **headerless, tab-separated** table that defines your samples. Modify this file or create a new project directory based on this template.

The `samples.tsv` file must contain **7 columns**:

| Column | Content         | Description                                        | Example                                          |
|--------|-----------------|----------------------------------------------------|--------------------------------------------------|
| 1      | `sample_id`     | Unique ID for the sample, used in all output files.| `Patient1_Tumor`                                 |
| 2      | `fastq1_path`   | Absolute path to the Read 1 FASTQ file.            | `/path/to/data/P1_T1_R1.fq.gz`                   |
| 3      | `fastq2_path`   | Absolute path to the Read 2 FASTQ file.            | `/path/to/data/P1_T1_R2.fq.gz`                   |
| 4      | `group`         | The group the sample belongs to (e.g., Tumor).     | `Tumor`                                          |
| 5      | `type`          | Analysis type, must be `wgs` or `wgbs`.            | `wgs`                                            |
| 6      | `fwd_adapter`   | Sequencing adapter for Read 1.                     | `AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA`                |
| 7      | `rev_adapter`   | Sequencing adapter for Read 2.                     | `AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG`     |

### Step 2: Run the Pipeline

Execute the following command from the project's root directory (`FragPipe/`) to start the analysis.
```bash
# Activate the environment
conda activate FragPipe-env

# Run Snakemake
snakemake \
    --snakefile workflow/snakefile \
    --cores 8 \
    --config project_name=Myproject1 \
    --use-conda
```

**Parameter Explanation**:
*   `--snakefile`: Path to the main Snakefile.
*   `--cores`: The maximum number of CPU cores to use.
*   `--config project_name=...`: **(Required)** The name of your project directory under `Projects/`.
*   `--use-conda`: Instructs Snakemake to manage environments using `environment.yaml`.

---

## Output Structure

Upon successful completion, all results will be neatly organized in the `Results/<YourProjectName>/` directory, with a dedicated folder for each sample.

```
Results/
└── Myproject1/
    └── sample_C_20k/
        ├── alignment/
        │   ├── sample_C_20k.markdup.bam      # Final alignment file (coordinate sorted)
        │   └── sample_C_20k.markdup.bam.bai  # BAM index
        ├── ctcf/
        │   └── sample_C_20k_..._mean.txt     # CTCF enrichment analysis results
        ├── motif/
        │   ├── sample_C_20k.4mer.motif.MDS # 4-mer motif diversity score
        │   └── sample_C_20k.size_ratio     # Short vs. long fragment ratio
        ├── qc/
        │   ├── sample_C_20k.fastp.html     # Fastp quality control report
        │   └── sample_C_20k.mkdup_metrics.txt # Duplication rate metrics
        ├── tss/
        │   └── sample_C_20k_TSS_mean.txt   # TSS enrichment analysis results
        └── wps/
            └── sample_C_20k.WPS_Score.txt  # Window Protection Score results
```

---

## License

This project is licensed under the [MIT License](LICENSE.md).

## Contact
lishaogang - [li.shaogang97@gmail.com](mailto:li.shaogang97@gmail.com
)

Project Link: [https://github.com/AlfredLSG/FragPipe](https://github.com/AlfredLSG/FragPipe)


-----------------------------中文版---------------------------------------------------



# FragPipe: cfDNA 片段组学分析流程

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/Made%20with-Snakemake-blue.svg)](https://snakemake.readthedocs.io)

**FragPipe** 是一个为**细胞游离DNA (cfDNA)** 的**全基因组测序 (WGS)** 或**全基因组亚硫酸氢盐测序 (WGBS)** 数据设计的、可复现的、端到端的 Snakemake 分析流程。它涵盖了从原始测序数据（FASTQ）到全面的**片段组学 (fragmentomics)** 分析，包括窗口保护分数 (WPS)、转录起始位点 (TSS) 富集、CTCF结合位点富集以及片段末端基序分析。

本流程旨在提供一个自动化、标准化的解决方案，以探索cfDNA片段模式在液体活检研究中的生物学意义。

---

## 核心功能

*   **自动化流程**: 从原始FASTQ文件到最终分析结果的全程自动化。
*   **可复现性**: 使用 Conda 进行严格的环境管理，确保在任何机器上都能得到一致的结果。
*   **灵活性**: 支持WGS和WGBS两种数据类型，并允许每个样本自定义测序接头。
*   **全面的片段组学分析**:
    *   **质量控制**: 使用 `fastp` 进行质控和接头去除。
    *   **序列比对**: 使用 `BWA-MEM` (WGS) 或 `bwameth` (WGBS) 进行比对。
    *   **BAM文件处理**: 排序、标记PCR重复 (`Picard MarkDuplicates`)。
    *   **窗口保护分数 (WPS)**: 计算基因组特定区域的片段保护情况。
    *   **TSS 富集分析**: 评估转录起始位点周围的片段富集模式。
    *   **CTCF 富集分析**: 评估CTCF结合位点周围的片段富集模式。
    *   **末端基序分析**: 分析cfDNA片段末端的核苷酸偏好性。
*   **清晰的结果结构**: 所有输出文件都按样本整齐地组织在各自的目录中。

---

## 流程图

![alt text](fragpipe_dag.svg)

---

## 部署与安装

请按照以下步骤部署和准备 FragPipe 运行环境。

### 第1步：克隆本仓库

首先，将本仓库克隆到您的本地机器。

```bash
git clone https://github.com/your_username/FragPipe.git
cd FragPipe
```

### 第2步：安装 Conda

本流程依赖 Conda 进行环境管理。如果您尚未安装 Conda，我们推荐安装 [Miniconda](https://docs.conda.io/en/latest/miniconda.html)（一个轻量级的 Conda 发行版）。

为了加快环境创建速度，强烈建议在安装完 Conda 后，安装 Mamba：

```bash
conda install -n base -c conda-forge mamba
```

### 第3步：创建 Conda 环境

本仓库提供了一个 `environment.yaml` 文件，其中定义了所有必需的软件依赖。

```bash
# 推荐使用 Mamba (速度更快)
mamba env create -f workflow/environment.yaml

# 或者使用 Conda (速度较慢)
# conda env create -f workflow/environment.yaml
```

这将会创建一个名为 `cfDNA` 的 Conda 环境。请在运行流程前始终激活此环境：

```bash
conda activate FragPipe-env
```

### 第4步：下载并准备参考基因组

我们提供了一个便捷的脚本来自动下载人类参考基因组 (hg38) 并为其建立 BWA 索引和bwa-meth。

```bash
# 1. 赋予脚本执行权限
chmod +x download_hg38_build_index.sh

# 2. 运行脚本 (请确保已激活 Conda 环境)
# 这个过程会下载约 1GB 的文件，并需要约 30-60 分钟来建立索引
./download_hg38_build_index.sh
```

---

## 使用方法

### 第1步：准备输入文件

1.  **将您的原始FASTQ文件** 放置在您选择的任意目录下。
2.  **编辑 `Projects/Myproject1/samples.tsv` 文件**。这是一个**无表头**的、用**Tab**分隔的表格，用于定义您的样本信息。请根据您的数据修改此文件，或创建一个新的项目目录。

`samples.tsv` 文件必须包含 **7** 列：

| 列号 | 内容          | 描述                                     | 示例                                     |
| ---- | ------------- | ---------------------------------------- | ---------------------------------------- |
| 1    | `sample_id`   | 样本的唯一ID，将用于所有输出文件名       | `Patient1_Tumor`                         |
| 2    | `fastq1_path` | Read 1 FASTQ 文件的绝对路径              | `/path/to/data/P1_T1_R1.fq.gz`           |
| 3    | `fastq2_path` | Read 2 FASTQ 文件的绝对路径              | `/path/to/data/P1_T1_R2.fq.gz`           |
| 4    | `group`       | 样本所属的分组（例如 Tumor, Normal）     | `Tumor`                                  |
| 5    | `type`        | 分析类型，必须是 `wgs` 或 `wgbs`         | `wgs`                                    |
| 6    | `fwd_adapter` | Read 1 的测序接头序列                    | `AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA`        |
| 7    | `rev_adapter` | Read 2 的测序接头序列                    | `AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG` |

### 第2步：运行流程

在项目根目录 (`FragPipe/`) 下，执行以下命令来启动分析流程。

```bash
# 激活环境
conda activate FragPipe-env

# 运行 Snakemake
snakemake \
    --snakefile workflow/snakefile \
    --cores 8 \
    --config project_name=Myproject1 \
    --use-conda
```

**参数说明**:
*   `--snakefile`: 指定流程定义文件。
*   `--cores`: 您希望使用的最大CPU核心数。
*   `--config project_name=...`: **(必需)** 指定您在 `Projects/` 目录下准备的项目名称。
*   `--use-conda`: 指示 Snakemake 使用 `environment.yaml` 来管理运行环境。

---

## 输出文件结构

分析完成后，所有的结果都会被整齐地存放在 `Results/<YourProjectName>/` 目录下，每个样本都有一个专属的文件夹。

```
Results/
└── Myproject1/
    └── sample_C_20k/
        ├── alignment/
        │   ├── sample_C_20k.markdup.bam      # 最终的比对文件 (坐标排序)
        │   └── sample_C_20k.markdup.bam.bai  # BAM 索引
        ├── ctcf/
        │   └── sample_C_20k_..._mean.txt     # CTCF 富集分析结果
        ├── motif/
        │   ├── sample_C_20k.4mer.motif.MDS # 4-mer 基序多样性得分
        │   └── sample_C_20k.size_ratio     # 短/长片段比例
        ├── qc/
        │   ├── sample_C_20k.fastp.html     # Fastp 质控报告
        │   └── sample_C_20k.mkdup_metrics.txt # 重复率统计
        ├── tss/
        │   └── sample_C_20k_TSS_mean.txt   # TSS 富集分析结果
        └── wps/
            └── sample_C_20k.WPS_Score.txt  # 窗口保护分数
```

---

## 许可证

本项目采用 [MIT License](LICENSE.md) 授权。

## 联系方式

lishaogang - [li.shaogang97@gmail.com
](mailto:li.shaogang97@gmail.com
)

项目链接: [https://github.com/AlfredLSG/FragPipe](https://github.com/AlfredLSG/FragPipe)

