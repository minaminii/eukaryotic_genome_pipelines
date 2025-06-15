# Eukaryotic Genome Assembly and Annotation Pipeline

This pipeline performs long read-based **Eukaryotic genome assembly** using either **Canu** or **wtdbg2** as the primary assembler, followed by multiple rounds of polishing with **Racon** (long-read based) and **Pilon** (short-read based). It also includes a modular script for genome annotation following the guide from [meiyang12/Genome-annotation-pipeline](https://github.com/meiyang12/Genome-annotation-pipeline).

---

## 🔬 Features

- **Primary Assembly**: Choose between `Canu` or `wtdbg2`
- **Long-read Polishing**: Performed with `Racon`
- **Short-read Polishing**: Final polishing step with `Pilon`
- **Annotation**: Script for genome annotation using tools defined in the annotation guide
- **Environment Management**: Conda YAML file included to reproduce the environment

---

## ⚙️ Dependencies

The pipeline requires the following tools (included in the conda environment):

- [Canu](https://github.com/marbl/canu)
- [wtdbg2](https://github.com/ruanjue/wtdbg2)
- [Racon](https://github.com/lbcb-sci/racon)
- [Pilon](https://github.com/broadinstitute/pilon)
- [Minimap2](https://github.com/lh3/minimap2)
- [SAMtools](http://www.htslib.org/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [Java](https://www.java.com/)
- Annotation dependencies from [meiyang12/Genome-annotation-pipeline](https://github.com/meiyang12/Genome-annotation-pipeline)

> [!IMPORTANT]
> This pipeline is a wrapper that orchestrates the execution of multiple external tools. 
>
> **It does not modify, package, or alter the behavior of those tools.**
>
> If you encounter issues related to the functionality, dependencies, or performance of any external tool (e.g., Canu, Racon, Pilon), please refer to the documentation or issue tracker of the respective project.


---

## Installation 

### Conda Environment Setup

Create and activate the conda environment:

```bash
conda env create -f environment.yaml
conda activate genome_assembly_env
```

## Usage

### 1. Genome Assembly

```bash
bash assembly.sh -l <long_reads.fq> -s <short_reads_R1.fq> <short_reads_R2.fq> -a canu -o output_directory
```

#### Parameters:

- `-l`: Long read file (ONT/PacBio FASTQ)
- `-s`: Short read pair files (Illumina FASTQ)
- `-a`: Assembler choice: canu or wtdbg2
- `-o`: Output directory

The script will:

1. Run assembly using selected assembler.
2. Polish the assembly using Racon with long reads.
3. Perform additional polishing using Pilon with short reads.

### 2. Genome Annotation

```bash
bash annotation.sh -i assembled_genome.fasta -s species_name -t threads
```

This script wraps the annotation process using the methodology outlined in [meiyang12/Genome-annotation-pipeline](https://github.com/meiyang12/Genome-annotation-pipeline).

## Output

output_directory/
- assembly.fasta: Final polished assembly
- intermediate/: All intermediate assemblies and polishing steps
- annotation/: GFF, protein, CDS files from annotation

## Contact

For questions or suggestions, feel free to open an issue or submit a pull request. 

Please note that this repository does not provide support for tool-specific bugs or installation problems.

## Third-Party Tools and Licenses

This project is licensed under the MIT License.

This project makes use of several external bioinformatics tools. Full credit goes to the original authors and maintainers of each tool. These tools are not bundled with this repository, and must be installed separately by the user. Each tool is governed by its own license, which remains fully applicable.

| Tool         | License     | Link                                   |
|--------------|-------------|----------------------------------------|
| Canu         | BSD-3-Clause| https://github.com/marbl/canu          |
| wtdbg2       | MIT         | https://github.com/ruanjue/wtdbg2      |
| Racon        | MIT         | https://github.com/lbcb-sci/racon      |
| Pilon        | GPLv3       | https://github.com/broadinstitute/pilon|
| Minimap2     | MIT         | https://github.com/lh3/minimap2        |
| BWA          | MIT         | http://bio-bwa.sourceforge.net/        |
| SAMtools     | MIT/BSD     | http://www.htslib.org/                 |
| Annotation tools (MAKER, etc.) | Various     | See respective documentation        |
