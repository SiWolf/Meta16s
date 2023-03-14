[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Meta16s
>An all-in-one solution for streamlined 16s analyses.

![Meta16s Workflow](https://github.com/SiWolf/Meta16s/blob/master/workflow.png)

Meta16s is an all-in-one pipeline for the analysis of 16s NGS data. Meta16s takes raw paired-end NGS reads (.fastq.gz) as an input in order to generate a cleaned taxonomic abundance file (biom and tsv).

Meta16s requires [Python2](https://www.python.org/), [Conda](https://docs.conda.io/en/latest/) and [USEARCH](https://drive5.com/usearch/) in order to run. All other dependencies will be downloaded automatically.

# Quick Start:
* Clone the Meta16s repository
* Extract the USEARCH binary (usearch11.0.667_i86linux32) into the Meta16s folder
* Copy read files into the "sequences" folder
* Make sure the paired-end reads are named "<unique_name>_R1.fastq.gz" and "<unique_name>_R2.fastq.gz" respectively
* Activate the Meta16s conda environment (Meta16s.yml)
* Run Meta16s (Reference_Analyzer_16s.py)

# Licence
This project is licensed under the GPLv3 Licence. See the [LICENSE](LICENSE) file for more details.