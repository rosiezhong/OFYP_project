## Data used
The data used in this study were collected, sequenced and assembled by the darwin tree of life project, taken from [Darwin Tree of Life](https://www.darwintreeoflife.org/data/). Hence, all codes and analysis described are designed for whole-genome assembly data, such as the ones provided by the site.

## Template.pipeline.slurm
This template script is used to run through all the steps of pre-alignment processing, genome alignment, variant calling and further processing to prepare the files for gene flow analysis. The script has been optimised where possible. For reference genomes ~200Mb (using whole-genome input files), it takes roughly 5 days to complete. For reference genomes ~60Mb, it takes 1-2 days to complete. Smaller references genomes can be completed within a day.

### Code dependencies
The tools needed to run this script are
minimap2/2.24
HTSlib/1.14
SAMtools/1.16.1
R/4.4.0
BCFtools/1.14

Further dependencies were installed using a conda environment. Additional packages needed(and corresponding dependencies auto downloaded by conda) are
seqkit/2.8.2-1
BBMap/39.09-0

The packages required by the R script used are
dplyr/1.1.4
VcfR/1.15.0
data.table/1.15.0

### Using the script
The script is designed to be run over different datasets with minimal editing. This script is designed for a SLURM scheduling system, so adjust accordingly. When running the script over a different species, only the following fields need to be edited(all fields are at the top of the script):

> species_dir=

This is the path of the directory where your raw data is stored, in the format of a gzipped fasta file. There should be 4 data files in total, corresponding to 2 haplotypes of each of the 2 species. The file names must be in the format <Species_name>1.1.fasta.gz(i.e. Haplotype 1 of Species 1), <Species_name>1.2.fasta.gz(i.e. Haplotype 2 of Species 1), <Species_name>2.1.fasta.gz, <Species_name>2.2.fasta.gz. The name of the directory should also be <Species_name>, case sensitive.

> ALN_DIR=

This is the directory path used to store your bam files after the alignment step.

> VAR_DIR=

This is the directory path used to store your vcf files after the variant calling step. It is also used to store the bedfiles used to further process the vcf files.

> SCR_DIR=

This is the directory you will submit your job from. This is primarily used to clean up intermediate bam files which take up a lot of storage space.

> WOR_DIR=

Additional scripts needed to run this script are split.sh, chr1.sh(if you are only using the first chromosome of the reference genome) and make_bedfiles.R. This is the folder that all these files are stored in.



