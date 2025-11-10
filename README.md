This is the code written for my OFYP project, <ins>Gene flow analysis using high-quality genome assemblies reveals porosity of species boundaries between insect sister-species</ins>. Paper can be found [here]([https://sites.google.com/view/barralab/home](https://dr.ntu.edu.sg/entities/publication/df41c785-677b-4c60-929d-430b05988696)).
The workflow is optimized for high-quality genome assemblies that applies Ancestral Recombination Graphs (ARGs) to detect gene flow.
Research was supervised by Professor Tim Barraclough, with lots of help from Tymek Pieszko[@Barralab](https://sites.google.com/view/barralab/home)


## Data used
The data used in this study were collected, sequenced and assembled by the darwin tree of life project, taken from [Darwin Tree of Life](https://www.darwintreeoflife.org/data/). Hence, all codes and analysis described are designed for whole-genome assembly data, such as the ones provided by the site.

## Template.pipeline.slurm
This template script is used to run through all the steps of pre-alignment processing, genome alignment, variant calling and further processing to prepare the files for gene flow analysis. The pipeline has been optimised where possible to minimise time taken and computational costs. For reference genomes ~200Mb (using whole-genome query files), it takes roughly 5 days to complete. For reference genomes ~60Mb, it takes 1-2 days to complete. Smaller references genomes can be completed within a day.

### Dependencies
The tools needed to run this script are:<br/>
minimap2/2.24<br/>
HTSlib/1.14<br/>
SAMtools/1.16.1<br/>
R/4.4.0<br/>
BCFtools/1.14<br/>

Further dependencies were installed using a conda environment. Additional packages needed(and corresponding dependencies auto downloaded by conda) are:<br/>
seqkit/2.8.2-1<br/>
BBMap/39.09-0<br/>
Since conda environment paths and names are created by you, remember to add the command to activate it in the script:
```
module load Anaconda3/2023.09-0
source activate <name of conda environment>
```

The packages required by the R script used are:<br/>
dplyr/1.1.4<br/>
VcfR/1.15.0<br/>
data.table/1.15.0<br/>

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

Additional scripts needed to run this script are split.sh, chr1.sh(if you are only using the first chromosome of the reference genome) and make_bedfiles.R. This is the folder that all these files are stored in. All additional scripts needed are provided in this repository.

### Output
A directory named <Species_name>.phased.vcfs will be created where you specified VAR_DIR, along with the associated bed files used to extract those vcfs stored in <Species_name>.bedfiles. The vcf files stored in <Species_name>.phased.vcfs are the vcf files required to run the next step.

## Template.singer.slurm
This script randomly samples a specified _vcf_num_ number of vcf files for gene flow analysis. It constructs Ancestral Recombination Graphs(ARGs), converts them into tree sequence files, and then classifies them into corresponding tree categories.

### Dependencies
A virtual python environment, containing the packages:<br/>
numPY<br/>
tskit/0.60<br/>
Please see how to create a virtual python environment using [Venv](https://docs.python.org/3/library/venv.html). Remember to add the command to activate it in the script:
```
module load Python/3.9.6
source <path_to_virtual_env>/bin/activate
```
This script generates ARGs from vcf files using the bioinformatics tool [SINGER](https://github.com/popgenmethods/SINGER). Please see how to install the tool on its github page. We have only tested the beta 1.17 version.

### Using the script
The script is designed to be run over different datasets with minimal editing. This script is designed for a SLURM scheduling system, so adjust accordingly. When running the script over a different species, only the following fields need to be edited(most fields are at the top of the script):

>VAR_DIR=

This is the directory path used to store your vcf files after the variant calling step. It is the same directory as the one used in Template.pipeline.slurm

>TREE_DIR=

This is the directory path used to store your tree sequence(.tskit) files.

>species=

This is the name of the species studied.

>ratio=

This is the ratio value used

>end=

This is the length of the reference genome

>vcf_num=

This is the number of randomly sampled vcf files to run SINGER over. We specified vcf_num as 100 for my project.

>SINGER_DIR=

This is the directory you installed SINGER in. The script currently only uses beta 0.1.7.

>mutation rate

This can be found in this line, after the -m flag, set as default 2.9e-9. It was not added to the editable field at the top of the script due to formatting issues. It is easiest specified using scientific notation.
```
"$SINGER_DIR/singer-0.1.7/singer_master" -m 2.9e-9 -ratio "$ratio" -vcf "$line" -output "$VAR_DIR/${species}.arg/${filename}.arg" -start 0 -end "$end" -n 100 -thin 20
```

> WOR_DIR=

An additional script needed to run this script is tree_type.py. This is the folder that this file is stored in. All additional scripts needed are provided in this repository. tree_type.py was provided by Tymek Pieszko

### Output
A file named <Species_name>.tree.txt will be created where you specified TREE_DIR. It is a tab delimited file containing information on the relevant trees and their relative proportions. 

## Running simulations
This study uses [msprime](https://tskit.dev/msprime/docs/stable/intro.html) to run simulations on Jupyter Notebook. The code used to generate simulations is
```
twopopmodel = msprime.Demography()
twopopmodel.add_population(name="A", initial_size=<Ne of A>)
twopopmodel.add_population(name="B", initial_size=<Ne of B>)
twopopmodel.add_population(name="C", initial_size=<Ne of ancestral pop. C>)
twopopmodel.add_population_split(time=<no. of generations>, derived=["A", "B"], ancestral="C")
twopopmodel.set_migration_rate(source='A', dest='B', rate=<gene flow rate per generation>)
Andrena_ts = msprime.sim_ancestry(
    recombination_rate=<recomb_rate>,
    sequence_length=<seq_length>,
    samples={"A": 1, "B": 1},
    demography=<name of demographic population>)
<species>_mts = msprime.sim_mutations(<species>_ts, rate=<mutation rate>)
```
The code can be reused for other species by adjusting parameter values. Simulated ancestries can then be used to generate a vcf file
```
with open("/path/simulated.vcf", "w") as vcf_file:
    <species>_mts.write_vcf(vcf_file, allow_position_zero=True)
```
A modified varation of Template.singer.slurm can then be used to run this vcf file through the gene flow analysis pipeline.

