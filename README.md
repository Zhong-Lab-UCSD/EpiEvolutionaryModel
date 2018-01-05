# EpiEvolutionaryModel
Development of a joint evolutionary model for the genome and the epigenome

## Introduction
Interspecies epigenome comparisons yielded functional information that cannot be revealed by genome comparison alone, begging for theoretical advances that enable principled analysis approaches. Whereas probabilistic genome evolution models provided theoretical foundation to comparative genomics studies, it remains challenging to extend DNA evolution models to epigenomes. Here we present an effort to develop ab initio evolution models for epigenomes, by explicitly expressing the joint probability of multispecies DNA sequences and histone modifications on homologous genomic regions. We developed four models representing four evolutionary hypotheses, namely dependence and independence of interspecies epigenomic variations to sequence mutations and to sequence insertions and deletions (indels). 

<center>

|               | Indel-independent          | Indel-dependent  |
| :-----------: |:--------------------------:| :---------------:|
| **Mutation-independent**     |  Model N    |       Model I    |
| **Mutation-dependent**       |  Model M    |       ModelB     |

</center>

For model fitting, we implemented a maximum likelihood method by coupling downhill simplex algorithm with dynamic programming. This repository contains codes for parameter estimation and likelihood calculation, as well as programs for downstream data analyses.

## Generating inputs for maximum likelihood estimation
The maximum likelihood estimation (MLE) programs take fastq-like files as inputs. The inputs can be generated with `InSilicoEvolution*.py` for simulation data or `InputToFastq_bed2.py` with actual ChIP-Seq peak calling results and a list of genomic regions.

### 1. Simulation data
Simulation datasets can be generated based on corresponding hypotheses of the four models and used to evaluate the performance of the MLE algorithm. `InSilicoEvolution.py` should be used for Model I and Model B, whereas `InSilicoEvolution_indel_v2.0.py` should be used for Model N and Model M. 

#### Options
**`input`**

Name of the paramter file. The first and second lines are values for parameters _s_ and _μ_, followed by H lines containing H values for parameters _κ_<sub>1</sub>, _κ_<sub>2</sub>, ..., _κ_<sub>H</sub>. Line (H+3) includes _π_<sub>A</sub>, _π_<sub>C</sub>, _π_<sub>G</sub>, _π_<sub>T</sub>, followed by H lines containing _π_<sub>0</sub> and _π_<sub>1</sub> for H histone modifications.

An example of the parameter file for one histone modificaion (H=1) can be found [here](https://github.com/Zhong-Lab-UCSD/EpiEvolutionaryModel/blob/master/Examples/parameter_ex.txt).

**`--seq_len`**

Length of the ancestral sequence. Default: 500.

**`--seq_num`**

Number of region pairs. Default: 1000.

**`--hypN`**

Use which hypothesis to generate the data. For `InSilicoEvolution.py`: hypN=1: Model I; hypN=2: Model B. For `InSilicoEvolution_indel_v2.0.py`: hypN=1: Model N; hypN=2: Model M.

**`-r\--ran_seed`**

Set random seed. The program will generate a random seed if not specified. You may consider setting a random seed when generating the simulation data so that you can reproduce it.

**`-o\--output`**

Output file name

**`-O\--refout`**

Optional. If specified, the program will output the reference alignment (correct answer).

#### Example

```
python InSiliconEvolution.py parameter.txt --seq_len=500 --seq_num=100 -o Simulation.out -O Simulation_ref.out
```

####Output

The output file is a fastq-like file. Each region pair consists of 4(H+1) lines: 2(H+1) for species 1 and 2(H+1) for species 2. For each species, the first line contains sequence name, and the second line contains genomic sequence. The 3 to 2(H+1) lines are histone modification signals for the H histone marks, seperated by "+" lines. An example of the output file can be found [here](https://github.com/Zhong-Lab-UCSD/EpiEvolutionaryModel/blob/master/Examples/Simulation.txt).


### 2. Actual data

`InputToFastq_bed2.py` can be used to generate the fastq-like files using actual ChIP-Seq peak calling results and a list of homologous region pairs. Peak regions need to be identified in advance with peak calling software, and homologous region pairs can be found using tools such as liftOver.

#### Options

**`q_region`**

Names of the two input files. The two files should contain coordinates of homologous query regions in two species in bed6 format. Note that region names of the same species should NOT be the same. The order of the two files should be the same as species specified in `--species`. 

**`-s\--species`**

Names of genome assemblies. Default: `["hg38","rheMac8"]`.

**`-f\--fasta_path`**

Path to the genome sequence files (.fa files). Note that the fasta file names in this folder should be the same as specified in `--species`.	

**`--histone`**

Name of histone modifications. Note that the list should have the same order as peak files. Default: `["H3K4me3"]`.

**`--bg`**

This file should contain paths to ChIP-Seq peak calling files (.bed files). An example can be found [here](https://github.com/Zhong-Lab-UCSD/EpiEvolutionaryModel/blob/master/Examples/hg38_rhe8_H3K4me3_peak.txt).

**`--s_path`**

Path to samtools. Default: samtools.

**`-p\--p_num`**

Number of processes to be used. Default: 5.

**`-o\--output`**

Output file name. 

#### Example

```
python InputToFastq_bed2.py hg38_region.bed rheMac8_region.bed --bg hg38_rhe8_H3K4me3_peak.txt -p 10 -o hg38_rheMac8_MLEInput.out
```

#### Output

The output file has the same format as that of `InSilicoEvolution*.py`

## Parameter estimation and likelihood calculation

After obtaining the inputs, model parameters can be estimated using `parameter_estimation_epi_MLE_simp*.py`. The programs have two modes: parameter estimation and likelihood calculation. For the first mode, initial guesses of the parameters need to be provided, and the MLE program will iterate to estimate parameters. For the second mode, estimated parameters have to be provided, and the program will not use the MLE algorithm, but simply compute the likelihoods of input regions. 

#### Options

**`Input`**

Input file name. The input here is the output of `InSilicoEvolution*.py` or `InputToFastq_bed2.py`. 

**`-n\--iter_non`**

If specified, the program will enter the likelihood calculation mode and only report the first likelihood. It can be used to calculate the likelihood when parameters are known. Otherwise the program will enter the parameter estimation mode and use the MLE algorithm for the purpose.

**`-e\--equil_file`**

If `--iter_non` is specified, the parameter file should contain estimated parameters and equilibrium probabilities. Otherwies, it should contain intial guesses for s, mu, k, and equilibrium probabilities estimated from the input data.

**`--hypN`**

Use which hypothesis to generate the data. For `parameter_estimation_epi_MLE_simp.py`: hypN=1: Model I; hypN=2: Model B. For `parameter_estimation_epi_MLE_simp_indel.py`: hypN=1: Model N; hypN=2: Model M.

**`-p\--p_num`**

Number of processes to be used. Default: 6.

**`-o\--output`**

Output file name. If -n is specified, the program will output region name and -log likelihood for each region. Otherwise it will output the -log likelihood in each iteration.

**`-O\--out_allvec`**

Output file name. The program will output parameters updated in each iteration. Please note that the algorithm may evaluate several sets parameters in each iteration to select the best one. These parameter sets will be printed to stdout.

#### Examples 

```
python parameter_estimation_epi_MLE_simp.py hg38_rheMac8_MLEInput.out -e initial_guesses.txt -p 50 --hypN=1 -o likelihood.txt -O allvec.txt > allparameters.txt
```

## Downstream analysis

A variety of downstream analyses can be performed after calculating likelihoods based on the four models. Here we represent codes used for two analyses: (1) likelihood comparison and (2) separating evolutionary impacts of sequence mutations and indels. Please see our manuscript for other analyses done in our work.

### Likelihood comparison

`Classification_parser_multiModel.R` can be used to do the likelihood comparison. It will assign each region to one of the four groups (Model I, Model B, Model N, Model M) based on the likelihoods.

#### Options

**`Input`**

The input file should contain 5 tab-separated columns. The first column is regionname, and the rest four are -log likelihoods yielded by the four model. The order of the columns can be specified by `--model_name`.

**`-n\--model_name`**

Order of the models. Default: `c("MI", "MB", "MN", "MM")`

**`-c\--split_name`**

If true, region names will be split based on the second underscore.

**`-s\--species`**

Species names. Default: `c("human","rhesus")`

**`-sp1`**

Region coordinates in species1. This is the input bed file for `InputToFastq_bed2.py`. 

**`-sp2`**

Region coordinates in species2. 

**`-o\--output`**

Output folder name.

#### Output

The program will create a folder containing five output files. 

**region_info.txt**: This file has 8 columns. Column 1: region name. Column 2-5: likelihoods yielded based on the four models. Column 6: trimmed region name. Column 7: group. Column 8: normalized -log likelihood difference between the largest and the smallest -log likelihoods.

**Type\*.txt**: Each file contains regions assigned to the corresponding group. Column 1: trimmed region name. Column 2-7: same as column 2-7 of region_info.txt. The last 6 columns are coordinates of the region in species 1 and 2.

### Separating evolutionary impacts of sequence mutations and indels

`Separate_SeqVar.R` can be used to separate impacts of sequence mutations and indels in each region.

#### Options

**`Input`**

The input file should contain at least 5 columns. The first column is regionname, and the rest four are -log likelihoods yielded by the four model. A header is required.

**`-o\--output`**

Output file name. 

#### Output

The output file has 9 columns. 

* Column 1: region name. 
* Column 2: negative log of the probability that epigenomic variation in the region is independent to mutations.
* Column 3: negative log of the probability that epigenomic variation in the region depends on mutations.
* Column 4: negative log of the probability that epigenomic variation in the region is independent on indels.
* Column 5: negative log of the probability that epigenomic variation in the region depends on indels.
* Column 6: group of the region
* Column 7: normalized difference between the largest and the smallest -log probabilities.
* Column 8: normliazed difference between Column 3 and 2. The extent of dependence on mutations increases as the value increases.
* Column 9: normalized difference between Column 5 and 4. The extent of dependence on indels increases as the value increases.