# EpiEvolutionaryModel
Development of a joint evolutionary model for the genome and the epigenome

## Introduction
Interspecies epigenome comparisons yielded functional information that cannot be revealed by genome comparison alone, begging for theoretical advances that enable principled analysis approaches. Whereas probabilistic genome evolution models provided theoretical foundation to comparative genomics studies, it remains challenging to extend DNA evolution models to epigenomes. 

## Generating inputs for maximum likelihood estimation
The maximum likelihood estimation (MLE) programs take fastq-like files as inputs. The inputs can be generated with `InSilicoEvolution*.py` for simulation data or `InputToFastq_bed2.py` with actual ChIP-Seq peak calling results and a list of genomic regions.

### Simulation data
Simulation datasets can be generated based on corresponding hypotheses of the four models and used to evaluate the performance of the MLE algorithm. `InSilicoEvolution.py` should be used for Model I and Model B, whereas `InSilicoEvolution_indel_v2.0.py` should be used for Model N and Model M. 

#### Options
**input**

Name of the paramter file. Required parameter. The first and second lines are values for parameters _s_ and _μ_, followed by H lines containing H values for parameters _κ_<sub>1</sub>, _κ_<sub>2</sub>, ..., _κ_<sub>H</sub>. 

An example of the parameter file for one histone modificaion (H=1) can be found here.


 
