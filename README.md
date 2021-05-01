# Function Optimization Benchmarking
AI method coursework function optimization algorithm MATLAB code and datasets

# Usage
## Benchmarking
- uncomment and comment different algorithms in [exampleexperiment](Algorithms/exampleexperiment.m)
- generate data with [2010 BBOB framework](https://coco.gforge.inria.fr/doku.php?id=bbob-2010-downloads) by running [exampleexperiment](Algorithms/exampleexperiment.m) or dedicated experiment script in MATLAB
	1. set path that points to fgeneric as "pwd" e.g. `addpath(pwd)`
	2. set fgeneric `delta precision = 1e-16`
	3. set `datapath = '../Datasets/ALGONAME'`
	4. set `maxfunevals = '5e3'`
	5. set random generator seed `rng('default')`
- post-process data with [numbbo framework](https://github.com/numbbo/coco/)

	```Sh
		python -m cocopp [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
	```
	Example command
	```Sh
		cd Datasets
		python -m cocopp -o "..\Post_Processing\PSO_ppdata" "PSO"
	```
	> **Note**   
	> Path to datasets will be used as algorithm name in ppdata, this is why we cd into the datasets

## Generate Excel Score
[see info2excel documentation](info2excel)

## Results
### Post Processing
<p align="center" float="left">
  <img src="screenshots/1.jpg" width="716"/>
</p>
### Excel Scores

All tested algorithm with their corresponding scores in 5d

| Algo       | Score 5d | In Report |
|------------|---------:|:---------:|
| CMAES      | 2344     | 0         |
| DEAE       | 2322     | 1         |
| JADEb      | 1941     | 0         |
| JADE       | 1868     | 1         |
| JADEctpb   | 1868     | 0         |
| DE         | 1083     | 1         |
| DEb        | 1083     | 0         |
| PSO_Bounds | 500      | 1         |
| PSO        | 234      | 1         |
| Bin_GA     | -14      | 1         |
| PSO_EDA    | -49      | 1         |


