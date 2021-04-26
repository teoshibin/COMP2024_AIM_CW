# Function Optimization Benchmarking
AI method coursework function optimization algorithm MATLAB code and datasets

# Usage
## Benchmarking
- uncomment and comment different algorithms in [exampleexperiment](Algorithms/exampleexperiment.m)
- generate data with [2010 BBOB framework](https://coco.gforge.inria.fr/doku.php?id=bbob-2010-downloads) by running [exampleexperiment](Algorithms/exampleexperiment.m) in MATLAB
	1. set path that points to fgeneric as "pwd" e.g. `addpath(pwd)`
	2. set `datapath = '../Datasets/ALGONAME'`
	3. set `maxfunevals = '5e4'`
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
[see info2excel documentation](info2excel_DOC.md)
