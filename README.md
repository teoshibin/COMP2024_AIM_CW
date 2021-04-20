# COMP2039_AIM_CW
AI method coursework function optimization algorithm MATLAB code and datasets

# Usage
- uncomment and comment different algorithms in [exampleexperiment](Algorithms/exampleexperiment.m)
- generate data with [2010 BBOB framework](https://coco.gforge.inria.fr/doku.php?id=bbob-2010-downloads) by running [exampleexperiment](Algorithms/exampleexperiment.m) in MATLAB
- post-process data with [numbbo framework](https://github.com/numbbo/coco/)

	```Sh
		python -m cocopp [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
	```

# Rules
1. set path that points to fgeneric as "pwd" e.g. `addpath(pwd)`
2. set `datapath = '../ALGONAME'`
3. set `maxfunevals = '5e4'`
