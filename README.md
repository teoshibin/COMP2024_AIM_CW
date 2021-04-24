# Function Optimization Benchmarking
AI method coursework function optimization algorithm MATLAB code and datasets

# Rules
1. set path that points to fgeneric as "pwd" e.g. `addpath(pwd)`
2. set `datapath = '../ALGONAME'`
3. set `maxfunevals = '5e4'`

# Usage
## Benchmarking
- uncomment and comment different algorithms in [exampleexperiment](Algorithms/exampleexperiment.m)
- generate data with [2010 BBOB framework](https://coco.gforge.inria.fr/doku.php?id=bbob-2010-downloads) by running [exampleexperiment](Algorithms/exampleexperiment.m) in MATLAB
- post-process data with [numbbo framework](https://github.com/numbbo/coco/)

	```Sh
		python -m cocopp [-o OUTPUT_FOLDERNAME] YOURDATAFOLDER [MORE_DATAFOLDERS]
	```

## Generate delta ftarget values Excel from info log
### Installation
1. navigate into the repository

	```Sh
		cd [REPO]
	```
2. create virtual environment
	```Sh
		python -m venv info2excel\venv
	```
3. activate virtual environment
	```Sh
		info2excel\venv\Scripts\activate
		# or
		info2excel\venv\Scripts\activate.bat
	```
4. install dependencies
	```Sh
		pip install -r info2excel\requirements.txt
	```
### Execution
1. activate virtual environment (step 3 of installation)

3. run python script
	```Sh
		info2excel\info2excel.py -i [DATASET]
		# or
		info2excel\info2excel.py -i [DATASET] -d[DIMENSION] -o [EXCELNAME]
		# optional d & o, default d = 5, o = 'output'
	```
	Example command
	```Sh
		info2excel\info2excel.py -i JADE
		info2excel\info2excel.py -i JADE -d 40
		info2excel\info2excel.py -i JADE -o JADE_FSMAP
	```

Output excel [here](output.xlsx)
