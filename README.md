# Function Optimization Benchmarking
AI method coursework function optimization algorithm MATLAB code and datasets

# Rules
1. set path that points to fgeneric as "pwd" e.g. `addpath(pwd)`
2. set `datapath = '../Datasets/ALGONAME'`
3. set `maxfunevals = '5e4'`

# Usage
## Benchmarking
- uncomment and comment different algorithms in [exampleexperiment](Algorithms/exampleexperiment.m)
- generate data with [2010 BBOB framework](https://coco.gforge.inria.fr/doku.php?id=bbob-2010-downloads) by running [exampleexperiment](Algorithms/exampleexperiment.m) in MATLAB
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
3. activate virtual environment (a prefix of "(venv)" will be displayed in your command prompt when succesfully activated)
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
		info2excel\info2excel.py -i [DATASET] 					# or
		info2excel\info2excel.py -i [DATASET] -d [DIMENSION] -o [EXCELNAME]
		# optional d & o by default is 5 & "[ALGONAME]_[DIMENSION]D"
	```
	Example command
	```Sh
		info2excel\info2excel.py -i JADE 		#or
		info2excel\info2excel.py -i JADE -d 40		#or
		info2excel\info2excel.py -i JADE -o JADE_FSMAP
	```
	Complete command flow
	```Sh
		cd info2excel
		info2excel.py -i PSO_Bounds
		# output : PSO_Bounds_5D.xlsx
	```
