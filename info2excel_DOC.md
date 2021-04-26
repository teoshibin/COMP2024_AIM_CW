# info2excel
A small script written for AI method coursework, to speed up the process of copying delta f target values

## some issue
retrieving negative and zeros values   
```
precision = 1e-8   
global_best = foptimal = 50   
ftarget = foptimal + precision   
delta fitness = fbest - ftarget   
```
Therefore, if fbest reached the global minimum, say at `50`, `ftarget = 50 + 1e-8`   
```
delta_fitness = 50 - 50 - 1e-8
delta fitness = -1e-8  (+1e-8)
```
The best way is to +1e-8 but intoducing 0 into excel causing log 0.

### Solution   
method 1: add precision, abs, add precision   
method 2: replace negative values with precision (default)

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
		info2excel\venv\Scripts\activate	# or
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
		info2excel\info2excel.py -i [DATASET] 								# or
		info2excel\info2excel.py -i [DATASET] -m [METHOD] -d [DIMENSION] -p [PRECISION] -o [EXCELNAME]
		# d by default 5
		# o by default "[ALGONAME]_[DIMENSION]D"
		# p by default 1e-8
		# m by default 2 (see SOME ISSUE)
	```
	Example command
	```Sh
		info2excel\info2excel.py -i CMAES 		#or
		info2excel\info2excel.py -i CMAES -m 2 -d 5 -p 1e-8 -o CMAES_1e-8_5D
		# both output : CMAES_1e-8_5D.xlsx
	```
