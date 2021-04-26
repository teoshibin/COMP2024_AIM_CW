
### Credits ###
# by Shi Bin Teo
# https://github.com/teoshibin

### Imports ###

import os
from numpy import double
import pandas as pd
import sys, getopt
from shutil import copyfile
from openpyxl import load_workbook

### Main Functions ###

# convert delta f opt from different files into 2d array
def infoTo2dlist(dataset_folder, name, extensions, number_benchmarks):
    # os.path.join(os.path.dirname(os.getcwd()), dataset_folder)
    paths = getPaths(dataset_folder, name, extensions, number_benchmarks)
    dimensions = getDimensions(paths[0]);
    line_to_read = [x*3+2 for x in range(0, dimensions)]
    data = []

    # get data
    for i in range(0, number_benchmarks):
        a_file = open(paths[i])
        for position, line in enumerate(a_file):
            if position in line_to_read:
                data.append(line.rstrip())

    # process data (remove unnecessary data)
    for i, element in enumerate(data):
        element = element.split(', ')
        element.pop(0)
        for j, el in enumerate(element):
            element[j] = el.split('|', 1)[1]
        data[i] = element

    return data

# convert 2d list to pandas dataframe
def createDataframe(datalist):
    number_instance = len(datalist[0])
    number_benchmarks = len(datalist)
    mycolumns = ['Instance ' + str(x) for x in range(1, number_instance + 1)]
    myindex = ['f' + str(x) for x in range(1, number_benchmarks + 1)]
    df = pd.DataFrame(datalist, columns=mycolumns, index=myindex)
    return df

# select certain dimension
def querySpecificDimension(data, number_benchmarks, selected_dimension_id):
    new_data = []
    number_dimensions = len(data) / number_benchmarks
    benchmark_count = 0
    for i, instances in enumerate(data):
        if i == selected_dimension_id + number_dimensions*benchmark_count:
            new_data.append(instances)
        if (i+1) % number_dimensions == 0:
            benchmark_count += 1
    return new_data

### Helper Functions ###

def getPaths(folder, name, extensions, number_benchmarks):
    paths = []
    for i in range(1, number_benchmarks + 1):
        paths.append(os.path.join(folder, name + str(i) + extensions))
    return paths


def getDimensions(path):
    return  int(countLine(path) / 3)


def countLine(path):
    file = open(path, "r")
    line_count = 0
    for line in file:
        if line != "\n":
            line_count += 1
    file.close()
    return line_count


def main(argv):

    ### Main Parameters ###
    ## benchmarks ##
    dimensions = [2, 3, 5, 10 , 20, 40]
    number_benchmarks = 24
    dimension = 5
    precision = 1e-8

    ## preset IO folder name ##
    dir = os.path.dirname(__file__)
    dataset_folder = os.path.join(dir, os.pardir, 'Datasets')
    output_folder = os.path.join(dir, os.pardir, 'ExcelScore')
    algorithm_name = ''
    excelname = ''

    ### Main script ###
    ## error msg ##
    wrong_syntax_msg = ('info2excel.py -i [ALGONAME] -d [DIMENSION] -p [precision] -o [EXCELNAME]\n'
                        'info2excel.py -i [ALGONAME]\n'
                        '-d & o by default is 5 & "[ALGONAME]_[DIMENSION]D"\n')
    wrong_dimensions_msg = '--dimension not found in' + str(dimensions)

    ## parse in args ##
    try:
        opts, args = getopt.getopt(argv,"i:d:p:o:",["ifile=", "dimension=", "precision=", "ofile="])
    except getopt.GetoptError:
        print(wrong_syntax_msg)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            algorithm_name = arg
        elif opt in ("-d", "--dimension"):
            dimension = abs(int(arg))
        elif opt in ("-p", "--precision"):
            precision = abs(double(arg))
        elif opt in ("-o", "--ofile"):
            excelname = arg
    
    ## validation ##
    if len(algorithm_name) == 0:
        print(wrong_syntax_msg)
        sys.exit(2)

    if not( dimension in dimensions):
        print(wrong_dimensions_msg)
        sys.exit(2)

    if len(excelname) == 0:
        excelname = algorithm_name + '_' + str(dimension) + 'D'

    ## data retrieval ##
    full_dataset_path = os.path.join(dataset_folder, algorithm_name)
    full_output_path = os.path.join(output_folder, excelname + '.xlsx')

    data = infoTo2dlist(full_dataset_path, 'bbobexp_f', '.info', number_benchmarks)
    data = querySpecificDimension(data, number_benchmarks, dimensions.index(dimension))

    data = double(data) # convert string to double
    df = createDataframe(data)
    df = df.add(precision*2) # shift delta to 0 as global optimum delta value and removing 0 by adding precision again
    # df = df.abs()

    print(df)

    df.where(df <= 1000, 1000, inplace=True)
    # when df <= 1000 then do not replace
    # when df > 1000 then replace with 1000

    ## output value into exisiting excel template ##
    copyfile(os.path.join(dir, 'template.xlsx'), full_output_path)

    book = load_workbook(full_output_path)
    writer = pd.ExcelWriter(full_output_path, engine='openpyxl') 
    writer.book = book

    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    df.to_excel(writer, 'Sheet1', index=False, header=False, startcol=1, startrow=1)
    writer.save()

if __name__ == "__main__":
    main(sys.argv[1:])