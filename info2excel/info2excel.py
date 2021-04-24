
import os
from numpy import double
import pandas as pd
import sys, getopt

# convert delta f opt from different files into 2d array
def infoTo2dlist(dataset_folder, name, extensions, number_benchmarks):
    os.path.join(os.path.dirname(os.getcwd()), dataset_folder)
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

# export 2d list to excel
def createDataframe(datalist, number_benchmarks):
    number_instance = len(datalist[0])
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

    dimensions = [2, 3, 5, 10 , 20, 40]
    number_benchmarks = 24
    datasetFile = ''
    excelname = 'output'
    dimension = 5

    try:
        opts, args = getopt.getopt(argv,"i:d:o:",["ifile=", "--dimension", "ofile="])
    except getopt.GetoptError:
        print('info2excel.py -i [FOLDERNAME] -d [DIMENSION] -o [EXCELFILENAME]\ninfo2excel.py -i [FOLDERNAME]')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--ifile"):
            datasetFile = arg
        elif opt in ("-d", "--dimension"):
            dimension = int(arg)
        elif opt in ("-o", "--ofile"):
            excelname = arg
    
    if len(datasetFile) == 0:
        print('info2excel.py -i [FOLDERNAME] -d [DIMENSION] -o [EXCELFILENAME]\ninfo2excel.py -i [FOLDERNAME]')
        sys.exit(2)

    data = infoTo2dlist(datasetFile, 'bbobexp_f', '.info', number_benchmarks)
    data = querySpecificDimension(data, number_benchmarks, dimensions.index(dimension))

    data = double(data) # convert string to double
    data = abs(data) # remove negative sign
    df = createDataframe(data, 24)
    print(df)
    df.where(df > 1000, 1000, inplace=True)

    df.to_excel(excelname + ".xlsx")

if __name__ == "__main__":
    main(sys.argv[1:])