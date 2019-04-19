import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D

def generateFile(fileToOpen, fileToGenerate):
    with open(fileToOpen) as csv_file:
        line_count = 0
        rows = []
        matrix_sizes = []
        precisions = []
        index = []
        errors = []
        times = []
        iterations = []
        while line_count <= 220:
            if line_count == 0:
                line = csv_file.readline()
                print('Column names are'.format(line_count, line.split(',')))
                line_count += 1
            else:
                line = csv_file.readline()
                row = line.split(',')
                print(row)
                rows.append(row)
                index.append(int(row[0]))
                matrix_sizes.append(int(row[1]))
                precisions.append(np.log10(float(row[2])))
                errors.append(float(row[3]))
                times.append(1e-6 * float(row[4]))
                iterations.append(int(row[5]))
                print('Pour la matrice ', row[0], 'L\'erreur est ', row[3])
                line_count += 1

        print('Processed', line_count, 'lines.')
        print(index, '\n', errors, '\n', matrix_sizes)  # , iterations)
        print(len(index))

    with open(fileToGenerate, mode='w') as data_traitees:
        data_writer = csv.writer(data_traitees, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        newrows = []
        print(len(index))
        for i in range(len(index)):
            newrow = [index[i], matrix_sizes[i], precisions[i], errors[i], times[i], iterations[i]]
            newrows.append(newrow)
            data_writer.writerow(newrows[i])

# with open('jacobi_traite2.txt', mode='w') as data_traitees2:
#     data_writer = csv.writer(data_traitees2, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
#     data_writer.writerow(index)
#     data_writer.writerow(matrix_sizes)
#     data_writer.writerow(errors)
#     data_writer.writerow(times)
#     data_writer.writerow(precisions)
#     data_writer.writerow(iterations)
