import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from mpl_toolkits.mplot3d import Axes3D

def plotGraph3D(fileGenerated, graphName, method_name):
    index = []
    matrix_sizes = []
    precisions = []
    errors = []
    times = []
    iterations = []
    with open(fileGenerated) as csv_file:
        for i, line in enumerate(csv_file):
            if (i % 2) == 0:
                line = line.split(',')
                index.append(int(line[0]))
                matrix_sizes.append(int(line[1]))
                precisions.append(float(line[2]))
                errors.append(float(line[3]))
                times.append(float(line[4]))
                iterations.append(float(line[5]))
                print(line)

    fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title(method_name)
    plt.xlabel('Taille de la matrice')
    plt.ylabel('Précision')
    ax.set_zlabel('Temps en ms / erreurs')
    ax.scatter(matrix_sizes, precisions, times, label='temps de résolution')
    ax.scatter(matrix_sizes, precisions, errors, label='erreur')
    ax.scatter(matrix_sizes, precisions, iterations, label='nb_itérations')
    plt.legend()
    plt.savefig(graphName)
    plt.show()

