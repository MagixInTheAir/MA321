import csv
import matplotlib.pyplot as plt


with open('bench1_gs_test.txt') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=';')
    line_count = 0
    rows = []
    index = []
    errors = []
    iterations = []
    for row in csv_reader:
        if line_count == 0:
            print('Column names are', row)
            line_count += 1
        else:
            rows.append(row)
            index.append(int(row[0]))
            errors.append(float(row[3]))
            iterations.append(int(row[7]))
            print('Pour la matrice ', row[0], 'L\'erreur est ', row[3])
            line_count += 1
    print('Processed', line_count, 'lines.')
    print(index, '\n', errors, '\n', iterations)
    print(len(index))
# print(rows)


import csv

with open('test1_traite', mode='w') as data_traitees:
    data_writer = csv.writer(data_traitees, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    newrows=[]
    for i in range(len(index)):
        newrow = [index[i], errors[i], iterations[i]]
        newrows.append(newrow)
        data_writer.writerow(newrows[i])

with open('test1_traite2', mode='w') as data_traitees2:
    data_writer = csv.writer(data_traitees2, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    data_writer.writerow(index)
    data_writer.writerow(errors)
    data_writer.writerow(iterations)



plt.subplot(211)
plt.title("Temps de résolution pour des matrices de taille 10")
plt.xlabel("Numéro de la matrice")
plt.ylabel("Temps de résolution (ns)")
plt.plot(index, errors)

plt.subplot(212)
plt.plot(index, iterations)
plt.xlabel("Numéro de la matrice")
plt.ylabel("Nombre d'itérations")

plt.show()
