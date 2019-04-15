from extractData import *
from analyzeData import *

generateFile('gauss-seidel_bench.txt', 'gauss-seidel_traite.txt')
plotGraph3D('gauss-seidel_traite.txt', 'gauss-seidel3D.png', 'Gauss-Seidel')

generateFile('jacobi_bench.txt', 'jacobi_traite.txt')
plotGraph3D('jacobi_traite.txt', 'jacobi3D.png', 'Jacobi')

generateFile('sor_bench.txt', 'sor_traite.txt')
plotGraph3D('sor_traite.txt', 'sor3D.png', 'SOR')
