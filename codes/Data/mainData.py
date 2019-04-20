from extractData import *
from analyzeData import *
from prepareBars import *
from prepare3D import *

# generateFile('gauss-seidel_bench_min.txt', 'gauss-seidel_traite_min.txt')
# plotGraph3D('gauss-seidel_traite_min.txt', 'gauss-seidel3D.png', 'Gauss-Seidel')
# prepareBars('gauss-seidel_traite_min.txt', 'bars-gauss-seidel.txt')
# prepare3D('gauss-seidel_traite_min.txt', '3D-gauss-seidel.txt')
#
# generateFile('jacobi_bench_min.txt', 'jacobi_traite_min.txt')
# plotGraph3D('jacobi_traite.txt', 'jacobi3D.png', 'Jacobi')
# prepareBars('jacobi_traite_min.txt', 'bars-jacobi.txt')
# prepare3D('jacobi_traite_min.txt', '3D-jacobi.txt')
#
#
# generateFile('sor_bench_min.txt', 'sor_traite_min.txt')
# plotGraph3D('sor_traite_min.txt', 'sor3D.png', 'SOR')
# prepareBars('sor_traite_min.txt', 'bars-sor.txt')
# prepare3D('sor_traite_min.txt', '3D-sor.txt')


generateFile('gmres_bench_min_2.txt', 'gmres_traite_min_2.txt')
plotGraph3D('gmres_traite_min_2.txt', 'gmres_2.png', 'GMRES')
prepareBars('gmres_traite_min_2.txt', 'bars-gmres_2.txt')
prepare3D('gmres_traite_min_2.txt', '3D-gmres_2.txt')
