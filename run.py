from sympy import Matrix
from func import diagonalize

# mat = Matrix([[1, 2],
#               [0, 3]])

mat = Matrix([[1, 2, 1],
              [6, -1, 0],
              [-1, -2, -1]])

diagonalize(mat)
