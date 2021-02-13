from sympy import *
import numpy as np

init_printing(use_unicode=True)

n = 6

u = symarray("u",n+1)
x = symarray("x",n+1)
A, B = symbols("A B")
g, m, t, h = symbols("g m t h")



exp = [x[1]]
x1_poly_list = []
u_poly_list = []
for i in range(1,n):
    exp.append(exp[i-1] + h*(A*exp[i-1] + B*u[i]))

    # factors of x1
    x1_poly = Poly(exp[-1], x[1]).coeffs()[0]
    coeff = x1_poly.args
    ls = [0]*n
    for item in coeff:
        if item == 1:
            ls[0] = item
            continue
        try:
            a = item.args[-1]
            ls[a.args[1]] = item
        except IndexError:
            ls[1] = item   
    x1_poly_list.append(ls)

    # factors of u
    u_poly = [Poly(exp[-1], u[j]).coeffs()[0] for j in range(1,i+1)]
    u_poly_list.append(u_poly)

# pad rows
for ls in [x1_poly_list, u_poly_list]:
    # max row size
    length = max([len(item) for item in ls])
    for item in ls:
        while len(item) < length:
            item.append(0)

test = Matrix([[1,1],[1,1]]) 
D = Matrix(x1_poly_list)
E = Matrix(u_poly_list)
pprint(D,wrap_line  =False)
pprint(E,wrap_line  =False)
# pprint(expand(X1_exp))
# pprint(expand(X2_exp))
# pprint(expand(X3_exp))
# pprint(expand(X4_exp))
# pprint(expand(X5_exp))