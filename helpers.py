import matplotlib.pyplot as plt
from matplotlib import cm
from sympy import *
import numpy as np
from copy import deepcopy

plt.style.use('_mpl-gallery')

#Translates Dictionaries to the Equivalent Cost Matrices
def dict2array(Q):
    temp_x = max([cord[0] for cord in Q.keys()])
    temp_y = max([cord[1] for cord in Q.keys()])
    res = [[0] * (temp_y + 1) for ele in range(temp_x + 1)]
 
    for (i, j), val in Q.items():
        res[i][j] = val
    return res

#Tryout Function: Helps to explore the gradients (for determining Langrage Multipliers later)
def gradient(expression, symbol):
    '''x = symbols('x')
    expr = x + 1
    expr.subs(x, 2)'''
    e_q = diff(expression, symbol[0]) 
    e_y = diff(expression, symbol[1]) 
    e_z = diff(expression, symbol[2]) 
    e_x = diff(expression, symbol[3])

    return np.array([e_q, e_y, e_z, e_x])

#Tryout Function: Evaluates the gradient at all possible values (for determining Langrage Multipliers later)
def gradient_evaluation(gradient, symbol, num_symbols, val):
    gradient_copy=deepcopy(gradient)
    for i in range(num_symbols):
        symbol[i] = symbols(str(symbol[i]))
        for j in range(len(gradient)):
            gradient_copy[j] = gradient_copy[j].subs(symbol[i], val[i])
    return gradient_copy

#Tryout Function: A dummy function which holds all the possible answers for 4 bit binary represantation (for determining Langrage Multipliers later)
def write_this(num_symbols):
    [int(bin(i).split("b")[1]) for i in range(16)]

    val_array=[[0,0,0,0],
               [0,0,0,1],
               [0,0,1,0],
               [0,0,1,1],
               [0,1,0,0],
               [0,1,0,1],
               [0,1,1,0],
               [0,1,1,1],
               [1,0,0,0],
               [1,0,0,1],
               [1,0,1,0],
               [1,0,1,1],
               [1,1,0,0],
               [1,1,0,1],
               [1,1,1,0],
               [1,1,1,1]]
    #bin(10).split("b")[1]
    #np.pad(a,(2,0), "constant", constant_values=(5))
    return val_array

#Tryout Function:
def get_langrage(potential_gradients_problem, potential_gradients_boundary, num_symbols):
    langrage_vector=[]
    for i in range(len(potential_gradients_boundary)):
        for j in range(num_symbols):
            if(potential_gradients_boundary[i][j]==0):
                langrage_vector.append(nan)
            continue
        langrage_vector.append((potential_gradients_problem[i]/potential_gradients_boundary[i])[0])
    return langrage_vector


#Tryout Function: Sketches the cost matrix in 3D space - For own considerations, no use in the code
def sketch(Q, num_symbols):
    # Make data
    X = np.arange(0,num_symbols,1)
    Y = np.arange(0,num_symbols,1)
    X, Y = np.meshgrid(X, Y)
    Z = Q

    # Plot the surface
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(X, Y, Z, vmin=Z.min() * 2, cmap=cm.Blues)

    ax.set(xticklabels=["x","x","x","x"],
        yticklabels=["y","y","y","y"],
        zticklabels=["z","z","z","z"])

    plt.show()

#Tryout Function: Idea was to somewhat Gauß-Filter the matrix around the mean value - Not complete
def matrix_relaxer(Q, num_symbols):
    mean = np.mean(Q)
    max = np.max(Q)
    min = np.min(Q)
    #sketch(Q, num_symbols)
    distance_to_min = abs(min - mean)
    distance_to_max = abs(max - mean)
    
    if (distance_to_max < distance_to_min): 
        distance = distance_to_min
    else:
        distance = distance_to_max

    factor = abs((distance + abs(mean))/mean)

    if (factor>4):
        scaling = (factor/4)**(-1)
        for i in range(num_symbols):
            for j in range(num_symbols):
                Q[i,j] = mean+(Q[i,j]-mean) * scaling

    #sketch(Q, num_symbols)
    return Q

#Tryout Function: Returns back if the given matrix is unitary - For own considerations, not complete
def is_unitary(matrix: np.ndarray) -> bool:

    unitary = True
    n = int(np.sqrt(matrix.size))
    error = np.linalg.norm(np.eye(n) - matrix.dot( matrix.transpose().conjugate()))

    if not(error < np.finfo(matrix.dtype).eps * 10.0 *n):
        unitary = False
    
    return unitary