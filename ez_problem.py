import helpers
import dwave
import dwave.inspector
#from dwave.system import LeapHybridSampler
from dwave.system import DWaveSampler
from dimod import BinaryQuadraticModel
import numpy as np
from collections import defaultdict
from sympy import symbols

### X=xyzq, binary ==> X=8*x+4*y+2*z+q
from sympy.abc import x, y, z, q 

REPS = 1

### Constant displacement
C = 6 

for i in range(REPS):
    symbol = [x,y,z,q]
    num_symbols=len(symbol)
    Q = np.zeros([num_symbols,num_symbols])
    #L = np.zeros([num_symbols,num_symbols]) #tryout for langrange multiplier
    expression_1 = (8*x+4*y+2*z+q)
    problem = (expression_1**2).expand() #.free_symbols ; len(problem.coeff(x*y, 1).free_symbols)
    expression_2= (8*x+4*y+2*z+q-C)
    boundary = (expression_2**2).expand()


    ### Saving coefficients of multiplicative variables -- PROBLEM
    for i in range(num_symbols):
        for j in range(num_symbols):
            coefficient = problem.coeff(symbol[i]*symbol[j],1)
            Q[i,j]+=coefficient  

    #Setting up the Langrange Multiplier : Using the average of the cost matrix
    '''langrage_hard = (np.sum(Q))/ (len(Q)**2) #np.max(Q)

    for i in range(num_symbols):
        for j in range(num_symbols):
            Q[i,j]+=2*langrage_hard 

    for i in range(num_symbols):
        Q[i,i]+= -3*langrage_hard'''

    #print (Q)

    ### Saving coefficients of multiplicative variables -- BOUNDARY CONDITION
    for i in range(num_symbols):
        for j in range(num_symbols):
            coefficient = boundary.coeff(symbol[i]*symbol[j],1)
            Q[i,j]+=coefficient
    #print (L)


    
    ### Saving coefficients of singular variables -- BOUNDARY CONDITION
    for i in range(num_symbols):
        coefficient = boundary.coeff(symbol[i],1)
        coefficient = coefficient.as_independent(x,y,z,q)[0]
        Q[i,i]+=coefficient
    #print (L)

    ### This part is for trying to find a suitable Langrange Multiplier - Not complete
    potential_gradients_problem=[]
    potential_gradients_boundary=[]
    val_array = helpers.write_this(num_symbols)

    
    gradient_problem = helpers.gradient(problem, symbol)
    gradient_boundary = helpers.gradient(boundary, symbol)
    for i in range(num_symbols**2):
        potential_gradients_problem.append(helpers.gradient_evaluation(gradient_problem, symbol, num_symbols,val_array[i]))
        potential_gradients_boundary.append(helpers.gradient_evaluation(gradient_boundary, symbol, num_symbols,val_array[i]))

    langrage_vector = helpers.get_langrage(potential_gradients_problem, potential_gradients_boundary, num_symbols)
    langrange_hard = langrage_vector[int(C/2)]

    #Setting up the Langrange Multiplier : Using the average of the cost matrix
    #langrage_hard = (np.sum(L))/ (len(L)**2) #np.max(Q)

    for i in range(num_symbols):
        for j in range(num_symbols):
            Q[i,j]+=2*langrange_hard 

    for i in range(num_symbols):
        Q[i,i]+= -3*langrange_hard

    #print (L)

    ### Offset Energy 
    e_offset=boundary.as_independent(x,y,z,q)[0]*langrange_hard

    #Q = Q + L
    #print(Q)

    #print(helpers.is_unitary(Q))
    #Q = helpers.matrix_relaxer(Q,num_symbols)
    #helpers.sketch(Q,num_symbols)

    bqm = BinaryQuadraticModel.from_qubo(Q, offset=e_offset)
    sampler = DWaveSampler()     
    #sampler = LeapHybridSampler()
    results = sampler.sample_qubo(bqm, num_reads=100)
    results = sampler.sample(bqm, label='ez_problem')
    smpl = results.first.sample
    #dwave.inspector.show(sampler)
    #Building decimal number from binary solution: MSB is smpl[0]
    answer = 8*smpl[0] + 4*smpl[1] + 2*smpl[2] + smpl[3]

    print(smpl)
    print("Answer is ", answer)
    dwave.inspector.show(smpl) 