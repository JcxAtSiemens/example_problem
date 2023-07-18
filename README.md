# example_problem

## A Basic Mathematical Optimisation Problem With A Basic Boundary Condition To Demonstrate QUBO Formulation

This script demonstrates how to formulate the optimisation problem

$$ min \text{ } x^2 $$

$$ s.t. \text{ } x^2-c \quad \text{for a constant c} $$

up to a positive integer constant of 15.

### Formulation of Minimization Problem and the Boundary Condition

The problem variable is expressed in binary manner to make it easier to express in the forms of 1's and 0's - Those are the only allowed values in the QUBO formulation.
In the code, this is expressed as
> expression_1 = (8x+4y+2z+q)

coefficients being the powers of 2. _expression_1_ is here just our variable -x- given in the minimization problem. The x in the expression itself is just a simpy symbol,
  and should not be confused with the one in the optimization problem.

Boundary condition is expressed in the same manner with a constant factor between 0 and 15 (4 bit's allowance span).

### Creation of the Cost Matrix

After the problem is created, coefficients of the variables are placed in the cost matrix. We don't differentiate between quadratic expressions and singular experrions(such as $x^2$ and $x$),
since $0^2=0$ and $1^2=1$. The coefficients of the singular variables are placed in the diagonal of the cost matrix, while the coefficients of the multiplicative variables are placed in the 
off-diagonal elements.

$$
Q_{m,n} = 
 \begin{pmatrix}
  q_{x,x} & q_{x,y} & q_{x,z} & q_{x,q} \\
  q_{y,x} & q_{y,y} & q_{y,z} & q_{y,q} \\
  q_{z,x} & q_{z,y}  & q_{z,z} & q_{z,q}  \\
  q_{q,x} & q_{q,y} & q_{q,z} & q_{q,q} 
 \end{pmatrix}
$$

### Choice of Lambda
_That's the big question!_

### Determination of the Offset Boundary
The constant part in the problem definition (that would be e.g. 36, for the case of c=6 or 144, for the case of c=12) is the energy which can not be decreased with the manipulation of the variables.
Hence, it stays in the problem. This energy is communicated with the quantum annealer 
>e_offset=boundary.as_independent(x,y,z,q)[0]*langrange_hard

### Reformulation of the Answer from the Annealer Back to the Integer
>answer = 8*smpl[0] + 4*smpl[1] + 2*smpl[2] + smpl[3]

### helpers.py
Here are some helping functions for the sake of keeping the main script clean. Hopefully they are well commented by the time!
