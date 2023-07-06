# Metropolis algorithm for reconstructing EAS parameters
> **Required gcc 6+** <br/>
**Start program** - ./main input_data config <br/>
**Description of data format:** <br/>
**Input data** - no. of particles in mip in each det. [0:no.of det.] + trigger time of det.[no. of det:2*no. of det.] + X, Y, $\theta, \phi$ (from model) <br/>
**Number of particles, times - int** <br/>
**Parameters - float** <br/>
The program has a choice of LDF approximation.  
Algorithm based on creation of a Markov chain (a sequence of random events with a finite or countable number of outcomes,
characterized by the fact that the next element of the sequence depends only on the previous one). 
After finding the initial approximation of the parameters, we will directly apply the Monte Carlo method. The essence of the method is as follows:

The coordinate number in the 4-dimensional space (**X**, **Y**, **Ne**, **s**) is played sequentially.
An increment is played for the played coordinate according to following formula. The remaining coordinates are fixed.

$\Delta =-1+2*\gamma*\varepsilon$,

where $\gamma$ is a uniformly distributed random variable on the interval [0;1],
$\varepsilon$ is an arbitrarily selected length that is different for each coordinate.
The result of these random events is compared with the possible movement to a new position.
Then the value of the maximum likelihood functional in this possible state is calculated.
If the value of the functional $\Im_{j}>\Im_{i}$, then we accept new coordinates and the system is considered to have passed into a new, $j-th$ state.
If $\Im_{j}\leq\Im_{i}$, then the drawing procedure is repeated.

$x^{(i)}\rightarrow x^{(j)}=\ x^{(i)}+\ ∆$ , where $x^{(i)}$ is the coordinate in phase space.

