# better-homomorphic-sine-evaluation

This program aims to find an approximate polynomial of a sine function for the bootstrapping for HEAAN (In fact, we approximate the function $\cos(2 \pi x)$, but there is no big difference). Homomorphic evaluation of a sine function is the key part of the bootstrapping for HEAAN, and we suggest the bootstrapping optimized way to approximate a sine function by a polynomial in the paper "Better Bootstrapping for Approximate Homomorphic Encryption" (https://eprint.iacr.org/2019/688.pdf). 



## Dependency

It is written in C++ and requires NTL library (with GMP).
You can download each library at 
 - NTL : https://www.shoup.net/ntl/doc/tour.html.
 - GMP : https://gmplib.org/

We checked that the program properly works on Ubuntu 18.04.1 LTS with the libraries NTL(ver. 11.3.2) and GMP(ver. 6.1.2)

## How to use it?

First, you can make "Makefile" by typing "cmake .".

Next, type "make all" and obtain the executable file "find_polynomial". 

You can use the executable file by typing "./find_polynomial 'Degree '
'Error' 'Scaling'". 

You should put integer for each argument, and each argument has the following meaning.

Degree : degree bound of interpolation polynomial 

Error : $- \log$ (base 2) of maximum deviation of test values from each $i-0.25$

Scaling :  the number of scaling
 
[Refer to "Better Bootstrapping for Approximate Homomorphic Encryption" (https://eprint.iacr.org/2019/688.pdf) for more details]

For example, you can type "./find_polynomial 30 10 2"

In summary, you should type the following commands on terminal.
 1. Type "cmake .". 
 2. Type "make all".
 3. Type "./find_polynomial 'degree' 'error' 'scaling'" 


## How to interpret results?
### 1. A result printed on terminal

[Example]

![](https://user-images.githubusercontent.com/30550389/63686143-1b1d2a80-c83c-11e9-9769-1ab1bcbb8ceb.png)

Degree of polynomial : degree of an interpolation polynomial of $\cos(2 \pi x)$.
 
Degree : Deg[i] is the number of nodes in $I_{i} = \left[i - 0.25 - e, i - 0.25 + e \right]$.
 
Max_Error : $\log$ of a maximum error of an interpolation polynomial computed through the BS-GS algorithm.


### 2. Results stored in files

#### (1) Errors
Errors between $\cos(2 \pi x)$ and an approximate polynomial on various tested values are stored in the directory "result/error/".

The file is of .csv format, and its name is the form of "DegxxErrxxScalingxx.csv", where xx's are argument values you put in.

[Example]

|Tested Values|Real Values|Approximate Values|Error|log of Error|
|--       |--         |--         |--         |--      |
|-0.250977|-0.00613588|-0.00613588|7.01431E-09|-27.0871|
|-0.250938|-0.00589045|-0.00589045|6.73471E-09|-27.1457|
|-0.250937|-0.00564502|-0.00564501|6.45502E-09|-27.2069|
|...      |...        |...        |...        |...     |


#### (2) Coefficients    
Coefficients for the Baby-step Giant-step algorithm introduced in the paper, which is needed for implementing the bootstrapping for HEAAN, are stored in the directory "result/coef/".

The file is of .csv format, and its name is the form of "DegxxErrxxScalingxx.csv", where xx's are argument values you put in.

[Example]

|||||||||
|--|--|--|--|--|--|--|--|
|0.0292182|-0.131917|-0.295158|-0.102132|0.857435|-0.0648289|0.384785|-0.0222439|
|0.0864300|-0.156425|-0.332357|-0.110424|-1.38072|-0.0652592|0.771865| -0.0215346|
|0.451401|-0.0301473|-0.906402|-0.0158559|0.394436|-0.00718066|-0.138797|-0.00200036|
|0.0281924|-0.00196015|-0.0162168|-0.000480076|0.000828514|-6.78823E-05| -0.000500756|0|

Let $q_{i} (x)$ be a polynomial corresponding to $i$-th row.

For example, $q_{i} (x) = 0.0292182 - 0.131917 x - 0.295158 x^{2} + \cdots - 0.0222439 x^{7}$. 

Then, the approximate polynomial $p(x)$ is given by $p(x) = q_{1}(x) + q_{2}(x) \cdot T_{8}(x) + (q_{3}(x) + q_{4}(x) \cdot T_{8}(x)) \cdot T_{16}(x)$, where $T_{i} (x)$ is the adjusted chebyshev polynomial of deg $i$.
