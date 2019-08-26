# better-homomorphic-sine-evaluation

This program aims to find an approximate polynomial of sine function for bootstrapping for HEAAN (In fact, we approximate the function ![](http://latex.codecogs.com/gif.latex?%5Ccos%282%20%5Cpi%20x%29), but there is no big difference). Homomorphic evaluation of sine function is the key part of bootstrapping for HEAAN, and we suggest the bootstrapping optimized way to approximate sine function by a polynomial in the paper "Better Bootstrapping for Approximate Homomorphic Encryption" (https://eprint.iacr.org/2019/688.pdf). 



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

You should put integer for each argument, and each argument has following meaning.

 - Degree : degree bound of interpolation polynomial 
 - Error : ![](http://latex.codecogs.com/gif.latex?-%5Clog_2) of maximum deviation of test values from each i-0.25
 - Scaling :  the number of scaling
 
[Refer to "Better Bootstrapping for Approximate Homomorphic Encryption" (https://eprint.iacr.org/2019/688.pdf) for more detail]

For example, you can type "./find_polynomial 30 10 2"

In summary, you should type following commands on terminal.
 1. Type "cmake .". 
 2. Type "make all".
 3. Type "./find_polynomial 'degree' 'error' 'scaling'" 


## How to interpret the result?
### 1. The result printed on terminal. 

[Example]

![](https://user-images.githubusercontent.com/30550389/63686143-1b1d2a80-c83c-11e9-9769-1ab1bcbb8ceb.png)

 - Degree of the polynomial : the degree of the interpolation polynomial of ![](http://latex.codecogs.com/gif.latex?%5Ccos%282%20%5Cpi%20x%29).
 - Degree : Deg[i] is the number of nodes in ![](http://latex.codecogs.com/gif.latex?I_i%20%3D%20%5Bi-0.25-e%2C%20i-0.25&plus;e%5D).
 - Max_Error : log_2 of maximum error of the interpolation polynomial computed through Bs-GS algorithm.


### 2. The result stored in files

#### (1) Errors
The errors between ![](http://latex.codecogs.com/gif.latex?%5Ccos%282%20%5Cpi%20x%29) and the approximate polynomial on various tested values are stored in the directory "result/error/".

The file is of .csv format, and its name is the form of "DegxxErrxxScalingxx.csv", where xx's are argument values you put.

[Example]

|Tested Values|Real Values|Approximate Values|Error|log_2 of Error|
|--       |--         |--         |--         |--      |
|-0.250977|-0.00613588|-0.00613588|7.01431E-09|-27.0871|
|-0.250938|-0.00589045|-0.00589045|6.73471E-09|-27.1457|
|-0.250937|-0.00564502|-0.00564501|6.45502E-09|-27.2069|
|...      |...        |...        |...        |...     |


#### (2) Coefficients    
The coefficients for Baby-step Giant-step algorithm introduced in the paper, which is needed for implementing bootstrapping for HEAAN, are stored in the directory "result/coef/".

The file is of .csv format, and its name is the form of "DegxxErrxxScalingxx.csv", where xx's are argument values you put.

[Example]

|||||||||
|--|--|--|--|--|--|--|--|
|0.0292182|-0.131917|-0.295158|-0.102132|0.857435|-0.0648289|0.384785|-0.0222439|
|0.0864300|-0.156425|-0.332357|-0.110424|-1.38072|-0.0652592|0.771865| -0.0215346|
|0.451401|-0.0301473|-0.906402|-0.0158559|0.394436|-0.00718066|-0.138797|-0.00200036|
|0.0281924|-0.00196015|-0.0162168|-0.000480076|0.000828514|-6.78823E-05| -0.000500756|0|

Let ![](http://latex.codecogs.com/gif.latex?q_i%28x%29) ![](http://latex.codecogs.com/gif.latex?%24%28i%3D1%2C2%2C3%2C4%29%24) be a polynomial corresponding to i-th row.

For example, ![](http://latex.codecogs.com/gif.latex?q_1%28x%29%20%3D%200.0292182%20-%200.131917%20x%20-0.295158%20x%5E2%20&plus;%20%5Ccdots%20-0.0222439%20x%5E7). 

Then, the approximate polynomial ![](http://latex.codecogs.com/gif.latex?p%28x%29) is given by ![](http://latex.codecogs.com/gif.latex?p%28x%29%20%3D%20q_1%28x%29%20&plus;%20q_2%28x%29%5Ccdot%20T_8%28x%29%20&plus;%20%28q_3%28x%29%20&plus;%20q_4%28x%29%20%5Ccdot%20T_8%28x%29%29%20%5Ccdot%20T_%7B16%7D%28x%29), where ![](http://latex.codecogs.com/gif.latex?T_i%28x%29) is the adjusted chebyshev polynomial of deg i.
