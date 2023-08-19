# Monte-Carlo-Error-Analysis
_This project was based on the information available in the book "Implementing Derivative Models" by Les Clewlow and Chris Strickland_
## Introduction
The Monte Carlo Method is a simulation method of estimating the value of probabilistic quantities that may follow a stochastic equation. This estimation is done by iterating over several possible cases using random variables then averaging their values. This averaging method is at the core of Monte Carlo Simulations. The standard error of this process can be estimated as follows,

```math
\frac{\sigma}{\sqrt{N}}
```
While technically this standard error is not exactly related to the accuracy of the final result corresponding to the physical quantity we are simulating, it does give us a general idea of the accuracy of the Monte Carlo formulation. One way to colloquially oversimplify this is to think of the process as throwing a large number of darts on a board. The efforts to increase your hopes of a bulls-eye can be approached in two ways - controlling the 'spread' of the darts and the other is increasing the number of darts that are being thrown. The standard error therefore contains the variance of the final data set - indicating spread - and a division by the total number of possible sample cases. Hence the standard error is a measure of the accuracy of the method employed and not that of the simulation that produces the final result.

The following project aims to achieve two things.

1. Analysis the magnitude of standard error of the simulation and find approaches to decrease it.
2. Improve the rate of convergence of the simulation.

## Setup

The test case use to analyze the use of the simulation in a real world scenario is that of option pricing. We use the Monte Carlo Simulation to simulate S, the stock price of the security and then find the payoff from that. The stochastic equation that the stock price follows is based on the Black Scholes Model which is as follows,
```math
\ln S_{i+1} = \ln S_i + (r- \frac{\sigma^2}{2})\Delta t + \sigma \sqrt{ \Delta t}N(0,1)
```

The payoff can be used to calculate the price of the option based on the no-arbitrage assumption. The code contains four functions:-

* calcvol - This uses the mibian python library to calculate the implied volatility of the market after the user has given it the parameters of an option contract whose price is already known. Usually the input parameters are the values taken from a similar contract with a different strike price.
* simulation - This function simulates the Geometric Brownian Motion of for a number of iterations i and time steps i_t. It returns the results in the form of an array containing all the spot prices for all time steps.
* payoff - This caclulates the payoff the call option and also finds out the standard error associated with it. It returns the price of the contract and its standard error. The inputs are the simulation list and the results from the simulation function.
* plotMCS - If required, the iterations can be plotted using this command. They are graphed against day-to-day progression to show how each spot prices fares as the time increases. The input for this function includes the simulation list and the results from the simulation function. 


![Blank diagram](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/2cfd075b-7fa4-4663-91e4-4f5815f9f5f5)

Each simulation is treated as a class defined by its input parameters and all the test cases have to be entered in a list which can then be passed to the functions. The use of numpy, matplotlib, datetime and mibian has been made. scipy has also been installed as a dependency for mibian. 

## Error Reduction

One way to reduce the error is to increase the number of possibilities being considered. This is a crude way to throw computational power at the problem. A more elegant approach is variance reduction through the use of certain techniques that help to reduce the variance of the entire simulation. We analyze three methods that could be used to accomplish the same

### Anti-Thetic Variates
This part is more along the lines of a statistical trick that employs the use of the total variance formula and the concept of the covariance between two variables. The total variance of a variable that is composed of two other variables is as follows
```math
\theta = \frac{N_1 + N_2}{2}
```
```math
var(\theta) = \frac{var(N_1) + var(N_2) + cov(N_1,N_2)}{4}
```
If the correlation between the two random variables is negative then it can lead to a total variance reduction. This is achieved through the help of taking the negative of the random variable we are simulating and taking the system through two oppositely mirrored states simulatneously. This allows for a perfectly negative correlation between the two variable and hence causes the variance and the standard error to drop. 

### Control Variates
This method employs the use of another random variable whose variance we can use to manipulate the total variance of the simulation. Let there be two random variables like follows
```math
X,Y
```
```math
E[X] = \mu
```
```math
E[Y] = \tau
```
Then the variable with the same mean as the original variable is as follows
```math
M = X + c(Y- \tau)
```
Therefore, the variance becomes 
```math
{\displaystyle {\textrm {Var}}(M)={\textrm {Var}}\left(X\right)+c^{2}\,{\textrm {Var}}\left(Y\right)+2c\,{\textrm {Cov}}\left(X,Y\right).}
```
Y is called the control variate and in our current case we take it to be the stock price. We can differentiate the above expression and get the optimal value of c that will lower the total variance. Doing so we get
```math
c = \frac{-cov(X,Y)}{Var(Y)}
```
Using the stock price, we find that the value of c corresponds to the physical value of delta of the option. Hence, this method of using control variates is equivalent to a very powerful concept in trading and finance - delta hedging. 

The final equation to calculate the pay off is as follows
```math
C_{initial} + \Sigma \frac{\partial C}{\partial S}(S_{i+1}-E[S_i]) = C_{final}
```
From the above equation we can calculate the initial value of the option price, the summation is over all iterations. 

## Error-Reduction Compiled Results

The simulaiton was run with 100 timesteps with 1000 possibilities per iteration. The test case was a call option of the Reliance group. 

|  Method Used  |Absolute Standard Error  |Relative Standard Error  |
|---------------|-------------------------|-------------------------|
|Simple MCS     |         3.53243         |         3.62544%        |       
|Anti-Thetic Var|         1.40395         |         1.42156%        |
|Control Var    |         0.78059         |         0.80146%        |
|ATV + Control  |         0.38696         |         0.39896%        |
____________________________________________________________________

![SIMPLE_MCS](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/21899bb7-8abd-437a-a9d6-90675834e722)
![timestep](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/e361b3fc-bb69-411d-a7a2-a99c9808fedf)
![Legend](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/f74a1de2-d4d7-41d9-9ec7-b6118cc5185b)

## Error Rate Convergence

Here we experiment with different number generators to see if the we can achieve a faster rate for convergence.

The standard generator is that of numpy's Pseudo Random Number Generator (PNRG). There is another method that is more suited for Monte Carlo Simulations - Quasi Random Number Generators. These generators don't actually simulate random numbers, they are more deterministic in nature. The use of primitive polynomials that prevent both excessive clustering and gaps in the number space are used. For the case of Monte Carlo, this is beneficial. 

The two Quasi Random Generators we use are from scipy's library - Halton and Sobol Distribution. The results for them are as follows

![Pseudo](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/e9bb3dd3-e6d5-44af-9fef-8cfd2cb65c8d)
![Halton](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/5be27db0-12de-41cd-b11a-2ea78b0e50b7)
![Sobol](https://github.com/areenraj/Monte-Carlo-Error-Analysis/assets/80944803/e3a653cd-47c1-46ea-910d-182ee4240454)
</br>
It can be clearly inferred that the Sobol generator converges faster. 
