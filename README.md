# Monte-Carlo-Error-Analysis

## Introduction
The Monte Carlo Method is a simulation method of estimating the value of probabilistic quantities that may follow a stochastic equation. This estimation is done by iterating over several possible cases using random variables then averaging their values. This averaging method is at the core of Monte Carlo Simulations. The standard error of this process can be estimated as follows,

```math
\frac{\sigma}{\sqrt{N}}
```
While technically this standard error is not exactly related to the accuracy of the final result corresponding to the physical quantity we are simulating, it does give us a general idea of the accuracy of the Monte Carlo formulation. One way to colloquially oversimplify this is to think of the process as throwing a large number of darts on a board. The efforts to increase your hopes of a bulls-eye can be approached in two ways - controlling the 'spread' of the darts and the other is increasing the number of darts that are being thrown. The standard therefore contains the variance of the final data set - indicating spread - and a division by the total number of possible sample cases. Hence the standard error is a measure of the accuracy of the method employed not the simulation that produces the final result.

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

Each simulation is treated as a class defined by its input parameters and all the test cases have to be entered in a list which can then be passed to the functions. The use of numpy, matplotlib, datetime and mibian has been made. scipy has also been installed as a dependency for mibian. 

## Error Reduction

One way to reduce the error is to increase the number of possibilities being considered. This is a crude way to throw computational power at the problem. A more elegant approach is variance reduction through the use of certain techniques that help to reduce the variance of the entire simulation. We analyze three methods 
