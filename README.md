# Monte-Carlo-Error-Analysis
The code in this repo takes input of certain parameters associated with a particular option and simulates the expected payoff using the Monte Carlo Method. This payoff can be used to calculate the price of the option based on the no-arbitrage assumption. The code contains four functions:-

- calcvol - This uses the mibian python library to calculate the implied volatility of the market after the user has given it the parameters of an option contract whose price is already known. Usually the input parameters are the values taken from a similar contract with a different strike price.
- simulation - This function simulates the Geometric Brownian Motion of for a number of iterations i and time steps i_t. It returns the results in the form of an array containing all the spot prices for all time steps.
- payoff - This caclulates the payoff the call option and also finds out the standard error associated with it. It returns the price of the contract and its standard error. The inputs are the simulation list and the results from the simulation function.
- plotMCS - If required, the iterations can be plotted using this command. They are graphed against day-to-day progression to show how each spot prices fares as the time increases. The input for this function includes the simulation list and the results from the simulation function. 

Each simulation is treated as a class defined by its input parameters and all the test cases have to be entered in a list which can then be passed to the functions. The use of numpy, matplotlib, datetime and mibian has been made. scipy has also been installed as a dependency for mibian. 
