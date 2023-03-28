import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import mibian
import scipy as spy

class Sim:   #We define a class where each member contains the input parameters of the simulation
    def __init__(self, S0, k, iv, rf, i_t, i, T):
        self.S0 = S0   #known spot price
        self.K = k      #strike price 
        self.iv = iv
        self.rf = rf      #risk-free rate
        self.i_t = i_t      #no of time step iterations
        self.i = i       #no. of iterationsp
        self.T = T       #time period from the past time to today or today to future time

def calcvol(S1, k, r, ttm, optprice, type):  #r is in % and ttm is in days, takes the input of already known parameters of a certain option
    vol = 0
    if type == 'call':
        c = mibian.BS([S1, k, r, ttm], callPrice=optprice)
        vol = c.impliedVolatility
    elif type == 'put':
        c = mibian.BS([S1, k, r, ttm], putPrice=optprice)
        vol = c.impliedVolatility

    return (vol/100)

#We simulate the GBM through the following equation
# delx = nu*delt + vol*sqrt(delt)*Z where x = lnS

def Simulation(simlist):  #returns a results array which holds the spot price values for all time steps for all iterations and for all sims in the simlist
    results = []    
    for sim in simlist:
        delt = sim.T/sim.i_t
        Z = np.random.normal(size=(sim.i_t+1,sim.i))  #This specifies the matrix that contains all the random variables
        delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z   #this calculates the incremental increase in the value of lnS or x
        delx[0] = 0
        lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0)  #this takes the cumulative sum along the rows, this means that it calculates the total values of increment at each time step
        St = np.exp(lnSt)
        results.append(St)
    
    return results


def PlotMCS(results, simlist):                   #takes in the results that you get from the simulation function and plots it 
    for i in range(len(results)):
        plt.plot(np.linspace(1, np.round(simlist[i].T*365),simlist[i].i_t+1), results[i])
        plt.title("Simulation Results")
        plt.xlabel("No. of days")
        plt.ylabel("Spot Price")
        plt.show()

def Payoff(results, simlist):  #takes in the list of the sim parameters and returns the price and the standard error 
    payoff = []
    deviation = []
    for i in range(len(simlist)): 
        temp = np.array(results[i])       
        P = np.exp(-simlist[i].rf*simlist[i].T)*(np.maximum(0, (temp[-1] - simlist[i].K)))    #calculating and discounting back all the payoff
        
        P0 = np.sum(P)/simlist[i].i                     #Mean of the payoffs
        A  = np.sqrt((np.sum((P-P0)**2))/(simlist[i].i-1)) #Standard Error
        
        payoff.append(P0)
        deviation.append(A/np.sqrt(simlist[i].i))
        
    return payoff, deviation

#Inputting Parameters
sim1 = Sim(101.15, 98.01, calcvol(101.15,97.01,2,61, 4.8, 'call'), 0.01, 10, 1000, (((dt.date(2023,5,27)-dt.date.today()).days + 1)/365.0))

sims = [sim1]   #this contains your total array of simulations, each member represents one case that is taken

#Simulating
S_t = Simulation(sims)

price, error = Payoff(S_t, sims)
print(price, error)

PlotMCS(S_t, sims)


    
        


