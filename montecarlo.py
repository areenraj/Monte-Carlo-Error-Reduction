import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import mibian
import scipy as spy
import scipy.stats as stats


class Sim:   #We define a class where each member contains the input parameters of the simulation
    def __init__(self, S0, k, iv, rf, i_t, i, T):
        self.S0 = S0    #known spot price
        self.K = k      #strike price 
        self.iv = iv
        self.rf = rf    #risk-free rate
        self.i_t = i_t  #no of time step iterations
        self.i = i      #no. of iterationsp
        self.T = T      #time period from the past time to today or today to future time

def calcvol_delta(S1, k, r, ttm, optprice, type):  #r is in % and ttm is in days, takes the input of already known parameters of a certain option
    vol = 0
    if type == 'call':
        c = mibian.BS([S1, k, r, ttm], callPrice=optprice)
        vol = c.impliedVolatility
    elif type == 'put':
        c = mibian.BS([S1, k, r, ttm], putPrice=optprice)
        vol = c.impliedVolatility

    return (vol/100)

#We simulate the GBM through the following equation
#delx = nu*delt + vol*sqrt(delt)*Z where x = lnS

def simulation(simlist,type):  #returns a results array which holds the spot price values for all time steps for all iterations and for all sims in the simlist
    results = [] 
    resultsminus = [] 
    for sim in simlist:
        if type == 'Sobol':
            sampler = stats.qmc.Sobol(d=sim.i_t+1, scramble = True)         #This if structure tree decides what random number generation to use
            Z = sampler.random_base2(m=np.log2(sim.i))
            Z = stats.norm.ppf(Z.T)
            delt = sim.T/sim.i_t                                             #This specifies the matrix that contains all the random variables
            delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z    #This calculates the incremental increase in the value of lnS or x
            delx[0] = 0
            dely = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*(-Z) #This calculates the incremental increase in the value of lnS or x
            dely[0] = 0
            lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0)                  #This takes the cumulative sum along the rows, this means that it calculates the total values of increment at each time step
            St = np.exp(lnSt)
            results.append(St)
            lnStminus = np.log(sim.S0) + np.cumsum(dely, axis=0)
            stminus = np.exp(lnStminus)
            resultsminus.append(stminus)
        elif type == 'Halton':
            sampler = stats.qmc.Halton(d=sim.i_t+1, scramble = True)
            Z = sampler.random(m=sim.i)
            Z = stats.norm.ppf(Z.T)
            delt = sim.T/sim.i_t  
            delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z   
            delx[0] = 0
            dely = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*(-Z)   
            dely[0] = 0
            lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0) 
            St = np.exp(lnSt)
            results.append(St)
            lnStminus = np.log(sim.S0) + np.cumsum(dely, axis=0)
            stminus = np.exp(lnStminus)
            resultsminus.append(stminus)
        else:
            Z = np.random.normal(size=(sim.i_t+1,sim.i))
            delt = sim.T/sim.i_t  
            delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z  
            delx[0] = 0
            dely = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*(-Z)   
            dely[0] = 0
            lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0)  
            St = np.exp(lnSt)
            results.append(St)
            lnStminus = np.log(sim.S0) + np.cumsum(dely, axis=0)
            stminus = np.exp(lnStminus)
            resultsminus.append(stminus)
    
    return results, resultsminus

def simulation_delta_gamma(simlist,type):  #This is similar to the original simulation function but uses Control Variates
    results = [] 
    control = []
    resultsminus = []
    controlminus = []
    for sim in simlist:
        if type == 'Sobol':
            k = np.log2(sim.i)
            delt = sim.T/sim.i_t

            sampler = stats.qmc.Sobol(d=sim.i_t+1, scramble=True)
            Z = sampler.random_base2(m=tick)
            Z =  stats.norm.ppf(Z.T)  
            
            delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z   
            delx[0] = 0
            dely = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*(-Z)  
            dely[0] = 0
            
            lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0)  
            St = np.exp(lnSt)
            lnStminus = np.log(sim.S0) + np.cumsum(dely, axis=0)
            stminus = np.exp(lnStminus)
            results.append(St)
            resultsminus.append(stminus)

            d1 = ((np.log(St[:-1].T/sim.K) + (sim.rf + (sim.iv**2)/2)*np.linspace(sim.T,0,sim.i_t))/(sim.iv*np.linspace(sim.T,0,sim.i_t))).T
            delta1 = stats.norm.cdf(d1,0,1)
            delcv = delta1*(St[1:]-(St[:-1]*(2.7183**(sim.rf*delt))))
            cv = np.cumsum(delcv, axis=0)
            control.append(cv)

            d1 = ((np.log(stminus[:-1].T/sim.K) + (sim.rf + (sim.iv**2)/2)*np.linspace(sim.T,0,sim.i_t))/(sim.iv*np.linspace(sim.T,0,sim.i_t))).T
            delta1 = stats.norm.cdf(d1,0,1)
            delcv = delta1*(stminus[1:]-(stminus[:-1]*(2.7183**(sim.rf*delt))))
            cvminus = np.cumsum(delcv, axis=0)
            controlminus.append(cvminus)

        elif type == 'Halton':
            delt = sim.T/sim.i_t

            sampler = stats.qmc.Halton(d=sim.i_t+1)
            Z = sampler.random(sim.i)
            Z =  stats.norm.ppf(Z.T)  
            
            delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z  
            delx[0] = 0
            dely = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*(-Z)   
            dely[0] = 0
            
            lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0)  
            St = np.exp(lnSt)
            lnStminus = np.log(sim.S0) + np.cumsum(dely, axis=0)
            stminus = np.exp(lnStminus)
            results.append(St)
            resultsminus.append(stminus)

            d1 = ((np.log(St[:-1].T/sim.K) + (sim.rf + (sim.iv**2)/2)*np.linspace(sim.T,0,sim.i_t))/(sim.iv*np.linspace(sim.T,0,sim.i_t))).T
            delta1 = stats.norm.cdf(d1,0,1)
            delcv = delta1*(St[1:]-(St[:-1]*(2.7183**(sim.rf*delt))))
            cv = np.cumsum(delcv, axis=0)
            control.append(cv)

            d1 = ((np.log(stminus[:-1].T/sim.K) + (sim.rf + (sim.iv**2)/2)*np.linspace(sim.T,0,sim.i_t))/(sim.iv*np.linspace(sim.T,0,sim.i_t))).T
            delta1 = stats.norm.cdf(d1,0,1)
            delcv = delta1*(stminus[1:]-(stminus[:-1]*(2.7183**(sim.rf*delt))))
            cvminus = np.cumsum(delcv, axis=0)
            controlminus.append(cvminus)
        else:
            delt = sim.T/sim.i_t

            Z = np.random.normal(size=(sim.i_t+1,sim.i))  
            
            delx = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*Z   
            delx[0] = 0
            dely = (sim.rf - 0.5*sim.iv**2)*delt + sim.iv*np.sqrt(delt)*(-Z)  
            dely[0] = 0
            
            lnSt = np.log(sim.S0) + np.cumsum(delx, axis=0)  
            St = np.exp(lnSt)
            lnStminus = np.log(sim.S0) + np.cumsum(dely, axis=0)
            stminus = np.exp(lnStminus)
            results.append(St)
            resultsminus.append(stminus)


            d1 = ((np.log(St[:-1].T/sim.K) + (sim.rf + (sim.iv**2)/2)*np.linspace(sim.T,0,sim.i_t))/(sim.iv*np.linspace(sim.T,0,sim.i_t))).T
            delta1 = stats.norm.cdf(d1,0,1)
            delcv = delta1*(St[1:]-(St[:-1]*(2.7183**(sim.rf*delt))))
            cv = np.cumsum(delcv, axis=0)
            control.append(cv)
            d1 = ((np.log(stminus[:-1].T/sim.K) + (sim.rf + (sim.iv**2)/2)*np.linspace(sim.T,0,sim.i_t))/(sim.iv*np.linspace(sim.T,0,sim.i_t))).T
            delta1 = stats.norm.cdf(d1,0,1)
            delcv = delta1*(stminus[1:]-(stminus[:-1]*(2.7183**(sim.rf*delt))))
            cvminus = np.cumsum(delcv, axis=0)
            controlminus.append(cvminus)
    
    return results, resultsminus, control, controlminus

def atv_payoff(results, resultsminus, simlist):    #This is used to calculate the Anti Thetic Variates Payoff
    payoff = [] 
    deviation = []
    for i in range(len(simlist)): 
        temp = np.array(results[i])  
        tempminus = np.array(resultsminus[i])     
        P = np.exp(-simlist[i].rf*simlist[i].T)*0.5*(np.maximum(0, (temp[-1] - simlist[i].K))+ np.maximum(0, (tempminus[-1] - simlist[i].K)))    #calculating and discounting back all the payoff 
        P0 = np.sum(P)/simlist[i].i                     #Mean of the payoffs
        A  = np.sqrt((np.sum((P-P0)**2))/(simlist[i].i-1)) #Standard Error
        
        payoff.append(P0)
        deviation.append(A/np.sqrt(simlist[i].i))
        
    return payoff, deviation

def delta_gamma_payoff(results, control, simlist):  #This is used for the Control Variates case
    payoff = []
    deviation = []
    for i in range(len(simlist)): 
        temp = np.array(results[i])
        cv = np.array(control[i])   
        P = np.exp(-simlist[i].rf*simlist[i].T)*(np.maximum(0, (temp[-1] - simlist[i].K))-cv[-1])    #calculating and discounting back all the payoff 
        P0 = np.sum(P)/simlist[i].i                    
        A  = np.sqrt((np.sum((P-P0)**2))/(simlist[i].i-1)) 
        
        payoff.append(P0)
        deviation.append(A/np.sqrt(simlist[i].i))
        
    return payoff, deviation

def normal_payoff(results,  simlist):  #Takes in the list of the sim parameters and returns the price and the standard error 
    payoff = []
    deviation = []
    for i in range(len(simlist)): 
        temp = np.array(results[i])       
        P = np.exp(-simlist[i].rf*simlist[i].T)*(np.maximum(0, (temp[-1] - simlist[i].K)))    #calculating and discounting back all the payoff
        
        P0 = np.sum(P)/simlist[i].i                    
        A  = np.sqrt((np.sum((P-P0)**2))/(simlist[i].i-1)) 
        
        payoff.append(P0)
        deviation.append(A/np.sqrt(simlist[i].i))
        
    return payoff, deviation

def atv_delta_gamma_payoff(results,resultsminus, control, controlminus, simlist):
    payoff = []
    deviation = []
    for i in range(len(simlist)): 
        temp = np.array(results[i])  
        tempminus = np.array(resultsminus[i])  
        cv = np.array(control[i])    
        cvminus = np.array(controlminus[i])
        P = (np.exp(-simlist[i].rf*simlist[i].T))*0.5*(np.maximum(0, (temp[-1] - simlist[i].K)) + np.maximum(0, (tempminus[-1] - simlist[i].K)) - cv[-1] - cvminus[-1])    #calculating and discounting back all the payoff 
        P0 = np.sum(P)/simlist[i].i                     
        A  = np.sqrt((np.sum((P-P0)**2))/(simlist[i].i-1)) 
        payoff.append(P0)
        deviation.append(A/np.sqrt(simlist[i].i))
        
    return payoff, deviation


def plotMCS(results, simlist):                   #Takes in the results that you get from the simulation function and plots it 
    for i in range(len(results)):
        plt.plot(np.linspace(1, np.round(simlist[i].T*365),simlist[i].i_t+1), results[i])
        plt.title("Simulation Results")
        plt.xlabel("No. of days")
        plt.ylabel("Spot Price")
    
    plt.show()


simtemp = Sim(855.4, 850, calcvol_delta(855.4,830,6.26,12,30.3,'call'), 0.0626, 100, 1000, (((dt.date(2023,8,31)-dt.date.today()).days + 1)/365.0))

sims = [simtemp]
    s
S_t, S_tminus, control, controlminus = simulation_delta_gamma(sims,'PNRG')

price1, error1 = normal_payoff(S_t, sims)
price2, error2 = atv_payoff(S_t, S_tminus, sims)
price3, error3 = delta_gamma_payoff(S_t,control,sims)
price4, error4 = atv_delta_gamma_payoff(S_t,S_tminus,control, controlminus, sims)

print(calcvol_delta(855.4,830,6.26,12,30.3,'call'))

        
print(price1,error1)    #Normal 
print(price2,error2)    #Anti-Thetic    
print(price3,error3)    #Delta Gamma Control Variates
print(price4,error4)    #Combined
 
plotMCS(S_t,sims)
