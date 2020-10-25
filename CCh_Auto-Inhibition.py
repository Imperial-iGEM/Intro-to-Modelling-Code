#Example code for Introduction to Mathematical Modelling in Synthetic Biology Package
#Written by Gabriel Swallow, Imperial College London

#####
#Import the modules and submodules necessary for the code
#####

import scipy as sp                 #Fundamental module for science
import matplotlib.pyplot as plt    #Module for plotting data
from scipy.integrate import odeint #specify the method "odeint" for efficiency

#####
#Define the the constants
#####

#Transcription:
ktx = 1e-3        #M/s transcription constant
kd = 1e-8         #M  dissociation constant
delta_mRNA = 1e-3 #/s degredation consant of mRNA
n=1
#Translation:
ktl = 2e-3              #/s translation constant
delta_Protein = 1e-3  #/s degredation constant for Protein


#####
#Define the ODEs
#####


def ODEs(variables , t):
    #variables = list of concentrations, so here, [mRNA , Protein]. t = time
    mRNA = variables[0]    #indexes variables to specify we want the first element, mRNA
    Protein = variables[1] #indexes variables to specify we want the second element, Protein


    dmRNA_dt = ktx*(1-(Protein**n / (kd + Protein**n))) - delta_mRNA*mRNA #mRNA conc. ODE
    dProtein_dt = ktl*mRNA - delta_Protein*Protein             #Protein conc. ODE


    return [dmRNA_dt , dProtein_dt] #returns a list of the ODEs

#####
#Solving the ODEs
#####

t0 = 0              #Initial time
t1 = 3600*2           #Final time
total =  100000     #Number of time steps (larger the better)

initial_conditions = [0.0 , 0.0001]        #set the initial values for [mRNA] and [Protein]
t = sp.linspace(t0,t1,total)    #set the array of time values to integrate over


solution = odeint(ODEs , initial_conditions , t) #Produces an 2d array of solutions
                                                #for [mRNA] and [Protein]
mRNA = solution[:,0]    #Index all values in first column
Protein = solution[:,1] #Index all values in second column


#####
#Plot the data
#####

#Set the parameters for the figure   (arbitrary values, vary as you like)
params = {
    'axes.labelsize':10,
    'font.size':15,
    'legend.fontsize':10,
    'xtick.labelsize':8,
    'ytick.labelsize':8,
    'figure.figsize': [8,8],
}
plt.rcParams.update(params)

#Ploting
plt.plot(t/3600 , mRNA, label="mRNA conc.")
plt.plot(t/3600 , Protein , label="Protein conc.")


plt.title("Variation of concentration with time")
plt.xlabel("Time (hours)")
plt.ylabel("Concentration (M)")
plt.grid()
plt.legend()
plt.show()
