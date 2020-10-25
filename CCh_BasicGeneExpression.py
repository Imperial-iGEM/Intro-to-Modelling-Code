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
delta_mRNA = 1e-3 #/s degredation consant of mRNA
#Translation:
ktl = 2e-3           #/s translation constant
delta_Protein = 1e-3  #/s degredation constant for Protein

#####
#Define the ODE to solve
#####

def ODEs(variables , t):
    #variables = list of concentrations, so here, [mRNA , Protein]. t = time
    mRNA = variables[0]    #indexes variables to specify we want the first element, mRNA
    Protein = variables[1] #indexes variables to specify we want the second element, Protein

    dmRNA_dt = ktx - delta_mRNA*mRNA #mRNA conc. ODE
    dProtein_dt = ktl*mRNA - delta_Protein*Protein             #Protein conc. ODE

    return [dmRNA_dt , dProtein_dt] #returns a list of the ODEs


#####
#Solving the ODEs
#####
t0 = 0              #Initial time
t1 = 36000/4        #Final time
total =  100000     #Number of time steps (larger the better)

initial_conditions = [0.0 , 0.0]        #set the initial values for [mRNA] and [Protein]
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
    'axes.labelsize':16,
    'font.size':20,
    'legend.fontsize':13,
    'xtick.labelsize':10,
    'ytick.labelsize':10,
    'figure.figsize': [9,9],
}
plt.rcParams.update(params)

#Ploting
plt.plot(t/3600 , mRNA, label = "mRNA concentration (M)")
plt.plot(t/3600 , Protein, label = "Protein concentration (M)")
plt.title("Variation of concentration with time")
plt.xlabel("time (hours)")
plt.ylabel("concentration (M)")
plt.grid()
plt.legend()
plt.show()
