#Example code for Introduction to Mathematical Modelling in Synthetic Biology Package
#Written by Gabriel Swallow, Imperial College London

#####
#Import the modules and submodules necessary for the code
#####

import scipy as sp                 #Fundamental module for science
import matplotlib.pyplot as plt    #Module for plotting data
from scipy.integrate import odeint #specify the method "odeint"

#####
#Define the the constants
#####

#####
#Define the ODE to solve
#####

def ODEs(variables , t):
    #variables = list of concentrations, t = time

    #Define ODEs

    return [ODE1 , ODE2 , ...] #returns a list of the ODEs

#####
#Solving the ODEs
#####
t0 = 0              #Initial time
t1 = 36000          #Final time
total =  100000     #Number of time steps (larger the better)
t = sp.linspace(t0,t1,total)    #set the array of time values to integrate over

initial_conditions = []         #set the initial values for variable concentrations as floats

solution = odeint(ODEs , initial_conditions , t) #Produces a 2d array of solutions
                                                 #for variable concentrations
Variable1 = solution[:,0]    #Index all values in first column for variable 1
Variable2 = solution[:,1]    #Index all values in second column for variable 2


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
plt.plot(t/3600 , Variable1, label = "Protein concentration (M)")
plt.plot(t/3600 , Variable2, label = "Protein concentration (M)")
...

plt.title("This is the Title")
plt.xlabel("x axis label")
plt.ylabel("y axis label")
plt.grid()
plt.legend()
plt.show()
