#Example code for Introduction to Mathematical Modelling in Synthetic Biology Package
#Written by Gabriel Swallow, Imperial College London

#####
#Import the modules and submodules necessary for the code
#####

import scipy as sp                 #Fundamental module for science
import matplotlib.pyplot as plt    #Module for plotting data
from scipy.integrate import odeint #specify the method "odeint" for efficiency

#####
#Define the Functions
#####

#Constants
ka = 0.005      #/s/M   association rate constant
kd = 0.00001    #//s    dissociation rate constant
Kd = kd/ka      #no units
n = 1           #no units

A_0 = 1.0          #Initial A Concentration
B_0 = 5.0          #Initial B Concentration
C_0 = 0.0          #Initial C Concentration

#ODE to solve
def ODEs(variables , t):
        #variables = list of concentrations, so here, [A, B, C]. t = time
    A = variables[0]  #indexes variables to specify we want the first element, A conc.
    B = variables[1]   #indexes variables to specify we want the second element, B conc.
    C = variables[2]    #indexes variables to specify we want the third element, C conc.

    DE = ka*A*(B**n) - kd*C

    dA_dt = -DE            #Unbound conc. ODE
    dB_dt = -DE            #Ligand conc. ODE
    dC_dt =  DE            #Bound conc. ODE


    return [dA_dt , dB_dt , dC_dt] #returns a list of the ODEs


#Also define the hill equation approximation:
Hill = A_0*(B_0**n / (Kd + B_0**n))
Hill2 = A_0*( 1 - (B_0**n / (Kd + B_0**n)))


######
#Solutions
######

t0 = 0              #Initial time
t1 = 100             #Final time
total =  100000     #Number of time steps (larger the better)

####
#Solving
####

initial_conditions = [A_0 , B_0 , C_0]        #set the initial values for [mRNA] and [Protein]
t = sp.linspace(t0,t1,total)    #set the array of time values to integrate over

solution = odeint(ODEs , initial_conditions , t) #Produces an 2d array of solutions
                                                 #for [mRNA] and [Protein]

A = solution[:,0]  #Index all values in first column
B = solution[:,1]   #Index all values in second column
C = solution[:,2]    #Index all values in third column

#Now we need to define a straight line, with concentration equal to Hill
#and Hill2 at all points in time
Hill_array = [Hill for i in range(total)]
Hill2_array = [Hill2 for i in range(total)]

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
plt.plot(t , A , label="Bound conc.")
plt.plot(t , B , label="Unbound conc.")
plt.plot(t , C , label="Ligand conc.")
plt.plot(t , Hill_array , label="Hill Approx conc.")
plt.plot(t , Hill2_array , label="Hill Approx 2 conc.")


plt.title("Variation of concentration with time")
plt.xlabel("Time (s)")
plt.ylabel("Cocentrations")
plt.grid()
plt.legend()
plt.show()
