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

#Values given here are fairly arbitrary, with some taken from my 2017 model


#Transcription:
a1 , a2 , a3 = 1e-3 , 1e-3 , 1e-3      #M/s transcription constant
a01 , a02 , a03 = 1e-6 , 1e-6 , 1e-6      #M/s transcription constant
#a01 , a02 , a03 = 0,0,0      #M/s transcription constant

Kd1 , Kd2 , Kd3 = 1e-6 , 1e-6 , 1e-6
delta_mRNA1 , delta_mRNA2 , delta_mRNA3 = 1e-3 , 1e-3 , 1e-3 #/s degredation consant of mRNA
n1 , n2 , n3 = 2,2,2
#Translation:
ktl1 , ktl2 , ktl3 = 0.002 , 0.002, 0.002               #/s translation constant
delta_Protein1 , delta_Protein2 , delta_Protein3 = 1e-3 , 1e-3 , 1e-3  #/s degredation constant for Protein

#####
#Define the ODE to solve
#####

def ODEs(variables , t):
    #variables = list of concentrations, so here, [mRNA , Protein]. t = time
    mRNA1 = variables[0]    #indexes variables to specify we want the first element, mRNA
    Protein1 = variables[1] #indexes variables to specify we want the second element, Protein
    mRNA2 = variables[2]    #indexes variables to specify we want the first element, mRNA
    Protein2 = variables[3] #indexes variables to specify we want the second element, Protein
    mRNA3 = variables[4]    #indexes variables to specify we want the first element, mRNA
    Protein3 = variables[5] #indexes variables to specify we want the second element, Protein

    dmRNA1_dt       = a1/(1+(Protein3**n1)/Kd1) + a01 - delta_mRNA1*mRNA1 #mRNA conc. ODE
    dProtein1_dt    = ktl1*mRNA1 - delta_Protein1*Protein1            #Protein conc. ODE
    dmRNA2_dt       = a2/(1+(Protein1**n2)/Kd2) + a02 - delta_mRNA2*mRNA2 #mRNA conc. ODE
    dProtein2_dt    = ktl2*mRNA2 - delta_Protein2*Protein2            #Protein conc. ODE
    dmRNA3_dt       = a3/(1+(Protein2**n3)/Kd3) + a03 - delta_mRNA3*mRNA3 #mRNA conc. ODE
    dProtein3_dt    = ktl3*mRNA3 - delta_Protein3*Protein3            #Protein conc. ODE


    return [dmRNA1_dt , dProtein1_dt , dmRNA2_dt , dProtein2_dt , dmRNA3_dt , dProtein3_dt] #returns a list of the ODEs


#####
#Solving the ODEs
#####
t0 = 0              #Initial time
t1 = 36000           #Final time
total =  100000     #Number of time steps (larger the better)

initial_conditions = [0.1 , 0.2 , 0.0 , 0.0 , 0.0 , 0.0]        #set the initial values for [mRNA] and [Protein]
t = sp.linspace(t0,t1,total)    #set the array of time values to integrate over

solution = odeint(ODEs , initial_conditions , t) #Produces an 2d array of solutions
                                                 #for [mRNA] and [Protein]

mRNA1 = solution[:,0]    #Index all values in first column
Protein1 = solution[:,1] #Index all values in second column
mRNA2 = solution[:,2]    #Index all values in first column
Protein2 = solution[:,3] #Index all values in second column
mRNA3 = solution[:,4]    #Index all values in first column
Protein3 = solution[:,5] #Index all values in second column

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
#plt.plot(t/3600 , mRNA1, label = "mRNA concentration (M)")
plt.plot(t/3600 , Protein1, label = "Protein concentration (M)")
#plt.plot(t/3600 , mRNA2, label = "mRNA concentration (M)")
plt.plot(t/3600 , Protein2, label = "Protein concentration (M)")
#plt.plot(t/3600 , mRNA3, label = "mRNA concentration (M)")
plt.plot(t/3600 , Protein3, label = "Protein concentration (M)")

plt.title("Variation of concentration with time")
plt.xlabel("time (hours)")
plt.ylabel("concentration (M)")
plt.grid()
plt.legend()
plt.show()
