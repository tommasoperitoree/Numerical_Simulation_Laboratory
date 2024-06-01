import os
import sys
import numpy as np

'''
This program builds a set of `input.in` files so that we can cycle over the parameters and extract the best parameters needed for the thermalization.
Input: 
- Please provide the name of the phase, one choice among: ['Liquid', 'Solid', 'Gas'], this will grant access to the specific folder.
'''

acceptedArgs = ['Liquid', 'Solid', 'Gas']                                               # Admissable command line inputs
phase_name = sys.argv[1]                                                                # Retrieve the phase name from the command line argument

# Check if the provided phase name is one of the accepted values
if phase_name in acceptedArgs:
    dir_name ='../'+ phase_name +"/INPUT/"
else: 
    # Raise an error if the phase name is not accepted
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs)

# Check if the directory exists, if not, create it
if not os.path.exists(dir_name):
    print('Making a new directory')
    os.makedirs(dir_name)

# Define the cycling variable depending on phase type
if phase_name == 'Solid' :
   cycling_var = list(np.arange(1.50, 1.58, 0.01)) #0.4,1.3,0.1 gala, io prima 1.50, 1.60, 0.01
if phase_name == 'Liquid' :
   cycling_var = list(np.arange(1.95, 2.03, 0.01)) # 
if phase_name == 'Gas' :
   cycling_var = list(np.arange(0.92, 1.0, 0.01)) #0.5,1.4,0.1 gala, io prima 0.92, 1.03, 0.01

cycling_var = ["{:.2f}".format(num) for num in cycling_var]

with open(dir_name+'temp_values.in', 'w') as temp_values:
   for i, cv in enumerate(cycling_var):
      temp_values.writelines(str(cv)+'\n')