# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 23:10:29 2014
scripts for data analysis of properties of the statisitical mechanics for the MAS
with varying absolute velocity, mainly focusing on the phase transition according
to the order parameter. 

adding functions:
(20140818)
1. calculating mean and variance of order parameter by numpy function np.mean() and ny.var()
2.
@author: 9020
"""

import os
import numpy as np

# defining parameters and constants
PARTICLE_NUM = 100 # the number of particles
NSTEPS = 200000 # the number of total steps of simulation

# cumulative average of the order parameter
filename_ord_para = 'order_parameter_'+ str(NSTEPS) + '_steps_' + str(PARTICLE_NUM) +'_particles.txt' # reading data of order parameter from saved file
if os.path.isfile(filename_ord_para):
    f_op = open(filename_ord_para, 'r')
    order_para = []
    for d in f_op:
        order_para.append(float(d))
    f_op.close()
    print 'loading order parameter success'
else:
    print 'loading order parameter failure'

cum_ave_list = []
relaxation_time = 40000 # how to derermine the relaxation time?
#cum_ave = sum(order_para[relaxation_time:-1])/len(order_para[relaxation_time:-1]) # cumulative average of the order parameter
op_ave = np.mean(order_para[relaxation_time:])
op_var = np.var(order_para[relaxation_time:])
print relaxation_time
print op_ave
print op_var
print '-------------------------------------'
