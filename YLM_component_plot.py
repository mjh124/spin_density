# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:54:57 2015

@author: Michael
"""

import numpy as np
import matplotlib.pyplot as plt

# Extract data form YLM Summary file

MyData = open('../Documents/superatom/Au13_ylm_l7_r5.txt', 'r')  
lines = MyData.readlines()[16:]
MyData.close()

YourData = open('../Documents/superatom/ylm_Au12Co1_edge_summary.txt', 'r')  
lines1 = YourData.readlines()[16:]
YourData.close()

HisData = open('../Documents/superatom/ylm_Au12Co1_middle_summary.txt', 'r')  
lines2 = HisData.readlines()[16:]
HisData.close()

def ReadYLMSum(lines, ef):
    
#    lines = data.readlines()[16:]
    
    en_up = []
    en_dw = []
    
    for line in lines:
        tokens = line.split()
        spin = int(tokens[0])
        if spin == 0:
            en_up.append(float(tokens[2])-ef)
        else:
            en_dw.append(float(tokens[2])-ef)
            
    return en_up, en_dw
    
def GetWeight(lines, ang_mom):

#    lines = data.readlines()[16:]
    l = 7 + ang_mom
    
    component_up = []
    component_dw = []
    total_up = []
    total_dw = []
    
    for line in lines:
        tokens = line.split()
        spin = int(tokens[0])
        if spin == 0:
            total_up.append(float(tokens[4]))
            component_up.append(float(tokens[l]))
        else:
            total_dw.append(float(tokens[4]))
            component_dw.append(float(tokens[l]))
            
    weight_up = []
    weight_dw = []
    for i in range(len(component_up)):
        up = component_up[i] / total_up[i]
        dw = component_dw[i] / total_dw[i]
        
        weight_up.append(up)
        weight_dw.append(dw)
        
    return weight_up, weight_dw
    
def GaussianBroadener(sig, data_up, data_dw, weight_up=np.ones(150), weight_dw=np.ones(150)):
        
    x_axis = np.arange(-16.0, 5.05, 0.05)

    x_mat = len(data_up)
    y_mat = len(x_axis)

    data_gauss_up = np.zeros(y_mat)
    data_gauss_dw = np.zeros(y_mat)
    for i in range(y_mat):
        tmp_up = 0.0
        tmp_dw = 0.0
        for j in range(x_mat):
            gauss_up = weight_up[j]*np.exp(-((x_axis[i] - data_up[j])/sig)**2)
            tmp_up += gauss_up
            gauss_dw = weight_dw[j]*np.exp(-((x_axis[i] - data_dw[j])/sig)**2)
            tmp_dw += gauss_dw
            
        data_gauss_up[i] = data_gauss_up[i] + tmp_up
        data_gauss_dw[i] = data_gauss_dw[i] - tmp_dw
    
    return data_gauss_up, data_gauss_dw

plt.clf()      
x_axis = np.arange(-16.0, 5.05, 0.05)
en_up, en_dw = ReadYLMSum(lines, -5.09393)
en1_up, en1_dw = ReadYLMSum(lines1, -4.90410)
en2_up, en2_dw = ReadYLMSum(lines2, -4.91592)

#s_up, s_dw = GetWeight(lines, 0)
#p_up, p_dw = GetWeight(lines, 1)
#d_up, d_dw = GetWeight(lines, 2)
#f_up, f_dw = GetWeight(lines, 3)
#g_up, g_dw = GetWeight(lines, 4)
#h_up, h_dw = GetWeight(lines, 5)
#i_up, i_dw = GetWeight(lines, 6)

#s_up, s_dw = GetWeight(lines, 6)
#s1_up, s1_dw = GetWeight(lines1, 6)
#s2_up, s2_dw = GetWeight(lines2, 6)
#
#sweighted_up, sweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=s_up, weight_dw=s_dw)
#s1weighted_up, s1weighted_dw = GaussianBroadener(0.05, en1_up, en1_dw, weight_up=s1_up, weight_dw=s1_dw)
#s2weighted_up, s2weighted_dw = GaussianBroadener(0.05, en2_up, en2_dw, weight_up=s2_up, weight_dw=s2_dw)
#
#plt.plot(x_axis, sweighted_up, '-k', label='Co13')
#plt.plot(x_axis, s1weighted_up, '-r', label='Au12Co13')
#plt.plot(x_axis, s2weighted_up, '-b', label='Co13PH3')
#
#plt.plot(x_axis, sweighted_dw, '-k')
#plt.plot(x_axis, s1weighted_dw, '-r')
#plt.plot(x_axis, s2weighted_dw, '-b')

unweighted_up, unweighted_dw = GaussianBroadener(0.05, en_up, en_dw)
unweighted1_up, unweighted1_dw = GaussianBroadener(0.05, en1_up, en1_dw)
unweighted2_up, unweighted2_dw = GaussianBroadener(0.05, en2_up, en2_dw)

#unweighted_up, unweighted_dw = GaussianBroadener(0.05, en_up, en_dw)
#sweighted_up, sweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=s_up, weight_dw=s_dw)
#pweighted_up, pweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=p_up, weight_dw=p_dw)
#dweighted_up, dweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=d_up, weight_dw=d_dw)
#fweighted_up, fweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=f_up, weight_dw=f_dw)
#gweighted_up, gweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=g_up, weight_dw=g_dw)
#hweighted_up, hweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=h_up, weight_dw=h_dw)
#iweighted_up, iweighted_dw = GaussianBroadener(0.05, en_up, en_dw, weight_up=i_up, weight_dw=i_dw)

#plt.plot(x_axis, sweighted_up, '-k', label='s')
#plt.plot(x_axis, pweighted_up, '-r', label='p')
#plt.plot(x_axis, dweighted_up, '-b', label='d')
#plt.plot(x_axis, fweighted_up, '-g', label='f')
#plt.plot(x_axis, gweighted_up, '-c', label='g')
#plt.plot(x_axis, hweighted_up, '-m', label='h')
#plt.plot(x_axis, iweighted_up, '-y', label='i')

#plt.plot(x_axis, sweighted_dw, '-k')
#plt.plot(x_axis, pweighted_dw, '-r')
#plt.plot(x_axis, dweighted_dw, '-b')
#plt.plot(x_axis, fweighted_dw, '-g')
#plt.plot(x_axis, gweighted_dw, '-c')
#plt.plot(x_axis, hweighted_dw, '-m')
#plt.plot(x_axis, iweighted_dw, '-y')

# Fill in area under the curve
#plt.fill(x_axis, sweighted_up)
#plt.fill(x_axis, pweighted_up)
#plt.fill(x_axis, dweighted_up)
#plt.fill(x_axis, fweighted_up)
#plt.fill(x_axis, gweighted_up)
#plt.fill(x_axis, hweighted_up)
#plt.fill(x_axis, iweighted_up)
#plt.fill(x_axis, sweighted_dw)
#plt.fill(x_axis, pweighted_dw)
#plt.fill(x_axis, dweighted_dw)
#plt.fill(x_axis, fweighted_dw)
#plt.fill(x_axis, gweighted_dw)
#plt.fill(x_axis, hweighted_dw)
#plt.fill(x_axis, iweighted_dw)

# Plot as sums of all smaller angular momenta

#sp_up = sweighted_up + pweighted_up
#sp_dw = sweighted_dw + pweighted_dw
#spd_up = sweighted_up + pweighted_up + dweighted_up
#spd_dw = sweighted_dw + pweighted_dw + dweighted_dw
#spdf_up = sweighted_up + pweighted_up + dweighted_up + fweighted_up
#spdf_dw = sweighted_dw + pweighted_dw + dweighted_dw + fweighted_dw
#spdfg_up = sweighted_up + pweighted_up + dweighted_up + fweighted_up + gweighted_up
#spdfg_dw = sweighted_dw + pweighted_dw + dweighted_dw + fweighted_dw + gweighted_dw
#spdfgh_up = sweighted_up + pweighted_up + dweighted_up + fweighted_up + gweighted_up + hweighted_up
#spdfgh_dw = sweighted_dw + pweighted_dw + dweighted_dw + fweighted_dw + gweighted_dw + hweighted_dw
#spdfghi_up = sweighted_up + pweighted_up + dweighted_up + fweighted_up + gweighted_up + hweighted_up + iweighted_up
#spdfghi_dw = sweighted_dw + pweighted_dw + dweighted_dw + fweighted_dw + gweighted_dw + hweighted_dw + iweighted_dw

plt.plot(x_axis, unweighted_up, '-b', label=r'$Au_{13}$')
plt.plot(x_axis, unweighted1_up, '-r', label=r'$Au_{12}Co_1 E$')
plt.plot(x_axis, unweighted2_up, '-k', label=r'$Au_{12}Co_1 M$')
#plt.plot(x_axis, sweighted_up, '-k')
#plt.plot(x_axis, sp_up, '-r')
#plt.plot(x_axis, spd_up, '-b')
#plt.plot(x_axis, spdf_up, '-g')
#plt.plot(x_axis, spdfg_up, '-c')
#plt.plot(x_axis, spdfgh_up, '-m')
#plt.plot(x_axis, spdfghi_up, '-y')
#
plt.plot(x_axis, unweighted_dw, '-b')
plt.plot(x_axis, unweighted1_dw, '-r')
plt.plot(x_axis, unweighted2_dw, '-k')
#plt.plot(x_axis, sweighted_dw, '-k')
#plt.plot(x_axis, sp_dw, '-r')
#plt.plot(x_axis, spd_dw, '-b')
#plt.plot(x_axis, spdf_dw, '-g')
#plt.plot(x_axis, spdfg_dw, '-c')
#plt.plot(x_axis, spdfgh_dw, '-m')
#plt.plot(x_axis, spdfghi_dw, '-y')

font = {'family' : ['Helvetica'],
        'weight' : 'normal',
        'size'   : 14}

#ax1.set_xlabel('System', fontsize=20, family='Helvetica')
#ax1.set_ylabel('Orbital Energy (eV)', fontsize=20, family='Helvetica')

plt.axhline(color='k')
plt.axvline(color='k', ls='dotted')
plt.axis([-4, 4, -6.0, 7.0])
plt.xlabel('Energy, E-Ef', fontsize=20)
plt.ylabel('Density of States (A.U.)', fontsize=20)

plt.legend(bbox_to_anchor=(1.0, 1), loc=1, borderaxespad=0.)

plt.show()

#fig = plt.gcf()
#fig.set_size_inches(12, 6)
#fig.savefig('../Desktop/Co13PH3Cl_X_shells_dos.png',dpi=120)