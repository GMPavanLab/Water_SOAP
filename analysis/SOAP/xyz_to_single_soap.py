#!/usr/bin/env python3
# coding: utf-8

import numpy as np
from ase.io import read
from dscribe.descriptors import SOAP
#from quippy import descriptors
from sklearn.decomposition import PCA

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#get_ipython().run_line_magic('matplotlib', 'inline')
import sys

def get_axes(L, max_col=3):
    cols = L if L <= max_col else max_col
    rows = int(L / max_col) + int(L % max_col != 0)
    fig, ax = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4))
    ax =  ax.flatten()
    return fig, ax



#--------------------------------------------------------------------
### Kernel (unit normalized)

def SOAPkernel(A,B,n):
    return ( np.dot(A, B) / (np.dot(A, A) * np.dot(B, B))**0.5 )**n

def SOAPdistance(ABker):
    np.seterr(all='raise')
    try:
        (2.0 - 2.0*ABker)**0.5
    except FloatingPointError:
        return 0
    return (2.0 - 2.0*ABker)**0.5

#--------------------------------------------------------------------
### Distance matrix (single frame)

def DISTmatrix(soap_frame, n_ker=1):
    from itertools import combinations
    comb = list(combinations(range(soap_frame.shape[0]), 2))
    DistMatrix = np.zeros( (soap_frame.shape[0], soap_frame.shape[0]) )
    # diagonal elements
    for j in range(40):
        DistMatrix[j][j] = SOAPdistance(SOAPkernel(soap_frame[j], soap_frame[j], n_ker))
    # non-diagonal elemets
    for i, (c1, c2) in  enumerate(comb):
        DistMatrix[c1][c2] = SOAPdistance(SOAPkernel(soap_frame[c1], soap_frame[c2], n_ker))
        DistMatrix[c2][c1] = SOAPdistance(SOAPkernel(soap_frame[c2], soap_frame[c1], n_ker))
    return DistMatrix


# read trajs

print('Reading trajectories...')

# SPC 
spc = read("../trajs/SPC.xyz",
               index=':',
               format="xyz")
spcbox = np.loadtxt("../trajs/SPC.box")
print('SPC read.')

# SPC/E
spce = read("../trajs/SPCE.xyz",
               index=':',
               format="xyz")
spcebox = np.loadtxt("../trajs/SPCE.box")
print('SPC/E read.')

# SPC/Eb
spceb = read("../trajs/SPCEb.xyz",
               index=':',
               format="xyz")
spcebbox = np.loadtxt("../trajs/SPCEb.box")
print('SPC/Eb read.')

# TIP3P
tip3p = read("../trajs/TIP3P.xyz",
               index=':',
               format="xyz")
tip3pbox = np.loadtxt("../trajs/TIP3P.box")
print('TIP3P read.')

# TIP3P-FB
tip3pfb = read("../trajs/TIP3P-FB.xyz",
               index=':',
               format="xyz")
tip3pfbbox = np.loadtxt("../trajs/TIP3P-FB.box")
print('TIP3P-FB read.')

# OPC3
opc3 = read("../trajs/OPC3.xyz",
               index=':',
               format="xyz")
opc3box = np.loadtxt("../trajs/OPC3.box")
print('OPC3 read.')

# TIP4P
tip4p = read("../trajs/TIP4P.xyz",
               index=':',
               format="xyz")
tip4pbox = np.loadtxt("../trajs/TIP4P.box")
print('TIP4P read.')

# TIP4P-EW
tip4pew = read("../trajs/TIP4PEW.xyz",
               index=':',
               format="xyz")
tip4pewbox = np.loadtxt("../trajs/TIP4PEW.box")
print('TIP4P-EW read.')

# TIP4P/2005
tip4p2005 = read("../trajs/TIP4P-2005.xyz",
               index=':',
               format="xyz")
tip4p2005box = np.loadtxt("../trajs/TIP4P-2005.box")
print('TIP4P/2005 read.')

# TIP4P-ICE
tip4pice = read("../trajs/TIP4P-ICE.xyz",
               index=':',
               format="xyz")
tip4picebox = np.loadtxt("../trajs/TIP4P-ICE.box")
print('TIP4P-ICE read.')

# TIP4P/ε
tip4peps = read("../trajs/TIP4P-eps.xyz",
               index=':',
               format="xyz")
tip4pepsbox = np.loadtxt("../trajs/TIP4P-eps.box")
print('TIP4P/ε read.')

# TIP4P-FB
tip4pfb = read("../trajs/TIP4P-FB.xyz",
               index=':',
               format="xyz")
tip4pfbbox = np.loadtxt("../trajs/TIP4P-FB.box")
print('TIP4P-FB read.')

# OPC
opc = read("../trajs/OPC.xyz",
               index=':',
               format="xyz")
opcbox = np.loadtxt("../trajs/OPC.box")
print('OPC read.')

# TIP5P
tip5p = read("../trajs/TIP5P.xyz",
               index=':',
               format="xyz")
tip5pbox = np.loadtxt("../trajs/TIP5P.box")
print('TIP5P read.')

# TIP5P-EW
tip5pew = read("../trajs/TIP5PEW.xyz",
               index=':',
               format="xyz")
tip5pewbox = np.loadtxt("../trajs/TIP5PEW.box")
print('TIP5P-EW read.')

# TIP5P-2018
tip5p2018 = read("../trajs/TIP5P-2018.xyz",
               index=':',
               format="xyz")
tip5p2018box = np.loadtxt("../trajs/TIP5P-2018.box")
print('TIP5P-2018 read.')


print('Adding PBC information...')
# correction to use periodic
for i in range(len(spc)):
    spc[i].set_cell(spcbox[i])
    spc[i].set_pbc([1, 1, 1])

for i in range(len(spce)):
    spce[i].set_cell(spcebox[i])
    spce[i].set_pbc([1, 1, 1])

for i in range(len(spceb)):
    spceb[i].set_cell(spcebbox[i])
    spceb[i].set_pbc([1, 1, 1])

for i in range(len(tip3p)):
    tip3p[i].set_cell(tip3pbox[i])
    tip3p[i].set_pbc([1, 1, 1])

for i in range(len(tip3pfb)):
    tip3pfb[i].set_cell(tip3pfbbox[i])
    tip3pfb[i].set_pbc([1, 1, 1])

for i in range(len(opc3)):
    opc3[i].set_cell(opc3box[i])
    opc3[i].set_pbc([1, 1, 1])

for i in range(len(tip4p)):
    tip4p[i].set_cell(tip4pbox[i])
    tip4p[i].set_pbc([1, 1, 1])

for i in range(len(tip4pew)):
    tip4pew[i].set_cell(tip4pewbox[i])
    tip4pew[i].set_pbc([1, 1, 1])

for i in range(len(tip4p2005)):
    tip4p2005[i].set_cell(tip4p2005box[i])
    tip4p2005[i].set_pbc([1, 1, 1])

for i in range(len(tip4pice)):
    tip4pice[i].set_cell(tip4picebox[i])
    tip4pice[i].set_pbc([1, 1, 1])

for i in range(len(tip4peps)):
    tip4peps[i].set_cell(tip4pepsbox[i])
    tip4peps[i].set_pbc([1, 1, 1])

for i in range(len(tip4pfb)):
    tip4pfb[i].set_cell(tip4pfbbox[i])
    tip4pfb[i].set_pbc([1, 1, 1])

for i in range(len(opc)):
    opc[i].set_cell(opcbox[i])
    opc[i].set_pbc([1, 1, 1])

for i in range(len(tip5p)):
    tip5p[i].set_cell(tip5pbox[i])
    tip5p[i].set_pbc([1, 1, 1])

for i in range(len(tip5pew)):
    tip5pew[i].set_cell(tip5pewbox[i])
    tip5pew[i].set_pbc([1, 1, 1])

for i in range(len(tip5p2018)):
    tip5p2018[i].set_cell(tip5p2018box[i])
    tip5p2018[i].set_pbc([1, 1, 1])


# SOAP Dscribe
species = ["O","H"]
rcut = float(sys.argv[1])
nmax = 8
lmax = 8
#soapDSC = SOAP(
#    species=species,
#    periodic=False,
#    rcut=rcut,
#    nmax=nmax,
#    lmax=lmax,
#    )
soapDSCAVE = SOAP(
    species=species,
    periodic=True,
    rcut=rcut,
    nmax=nmax,
    lmax=lmax,
    average='off'
    )

#positionAA=[i for i in range(128,256)]
#positionCG=[i for i in range(1,512,4)]

wat_indexes = [i for i in range(0,3072,3)]


out_folder = 'results_single_molecule'


#print('SOAPing SPC...')
#soapspc       = soapDSCAVE.create(spc,positions=[wat_indexes for x in range(1001)])
#np.savez_compressed(out_folder + '/spc_soap.npz',soapspc)

print('SOAPing SPC/E...')
soapspce      = soapDSCAVE.create(spce,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/spce_soap.npz',soapspce)

print('SOAPing SPC/Eb...')
soapspceb     = soapDSCAVE.create(spceb,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/spceb_soap.npz',soapspceb)

print('SOAPing TIP3P...')
soaptip3p     = soapDSCAVE.create(tip3p,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip3p_soap.npz',soaptip3p)

print('SOAPing TIP3P-FB...')
soaptip3pfb   = soapDSCAVE.create(tip3pfb,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip3pfb_soap.npz',soaptip3pfb)

print('SOAPing OPC3...')
soapopc3       = soapDSCAVE.create(opc3,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/opc3_soap.npz',soapopc3)

print('SOAPing TIP4P...')
soaptip4p     = soapDSCAVE.create(tip4p,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip4p_soap.npz',soaptip4p)

print('SOAPing TIP4P-EW...')
soaptip4pew   = soapDSCAVE.create(tip4pew,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip4pew_soap.npz',soaptip4pew)

print('SOAPing TIP4P/2005...')
soaptip4p2005 = soapDSCAVE.create(tip4p2005,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip4p2005_soap.npz',soaptip4p2005)

print('SOAPing TIP4P-ICE...')
soaptip4pice = soapDSCAVE.create(tip4pice,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip4pice_soap.npz',soaptip4pice)

print('SOAPing TIP4P/ε...')
soaptip4peps = soapDSCAVE.create(tip4peps,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip4peps_soap.npz',soaptip4peps)

print('SOAPing TIP4P-FB...')
soaptip4pfb   = soapDSCAVE.create(tip4pfb,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip4pfb_soap.npz',soaptip4pfb)

print('SOAPing OPC...')
soapopc       = soapDSCAVE.create(opc,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/opc_soap.npz',soapopc)

print('SOAPing TIP5P...')
soaptip5p     = soapDSCAVE.create(tip5p,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip5p_soap.npz',soaptip5p)

print('SOAPing TIP5P-EW...')
soaptip5pew   = soapDSCAVE.create(tip5pew,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip5pew_soap.npz',soaptip5pew)

print('SOAPing TIP5P-2018...')
soaptip5p2018   = soapDSCAVE.create(tip5p2018,positions=[wat_indexes for x in range(1001)])
np.savez_compressed(out_folder + '/tip5p2018_soap.npz',soaptip5p2018)

#print(soapspc.shape)




#print('Saving SOAP representations...')

#We close here
sys.exit()

# In[45]:

def comp_evo(soapvec,n=1):
    step_evolution = []
    step_comparison = []
    step_evolution.append(SOAPdistance(SOAPkernel(soapvec[0],soapvec[0],n)))
    step_comparison.append(SOAPdistance(SOAPkernel(soapvec[0],soapvec[0],n)))
    for t in range(soapvec.shape[0]-1):
        step_evolution.append(SOAPdistance(SOAPkernel(soapvec[t],soapvec[t+1],n)))
        step_comparison.append(SOAPdistance(SOAPkernel(soapvec[0],soapvec[t+1],n)))
    return step_evolution, step_comparison

# Average Kernel Investigation
#----------------------------------------------------------------

# Markovian-like investigation

evo_c36m, comp_c36m = comp_evo(soapc36m)
evo_slip, comp_slip = comp_evo(soapslip)
evo_lip17, comp_lip17 = comp_evo(soaplip17)
evo_22, comp_22 = comp_evo(soap22)
evo_30b2, comp_30b2 = comp_evo(soap30b2)
evo_22p, comp_22p = comp_evo(soap22p)
evo_23p, comp_23p = comp_evo(soap23p)
evo_dry, comp_dry = comp_evo(soapdry)
evo_dry16, comp_dry16 = comp_evo(soapdry16)


# In[46]:


Xrange = np.arange(0,200.1,0.1)
Xrange


rows=9
cols=2
fig, ax = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), dpi=200)
ax.flatten()

# charmm36m
ax[0,0].plot(Xrange,evo_c36m, c="Tab:blue", alpha=0.5, label='Charmm36m')
ax[0,0].set_title("Step Evolution")
ax[0,0].set_xlabel('time (ns)')
ax[0,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[0,1].plot(Xrange,comp_c36m, c="Tab:blue", alpha=0.5, label='Charmm36m')
ax[0,1].set_title("Step Comparison (with origin)")
ax[0,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[0,1].set_xlabel('time (ns)')
ax[0,1].legend()

# slipids
ax[1,0].plot(Xrange,evo_slip, c="Tab:blue", alpha=0.5, label='Slipids')
ax[1,0].set_title("Step Evolution")
ax[1,0].set_xlabel('time (ns)')
ax[1,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[1,1].plot(Xrange,comp_slip, c="Tab:blue", alpha=0.5, label='Slipids')
ax[1,1].set_title("Step Comparison (with origin)")
ax[1,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[1,1].set_xlabel('time (ns)')
ax[1,1].legend()

# lipid17
ax[2,0].plot(Xrange,evo_lip17, c="Tab:blue", alpha=0.5, label='Lipid17')
ax[2,0].set_title("Step Evolution")
ax[2,0].set_xlabel('time (ns)')
ax[2,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[2,1].plot(Xrange,comp_lip17, c="Tab:blue", alpha=0.5, label='Lipid17')
ax[2,1].set_title("Step Comparison (with origin)")
ax[2,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[2,1].set_xlabel('time (ns)')
ax[2,1].legend()

# Martini 22
ax[3,0].plot(Xrange,evo_22, c="Tab:blue", alpha=0.5, label='Martini22')
ax[3,0].set_title("Step Evolution")
ax[3,0].set_xlabel('time (ns)')
ax[3,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[3,1].plot(Xrange,comp_22, c="Tab:blue", alpha=0.5, label='Martini22')
ax[3,1].set_title("Step Comparison (with origin)")
ax[3,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[3,1].set_xlabel('time (ns)')
ax[3,1].legend()

# Martini 3.0 beta
ax[4,0].plot(Xrange,evo_30b2, c="Tab:blue", alpha=0.5, label='Martini30b2')
ax[4,0].set_title("Step Evolution")
ax[4,0].set_xlabel('time (ns)')
ax[4,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[4,1].plot(Xrange,comp_30b2, c="Tab:blue", alpha=0.5, label='Martini23')
ax[4,1].set_title("Step Comparison (with origin)")
ax[4,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[4,1].set_xlabel('time (ns)')
ax[4,1].legend()

# Martini 22p
ax[5,0].plot(Xrange,evo_22p, c="Tab:blue", alpha=0.5, label='Martini22p')
ax[5,0].set_title("Step Evolution")
ax[5,0].set_xlabel('time (ns)')
ax[5,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[5,1].plot(Xrange,comp_22p, c="Tab:blue", alpha=0.5, label='Martini22p')
ax[5,1].set_title("Step Comparison (with origin)")
ax[5,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[5,1].set_xlabel('time (ns)')
ax[5,1].legend()

# Martini 23p
ax[6,0].plot(Xrange,evo_23p, c="Tab:blue", alpha=0.5, label='Martini23p')
ax[6,0].set_title("Step Evolution")
ax[6,0].set_xlabel('time (ns)')
ax[6,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[6,1].plot(Xrange,comp_23p, c="Tab:blue", alpha=0.5, label='Martini23p')
ax[6,1].set_title("Step Comparison (with origin)")
ax[6,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[6,1].set_xlabel('time (ns)')
ax[6,1].legend()

# Dry Martini (pre-2016)
ax[7,0].plot(Xrange,evo_dry, c="Tab:blue", alpha=0.5, label='Dry Martini')
ax[7,0].set_title("Step Evolution")
ax[7,0].set_xlabel('time (ns)')
ax[7,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[7,1].plot(Xrange,comp_dry, c="Tab:blue", alpha=0.5, label='Dry Martini')
ax[7,1].set_title("Step Comparison (with origin)")
ax[7,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[7,1].set_xlabel('time (ns)')
ax[7,1].legend()

# Dry Martini (2016)
ax[8,0].plot(Xrange,evo_dry16, c="Tab:blue", alpha=0.5, label='Dry Martini 2016')
ax[8,0].set_title("Step Evolution")
ax[8,0].set_xlabel('time (ns)')
ax[8,0].set_ylabel(r'$d^{SOAP}(p_t,p_{t+1})$')

ax[8,1].plot(Xrange,comp_dry16, c="Tab:blue", alpha=0.5, label='Dry Martini 2016')
ax[8,1].set_title("Step Comparison (with origin)")
ax[8,1].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
ax[8,1].set_xlabel('time (ns)')
ax[8,1].legend()

#ax[2].scatter(step_evolution,step_comparison, c="Tab:blue", alpha=0.5, label='Charmm36m')
#ax[2].scatter(step_evolution1,step_comparison1, c="Tab:red", alpha=0.5, label='Slipids')
#ax[2].set_xlabel(r'$d^{SOAP}(p_t,p_{t+1})$')
#ax[2].set_ylabel(r'$d^{SOAP}(p_{t=0},p_{t})$')
#ax[2].legend()

fig.tight_layout()
fig.savefig(out_folder + "/averages_comparison.png")


# In[ ]:




