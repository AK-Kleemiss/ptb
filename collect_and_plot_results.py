import os
import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator) 
import traceback

def truncstring(thres):
    if float(thres)==0:
        return "no truncation"
    else:
        return r'truncation $10^{'+str(int(math.log10(float(thres))))+'}$'

rundir=sys.argv[1].replace("/","")

thres=["0.001","0.0001","0.00001","0.000001"]
modes=["full_matrix","kmeans_r_heuristic"]
systems=[]
charges=[]
dims=[]

sygvd={}

for s in os.listdir("systems/"):
    if s.endswith(".xyz"):
        s2=s.replace(".xyz","")
        if s2!="chrysophanol" and s2!="ethane":
            systems.append(s2)
            f = open("systems/"+s2+".charge")
            charge=int(f.read())
            f.close()
            charges.append(charge)
print(systems)
print(charges)

data={}

for m in modes:
    for f in thres:
        for i in range(len(systems)):
            tsolve2=0
            tsm=0
            dim=0
            nsm=0
            de=0
            mse=0
            inf=0
            try:
                fi=open(rundir+"/"+systems[i]+"_"+str(charges[i])+"/filter_"+f+"/"+m+"/out")
                print(rundir+"/"+systems[i]+"_"+str(charges[i])+"/filter_"+f+"/"+m+"/out")
                pop=False
                bonds=False
                dipole=False
                for l in fi.read().splitlines():
                    if l.startswith("Total solve2"):
                        tsolve2=float(l.split()[2])
                    if l.startswith("Submatrix "):
                        tsm=float(l.split()[1])
                    if l.startswith(" basis setup done. Ndim"):
                        dim=int(l.split()[4])
                    if l.startswith(" number of submatrices="):
                        nsm=int(l.split()[3])
                    if l.startswith("* sygvd solver"):
                        sygvd[dim]=float(l.split()[3])
                    if l.startswith("PTB H matrix iteration                      2"):
                        break
                    if l.startswith(" energy deviation per atom"):
                        de=abs(float(l.split()[7].replace("E","e")))
                    if l.startswith(" MSE of density matrix"):
                        mse=abs(float(l.split()[4].replace("E","e")))
                    if l.startswith(" Inf-norm of density matrix"):
                        inf=abs(float(l.split()[4].replace("E","e")))
                    if l.startswith(" largest Wiberg (=(PS)*(SP)) bond orders for each atom"):
                        pop=False
                        bonds=False
                        dipole=False
                    if l.startswith("dipole moment  X         Y          Z"):
                        pop=False
                        bonds=False
                        dipole=False
                    if l.startswith(" CT cal"):
                        pop=False
                        bonds=False
                        dipole=False
                    if l.startswith("  # element  Z   pop."):
                        pop=True
                        bonds=False
                        dipole=False
                        data[m+"_"+str(f)+"_"+str(systems[i])+"_atom"]=[]
                        data[m+"_"+str(f)+"_"+str(systems[i])+"_pop"]=[]
                        data[m+"_"+str(f)+"_"+str(systems[i])+"_shells"]=[]
                    else:
                        if pop:
                            ls=l.split()
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_atom"].append(float(ls[0]))
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_pop"].append(float(ls[3]))
                            x=[]
                            for j in range(4,len(ls)):
                                x.append(float(ls[j]))
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_shells"].append(x)
                    if l.startswith("           total WBO             WBO to atom > 0.05 ..."):
                        pop=False
                        bonds=True
                        dipole=False
                        data[m+"_"+str(f)+"_"+str(systems[i])+"_bond_atom"]=[]
                        data[m+"_"+str(f)+"_"+str(systems[i])+"_bond_tot"]=[]
                        data[m+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"]=[]
                    else:
                        if bonds:
                            ls=l.split()
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_bond_atom"].append(float(ls[0]))
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_bond_tot"].append(float(ls[2]))
                            x=[]
                            for j in range(100):
                                try:
                                    x.append(float(ls[3+j*3+2]))    
                                except: 
                                    break
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"].append(x)
                    if l.startswith("dipole moment  X         Y          Z"):
                        pop=False
                        bonds=False
                        dipole=True
                    else:
                        if dipole:
                            x=float(l.split()[0])
                            y=float(l.split()[1])
                            z=float(l.split()[2])
                            tot=float(l.split()[6])
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_dipole_x"]=x
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_dipole_y"]=y
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_dipole_z"]=z
                            data[m+"_"+str(f)+"_"+str(systems[i])+"_dipole_tot"]=tot
                            dipole=False
                fi.close()
                data[m+"_"+str(f)+"_"+str(systems[i])+"_nsm"]=nsm
                data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"]=dim
                data[m+"_"+str(f)+"_"+str(systems[i])+"_tsolve2"]=tsolve2
                data[m+"_"+str(f)+"_"+str(systems[i])+"_tsm"]=tsm
                data[m+"_"+str(f)+"_"+str(systems[i])+"_de"]=de
                data[m+"_"+str(f)+"_"+str(systems[i])+"_mse"]=mse
                data[m+"_"+str(f)+"_"+str(systems[i])+"_inf"]=inf
            except Exception as error:
                print(traceback.format_exc())
                print("An exception occurred:", type(error).__name__,error)                
                pass
for d in data:
    print(d,data[d])
#from scipy.optimize import curve_fit
#import numpy as np
#def func(x, a, b):
#    return a * x**b
#
#x=np.zeros(len(sygvd))
#y=np.zeros(len(sygvd))
#
#x=sorted(sygvd)
#i=0
#for xi in x:
#    y[i]=sygvd[xi]
#    i=i+1
#print(x,y)
#
#popt, pcov = curve_fit(func, x[8:], y[8:])
#print(popt)

#plot
col=[]
for i in mcolors.TABLEAU_COLORS:
    col.append(i)
ph=4
pw=ph*0.6666
bb=0.02

fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'Speedup')
ax.set_xscale("log")
#ax.set_yscale("log")

c=0
f=0.0001
m="full_matrix"
x0=[]
y0=[]
for i in range(len(systems)):
    x0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
    y0.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_tsolve2"])
x0, y0 = zip(*sorted(zip(x0, y0)))
#ax.plot(x0,y0,linewidth=1,label="truncation "+str(f),color=col[c])
#c=c+1

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_tsm"])
    x, y = zip(*sorted(zip(x, y)))
    y=list(y)
    for i in range(len(systems)):
        y[i]=y0[i]/y[i]
    ax.plot(x,y,linewidth=1,label=truncstring(f),color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_speedup_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_speedup_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'Number of Submatrices')
ax.set_xscale("log")
#ax.set_yscale("log")
c=0

f="0.0001"
m="full_matrix"

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_nsm"])
    x, y = zip(*sorted(zip(x, y)))
    ax.plot(x,y,linewidth=1,label=truncstring(f),color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_nsm_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_nsm_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)


#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph,pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'$|\Delta E_\mathrm{tr}$| / a.u.')
ax.set_xscale("log")
ax.set_yscale("log")
plt.ylim(1e-7,1e-4)
c=0

f="0.0001"
m="full_matrix"

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_de"])
    x, y = zip(*sorted(zip(x, y)))
    ax.plot(x,y,linewidth=1,label=truncstring(f),color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.tight_layout()
fig.savefig("submatrix_error_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_error_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'MSE of Density Matrix')
ax.set_xscale("log")
ax.set_yscale("log")
c=0

f="0.0001"
m="full_matrix"

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_mse"])
    x, y = zip(*sorted(zip(x, y)))
    ax.plot(x,y,linewidth=1,label=truncstring(f),color=col[c],marker=".",markersize=4)
    c=c+1


plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.ylim(1e-12,1e-6)
plt.tight_layout()
fig.savefig("submatrix_errorP_MSE_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_errorP_MSE_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)
#
#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'Number of Basis Functions')
ax.set_ylabel(r'Inf-Norm of Density Matrix')
ax.set_xscale("log")
ax.set_yscale("log")
c=0

for f in thres:
    m="kmeans_r_heuristic"
    x=[]
    y=[]
    for i in range(len(systems)):
        x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
        y.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_inf"])
    x, y = zip(*sorted(zip(x, y)))
    ax.plot(x,y,linewidth=1,label=truncstring(f),color=col[c],marker=".",markersize=4)
    c=c+1

plt.legend(loc='upper left',ncol=1,fontsize="7")
plt.ylim(1e-6,1e-1)
plt.tight_layout()
fig.savefig("submatrix_errorP_inf_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_errorP_inf_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'pop log10(|ref-sm|)')
ax.set_ylabel(r'Counts')
#ax.set_xscale("log")
ax.set_yscale("log")
c=0


for f in thres:
    m1="kmeans_r_heuristic"
    m2="full_matrix"
    x=[]
    y=[]
    for i in range(len(systems)):
        try:
            for j in range(len(data[m1+"_"+str(f)+"_"+str(systems[i])+"_pop"])):
                #x.append(data[m1+"_"+str(f)+"_"+str(systems[i])+"_pop"][j])
                y.append(math.log10(1e-20+abs(data[m1+"_"+str(f)+"_"+str(systems[i])+"_pop"][j]-data[m2+"_"+str(f)+"_"+str(systems[i])+"_pop"][j])))
        except:
            pass
    ax.hist(y,linewidth=1,label=truncstring(f),color=col[c],range=(-6,-1),bins=100,alpha=0.5) #,marker=".")#,markersize=1)
    c=c+1

plt.legend(loc='upper right',ncol=1,fontsize="7")
plt.xlim(-6,-1)
#plt.ylim(1e-6,1e-1)
plt.tight_layout()
fig.savefig("submatrix_pop_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_pop_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
fig = plt.figure()
fig, ax = plt.subplots(figsize=(ph, pw))
ax.set_xlabel(r'bonds log10(|ref-sm|)')
ax.set_ylabel(r'Counts')
#ax.set_xscale("log")
ax.set_yscale("log")
c=0


for f in thres:
    m1="kmeans_r_heuristic"
    m2="full_matrix"
    x=[]
    y=[]
    for i in range(len(systems)):
        try:
            for j in range(len(data[m1+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"])):
                if len(data[m1+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j])!=len(data[m2+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j]):
                    print(data[m1+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j],data[m2+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j])
                else:
                    for k in range(len(data[m1+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j])):
                        y.append(math.log10(1e-20+abs(data[m1+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j][k]-data[m2+"_"+str(f)+"_"+str(systems[i])+"_bond_atoms"][j][k])))
        except:
            pass
    ax.hist(y,linewidth=1,label=truncstring(f),color=col[c],range=(-6,-1),bins=100,alpha=0.5) #,marker=".")#,markersize=1)
    c=c+1

plt.legend(loc='upper right',ncol=1,fontsize="7")
plt.xlim(-6,-1)
#plt.ylim(1e-6,1e-1)
plt.tight_layout()
fig.savefig("submatrix_bonds_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
fig.savefig("submatrix_bonds_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
plt.close(fig)

#######################################################################################
for d in ["x","y","z","tot"]:
    fig = plt.figure()
    fig, ax = plt.subplots(figsize=(ph, pw))
    ax.set_xlabel(r'Number of Basis Functions')
    ax.set_ylabel(r'dipole moment '+d+' |ref-sm|/a.u. per atom')
    ax.set_xscale("log")
    ax.set_yscale("log")
    c=0

    for f in thres:
        m1="kmeans_r_heuristic"
        m2="full_matrix"
        x=[]
        y=[]
        for i in range(len(systems)):
            x.append(data[m+"_"+str(f)+"_"+str(systems[i])+"_dim"])
            y.append(abs(data[m1+"_"+str(f)+"_"+str(systems[i])+"_dipole_"+d]-data[m2+"_"+str(f)+"_"+str(systems[i])+"_dipole_"+d])/len(data[m2+"_"+str(f)+"_"+str(systems[i])+"_atom"]))
#        ax.hist(y,linewidth=1,label=truncstring(f),color=col[c],range=(-6,-1),bins=100,alpha=0.5) #,marker=".")#,markersize=1)
        x, y = zip(*sorted(zip(x, y)))
        ax.plot(x,y,linewidth=1,label=truncstring(f),color=col[c],marker=".",markersize=4)
        c=c+1

    plt.legend(loc='upper left',ncol=1,fontsize="7")
    #plt.xlim(-6,-1)
    plt.ylim(1e-8,1e-2)
    plt.tight_layout()
    fig.savefig("submatrix_dipole_"+d+"_"+rundir+".pdf",dpi=600,bbox_inches = 'tight',pad_inches = bb)
    fig.savefig("submatrix_dipole_"+d+"_"+rundir+".png",dpi=600,bbox_inches = 'tight',pad_inches = bb)
    plt.close(fig)
