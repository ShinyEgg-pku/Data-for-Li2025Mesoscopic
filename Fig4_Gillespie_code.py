import numpy as np
import matplotlib.pyplot as plt
import math
def Exp(b):
    F=np.random.uniform(0,1,1)[0]
    return -1/b*np.log(1-F)
def R(E,b):
    return np.exp(-b*E)
def Puga1(gamma0,gamma1,cc,b0,b,N,lock):
    t=np.zeros(N)
    E=np.zeros(N)
    Et=np.zeros((len(cc),N), dtype=np.float32)
    De=np.zeros(N)
    Pu=np.zeros(len(cc))
    for j in range(0,len(cc)):
        for i in range(0,N):
            if t[i]<10**(cc[j]) :
                if E[i]==0:
                    E[i]+=De[i]
                elif E[i]>0 and E[i]!=lock:
                    E[i]=0
            while t[i]<=10**(cc[j]):
                if E[i]==0:
                    r=np.random.rand()
                    if r<=gamma0/(gamma0+gamma1):
                        de=Exp(b0)
                        De[i]=de
                        r1 = np.random.rand()
                        dt=1/(gamma0*R(0,b))*np.log(1/r1)
                        if t[i]+dt<10**(cc[j]):
                            t[i]+=dt
                            E[i]+=de
                        else:
                            t[i]+=dt
                    elif r>gamma0/(gamma0+gamma1):
                        de=lock
                        De[i]=de
                        r1 = np.random.rand()
                        dt=1/(gamma1*R(0,b))*np.log(1/r1)
                        if t[i]+dt<10**(cc[j]):
                            t[i]+=dt
                            E[i]+=de
                        else:
                            t[i]+=dt
                else:
                    r1 = np.random.rand()
                    dt=1/(gamma0*R(E[i],b))*np.log(1/r1)
                    if t[i]+dt<10**(cc[j]):
                        t[i]+=dt
                        E[i]=0
                    else:
                        t[i]+=dt
        Pu[j]=1-np.count_nonzero(E)/len(E)
        Et[j,:]=E
    return Pu, Et
N=10000
gamma0=1
lock=30
b0=2
k=1
b=1
a=b0/b
cc=np.linspace(-2,6,1000)
gamma1=[0.1,0.01,0.005,0.001]
Put=np.zeros((len(cc),len(gamma1)))
plt.figure(figsize=[10, 8])
plt.rcParams['figure.dpi'] = 100

for j in range(len(gamma1)):
    Pu, Et = Puga1(gamma0,gamma1[j],cc,b0,b,N,lock)
    Put[:,j] = Pu  
    plt.loglog(10**cc, Put[:,j], linewidth=4, label='α={},Γ1={}'.format(a,gamma1[j]))

plt.xlabel('t', fontsize=30)
plt.ylabel('$P(t)$', fontsize=30)
plt.tick_params(axis='both', which='major', labelsize=20)  
plt.title('α={}'.format(a), fontsize=30)
plt.grid()
plt.legend(fontsize=18, frameon=True) 
plt.show()

np.savetxt('Pt_a2_E30_N10000.txt', Put, delimiter=',')