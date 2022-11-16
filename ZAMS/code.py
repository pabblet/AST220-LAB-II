import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.io import ascii
from scipy.optimize import curve_fit

plt.style.use('ggplot')

#1
data= ascii.read('mass.txt', data_start=9)
data['r']= (data['r']*u.m).to(u.Rsun)
Lstar= list(data['L_r'])[0]

#b)
ite=0
k=0
for i in list(data['L_r']):
    if i <= 0.99*Lstar:
        ite= list(data['L_r'])[k]
        break
    k+=1
T99= list(data['T'])[k]
print('T99 =', T99, '[K]')

ite=0
k=0
for i in list(data['L_r']):
    if i <= 0.50*Lstar:
        ite= list(data['L_r'])[k]
        break
    k+=1
T50= list(data['T'])[k]
print('T50 =', T50, '[K]')

data['L_r']= (data['L_r']*u.W).to(u.Lsun)

#c)
Z= 1
e= 1.602*(10**(-19))
u_= 1.67*(10**(-27))/2
eo= 8.85*(10**(-12))
h= 6.63*(10**(-34))
k= 1.38*(10**(-23))
Tq= ((Z**2)*(Z**2)*(e**4)*(u_))/((12*(np.pi)**2)*(eo**2)*(h**2)*k)
print('Ecuacion', Tq, '[K]')

#d)
k=0
s=True
flag=True
while k<=len(data['tau'])-1 and s==True:
    if list(data['tau'])[k]>1:
        ite=k #fila donde empieza la photosfera
        s=False
    k=k+1
print('% Radio donde comienza la photosfera:', 100*(data['r'][ite]/data['r'][0]))
print('% Masa encerrada:', 100*1-(data['1-M_r/Ms'])[ite], '[M*]')

#a)
#P, r
plt.scatter(data['r'][:-3], data['P'][:-3], s=15, color='#444444', label='Data')
plt.xlabel('Star radius [Rsun]')
plt.ylabel('Pressure [Pa]')
plt.legend(loc='best')
plt.show()

#Mr, r
plt.scatter(data['r'], 1-data['1-M_r/Ms'], s=15, color='#444444', label='Data')
plt.xlabel('Star Sadius [Rsun]')
plt.ylabel('Star Mass [Msun]')
plt.legend(loc='best')
plt.show()

#Lr, r
plt.scatter(data['r'][:-3], data['L_r'][:-3], s=15, color='#444444', label='Data')
plt.xlabel('Star Sadius [Rsun]')
plt.ylabel('Star Luminosity [Lsun]')
plt.legend(loc='best')
plt.show()

#T, r
plt.scatter(data['r'][:-3], data['T'][:-3], s=15, color='#444444', label='Data')
plt.xlabel('Star Sadius [Rsun]')
plt.ylabel('Temperature [K]')
plt.legend(loc='best')
plt.show()

#2
data= ascii.read('data.txt')
data['M*']= data['M*']*u.Msun
data['L*']= data['L*']*u.Lsun
data['R*']= data['R*']*u.Rsun
data['Rc']= (data['Rc']*u.m).to(u.Rsun)
data['Tc']= data['Tc']*u.K
data['ec']= data['ec']*u.W/u.kg
data['rhoc']= data['rhoc']*u.kg/u.m**3

#c)
Mc= data['M*']*(1-data['m'])
M= Mc/data['M*']
R= data['Rc']/data['R*']

grade= 2
fit= np.polyfit(data['M*'], M, grade)
pol= np.poly1d(fit)
n= 300
StarMasses= np.linspace(data['M*'][0], data['M*'][-1], n)

plt.plot(StarMasses, pol(StarMasses), color='#2FA5E9', label='Ajuste')
plt.scatter(data['M*'], M, s=15, color='#444444', label='Data')
plt.xlabel('Star Mass [Msun]')
plt.ylabel('Fraction Mass')
plt.legend(loc='best')
plt.show()

grade= 2
fit= np.polyfit(data['M*'], R, grade)
pol= np.poly1d(fit)
n= 300
StarMasses= np.linspace(data['M*'][0], data['M*'][-1], n)

plt.plot(StarMasses, pol(StarMasses), color='#2FA5E9', label='Ajuste')
plt.scatter(data['M*'], R, s=15, color='#444444', label='Data')
plt.xlabel('Star Mass [Msun]')
plt.ylabel('Fraction Radius')
plt.legend(loc='best')
plt.show()

#d)
m4= [0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0]
e4= [2.63e-05, 0.000175, 0.000297, 0.00165, 0.00174, 0.0015, 0.00361]
T4= [6800000.0, 9460000.0, 11100000.0, 15500000.0, 16500000.0, 17600000.0, 19700000.0]

fit= np.polyfit(m4, e4, grade)
pol= np.poly1d(fit)
n= 300
StarMasses= np.linspace(m4[0], m4[-1], n)

plt.plot(StarMasses, pol(StarMasses), color='#2FA5E9', label='Ajuste')
plt.scatter(m4, e4, s=15, color='#444444', label='data')
plt.xlabel('Star Mass [Msun]')
plt.ylabel('Energy production rate [W/kg]')
plt.legend(loc='best')
plt.show()


fit= np.polyfit(m4, T4, grade)
pol= np.poly1d(fit)
n= 300
StarMasses= np.linspace(m4[0], m4[-1], n)

plt.plot(StarMasses, pol(StarMasses), color='#2FA5E9', label='Ajuste')
plt.scatter(m4, T4, s=15, color='#444444', label='data')
plt.xlabel('Star Mass [Msun]')
plt.ylabel('Temperature [K]')
plt.legend(loc='best')
plt.show()


fit= np.polyfit(T4, e4, grade)
pol= np.poly1d(fit)
n= 300
StarMasses= np.linspace(T4[0], T4[-1], n)

plt.plot(StarMasses, pol(StarMasses), color='#2FA5E9', label='Ajuste')
plt.scatter(T4, e4, s=15, color='#444444', label='data')
plt.xlabel('Temperature [K]')
plt.ylabel('Energy production rate [W/kg]')
plt.legend(loc='best')
plt.show()

#e)
#fit m<1.5
M=[0.5, 0.7, 1.0]
L=[0.02, 0.13, 0.919]
def modelo(Mstar, a):
    return a*Mstar
fit, cov= curve_fit(modelo, np.log10(M), np.log10(L))
ans= fit[0]
print('alpha 1:', ans)
t= np.linspace(np.log10(M[0]), np.log10(M[-1]), 200)

plt.plot(t, modelo(t, ans), c='#4DADC0', label='m<1.5')
plt.scatter(np.log10(M), np.log10(L), s=15, color='#444444')

#fit 1.5<m<4.5
M=[2.0, 3.0, 4.0]
L=[22, 112, 330]
def modelo(Mstar, a):
    return a*Mstar
fit, cov= curve_fit(modelo, np.log10(M), np.log10(L))
ans= fit[0]
print('alpha 2:', ans)
t= np.linspace(np.log10(M[0]), np.log10(M[-1]), 200)

plt.plot(t, modelo(t, ans), c='#EA5656', label='1.5<m<4.5')
plt.scatter(np.log10(M), np.log10(L), s=15, color='#444444')

#fit m>4.5
M=[8.0, 10.0, 12.0]
L=[3333, 6800, 12500]
def modelo(Mstar, a):
    return a*Mstar
fit, cov= curve_fit(modelo, np.log10(M), np.log10(L))
ans= fit[0]
print('alpha 3:', ans)
t= np.linspace(np.log10(M[0]), np.log10(M[-1]), 200)

plt.plot(t, modelo(t, ans), c='#7AB86B', label='m>4.5')
plt.scatter(np.log10(M), np.log10(L), s=15, color='#444444', label='data')

plt.xlabel('log(M*/Msun)')
plt.ylabel('log(L*/Lsun)')
plt.legend(loc='best')
plt.show()