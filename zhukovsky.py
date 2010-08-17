#!/usr/bin/env python

#IMPORTS---------------------------------------
from pylab import *
from scipy.integrate import quad
#----------------------------------------------

#CONFIGS N: number of panels ------------------
N=50
twopi = 2*pi
#----------------------------------------------

#CONSTANT a, alpha, beta, uinf
a = 0.75
alphadg = 5
alpha = radians(alphadg)
beta=-0.05
Uinf = 1
#----------------------------------------------

#FUNCTIONS-------------------------------------

def Zeta(theta):
	return (cos(theta)+beta)+(sin(theta)+alpha)*1j;

def Zhukovsky(zeta):
	return zeta+a*pow(zeta,-1)

#----------------------------------------------

#DEFINE theta, zeta, z-------------------------
theta = arange(0, twopi, twopi/N)
zeta = Zeta(theta)
z = Zhukovsky(zeta)
x = z.real; y = z.imag
#----------------------------------------------

#INITIAL VISUALIZATION-------------------------

figure()
plot(zeta.real, zeta.imag,'g',label=r'$\zeta(\theta) = \xi + i\eta$')
plot(x, y,'r-o',label=r'$z(\theta) = x + iy$')
plot(beta, alpha, 'g+')                    
grid()
legend()
axis('scaled')
title('Aerofolio fino')
xlabel(r'$a = {0},\, \alpha = {1} graus, \beta = {2}$'.format(a, alphadg, beta)) 
plot(0, 0, 'k+', ms=10)
#-----------------------------------------------

#GET COEFFICIENTS INFLUENCE MATRIX AND APPLY KUTTA CONDITION
A = mat(zeros((N+1,N+1)))
A[N,0]  =  1
A[N,N-1]= -1

for i in range(N):
	A[i,N]=1
	for j in range(N):
		if i==j :
			A[i,i]=0.5
		else:
			integrand = lambda t: log(abs((z[i]+z[(i+1)%N])/2-Zhukovsky(Zeta(t))))
			A[i,j] = quad(integrand, theta[j], theta[(j+1)%N])[0]/twopi


b = Uinf*(concatenate((y,[0])))
#-----------------------------------------------

#SOLVE LINEAR SYSTEM TO OBTAIN gamma AND Cp-----
sol = solve(A, b)

gamma = sol[0:N]
C = sol[N]

print 'C={0}'.format(C)

Cp = 1-pow(gamma,2)/(Uinf**2)
#-----------------------------------------------

#PLOT EVERYTHING, SHOW AND SAVE-----------------

figure()
xn = x - min(x)
plot(xn/max(xn), -Cp, 'b-o', label = r'$-C_p(x)$')
legend()

show()
#-----------------------------------------------
