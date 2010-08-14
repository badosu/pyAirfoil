#!/usr/bin/env python

#IMPORTACAO------------------------------------
from pylab import *
from scipy.integrate import quad
from scipy.io.numpyio import fwrite, fread
#----------------------------------------------

#CONFIGURACOES N: numero de paineis -----------
N = 50
twopi = 2*pi
#----------------------------------------------

#CONSTANTES a, alpha, beta, uinf
a=0.75
alphadg = 7
alpha=radians(alphadg)
beta=-0.07
uinf = 1
#----------------------------------------------

#FUNCOES --------------------------------------

def Zeta(theta):
    return (cos(theta)+beta)+(sin(theta)+alpha)*1j;

def Zhukovsky(zeta):
    return zeta+a*pow(zeta,-1)

#----------------------------------------------

#DEFINE theta, zeta, z-------------------------
theta = arange(0, twopi, twopi/N)
print 'len(theta)={0}'.format(len(theta))
print 'N = {0}, 2pi/N={1}'.format(N, twopi/N)
zeta = Zeta(theta)
z = Zhukovsky(zeta)
#----------------------------------------------

#EXIBICAO INICIAL------------------------------
figure(1)
plot(zeta.real, zeta.imag,'g',label=r'$\zeta(\theta) = \xi + i\eta$')
plot(z.real, z.imag,'r-o',label=r'$z(\theta) = x + iy$')
plot(beta, alpha, 'g+')                    
grid()
legend()
axis('scaled')
title('Aerofolio fino')
xlabel(r'$a = {0},\, \alpha = {1} graus, \beta = {2}$'.format(a, alphadg, beta)) 
plot(0,0,'k+',ms=10)
#-----------------------------------------------

#CALCULA MATRIX DE COEF DE INFLUENCIA-----------
A = mat(zeros((N, N)))
for i in range(N):
	for j in range(N):
		A[i,j] = (quad(lambda t: log(abs((z[i]+(z[(i+1)%N])/2-Zhukovsky(Zeta(t))))), theta[j], theta[(j+1)%N])[0])/twopi
#-----------------------------------------------

#RESOLVE SISTEMA LINEAR-------------------------
gamma = zeros(N); y = z.imag
#-----------------------------------------------

#EXIBE TUDO-------------------------------------
show()
#-----------------------------------------------
