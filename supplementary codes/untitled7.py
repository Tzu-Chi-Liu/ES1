import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Discrete Grid Points
# =============================================================================
ng            = 10
NG            = 1000
L             = 2.*np.pi
dx            = L/ng
dX            = L/NG
x             = np.arange(0,L,dx)
X             = np.arange(0,L,dX)
print('dx = %.3f, minimum wavelength = 2*dx = %.3f, maximum wave number = %.3f'
      %(dx,2.*dx,np.pi/dx))

# =============================================================================
# Signal
# =============================================================================
k             = 1
y             = np.sin(k*x)
Y             = np.sin(k*X)

# =============================================================================
# Discrete Gradient operator
# =============================================================================
gradient      = (np.diagflat(np.ones(ng-1),1)-np.diagflat(np.ones(ng-1),-1))*1/(2.*dx)
gradient[0,ng-1]=-1./(2.*dx)
gradient[ng-1,0]=1./(2.*dx)
    
Gradient      = (np.diagflat(np.ones(NG-1),1)-np.diagflat(np.ones(NG-1),-1))*1/(2.*dX)
Gradient[0,NG-1]=-1./(2.*dX)
Gradient[NG-1,0]=1./(2.*dX)

# =============================================================================
# Numerical Differentiation
# =============================================================================
dy            = gradient@y
dY            = Gradient@Y

dy_ref        = k*np.cos(k*x)
dY_ref        = k*np.cos(k*X)

# =============================================================================
# Noise
# =============================================================================
k_noise       = 4
lambda_noise  = 2.*np.pi/k_noise

y_noise       = 0.005*np.cos(k_noise*x)
Y_noise       = 0.005*np.cos(k_noise*X)
y_tot         = y+y_noise
Y_tot         = Y+Y_noise

dy_totnum     = gradient@y_tot
dY_totnum     = Gradient@Y_tot
dy_totanly    = np.cos(k*x)-k_noise*0.005*np.sin(k_noise*x)
dY_totanly    = np.cos(k*X)-k_noise*0.005*np.sin(k_noise*X)
print('noise wavelength lambda = %.3f, noise wave number k = %.3f'
      %(lambda_noise,k_noise))

# =============================================================================
# Time Domain Sampling
# =============================================================================
plt.figure()
# plt.plot(x,y,'-o',label='y')
# plt.plot(X,Y,label='Y')
plt.plot(x,y_tot,'-o',label='y_tot')
plt.plot(X,Y_tot,label='Y_tot')
plt.plot(x,dy_totnum,'-o',label='dy_tot_num')
plt.plot(x,dy_totanly,'-o',label='dy_totref')
# plt.plot(x,dy_tot-dy_totref,'--',label='dy_tot-dy_totref')
# plt.plot(X,dY,label='dY')
# plt.plot(X,dY_totref,label='dY_totref')
# plt.plot(X,dY-dY_totref,'--')
plt.legend(loc='lower left')
plt.title(r'dx = L/%i = %.4f, $k_c$ = %.4f, k = %.4f, k_noise = %.4f'
          %(ng,dx,np.pi/dx,k,k_noise))

# =============================================================================
# FFT
# =============================================================================
k_FFT=2*np.pi*np.fft.rfftfreq(ng,dx)
K_FFT=2*np.pi*np.fft.rfftfreq(NG,dX)

y_FFT=np.fft.rfft(y)
Y_FFT=np.fft.rfft(Y)

y_tot_FFT=np.fft.rfft(y_tot)
Y_tot_FFT=np.fft.rfft(Y_tot)

# =============================================================================
# Frequency Domain
# =============================================================================
# plt.figure()
# plt.plot(k_FFT,np.abs(y_FFT))
# plt.plot(K_FFT[0:20],np.abs(Y_FFT)[0:20])
# plt.plot(k_FFT,np.abs(y_tot_FFT))
# plt.plot(K_FFT[0:20],np.abs(Y_tot_FFT)[0:20])