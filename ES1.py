import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# physical constants
# =============================================================================
e_e                 = 5e-2       # elementary charge
epsilon0            = 1.0        # vacuum permittivity
m_e                 = 1.0        # electron mass
m_p                 = 2000.0     # proton mass

# =============================================================================
# simulation parameters
# =============================================================================
N                   = 40000      # number of particles
L                   = 4.0        # length of box
NG                  = 200         # number of grids
dx=L/NG
x=np.arange(0,L,dx)

T_end               = 20.0      # total simulation time
dt                  = 0.01       # time step
T=np.arange(0,T_end,dt)

weight              = 'CIC'
scheme              = 'leapfrog'
solver              = 'FFT'

InitialCondition    = 'Two Stream Instability'

Tracker             = 1          # number of tracker particles
plot_animation      = True       # plot animation of snapshots
plot_history        = True       # plot history
plot_omegak         = True       # plot dispersion relation
plot_trajectory     = False      # plot trajectory of tracker particle(s)

omega_plasma        = np.sqrt(((N/(2*L))*e_e**2)/(epsilon0*m_e))
print('Assuming half of N are electrons, other half are protons,\n'+
      'plasma frequency = %.4f (s^{-1}),\nomega_plasma*dt = %.4f, total steps = %i'
      %(omega_plasma,omega_plasma*dt,int(T_end/dt)))

# =============================================================================
# initial conditions
# =============================================================================
if InitialCondition=='Particle Pair Oscillation':
    N=2
    m=np.array([m_e,m_e])
    q=np.array([e_e,-e_e])
    r=np.array([L/4.,2*L/4.])
    v=np.array([0.,0.])
    
if InitialCondition=='Plasma Oscillation':
    A          = 1e-5
    v0         = 1.         # velocity
    k          = 2.*np.pi/L # velocity wavenumber 
    n0         = 0.5*N/L
    v0_sigma   = 0.0
    
    # first half electrons, second half protons
    m=np.ones(N)*m_e
    m[N//2:]=m_p
    q=-1.*np.ones(N)*e_e
    q[N//2:]*=-1.
    r=np.append(np.linspace(0,L,N//2),np.linspace(0,L,N//2))
    v=np.append(np.random.normal(0,v0_sigma,size=(1,N//2)),np.zeros(N//2))
    
    # sinusoidal perturbation
    r[:N//2]+=A/(n0*k)*np.sin(k*r[:N//2])
    r=r%L
    
    # excite all modes
    # for i in range(1,NG//2+1):
    #     r[:N//2]+=A*i/(n0)*np.sin(i*k*r[:N//2])
    # r=r%L
    
    # excite first few modes
    # for i in range(1,51):
    #     r[:N//2]+=A/(n0)*np.sin(i*k*r[:N//2])
    # r=r%L
    
    # r[15*N//64:17*N//64]+=A/(n0)

if InitialCondition=='Two Stream Instability':
    v0         = 1.0        # velocity
    v0_sigma   = 0.1        # velocity width
    A          = 0.01       # sinusoidal perturbation amplitude
    k          = 2.*np.pi/L # sinusoidal perturbation wavevector 
    
    m=np.ones(N)*m_e
    q=-1*np.ones(N)*e_e
    # q[N//2:]*=-1.
    r=np.random.random(N)*L
    # v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
    # v[N//2:]*=-1.
    v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
                np.random.normal(-v0,v0_sigma,size=(1,N//2)))
    # v+=A*np.sin(k*r) # add perturbation

if InitialCondition=='Single Beam':
    v0         = 1.0        # velocity
    v0_sigma   = 0.0        # velocity width
    A          = 0.01       # sinusoidal perturbation amplitude
    k          = 2.*np.pi/L # sinusoidal perturbation wavevector 
    
    m=np.ones(N)*m_e
    q=-1*np.ones(N)*e_e
    # q[N//2:]*=-1.
    r=np.random.random(N)*L
    # v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
    # v[N//2:]*=-1.
    v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
    
# =============================================================================
# diagnostics physical quantites
# =============================================================================
R     = np.zeros((1,len(T))) # tracker particle phase space trajectory
V     = np.zeros((1,len(T)))

P     = np.zeros(len(T)) # total momentum
E_D   = np.zeros(len(T)) # drift kinetic energy
E_T   = np.zeros(len(T)) # thermal kinetic energy
E_F   = np.zeros(len(T)) # field energy

grid_history = np.zeros((NG,len(T)))

# =============================================================================
# external fields
# =============================================================================
def external_field(x,t):
    phiX=np.zeros(NG)
    EX=np.zeros(NG)
    
    return phiX,EX

# =============================================================================
# source weighting
# =============================================================================
def weightrho(r,e,weight):
    # Weighting particle positions r to grid density rho_grid
    rho_grid=np.zeros(NG)
    for i in range(N):
        x=int(r[i]//dx)
        if weight=='NGP' or weight=='0th':
            if np.abs(r[i]%dx)<0.5*dx:
                rho_grid[x%NG]+=e[i]/dx
            else:
                rho_grid[(x+1)%NG]+=e[i]/dx
        if weight=='CIC' or weight=='1st' or weight=='energy conserving':
            rho_grid[x%NG]+=e[i]/dx*(dx-r[i]%dx)/dx
            rho_grid[(x+1)%NG]+=e[i]/dx*(r[i]%dx)/dx
            
    return rho_grid

# =============================================================================
# Poisson solver
# =============================================================================
def solve_poisson(r,e,weight,solver='FFT',operator='local'):
    rho_grid=weightrho(r,e,weight)
    
    if solver=='FFT':
        # Use FFT to solve Poisson equation with periodic boundary condition
        rho_k=np.fft.rfft(rho_grid)
        k=2*np.pi*np.fft.rfftfreq(NG,dx)
        k[0]=1.     # to avoid division by 0
        rho_k[0]=0. # set DC component (average) of rho to 0 (average charge density = 0)
        
        # laplacian
        if operator=='nonlocal':
            phi_k=rho_k/(epsilon0*k**2)
        if operator=='local':
            K=k*(np.sin(0.5*k*dx)/(0.5*k*dx)) # Fourier transform of discrete laplacian
            phi_k=rho_k/(epsilon0*K**2)
        phi_grid=np.fft.irfft(phi_k)
    
    # Differentiate potential phi to obtain field E 
    gradient=(np.diagflat(np.ones(NG-1),1)-np.diagflat(np.ones(NG-1),-1))*1/(2.*dx)
    # gradient operator also has periodic BC
    gradient[0,NG-1]=-1./(2.*dx)
    gradient[NG-1,0]=1./(2.*dx)
    E_grid=-gradient@phi_grid
    
    return phi_grid,rho_grid,E_grid,k,rho_k,phi_k

# =============================================================================
# force weighting
# =============================================================================
def weightField(m,q,r,v,E_grid,weight):
    # Weighting grid field E to obtain field E_p at particle position r
    X=int(r//dx)
    if weight=='0th' or weight=='NGP' or weight=='energy conserving':
        if np.abs(r%dx)<(0.5)*dx:    
            E_par=E_grid[X%NG]
        else:
            E_par=E_grid[(X+1)%NG]
    if weight=='1st' or weight=='CIC':
        E_par=E_grid[X%NG]*(dx-r%dx)/dx+E_grid[(X+1)%NG]*(r%dx)/dx
        
    return E_par

# =============================================================================
# integration of EOM
# =============================================================================
def update_motion(m,q,r,v,E_grid,scheme,weight):
    global t 
    E_par=weightField(m,q,r,v,E_grid,weight)
    a=q*E_par/m
    
    if scheme=='leapfrog':
        if t==0:
            v-=a*0.5*dt
            
        v+=a*dt
        
        r+=v*dt
        r=r%L
        
    if scheme=='KDK':
        v+=a*0.5*dt
        
        r+=v*dt
        r=r%L
        
        E_par=weightField(m,q,r,v,E_grid,weight)
        a=q*E_par/m
        v+=a*0.5*dt
        
    if scheme=='DKD':
        r+=v*0.5*dt
        r=r%L
        
        v+=a*dt
        
        r+=v*0.5*dt
        r=r%L
        
    return r,v

# =============================================================================
# diagnostics
# =============================================================================
def diagnostics_animation(r,v,phi_grid,E_grid,rho_grid):
    global fig
    if plot_animation:    
        if t==0:
            fig=plt.figure(figsize=(15,12))
        
        # 
        ax1=fig.add_subplot(221)
        ax3=fig.add_subplot(222)
        ax4=fig.add_subplot(223)
        ax2=fig.add_subplot(224)
        ax8=ax2.twinx()
        
        # labels and title
        ax1.set_xlabel('position x (m)')
        ax1.set_ylabel('velocity v (m/s)')
        ax1.set_title('Phase Space')
        # ax1.set_xlim(0.,L)
        ax1.set_ylim(-1.1*np.max(np.abs(v)),1.1*np.max(np.abs(v)))
        # ax1.set_ylim(-7e-10,7e-10)
        ax2.set_xlabel('position x (m)')
        ax2.set_ylabel(r'potential $\phi$ (V)',color='blue')
        ax8.set_ylabel('field E (V/m)',color='orange')
        ax2.set_title('Potential and Fields')
        ax2.set_ylim(-1.1*np.max(np.abs(phi_grid)),1.1*np.max(np.abs(phi_grid)))
        ax8.set_ylim(-1.1*np.max(np.abs(E_grid)),1.1*np.max(np.abs(E_grid)))
        # ax2.set_ylim(-1.5e-8,1.5e-8)
        # ax8.set_ylim(-1e-7,1e-7)
        # ax2.set_title('Potential')
        ax3.set_xlabel('number of particles f(v)dv')
        ax3.set_ylabel('velocity v (m/s)')
        ax3.set_title('Velocity Distribution')
        ax3.set_ylim(-1.1*np.max(np.abs(v)),1.1*np.max(np.abs(v)))
        # ax3.set_ylim(-7e-10,7e-10)
        ax4.set_xlabel('position x (m)')
        ax4.set_ylabel(r'charge density $\rho$ $(m^{-1})$')
        ax4.set_title('Charge Density')
        ax4.set_ylim(-1.1*np.max(np.abs(rho_grid-np.mean(rho_grid))),
                      1.1*np.max(np.abs(rho_grid-np.mean(rho_grid))))
        # ax4.set_ylim(-6e-7,6e-7)
        # ax4.set_title('Density')
        fig.suptitle(InitialCondition+'\n'+'Snapshot at t = %.3f (s)'%T[t])
        
        # phase space
        ax1.scatter(r[:N//2],v[:N//2],fc=(1,0,0,0.3),s=1)
        ax1.scatter(r[N//2:],v[N//2:],fc=(0,0,1,0.3),s=1)
        # ax1.hexbin(r,v,gridsize=200)
        # tracker particles 
        ax1.scatter(r[N//4],v[N//4],c='black',s=30,lw=1,edgecolor='white') # tracker 1
        ax1.scatter(r[-1],v[-1],c='white',s=30,lw=1,edgecolor='black') # tracker 2
        
        # potential and field at grid positions
        ax2.plot(x,phi_grid,'blue')
        ax8.plot(x,E_grid,'orange')
        ax2.plot(x,np.zeros(NG),'black',lw=1)
        
        # velocity distribution
        ax3.hist(v[:N//2],bins=100,range=((np.min(v),np.max(v))),
                  color='r',orientation='horizontal')
        ax3.hist(v[N//2:],bins=100,range=((np.min(v),np.max(v))),
                  color='b',alpha=0.5,orientation='horizontal')
        ax3.hist(v,bins=100,range=((np.min(v),np.max(v))),color='white',
                  edgecolor='black',lw=1.,orientation='horizontal',histtype='step')
        
        # density at grid positions
        ax4.plot(x,rho_grid-np.mean(rho_grid))
        # ax4.hist(r,bins=NG)
        
        plt.tight_layout()
        plt.pause(0.1)
        plt.clf()
        
        return t
    
# =============================================================================
# history
# =============================================================================
def history(T,E_D,E_T,E_F,P):
    if plot_history:
        fig2,ax5=plt.subplots()
        ax5.plot(T,E_D,label='drift')
        ax5.plot(T,E_T,label='thermal')
        ax5.plot(T,E_F,label='field')
        ax5.plot(T,E_D+E_T+E_F,label='total')
        ax5.set_xlabel('time t (s)')
        ax5.set_ylabel('energy E (J)')
        ax5.legend(bbox_to_anchor=(0,1.01,1,0.2),loc='lower left'
                   ,mode='expand',borderaxespad=0.,ncol=4)
        ax5.set_title('Energy History',y=1.07)
        fig2.canvas.set_window_title('energy')
        
        fig3,ax6=plt.subplots()
        ax6.plot(T,P-P[0])
        ax6.set_xlabel('time t (s)')
        ax6.set_ylabel('total momentum P (m/s)')
        ax6.set_title('Total Momentum Change History')
        fig3.canvas.set_window_title('momentum')
    
    return True

# =============================================================================
# dispersion relation
# =============================================================================
def dispersion_relation(grid_history):
    global x,T
    if plot_omegak:
        # plot grid_history(x,t)
        fig4=plt.figure()
        ax7=fig4.gca(projection='3d')
        x_mesh,T_mesh=np.meshgrid(x,T)
        grid_history=np.transpose(grid_history)
        ax7.plot_surface(T_mesh,x_mesh,grid_history)
        ax7.set_xlabel('time t (s)')
        ax7.set_ylabel('position x (m)')
        ax7.set_zlabel(r'potential $\phi$ (V)')
        ax7.set_title('Field history')
        fig4.canvas.set_window_title('field history')
        
        # FFT grid_history(x,t) to obtain grid_omegak(k,omega)
        grid_omegak=np.fft.rfft2(grid_history)
        k=2*np.pi*np.fft.rfftfreq(NG,dx)
        dk=k[1]-k[0]
        lenk=len(k)
        omega=2*np.pi*np.fft.rfftfreq(len(T),dt)
        domega=omega[1]-omega[0]
        lenomega=len(omega)
        
        # plot grid_omegak(k,omega)
        grid_omegak=grid_omegak[:grid_omegak.shape[0]//2,:]
        grid_omegak=grid_omegak[::-1,:]
        fig5,ax9=plt.subplots()
        ax9.set_xlabel(r'wave number $k$ (1/m)')
        ax9.set_ylabel(r'angular frequency $\omega$ (1/s)')
        ax9.set_title(InitialCondition+' Dispersion Relation')
        fig5.canvas.set_window_title('dispersion relation')
        ax9.imshow(np.abs(grid_omegak)[9*grid_omegak.shape[0]//10:,:]
                   ,extent=(k[0]-0.5*dk,k[-1]-0.5*dk
                            ,omega[0]-0.5*domega,omega[lenomega//10]-0.5*domega))
        
        # theoretical dispersion relation
        if InitialCondition=='Plasma Oscillation':
            ax9.plot(k,omega_plasma*np.cos(0.5*dx*k),'red',lw=0.5,label='theoretical')
            ax9.legend(loc='upper right')
        
        # plt.tight_layout()
    
    return omega,k,grid_omegak

# =============================================================================
# main simulation loop
# =============================================================================
for t in range(len(T)):
    #solve field
    phi_grid,rho_grid,E_grid=solve_poisson(r,q,weight)[0:3]
    k,rho_k,phi_k=solve_poisson(r,q,weight)[3:6]
    
    v_old=np.copy(v)
    diagnostics_animation(r,0.5*(v_old+v),phi_grid,E_grid,rho_grid)
    
    # update particle i velocity v[i] and position r[i] at time t
    for i in range(N):
        r[i],v[i]=update_motion(m[i],q[i],r[i],v[i],E_grid,scheme,weight)
        # tracker particle
        if i==N//4:
            print('r = %.6f(m), v = %.6f(m/s), i = %i, t = %.3f(s)'%(r[i],v[i],t,T[t]))
            R[0,t]=r[i]

    # history
    E_F[t]=0.5*np.dot(phi_grid*dx,rho_grid)
    P[t]=np.sum(m*0.5*(v_old+v))
    E_D[t]=np.sum(0.5*m*v*v_old)
    grid_history[:,t]=E_grid
        
# history
history(T,E_D,E_T,E_F,P)

# dispersion relation
omega,k,grid_omegak=dispersion_relation(grid_history)