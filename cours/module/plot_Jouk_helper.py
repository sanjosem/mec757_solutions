import matplotlib.pyplot as plt
import module.banque_ecoulements as bq
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable    

def Jouk(z):
    return (z+1.0/z)
    
def dJouk(z):
    deriv = (1.0-(1.0/z)**2)
    return deriv    
    
def define_mplt_par():
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['font.size'] = 10
    
def create_cyltrJouk(Vinf=1.0,center=[0.,0.],alpha=0.0,Gamma=0.0,nr=100,nt=720,apply_Kutta=False):

    Zcenter = (center[0] + 1j*center[1]) 
    ZTE = (1.0 + 1j * 0.0) 
    R0 = np.abs(ZTE-Zcenter) 
    beta = np.angle(ZTE-Zcenter)
    Kappa = 2*np.pi*R0**2*Vinf

    if apply_Kutta:
        Gamma = -4.0 * np.pi * Vinf * R0 * np.sin(beta-alpha)

    grid = bq.create_2Dgrid_cyl(rbounds=[R0-1.0e-5,10*R0],center=center,
                                alpha=alpha,nr=nr,nt=nt)
    
    cercle = bq.create_2Dgrid_cyl(rbounds=[R0,R0],center=center,
                                  alpha=alpha,nr=1,nt=nt)

    ecoul = bq.create_cylinder_flow(grid,Vinf=Vinf,Kappa=Kappa,Gamma=Gamma,R0=R0,center=center)
    cercle_ecoul = bq.create_cylinder_flow(cercle,Vinf=Vinf,Kappa=Kappa,Gamma=Gamma,R0=R0,center=center)
    
    return grid,ecoul,cercle,cercle_ecoul
    
def trJecoul(grid,ecoul,Vinf=1.0):
    
    Zgrid = grid['x'] + 1j*grid['y']
    W = ecoul['u'] - 1j * ecoul['v']
    w = W/(dJouk(Zgrid))
    w_prof = np.abs(w)
    Cp_prof = 1-(w_prof/Vinf)**2
    return w_prof,Cp_prof
# 
# 
# Zcercle = (cercle['x'] + 1j * cercle['y']) 
# Prof = Jouk(Zcercle)
# LocTE = Jouk(ZTE)
# 
# W_cercle = cercle_ecoul['u'] - 1j * cercle_ecoul['v']
# w_cercle = W_cercle/(dJouk(Zcercle))
# w_cercle_prof = np.abs(w_cercle)
# Cp_cercle_prof = 1-(w_cercle_prof/Vinf)**2
    

def plot_phi_psi_cylinder(Vinf=1.0,center=[0.,0.],alpha=0.0,Gamma=0.0,nr=100,nt=720,apply_Kutta=False):
    
    Ncontours = 100
    
    grid,ecoul,_,_ = create_cyltrJouk(Vinf=Vinf,center=center,alpha=alpha,Gamma=Gamma,nr=nr,nt=nt,apply_Kutta=apply_Kutta)

    Zcenter = (center[0] + 1j*center[1]) 
    Zgrid = grid['x'] + 1j * grid['y']
    zgrid = Jouk(Zgrid)
    ZTE = 1.0 + 0.0*1j
    LocTE = Jouk(ZTE)

    # contour of the cylinder 
    psi0 = 0. 

    fig,axs = plt.subplots(2,2,figsize=(12,7))
    axs[0,0].plot(Zcenter.real,Zcenter.imag,'x',color='red')
    axs[0,0].contour(grid['x'],grid['y'],ecoul['psi'],Ncontours,cmap=plt.cm.Reds)
    axs[0,0].contour(grid['x'],grid['y'],ecoul['phi'],Ncontours,linestyles='--',cmap=plt.cm.Blues)
    axs[0,0].contour(grid['x'],grid['y'],ecoul['psi'],[psi0,],colors='black')
    axs[1,0].plot(Zcenter.real,Zcenter.imag,'x',color='red')
    axs[1,0].contour(grid['x'],grid['y'],ecoul['psi'],Ncontours,cmap=plt.cm.Reds)
    axs[1,0].contour(grid['x'],grid['y'],ecoul['phi'],Ncontours,linestyles='--',cmap=plt.cm.Blues)
    axs[1,0].contour(grid['x'],grid['y'],ecoul['psi'],[psi0,],colors='black')

    # axs[1].plot([-np.pi,np.pi],[0,0],color='grey')
    axs[0,1].contour(zgrid.real,zgrid.imag,ecoul['psi'],Ncontours,cmap=plt.cm.Reds)
    axs[0,1].contour(zgrid.real,zgrid.imag,ecoul['phi'],Ncontours,linestyles='--',cmap=plt.cm.Blues)
    axs[0,1].contour(zgrid.real,zgrid.imag,ecoul['psi'],[psi0,],colors='black')
    axs[1,1].contour(zgrid.real,zgrid.imag,ecoul['psi'],Ncontours,cmap=plt.cm.Reds)
    axs[1,1].contour(zgrid.real,zgrid.imag,ecoul['phi'],Ncontours,linestyles='--',cmap=plt.cm.Blues)
    axs[1,1].contour(zgrid.real,zgrid.imag,ecoul['psi'],[psi0,],colors='black')

    axs[0,0].plot(ZTE.real,ZTE.imag,'+',color='green')
    axs[0,1].plot(LocTE.real,LocTE.imag,'+',color='green')
    axs[1,0].plot(ZTE.real,ZTE.imag,'+',color='green')
    axs[1,1].plot(LocTE.real,LocTE.imag,'+',color='green')

    for ax in axs[:,0]:
        ax.axis('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(-2,2)    
        ax.set_ylim(-2,2)
        ax.axvline(Zcenter.real,color='red',linestyle='dotted')
        ax.axhline(Zcenter.imag,color='red',linestyle='dotted')
        ax.axvline(ZTE.real,color='green',linestyle='dotted')
        ax.axhline(ZTE.imag,color='green',linestyle='dotted')
    for ax in axs[:,1]:
        ax.axis('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim(-4,4)    
        ax.set_ylim(-1,1)
    axs[1,0].set_xlim(0.8,1.2)    
    axs[1,0].set_ylim(-0.2,0.2)
    axs[1,1].set_xlim(1.8,2.2)    
    axs[1,1].set_ylim(-0.5,0.5)
    plt.show()



def plot_vel_cylinder(Vinf=1.0,center=[0.,0.],alpha=0.0,Gamma=0.0,nr=100,nt=720):
    
    grid,ecoul,_,_ = create_cyltrJouk(Vinf=Vinf,center=center,alpha=alpha,Gamma=Gamma,nr=nr,nt=nt,apply_Kutta=True)
    
    w_prof,_ = trJecoul(grid,ecoul,Vinf=Vinf)

    Zcenter = (center[0] + 1j*center[1]) 
    Zgrid = grid['x'] + 1j * grid['y']
    zgrid = Jouk(Zgrid)
    ZTE = 1.0 + 0.0*1j
    LocTE = Jouk(ZTE)

    # contour of the cylinder 
    psi0 = 0. 


    fig,axs = plt.subplots(1,2,figsize=(12,3))
    axs[0].plot(Zcenter.real,Zcenter.imag,'x',color='red')
    cf = axs[0].contourf(grid['x'],grid['y'],ecoul['V'],cmap=plt.cm.jet)
    axs[0].contour(grid['x'],grid['y'],ecoul['psi'],[psi0,],colors='black',
                   linewidths=0.5)
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cf, cax=cax)

    cf2 = axs[1].contourf(zgrid.real,zgrid.imag,w_prof,linestyles='--',cmap=plt.cm.jet)
    axs[1].contour(zgrid.real,zgrid.imag,ecoul['psi'],[psi0,],colors='black',
                   linewidths=0.5)
    divider = make_axes_locatable(axs[1])
    cax2 = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cf2, cax=cax2)

    axs[0].plot(ZTE.real,ZTE.imag,'+',color='red')
    axs[1].plot(LocTE.real,LocTE.imag,'+',color='red')

    for ax in axs:
        ax.axis('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    axs[0].set_xlim(-3,3)
    axs[0].set_ylim(-2,2)
    axs[1].set_ylim(-2,2)
    axs[1].set_xlim(-4,4)
    plt.suptitle(r'Norme de la vitesse $V$')
    plt.show()

def plot_press_cylinder(Vinf=1.0,center=[0.,0.],alpha=0.0,Gamma=0.0,nr=100,nt=720):
    
    grid,ecoul,_,_ = create_cyltrJouk(Vinf=Vinf,center=center,alpha=alpha,Gamma=Gamma,nr=nr,nt=nt,apply_Kutta=True)
    
    _,Cp_prof = trJecoul(grid,ecoul,Vinf=Vinf)

    Zcenter = (center[0] + 1j*center[1]) 
    Zgrid = grid['x'] + 1j * grid['y']
    zgrid = Jouk(Zgrid)
    ZTE = 1.0 + 0.0*1j
    LocTE = Jouk(ZTE)

    # contour of the cylinder 
    psi0 = 0. 


    fig,axs = plt.subplots(1,2,figsize=(12,3))
    axs[0].plot(Zcenter.real,Zcenter.imag,'x',color='red')
    cf = axs[0].contourf(grid['x'],grid['y'],ecoul['Cp'],cmap=plt.cm.jet)
    axs[0].contour(grid['x'],grid['y'],ecoul['psi'],[psi0,],colors='black',
                   linewidths=0.5)
    divider = make_axes_locatable(axs[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cf, cax=cax)

    cf2 = axs[1].contourf(zgrid.real,zgrid.imag,Cp_prof,linestyles='--',cmap=plt.cm.jet)
    axs[1].contour(zgrid.real,zgrid.imag,ecoul['psi'],[psi0,],colors='black',
                   linewidths=0.5)
    divider = make_axes_locatable(axs[1])
    cax2 = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(cf2, cax=cax2)

    axs[0].plot(ZTE.real,ZTE.imag,'+',color='red')
    axs[1].plot(LocTE.real,LocTE.imag,'+',color='red')

    for ax in axs:
        ax.axis('equal')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    axs[0].set_xlim(-3,3)
    axs[0].set_ylim(-2,2)
    axs[1].set_ylim(-2,2)
    axs[1].set_xlim(-4,4)
    plt.suptitle(r'Coefficient Pression $C_p$')
    plt.show()



def plot_Cp_peau(Vinf=1.0,center=[0.,0.],alpha=0.0,Gamma=0.0,nr=100,nt=720):
    
    _,_,cercle,cercle_ecoul = create_cyltrJouk(Vinf=Vinf,center=center,alpha=alpha,Gamma=Gamma,nr=nr,nt=nt,apply_Kutta=True)
    
    _,Cp_cercle_prof = trJecoul(cercle,cercle_ecoul,Vinf=Vinf)

    Zcercle = (cercle['x'] + 1j * cercle['y']) 
    Prof = Jouk(Zcercle)
    
    fig,axs = plt.subplots(1,2,figsize=(12,3))
    axs[0].plot(Prof.real[:-1]/2.,Prof.imag[:-1]/2.) #,marker='o')
    axs[0].set_xlabel(r'x/2')
    axs[0].set_ylabel(r'y/2')
    axs[0].axis('equal')
    axs[1].plot(Prof.real[:-1]/2.,Cp_cercle_prof[:-1]) #,marker='o')
    # plt.ylim(-2.0,1.2)
    axs[1].set_xlabel(r'x/2')
    axs[1].set_ylabel(r'$C_p$')
    axs[1].invert_yaxis()
    axs[1].grid()
    plt.show()
