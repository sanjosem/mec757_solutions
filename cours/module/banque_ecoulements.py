# Banque d'ecoulements elementaires
import numpy as np

# Fonctions pour creer des grilles du plan d'ecoulement
# Tous les parametres sont optionnels
# Ils permettent de definir les bornes du plan et le nombre de points de discretisation
# Les objets grid sont des dictionnaires contenant les tableaux de coordonnees cartesiennes

def create_2Dgrid_cart(xbounds=[-5.,5.],ybounds=[-3.,3],nx=500,ny=300):
    grid = dict()
    x = np.linspace(xbounds[0],xbounds[1],nx)
    y = np.linspace(ybounds[0],ybounds[1],ny)
    grid['x'],grid['y'] = np.meshgrid(x,y)
    return grid
    
def create_2Dgrid_cyl(center=[0.,0.],alpha=0.0,rbounds=[1.0e-3,6.],nr=500,nt=360):
    grid = dict()
    rad = np.linspace(rbounds[0],rbounds[1],nr)
    theta = np.linspace(0,2*np.pi,nt)
    R,T = np.meshgrid(rad,theta)
    grid['x'] = R * np.cos(T + alpha) + center[0]
    grid['y'] = R * np.sin(T + alpha) + center[1]
    grid['alpha'] = alpha
    return grid

# Fonction pour calculer les coordonneers polaire a partir des coordonnees cartesienne
# Il est possible de redefinir le centre du repere polaire a l'aide du parametre
# optionnel center. Par defaut il est en [0,0]
    
def cart2cyl(grid,center=[0,0]):
    if 'alpha' in grid.keys():
        alpha = grid['alpha']
    else:
        alpha = 0.
    # print(f"Using AoA: {alpha} rad")
    rad = np.sqrt((grid['x'] - center[0])**2 + (grid['y'] - center[1])**2)
    theta = np.arctan2(grid['y']-center[1],grid['x']-center[0]) - alpha
    return rad,theta



# Les objets ecoul sont des dictionnaires contenant les tableaux permettant 
# de definir les ecoulements elementaires sur la grille du plan fournie. 
# Ils contiennent en particulier les valeurs sur la grille de 
# la fonction courant psi, la fonction potentielle phi, les composantes de vitesse


# Ecoulement uniforme d'incidence alpha en RADIANS ! 

def uniform(grid,Vinf,center=[0.,0.]):
    ecoulement = dict()

    if 'alpha' in grid.keys():
        alpha = grid['alpha']
    else:
        alpha = 0.
    rad,theta = cart2cyl(grid,center)

    ecoulement['phi'] = Vinf*np.cos(alpha)*(grid['x']-center[0]) + Vinf*np.sin(alpha)*(grid['y']-center[1])
    ecoulement['psi'] = -Vinf*np.sin(alpha)*(grid['x']-center[0]) + Vinf*np.cos(alpha)*(grid['y']-center[1])
    ecoulement['u'] = Vinf * np.ones_like(rad)
    ecoulement['v'] = np.zeros_like(rad)
    ecoulement['ur'] = ecoulement['u'] * np.cos(theta) + ecoulement['v'] * np.sin(theta) 
    ecoulement['ut'] = -ecoulement['u'] * np.sin(theta) + ecoulement['v'] * np.cos(theta) 

    return ecoulement


# Source ponctuelle

def source(grid,La,center=[0,0]):
    ecoulement = dict()

    rad,theta = cart2cyl(grid,center)

    ecoulement['phi'] = La/(2*np.pi) * np.log(rad)
    ecoulement['psi'] = La/(2*np.pi) * theta
    ecoulement['ur'] = La/(2*np.pi*rad) 
    ecoulement['ut'] = np.zeros_like(rad)
    ecoulement['u'] = ecoulement['ur'] * np.cos(theta) - ecoulement['ut'] * np.sin(theta)
    ecoulement['v'] = ecoulement['ur'] * np.sin(theta) + ecoulement['ut'] * np.cos(theta)

    return ecoulement


# Dipole ponctuel

def dipole(grid,Ka,center=[0,0]):
    ecoulement = dict()

    rad,theta = cart2cyl(grid,center)

    ecoulement['phi'] = Ka/(2*np.pi) * np.cos(theta)/rad
    ecoulement['psi'] = -Ka/(2*np.pi) * np.sin(theta)/rad
    ecoulement['ur'] = -Ka/(2*np.pi) * np.cos(theta)/rad**2
    ecoulement['ut'] = -Ka/(2*np.pi) * np.sin(theta)/rad**2
    ecoulement['u'] = ecoulement['ur'] * np.cos(theta) - ecoulement['ut'] * np.sin(theta)
    ecoulement['v'] = ecoulement['ur'] * np.sin(theta) + ecoulement['ut'] * np.cos(theta)

    return ecoulement


# Tourbillon ponctuel

def tourbillon(grid,Ga,center=[0,0],R0=1.0):
    ecoulement = dict()

    rad,theta = cart2cyl(grid,center)

    ecoulement['phi'] = -Ga/(2*np.pi) * theta
    ecoulement['psi'] = Ga/(2*np.pi) * np.log(rad/R0)
    ecoulement['ur'] = np.zeros_like(rad)
    ecoulement['ut'] = - Ga/(2*np.pi * rad)
    ecoulement['u'] = ecoulement['ur'] * np.cos(theta) - ecoulement['ut'] * np.sin(theta)
    ecoulement['v'] = ecoulement['ur'] * np.sin(theta) + ecoulement['ut'] * np.cos(theta)

    return ecoulement


# Fonction pour calculer les coordonneers polaire a partir des coordonnees cartesienne
# Il est possible de redefinir le centre du repere polaire a l'aide du parametre
# optionnel center. Par defaut il est en [0,0]
    
def superpose_ecoulement(ecoul1,ecoul2,grid):
    ecoulement = dict()    
    rad,theta = cart2cyl(grid,center=[0.,0.])    
    for var in ['phi','psi','u','v']:
        ecoulement[var] = ecoul1[var] + ecoul2[var]
    ecoulement['ur'] = ecoulement['u'] * np.cos(theta) + ecoulement['v'] * np.sin(theta) 
    ecoulement['ut'] = - ecoulement['u'] * np.sin(theta) + ecoulement['v'] * np.cos(theta) 
    return ecoulement

def create_cylinder_flow(grid,Vinf=1.0,Kappa=1.0,Gamma=0.0,
                         center=[0.,0.],R0=1.0):
    
    # Alpha must be provided in the grid
    if 'alpha' not in grid.keys():
        raise RuntimeError('No alpha set with the grid')
    else: 
        alpha = grid['alpha']
    print(f"Uniform flow around cylinder with :")
    print(f"     -> Vinf {Vinf:.2f} m/s")
    print(f"     -> AoA {alpha:.2f} rad")
    print(f"Circulation around cylinder:")
    print(f"     -> Gamma {Gamma:.2f} m2/s")
    unif = uniform(grid,Vinf,center=center)
    dipol = dipole(grid,Kappa,center=center)
    tourb = tourbillon(grid,Gamma,R0=R0,center=center)
    ecoul = superpose_ecoulement(unif,dipol,grid)
    ecoul = superpose_ecoulement(ecoul,tourb,grid)
    ecoul['V'] = (ecoul['u']**2 + ecoul['v']**2)**0.5
    ecoul['Cp'] = 1 - (ecoul['V']/Vinf)**2
    return ecoul
