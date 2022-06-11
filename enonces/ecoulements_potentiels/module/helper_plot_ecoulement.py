import matplotlib.pyplot as plt

def set_nice_plot_params():
  plt.rcParams['figure.dpi'] = 100
  plt.rcParams['font.size'] = 12


def plot_contours_phi_psi(grid,ecoulement,titre=None,Ncontours=25,isoPsi=None,xrange=None,yrange=None):
    """
    X et Y sont des tableaux de coordonnées du plan
    PHI et PSI sont les fonctions potentiel et courant évaluées sur cette grille
    Vous pouvez ajouter un titre du graphique sous la forme d'une chaine de caractères.
    Vous pouvez changez le nombres de contours avec la variable Ncontours
    Vous pouvez ajouter des iso-courant particuliers (en noir) dans la liste isoPsi
    """
    from numpy import asarray
    plt.contour(grid['x'],grid['y'],ecoulement['phi'],Ncontours,cmap=plt.cm.Blues)
    plt.contour(grid['x'],grid['y'],ecoulement['psi'],Ncontours,linestyles='dashed',cmap=plt.cm.Reds)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('image')
    if not titre is None:
        plt.title(titre)
    if not isoPsi is None:
        plt.contour(grid['x'],grid['y'],ecoulement['psi'],asarray(isoPsi),colors='black')
    if not xrange is None:
        plt.xlim(xrange[0],xrange[1])
    if not yrange is None:
        plt.ylim(yrange[0],yrange[1])
    plt.show()

def plot_lignes_courant(grid,ecoulement,titre=None,isoPsi=None,xrange=None,yrange=None):
    """
    X et Y sont des tableaux de coordonnées du plan
    U et V sont les composantes cartésiennes de vitesse évaluées sur cette grille
    Vous pouvez ajouter un titre du graphique sous la forme d'une chaine de caractères.
    Vous pouvez ajouter des iso-courant particuliers (en noir) dans la liste isoPsi
    """
    from numpy import asarray
    plt.streamplot(grid['x'],grid['y'],ecoulement['u'],ecoulement['v'],)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('image')
    if not titre is None:
        plt.title(titre)
    if not isoPsi is None:
        plt.contour(grid['x'],grid['y'],ecoulement['psi'],asarray(isoPsi),colors='black')
    if not xrange is None:
        plt.xlim(xrange[0],xrange[1])
    if not yrange is None:
        plt.ylim(yrange[0],yrange[1])
    plt.show()

def plot_champs(grid,ecoulement,FIELD,titre=None,isoPsi=None,legend=None,range=None,
                xrange=None,yrange=None):
    """
    X et Y sont des tableaux de coordonnées du plan
    FIELD est le champs à afficher sur cette grille
      * soit le nom d'un champ dans ecoulement
      * soit un tableau d'une variable
    Vous pouvez ajouter un titre du graphique sous la forme d'une chaine de caractères.
    Vous pouvez ajouter des iso-courant particuliers (en noir) dans la liste isoPsi
    Vous pouvez ajouter les bornes de la colormap par exemple range=[0.,5.0]
    """
    import numpy as np
    if type(FIELD) is str:
        if FIELD in ecoulement.keys():
            if legend is None:
                legend = FIELD
            FIELD = ecoulement[FIELD]
        else:
            raise RuntimeError("I don't know this field !")
    elif type(FIELD) is np.ndarray:
        if FIELD.shape != grid['x'].shape:
            raise RuntimeError('Not a valid field! Check dimensions')
    else:
        raise RuntimeError('Check what you provided as FIELD input')
    if not range is None:
        levels = np.linspace(range[0],range[1],11)
        cs = plt.contourf(grid['x'],grid['y'],FIELD,levels,cmap=plt.cm.jet)
    else:
        cs = plt.contourf(grid['x'],grid['y'],FIELD)
    ax = plt.gca()
    cbar = plt.colorbar(cs,ax=ax)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('image')
    if not titre is None:
        plt.title(titre)
    if not isoPsi is None:
        plt.contour(grid['x'],grid['y'],ecoulement['psi'],np.asarray(isoPsi),colors='black')
    if not legend is None:
        cbar.set_label(legend) #, rotation=270)
    if not xrange is None:
        plt.xlim(xrange[0],xrange[1])
    if not yrange is None:
        plt.ylim(yrange[0],yrange[1])
    plt.show()
