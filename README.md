# MEC757 Students

Ce dépôt s'adresse aux étudiants du cours de MEC757 _Introduction à l'aérodynamique_ de l'École de Technologie Supérieure (Montréal, Canada, QC).

Il réunit les énoncés et les sources de programmes qui accompagnent les étudiants lors de leurs TP (_laboratoires informatiques_) et leurs devoirs dans leur apprentissage de l'aérodynamique, en particulier de la méthode des panneaux.

Il est créé pour la première session à l'Automne 2021.

Une partie des énoncés s'inspire librement de [barbagroup/AeroPython](https://github.com/barbagroup/AeroPython) : 
Barba, Lorena A., and Mesnard, Olivier (2019). Aero Python: classical aerodynamics of potential flow using Python. _Journal of Open Source Education_, 2(15), 45, https://doi.org/10.21105/jose.00045

## Installation de python 

Pour une installation sur **votre ordinateur personnel** des outils python, je vous conseille d'installer [Anaconda | Individual Edition](https://www.anaconda.com/products/individual-d) ou sa version minimaliste [Miniconda](https://docs.conda.io/en/latest/miniconda.html) avec module `numpy`, `matplotlib` et `jupyter-lab`. 

Nous utiliserons des `notebook` la plupart du temps pour l'interactivité. Ce sont des fichiers `*.pynb` qui contiennent du texte (_markdown_), du code et son execution (affichage). Ils peuvent facilement être exportés (pdf, html,...). Les fonctions et classes réutilisées plus complexes seront implémentées dans des modules séparés (fichier texte`*.py`) que nous pourront facilement charger.

Lors des TP numériques en classe, nous préférerons une version _cloud_ propulsée par The Pacific Institute for the Mathematical Sciences (PIMS) en collaboration et grace à l'infrastructure de Calcul Canada et Cybera. 

[PIMS Jupyter Lab](https://pims.syzygy.ca/)

Pour y accéder, cela vous prend simplement une authentification avec un compte google (personnel ou `etsmtl.net` si vous l'avez activé). Les modules nécessaires sont déjà tous disponibles.


## Liste des TP numériques

### 1. Introduction à Python

Nous balayons les commandes de base de python par analogie à ce que nous connaissons sous matlab. 

[enonces/introduction/IntroNum.ipynb](./enonces/introduction/IntroNum.ipynb)

Vous pouvez utilisez directement l'adresse suivante pour obtenir une copie de travail sur le serveur PIMS: 
https://pims.syzygy.ca/jupyter/user-redirect/git-pull?repo=https://github.com/sanjosem/mec757_students&subPath=enonces/introduction/IntroNum.ipynb&branch=main
