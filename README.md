Petite note : La mise en place du répertoire Git  (dépot du projet sur git+ README.md) a été en retard, je suis vraiment désolé. Le script python (projet_transmb.py) n'a pas été modifié après 14h.

Le programme principale se nomme : projet_transmb.py 

Prérequis : 
-dssp 3.0.0

-biopython 1.79

-matplotlib 3.3.4

-mayavi 4.7.2

-pandas 1.3.2

-numpy 1.20.3

Usage : 
```
$ python projet_transmb.py -i chemin_du_fichier_pdb 
$ python projet_transmb.py -i /home/user/Projet/2n90.pdb
$ python projet_transmb.py -i 6cmo.pdb  
```

Exemple de sortie (pour PDB id : 6cmo) :  
Sur le longueur du vecteur, le meilleure score d'hydrophobicité est de : 1.0, la position de la membrane se trouve à : -0.38035145642583484*x+0.870005222050981*y+0.3137254901960783*z+0



Lorsqu'on est orthogonal au vecteur, le meilleure score d'hydrophobicité est de : 1.0, la position de la membrane se trouve à : -0.38035145642583484*x+0.870005222050981*y+0.3137254901960783*z+0

La visualisation s'ouvre dans une nouvelle fenêtre
