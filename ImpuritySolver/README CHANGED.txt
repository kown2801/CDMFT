Changement d'IS.C, pour ne charger que des ficheirs .json
Rajout du fichier Montecarlo.h pour la boucle, pas de modifications sur ce fichier
Copie du fichier MPIUtilities sans modifications
Lecture des données d'une ancienne configuration, pour pouvoir reprendre la simulation au milieu si besoin. (config<numero>.json, a modifier dans Trace.h pour la lecture et l'ecriture)
Remplacement de la variable measurements
		(dans Link, elle ne fait que passer)
		(dans Trace, elle permet d'enregistrer les variables sur les differents sites)
		(dans Green, elle permet d'enregistrer les fonctions de green (les trois composantes))
Ajout de la class Ut::Measurements et Ut::Simulation dans le fichier Utilities.h
Avec Measurements, même hangement effectué
Dans Green.h probleme de lecture de fichiers (green's function)


Il faut changer Ut::Simulation pour avoir les bons parametres (params())