/***************************************************************************************************************************************************************************************************************************
 *
 * MODULE:       r.waterbalance
 *
 * AUTHOR(S):    Ian Ondo
 *              
 * PURPOSE:      Ce programme permet de modeliser la redistribution des eaux de ruissellement pour chaque pixels de la zone d'étude à partir de l'equation d'onde diffusive.
 *				 Cette approche consiste a determiner le temps de trajet d'un point de depart vers un point d'arrivee quelconque situe en aval en suivant un chemin d'ecoulement.
 *               Une fonction de reponse basee sur la moyenne et la variance du temps d'ecoulement, est modelisee par la fonction de densite du premier temps de passage.
 *               Elle permet de determiner pour chaque point du paysage la quantite de ruissellement reçu a chaque instant t donne.
 *               Le module calcule pour un pas de temps donne la quantite d'eau drainant depuis chaque pixel vers chaque point situe en aval le long d'un chemin d'ecoulement.				 
 *               La sortie du modele est donc une carte raster representant a un instant t la redistribution laterale d'un flux d'eau le long d'un versant.
 *
 ************************************************************************************************************************************************************************************************************************/

	
	/*
	A FAIRE : ??

	*/

#define _GNU_SOURCE 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <grass/config.h>
#include <grass/gis.h>
//#include <grass/defs/site.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <grass/segment.h>
#include <fcntl.h>

#include "head.h"
#include "Queue.h"
#include "utils.h"

#define _USE_MATH_DEFINES
#define EPS 0.01
#define UNDEF -1 
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define NUM_ELEM(x) (sizeof (x) / sizeof (*(x)))
#define M_DEG_TO_RAD  (M_PI/180.0)
#define M_RAD_TO_DEG  (180.0/M_PI)
#define M_PI_045 (M_PI/4.0)
#define M_PI_090 (M_PI/2.0)
#define M_PI_270 (3.0*(M_PI/2.0))
#define M_PI_360 (2.0*M_PI)
#define DAY_TO_HOUR 24.0
#define SECOND_TO_HOUR 3600.0
#define SEGCOLSIZE 64
#define FREE(x) { if (x) free(x); x = NULL; }

	/* ****************************************************************************************/
	/* Fonction main, execute le programme du module et creer la carte d'accumulation de flux */
	/* ****************************************************************************************/
	
	int main(int argc, char *argv[]){
		parseOptions(argc, argv);
		createSEGMENT();
		Init();
		Process();
		exit(EXIT_SUCCESS);
	}

	/* ****************************************************************/
	/* Cette fonction gere les options et les parametres du programme */
	/* ****************************************************************/
	
	void parseOptions(int argc, char *argv[]){

	/* Initalise l'environnement GIS */
	G_gisinit(argv[0]);													/* Lit l'environnement GRASS et stocke le nom du programme dans G_program_name() */

	/* Initialise le module */
    module = G_define_module();
	G_add_keyword(_("raster"));
	G_add_keyword(_("hydrologie"));
	G_add_keyword(_("bilan hydrique"));
	G_add_keyword(_("flux lateraux"));
	G_add_keyword(_("redistribution"));
	G_add_keyword(_("eau du sol"));
    module->description =
        _("Calcule le bilan hydrique du sol sur chaque pixel de la carte, "
          "en tenant compte ou non des flux latéraux, et selon un pas de temps definit ");
				
	/* Definit les differentes options comme definies dans gis.h */
	
	// ATMOSPHERIC INPUTS
	parm.prec = G_define_standard_option(G_OPT_R_INPUTS);
	parm.prec->key = "Precipitation [L]";
	parm.prec->required = NO;
	parm.prec->guisection = _("Atmosheric maps");
	
	parm.etp = G_define_standard_option(G_OPT_R_INPUTS);
	parm.etp->key = "Evapotranspiration [L]";
	parm.etp->required = NO;	
	parm.etp->guisection = _("Atmosheric maps");
	
	// TOPOGRAPHIC INPUTS	
    parm.altitude = G_define_standard_option(G_OPT_R_INPUT);
	parm.altitude->key = "altitude[L]";
    parm.altitude->description =
        _("Nom de la couche raster d'entree dont les valeurs entieres representent l'altitude en m.");
	parm.altitude->type = TYPE_STRING;
	parm.altitude->required = NO;
	parm.altitude->multiple = NO;
	parm.altitude->guisection = _("Topographic parameters");
	
	parm.slope = G_define_standard_option(G_OPT_R_INPUT);
	parm.slope->key = "slope[%]";
    parm.slope->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la pente en pourcentage.");
	parm.slope->type = TYPE_STRING;
	parm.slope->required = YES;
	parm.slope->multiple = NO;
	parm.slope->guisection = _("Topographic parameters");
	
	// SUBSURFACE INPUTS
	parm.sat = G_define_standard_option(G_OPT_R_INPUT);
	parm.sat->key = "saturation[L^3.L^-3]";
    parm.sat->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la teneur en eau du sol à saturation en m3.m-3.");
	parm.sat->type = TYPE_STRING;
	parm.sat->required = YES;
	parm.sat->multiple = NO;
	parm.sat->guisection = _("Soil hydraulic parameters");
	
	parm.fc = G_define_standard_option(G_OPT_R_INPUT);
	parm.fc->key = "field_capacity[L^3.L^-3]";
    parm.fc->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la teneur en eau du sol à la capacité au champ en m3.m-3.");
	parm.fc->type = TYPE_STRING;
	parm.fc->required = YES;
	parm.fc->multiple = NO;
	parm.fc->guisection = _("Soil hydraulic parameters");
	
	parm.pwp = G_define_standard_option(G_OPT_R_INPUT);
	parm.pwp->key = "permanent_wilting_point[L^3.L^-3]";
    parm.pwp->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la teneur en eau du sol au point de fletrissement permanent en m3.m-3.");
	parm.pwp->type = TYPE_STRING;
	parm.pwp->required = YES;
	parm.pwp->multiple = NO;
	parm.pwp->guisection = _("Soil hydraulic parameters");
	
	parm.rum = G_define_standard_option(G_OPT_R_INPUT);
	parm.rum->key = "maximum_plant_available_water[L]";
    parm.rum->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la quantite maximale d eau disponible pour les plantes en mm.");
	parm.rum->type = TYPE_STRING;
	parm.rum->required = YES;
	parm.rum->multiple = NO;
	parm.rum->guisection = _("Soil hydraulic parameters");
	
	parm.ksat = G_define_standard_option(G_OPT_R_INPUT);
	parm.ksat->key = "saturated_hydraulic_conductivity[L.T-1]";
    parm.ksat->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la conductivité hydraulique du sol à saturation en m.j-1.");
	parm.ksat->type = TYPE_STRING;
	parm.ksat->required = YES;
	parm.ksat->multiple = NO;
	parm.ksat->guisection = _("Soil hydraulic parameters");
	
	parm.depth = G_define_standard_option(G_OPT_R_INPUT);
	parm.depth->key = "soil depth[L]";
    parm.depth->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la profondeur du sol en m");
	parm.depth->type = TYPE_STRING;
	parm.depth->required = YES;
	parm.depth->multiple = NO;
	parm.depth->guisection = _("Soil hydraulic parameters");
	
	//SURFACE/SUBSURFACE INPUTS	
	parm.flow_speeds = G_define_standard_option(G_OPT_R_INPUT);
	parm.flow_speeds->key = "v[L.T-1]";
    parm.flow_speeds->description =
        _("Nom de(s) la couche(s) raster(s) dont les valeurs "
		  "representent la vitesse a laquelle se propage "
		  "une perturbation le long d'un trajet d'ecoulement en m.h-1");
	parm.flow_speeds->type = TYPE_STRING;
	parm.flow_speeds->required = NO;
	parm.flow_speeds->multiple = YES;
	parm.flow_speeds->guisection = _("Soil hydraulic parameters");
	
	parm.flow_disps = G_define_standard_option(G_OPT_R_INPUT);
	parm.flow_disps->key = "D[L².T-1]";
    parm.flow_disps->description =
        _("Nom de(s) la couche(s) raster(s) dont les valeurs "
		  "representent la tendance d'une perturbation a se disperser "
		  "longitudinalement le long d'un trajet d'ecoulement.");
	parm.flow_disps->type = TYPE_STRING;
	parm.flow_disps->required = NO;
	parm.flow_disps->multiple = YES;
	parm.flow_disps->guisection = _("Soil hydraulic parameters");
	
	parm.smax = G_define_standard_option(G_OPT_R_INPUT);
	parm.smax->key = "smax[L]";
    parm.smax->description =
        _("Nom de la couche raster dont les valeurs "
		  "representent la vitesse a laquelle se propage "
		  "une perturbation le long d'un trajet d'ecoulement.");
	parm.smax->type = TYPE_STRING;
	parm.smax->required = NO;
	parm.smax->multiple = NO;
	parm.smax->guisection = _("Soil hydraulic parameters");
	
	parm.w = G_define_standard_option(G_OPT_R_INPUT);
	parm.w->key = "w[-]";
    parm.w->description =_(".");
	parm.w->type = TYPE_STRING;
	parm.w->required = NO;
	parm.w->multiple = YES;	
	parm.w->guisection = _("Soil hydraulic parameters");
	
	//LAND COVER
	parm.waterbodies = G_define_standard_option(G_OPT_R_INPUT);
	parm.waterbodies->key = "waterbodies[-]";
    parm.waterbodies->description =_(".");
	parm.waterbodies->type = TYPE_STRING;
	parm.waterbodies->required = NO;
	parm.waterbodies->multiple = NO;
	parm.waterbodies->guisection = _("Landscape features");	
	
	parm.riparian = G_define_standard_option(G_OPT_R_INPUT);
	parm.riparian->key = "riparian[-]";
    parm.riparian->description =_(".");
	parm.riparian->type = TYPE_STRING;
	parm.riparian->required = NO;
	parm.riparian->multiple = NO;
	parm.riparian->guisection = _("Landscape features");	
	
	//OUTPUT(S)
	parm.outputs = G_define_standard_option(G_OPT_R_OUTPUT);
	parm.outputs->key = "output(s)";
    parm.outputs->description =
    _("Nom de(s) la couche(s) raster(s) de sortie.");
	parm.outputs->type = TYPE_STRING;
	parm.outputs->options = build_outputs_list();
	parm.outputs->required = NO;
	parm.outputs->multiple = YES;
	parm.outputs->answers[0] ="PAW";
	parm.outputs->guisection = _("Outputs maps");
	
	//OPTION(S)
	parm.start = G_define_option();
	parm.start->key = "start";
    parm.start->description =
        _("Numero du jour/mois par lequel le calcul doit commencer.");
	parm.start->type = TYPE_INTEGER;
	parm.start->required = NO;
	parm.start->multiple = NO;
	parm.start->answer = "1";
	parm.start->guisection = _("Settings");
	
	parm.outiter = G_define_option();
	parm.outiter->key = "outiter[-]";
    parm.outiter->description =
        _("Frequence de calcul du bilan hydrique");
	parm.outiter->type = TYPE_INTEGER;
	parm.outiter->required = NO;
	parm.outiter->multiple = NO;
	parm.outiter->answer = "1";
	parm.outiter->guisection = _("Settings");
	
	parm.mem = G_define_option();
	parm.mem->key = "memory";
	parm.mem->type = TYPE_INTEGER;
	parm.mem->key_desc = "value";
    parm.mem->required = NO;
    parm.mem->multiple = NO;
    parm.mem->answer = "300";
    parm.mem->description = _("Memoire maximale a utiliser en MB");
	parm.mem->guisection = _("Settings");
	
	parm.method = G_define_option();
	parm.method->key = "method(s)";
	parm.method->type = TYPE_STRING;
	parm.method->options = build_method_list();
	parm.method->description = _("Methode(s) de calcul du bilan hydrique");
	parm.method->required = YES;
	parm.method->multiple = NO;	
	parm.method->guisection = _("Settings");
	
	parm.algorithm = G_define_option();
	parm.algorithm->key = "algorithm(s)";
	parm.algorithm->type = TYPE_STRING;
	parm.algorithm->options = build_algorithm_list();
	parm.algorithm->description = _("Algorithme(s) de calcul de l'aire de drainage amont");
	parm.algorithm->required = NO;
	parm.algorithm->multiple = NO;
	parm.algorithm->answer = "MFDmd";	
	parm.algorithm->guisection = _("Settings");
	
	parm.init_abs = G_define_option();
	parm.init_abs->key = "Ia[-]";
	parm.init_abs->type = TYPE_STRING;
	parm.init_abs->description = _("Methode(s) de calcul de la quantite d eau evapotranspiree, infiltree,"
								   "stockee a la surface...etc, avant que le ruisselement de surface ne se produise ");
	parm.init_abs->answer = "alternative";
	parm.init_abs->required = NO;
	parm.init_abs->multiple = NO;
	parm.init_abs->options = "traditional, alternative";
	parm.init_abs->guisection = _("Settings");
	
	parm.drainage_times = G_define_option();
    parm.drainage_times->key = "drainage times[T]";
    parm.drainage_times->type = TYPE_DOUBLE;
    parm.drainage_times->key_desc = "temps[h]";
    parm.drainage_times->required = NO;
    parm.drainage_times->multiple = YES;
    parm.drainage_times->answers[0] ="1.0";
	parm.drainage_times->options = "1.0-24.0";
    parm.drainage_times->description = _("Temps de redistribution de l eau du bassin versant en heures");	
	parm.drainage_times->guisection = _("Settings");
	
	parm.convergence = G_define_option();
    parm.convergence->key = "convergence[-]";
    parm.convergence->type = TYPE_DOUBLE;
    parm.convergence->required = NO;
    parm.convergence->multiple = NO;
    parm.convergence->answer = "5.0";
	parm.convergence->options = "1.0-10.0";
	parm.convergence->guisection = _("Settings");
    parm.convergence->description = _("Exposant controlant la distribution (divergence) du flux vers l aval"
									  "Holmgren,  P.  1994.  Multiple  flow  direction  algorithms  for runoff"
									  " modelling   in   grid   based   elevation   models: An empirical  evaluation."  
									  " Hydrological  Processes,  Vol.  8(4),  pp. 327-334");
	
	// parm.lower = G_define_standard_option(G_OPT_R_OUTPUT);
    // parm.lower->key = "lower";
    // parm.lower->type = TYPE_DOUBLE;
    // parm.lower->required = NO;
    // parm.lower->multiple = NO;
    // parm.lower->answer = "0.0";
	// parm.lower->guisection = "t1[T]";
    // parm.lower->description = _("Limite inferieure de la periode de temps requise");
	
	// parm.upper = G_define_standard_option(G_OPT_R_OUTPUT);
    // parm.upper->key = "upper";
    // parm.upper->type = TYPE_DOUBLE;
    // parm.upper->required = NO;
    // parm.upper->multiple = NO;
    // parm.upper->answer = "86400.0";
	// parm.upper->guisection = "t2[T]";
    // parm.upper->description = _("Limite superieure de la periode de temps requise");
	
	flag = G_define_flag();
	flag->key = 'i';
	flag->description = _("Afficher les informations sur le calcul du bilan hydrique" 
	" ainsi que l'espace disque et les besoins en memoire");
	
	flag2 = G_define_flag();
	flag2->key = 'n';
	flag2->description = _("Garder les valeurs nulles dans la carte raster de sortie");
	
	flag3 = G_define_flag();
	flag3->key = 'j';
	flag3->description = _("Calculer le bilan hydrique avec un pas de temps journalier");
	
	flag4 = G_define_flag();
	flag4->key = 'm';
	flag4->description = _("Calculer le bilan hydrique avec un pas de temps mensuel");
	
	flag5 = G_define_flag();
	flag5->key = 'f';
	flag5->description = _("Calculer le bilan hydrique avec un pas de temps spécifique");
	
	flag6 = G_define_flag();
	flag6->key = 'z';
	flag6->description = _("Calculer le bilan hydrique avec prise en compte de la zone alluviale");	
	
    /*  Analyse la ligne de commande */
    if (G_parser(argc, argv))
	{
		exit(EXIT_FAILURE);
	}
	
    /*  Obtenir les parametres d'en-tête de la base de donnees */
	Rast_get_window(&window);
   /* if (Rast_get_window(&window) < 0){
    G_fatal_error(_("Impossible de lire les parametres d'en-tête actuels"));
    }*/
	
	/* Vérifie que les données d'entrée ont le même nombre de valeurs */
	for (i = 0; parm.prec->answers[i]; i++)
	    ;
	num_inputs = i;
		if (num_inputs < 1) 
			G_fatal_error(_("Carte(s) raster(s) non trouvee(s)"));
		
	for (i = 0; parm.etp->answers[i]; i++)
	        ;
	if (num_inputs != i)
		G_fatal_error(_("Les listes des rasters d entree prec= et etp= doivent avoir la meme longueur."));
	
	/* Vérifie les options et calcule le nombre de couches de sortie */
	if(!flag3 && !flag4 && !flag5)
		G_fatal_error(_("Choisissez au moins une des options suivantes:"
						" flag -j : Calcul du bilan hydrique journalier\n"
						" flag -m : Calcul du bilan hydrique mensuel\n"
						" flag -f : Calcul du bilan hydrique avec un pas de temps specifique\n"));
	if(flag3 && flag4 && flag5)
		G_fatal_error(_("Impossible de choisir plus de deux options simultanement."));				

	num_outputs_names = 0;
	int count = 0;
	sum_days = 0;
	for (i = 0; parm.outputs->answers[i]; i++)
		num_outputs = num_outputs_names = i;
		
	num_outputs_names =	(num_outputs_names==0) ? 1 : num_outputs_names;
	
	if(flag3 && !flag4 && !flag5 || flag4 && !flag3 && !flag5){
	options = 1;
	num_outputs *= num_inputs;
	}
	if(flag3 && flag4){
	options = 2;
	i = 0;
	while (num_inputs>sum_days){
		num_outputs = num_outputs_names * ++count;
		sum_days += num_days[i];
		i++;
	}
	if(sum_days<num_inputs || num_inputs<28)
		G_warning(_("Attention !!! La liste des rasters d entree prec= et etp= est inférieure a la liste de rasters attendue en sortie..."));
	else if(sum_days>num_inputs)
		G_warning(_("Attention !!! La liste des rasters d entree prec= et etp= est superieure a la liste de rasters attendue en sortie..."));
	}
	else if (flag3 && flag5 || flag4 && flag5){
	options = 3;
		int outiter = atoi(parm.outiter->answer);					  
		if (sscanf(parm.outiter->answer, "%i", &outiter) != 1 || outiter < 1)
			G_fatal_error(_("Frequence de calcul du bilan hydrique (outiter=) inappropriee : %i"), outiter);

		for(i=0;i<num_inputs;i++){
			if((i+1)%outiter==0){
			num_outputs = num_outputs_names * ++count;
			}
		}			
	}
	/* Verifie si le(s) nom de raster(s) de sortie specifié(s) est(sont) compatible(s) avec l'option de calcul */
	for(i=0; i<num_outputs_names; i++){
		int indice = find_method(parm.outputs->answers[i]);
		if( (options==1 && indice>2) || (options==2 && indice>2 && indice<5) || (options==3 && indice>4) )
			G_fatal_error(_("L option %i ne permet pas de calculer : %s en sortie"), options, menu_outputs[indice].text);
	}
	
	/* Vérifie que les cartes identifiant les surfaces en eau et la zone alluviale sont renseignées
	pour la prise en compte de la zone alluviale */
	if(flag6 && !parm.waterbodies->answer && !parm.riparian->answer)
		G_fatal_error(_("Pour la prise en compte de la zone alluviale (flag: '-z') veuillez remplir les champs waterbodies= et riparian="));
	
	/* Vérifie le montant de la mémoire spécifié */
	if (sscanf(parm.mem->answer, "%d", &maxmem) != 1 || maxmem <= 0)
		G_fatal_error(_("Quantite de memoire inappropriee: %d"), maxmem);

	/* Récupère les paramètres renseignés */
	method 			= find_method(parm.method->answer);
	method_ia		= find_ia_method(parm.init_abs->answer);
	if(method){
		algorithm	= find_algorithm_method(parm.algorithm->answer);
			if(algorithm==2||algorithm==4)
				mfd_converge = atof(parm.convergence->answer);
	}
	
	for (i = 0; parm.flow_speeds->answers[i]; i++)
		;
	if(method>1 && i<1)
		G_fatal_error(_("Le champ flow_speeds= doit etre renseigner pour le calcul du ruissellement"));
	if(method==4 && i<2)
		G_fatal_error(_("Le champ flow_speeds= doit comporter le nom de 2 cartes rasters de vitesse d ecoulement de l eau: a la surface et en subsurface"));
	if(method==4 && i==2)	
		G_verbose_message(_("Le module suppose que les vitesses d ecoulement de l eau sont dans l ordre suivant : vitesse d ecoulement surfacique, vitesse d ecoulement subsurfacique"));
	
	for (i = 0; parm.flow_disps->answers[i]; i++)
		;
	if(method>0 && i<1)
		G_fatal_error(_("Le champ flow_disps= doit etre renseigner pour le calcul du ruissellement"));
	if(method==3 && i<2)
		G_fatal_error(_("Le champ flow_disps= doit comporter le nom de 2 cartes rasters de coefficient de diffusion de l eau: a la surface et en subsurface"));
	if(method==3 && i==2)	
		G_verbose_message(_("Le module suppose que les coefficients de diffusion de l eau sont dans l ordre suivant : diffusion de l ecoulement surfacique, diffusion de l ecoulement subsurfacique"));
				
	for(i=0;parm.drainage_times->answers[i]; i++)
		;
	double drainage_times[2] ={ ( (flag3 && !flag4 && !flag5) || (flag3 && flag5) ) ? atof(parm.drainage_times->answers[0]) : atof(parm.drainage_times->answers[0]) * 30., ( (flag3 && !flag4 && !flag5) || (flag3 && flag5) ) ? atof(parm.drainage_times->answers[(i>1)?1:0]) : atof(parm.drainage_times->answers[(i>1)?1:0]) * 30. };

	// dimensions et résolution de la fenêtre d'etude
	RES 	= (double) window.ew_res;
	nrows 	= Rast_window_rows();
    ncols 	= Rast_window_cols();

	Rast_set_d_null_value(&null_val, 1);
		
    if ((double) nrows * ncols > 200000000) srows = scols = SEGCOLSIZE / 2;
	else   srows = scols = SEGCOLSIZE;
	
	/* Calcule le nombre total de segments */
	nseg = ((nrows + srows - 1) / srows) * ((ncols + scols - 1) / scols);
	
	/* Calcule l'espace disque et les besoins en mémoire */
	/* (nrows + ncols) * 8. * 20.0 / 1048576. for Dijkstra search */	
	pq_mb = ((double)nrows + ncols) * 8. * 20.0 / 1048576.;
    G_debug(1, "pq MB: %g", pq_mb);
    maxmem -= pq_mb;
    if (maxmem < 10) maxmem = 10;
    disk_mb = (double) nrows * ncols * 24. / 1048576.;
    segments_in_memory = maxmem / ((double) srows * scols * (24. / 1048576.));
    if (segments_in_memory < 4) segments_in_memory = 4;
    if (segments_in_memory > nseg) segments_in_memory = nseg;
    mem_mb = (double) srows * scols * (24. / 1048576.) * segments_in_memory;

	// options 1

	//options2
	//	G_message(_("Vous avez selectionne simultanément les options -j et -m;"
	//				"cela signifie que le module va calculer un bilan hydrique mensuel"
	//				"a partir de cartes rasters journalieres."));
	//options3		
	//	G_message(_("Vous avez selectionne simultanément les options %s et -f;"
	//			  "cela signifie que le module va calculer un bilan hydrique %s"
	//			  "tous les %i%s."),(flag3)?"-j":"-m",(flag3)?"journalier":"mensuel",outiter,(flag3)?"jours":"mois");
	
	if (flag->answer) {
	fprintf(stdout, _("Bilan hydrique :%s\n"),
	(options==1 && flag3)?"journalier":(options==1 && flag4)?"mensuel":(options==2)?"mensuel (avec donnees journalieres)":(options==3 && flag3)?"journalier":"mensuel");
	fprintf(stdout, "\n");
	if(options==3)
		fprintf(stdout, _(" Frequence des calculs :tous les %i%s\n"),outiter,(flag3)?"jours":"mois");	
    fprintf(stdout, _("Methode de calcul:%s -%s-"), menu[method].name, menu[method].text);
    fprintf(stdout, "\n");
		if(method){
			fprintf(stdout, _("Algorithme de redistribution du ruissellement:%s"), menu_algorithm[algorithm].name);
			fprintf(stdout, "\n");
			fprintf(stdout, _("%s"), menu_algorithm[algorithm].text);
			if(algorithm==2||algorithm==4)
				fprintf(stdout, _("Exposant de partitionnement du ruissellement:%.1f"), mfd_converge);
			if(method==1||method==3)
				fprintf(stdout, _("Temps de drainage du bassin versant par ruissellement de surface:%.1f h"), drainage_times[1]);			
			if(method>1)
				fprintf(stdout, _("Temps de drainage du bassin versant par ruissellement de subsurface:%.1f h"), drainage_times[2]);			
			if(method==1||method==3)
				fprintf(stdout, _("Technique de calcul de l abstraction initiale:%s -%s-"), menu_ia[method_ia].name,menu_ia[method_ia].text);
			fprintf(stdout, "\n");			
		}
    fprintf(stdout, _("Dimensions de la carte :\nNombre de lignes:%.1f\nNombre de colonnes:%.1f\nResolution des pixels:%.1fmx%.1fm"),nrows,ncols,RES,RES);
    fprintf(stdout, "\n");	
	fprintf(stdout, _("Vous aurez besoin d'au moins %.2f MB d'espace de disque"), disk_mb);
	fprintf(stdout, "\n");
    fprintf(stdout, _("Vous aurez besoin d'au moins %.2f MB de memoire"), mem_mb);
    fprintf(stdout, "\n");
    fprintf(stdout, _("%d des %d segments sont gardes en memoire"), segments_in_memory, nseg);
    fprintf(stdout, "\n");
    
    exit(EXIT_SUCCESS);
    }
		
	//if(flag3 && !flag4 && !flag5 && num_outputs>0 && num_outputs!=num_inputs || flag4 && !flag4 && !flag5 && num_outputs>0 && num_outputs!=num_inputs)
	//G_fatal_error(_("Les listes des rasters d entree prec=, etp= et de sortie outputs= doivent avoir la meme longueur lorsque les options -j ou -m sont specifiees seules."));

}

	/* ******************************************** */
	/* Crée et écrit un format de fichier segmenté  */
	/* ******************************************** */

	void createSEGMENT(){

    G_verbose_message(_("Cree un fichier temporaire..."));
	
    if (Segment_open(&parms_seg, G_tempfile(), nrows, ncols, srows, scols, sizeof(struct Parm), segments_in_memory) != 1)
		G_fatal_error(_("Ne peux pas creer le fichier temporaire"));


	/* DECLARE */
	int skip_nulls;
	int sat_fd, fc_fd, rum_fd, depth_fd;
	int sat_dsize, fc_dsize, rum_dsize, depth_dsize;
	double p_sat, p_fc, p_rum, p_depth;
	void *sat_cell, *fc_cell, *rum_cell, *depth_cell;

	// method>0	
	int alt_fd, pwp_fd, slope_fd;
	int alt_dsize,  pwp_dsize, slope_dsize;
	double  p_pwp, p_slope;
	float p_alt;
	void *alt_cell, *pwp_cell, *slope_cell;
	// method==1||method==3
	int speed_sf_fd, disp_sf_fd;
	int speed_sf_dsize, disp_sf_dsize;
	double p_speed_sf, p_disp_sf;
	void *speed_sf_cell, *disp_sf_cell;
	//method>1
	int speed_ssf_fd, disp_ssf_fd, ksat_fd;
	int speed_ssf_dsize, disp_ssf_dsize, ksat_dsize;
	double p_speed_ssf, p_disp_ssf, p_ksat;
	void *speed_ssf_cell, *disp_ssf_cell, *ksat_cell;	// pointeurs vers un bloc de mémoire				

	
	/* INITIALISE */
	
	// Fichiers descripteurs des couches raster : (1)
	sat_fd 		= openLayer(parm.sat->answer);
	fc_fd 		= openLayer(parm.fc->answer);
	rum_fd 		= openLayer(parm.rum->answer);
	depth_fd 	= openLayer(parm.depth->answer);	
	// Types de données : (2)
	sat_data_type 	= Rast_get_map_type(sat_fd);
	fc_data_type 	= Rast_get_map_type(fc_fd);
	rum_data_type 	= Rast_get_map_type(rum_fd);
	depth_data_type = Rast_get_map_type(depth_fd);	
	// Taille des données : (3)
	sat_dsize 		= Rast_cell_size(sat_data_type);
	fc_dsize 		= Rast_cell_size(fc_data_type);
	rum_dsize 		= Rast_cell_size(rum_data_type);
	depth_dsize 	= Rast_cell_size(depth_data_type);	
	// Alloue de la mémoire et initialise l'adresse des pointeurs : (4)	
	sat_cell 		= Rast_allocate_buf(sat_data_type);
	fc_cell 		= Rast_allocate_buf(fc_data_type);
	rum_cell 		= Rast_allocate_buf(rum_data_type);
	depth_cell 		= Rast_allocate_buf(depth_data_type);	
	
	p_sat = 0.0; p_fc = 0.0; p_rum = 0.0; p_depth = 0.0;
	
	
	if(method>0){
		
		// (1)
		alt_fd 			= openLayer(parm.altitude->answer);
		pwp_fd 			= openLayer(parm.pwp->answer);
		slope_fd 		= openLayer(parm.slope->answer);
		// (2)
		alt_data_type 	= Rast_get_map_type(alt_fd);
		pwp_data_type 	= Rast_get_map_type(pwp_fd);
		slope_data_type = Rast_get_map_type(slope_fd);
		// (3)
		alt_dsize 		= Rast_cell_size(alt_data_type);
		pwp_dsize 		= Rast_cell_size(pwp_data_type);
		slope_dsize 	= Rast_cell_size(slope_data_type);
		// (4)		
		alt_cell		= Rast_allocate_buf(alt_data_type);	
		pwp_cell 		= Rast_allocate_buf(pwp_data_type);
		slope_cell 		= Rast_allocate_buf(slope_data_type);

		p_alt =  0.0; p_pwp = 0.0; p_slope = 0.0;
	}	
		if(method==1||method==3){

		// (1)			
			speed_sf_fd = openLayer(parm.flow_speeds->answers[1]);
			disp_sf_fd 	= openLayer(parm.flow_disps->answers[1]);
		// (2)			
			speed_sf_data_type 	= Rast_get_map_type(speed_sf_fd);
			disp_sf_data_type 	= Rast_get_map_type(disp_sf_fd);
		// (3)			
			speed_sf_dsize 	= Rast_get_map_type(speed_sf_data_type);
			disp_sf_dsize 	= Rast_get_map_type(disp_sf_data_type);
		// (4)			
			speed_sf_cell 	= Rast_allocate_buf(speed_sf_data_type);
			disp_sf_cell 	= Rast_allocate_buf(disp_sf_data_type);
			
			p_speed_sf = 0.0; p_disp_sf = 0.0;
		}
			if(method>1){

		// (1)				
				speed_ssf_fd = openLayer(parm.flow_speeds->answers[2]);
				disp_ssf_fd = openLayer(parm.flow_disps->answers[2]);	
				ksat_fd 	= openLayer(parm.ksat->answer);
		// (2)
				speed_ssf_data_type = Rast_get_map_type(speed_ssf_fd);
				disp_ssf_data_type 	= Rast_get_map_type(disp_ssf_fd);	
				ksat_data_type 		= Rast_get_map_type(ksat_fd);
		// (3)				
				speed_ssf_dsize = Rast_cell_size(speed_ssf_data_type);
				disp_ssf_dsize 	= Rast_cell_size(disp_ssf_data_type);
				ksat_dsize		= Rast_cell_size(ksat_data_type);
		// (4)				
				speed_ssf_cell 	= Rast_allocate_buf(speed_ssf_data_type);
				disp_ssf_cell 	= Rast_allocate_buf(disp_ssf_data_type);	
				ksat_cell 		= Rast_allocate_buf(ksat_data_type);
				
				p_speed_ssf = 0.0; p_disp_ssf = 0.0; p_ksat = 0.0;
			}
	// Fixe les valeurs nulles
	Rast_set_d_null_value(&dnullval, 1);
	//Esquive les valeurs nulles ??
	skip_nulls = Rast_is_d_null_value(&null_val);

	total_cells = nrows * ncols;

	    for (row = 0; row < nrows; row++) {
		
	    G_percent(row, nrows, 2);
		
			Rast_get_row(sat_fd, sat_cell, row, sat_data_type);
		//if( Rast_get_row(sat_fd, sat_cell, row, sat_data_type) < 0 )
			//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.sat->answer, row);
			Rast_get_row(fc_fd, fc_cell, row, fc_data_type);
		//if( Rast_get_row(fc_fd, fc_cell, row, fc_data_type) < 0 )
			//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.fc->answer, row);
			Rast_get_row(rum_fd, rum_cell, row, rum_data_type);
		//if( Rast_get_row(rum_fd, rum_cell, row, rum_data_type) < 0 )
			//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.rum->answer, row);
			Rast_get_row(depth_fd, depth_cell, row, depth_data_type);
		//if( Rast_get_row(depth_fd, depth_cell, row, depth_data_type) < 0 )
			//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.depth->answer, row);
	    	
		if(method>0){
			Rast_get_row(alt_fd, alt_cell, row, alt_data_type);
			//if( Rast_get_row(alt_fd, alt_cell, row, alt_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.altitude->answer, row);
				Rast_get_row(pwp_fd, pwp_cell, row, pwp_data_type);
			//if( Rast_get_row(pwp_fd, pwp_cell, row, pwp_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.pwp->answer, row);
				Rast_get_row(slope_fd, slope_cell, row, slope_data_type);
			//if( Rast_get_row(slope_fd, slope_cell, row, slope_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.slope->answer, row);
		}
		
		if(method==1||method==3){
			Rast_get_row(speed_sf_fd, speed_sf_cell, row, speed_sf_data_type);
			//if( Rast_get_row(speed_sf_fd, speed_sf_cell, row, speed_sf_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.flow_speeds->answers[0], row);
			Rast_get_row(disp_sf_fd, disp_sf_cell, row, disp_sf_data_type);
			//if( Rast_get_row(disp_sf_fd, disp_sf_cell, row, disp_sf_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.flow_disps->answers[0], row);
		}
		if(method>1){
			Rast_get_row(speed_ssf_fd, speed_ssf_cell, row, speed_ssf_data_type);
			//if( Rast_get_row(speed_ssf_fd, speed_ssf_cell, row, speed_ssf_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.flow_speeds->answers[1], row);
				Rast_get_row(disp_ssf_fd, disp_ssf_cell, row, disp_ssf_data_type);
			//if( Rast_get_row(disp_ssf_fd, disp_ssf_cell, row, disp_ssf_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.flow_disps->answers[1], row);
				Rast_get_row(ksat_fd, ksat_cell, row, ksat_data_type);
			//if( Rast_get_row(ksat_fd, ksat_cell, row, ksat_data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.ksat->answer, row);
		}
		
	    // INPUT NULL VALUES: ???

		ptr2 = sat_cell;
		ptr3 = fc_cell;
		ptr4 = rum_cell;
	    ptr7 = depth_cell;

		if(method>0){
			ptr  = alt_cell;		
			ptr5 = pwp_cell;
			ptr6 = slope_cell;
		}
		if(method==1||method==3){
			ptr8 = speed_sf_cell;
			ptr9 = disp_sf_cell;
		}
		if(method>1){
			ptr10= speed_ssf_cell;
			ptr11= disp_ssf_cell;
			ptr12= ksat_cell;
		}
	        for (col = 0; col < ncols; col++){
  
	            if (Rast_is_null_value(ptr2, sat_data_type)){
	                p_sat = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr2, sat_data_type)) total_cells--;
	            }
	            else{
						switch (sat_data_type){
					  
						case CELL_TYPE:
	                        p_sat = *(CELL *)ptr2;
	                        break;
						case FCELL_TYPE:
	                        p_sat = *(FCELL *)ptr2;
	                        break;
						case DCELL_TYPE:
	                        p_sat = *(DCELL *)ptr2;
	                        break;
	                    }
	                }
	                parms.sat = p_sat;
					
				if (Rast_is_null_value(ptr3, fc_data_type)){
	                p_fc = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr3, fc_data_type)) total_cells--;
	            }
	            else{
				    	switch (fc_data_type){
						
						case CELL_TYPE:
	                        p_fc = *(CELL *)ptr3;
	                        break;
	                    case FCELL_TYPE:
	                        p_fc = *(FCELL *)ptr3;
	                        break;
	                    case DCELL_TYPE:
	                        p_fc = *(DCELL *)ptr3;
	                        break;
	                    }
	                }
	                parms.fc = p_fc;
					
				if (Rast_is_null_value(ptr4, rum_data_type)){
	                p_rum = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr4, rum_data_type)) total_cells--;
				}
				else{
				    	switch (rum_data_type){
						
						case CELL_TYPE:
	                        p_rum = *(CELL *)ptr4;
	                        break;
	                    case FCELL_TYPE:
	                        p_rum = *(FCELL *)ptr4;
	                        break;
	                    case DCELL_TYPE:
	                        p_rum = *(DCELL *)ptr4;
	                        break;
	                    }
	                }
	                parms.rum = p_rum;				

				if (Rast_is_null_value(ptr7, depth_data_type)){
	                p_depth = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr7, depth_data_type)) total_cells--;
				}
				else{
				    	switch (depth_data_type){
						
						case CELL_TYPE:
	                        p_depth = *(CELL *)ptr7;
	                        break;
	                    case FCELL_TYPE:
	                        p_depth = *(FCELL *)ptr7;
	                        break;
	                    case DCELL_TYPE:
	                        p_depth = *(DCELL *)ptr7;
	                        break;
	                    }
	                }
	                parms.depth = p_depth;					
	
	if(method>0){
					
	            if (Rast_is_null_value(ptr, alt_data_type)){
	                p_alt = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr, alt_data_type)) total_cells--;
	            }
	            else{
						switch (alt_data_type){
					  
						case CELL_TYPE:
	                        p_alt = *(CELL *)ptr;
	                        break;
						case FCELL_TYPE:
							p_alt = *(FCELL *)ptr;
							break;
						case DCELL_TYPE:
	                        p_alt = *(DCELL *)ptr;
	                        break;
	                    }
					}
	                parms.altitude = p_alt;					

				if (Rast_is_null_value(ptr5, pwp_data_type)){
	                p_pwp = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr5, pwp_data_type)) total_cells--;
				}
				else{
				    	switch (pwp_data_type){
						
						case CELL_TYPE:
	                        p_pwp = *(CELL *)ptr5;
	                        break;
	                    case FCELL_TYPE:
	                        p_pwp = *(FCELL *)ptr5;
	                        break;
	                    case DCELL_TYPE:
	                        p_pwp = *(DCELL *)ptr5;
	                        break;
	                    }
	                }
	                parms.pwp = p_pwp;

				if (Rast_is_null_value(ptr6, slope_data_type)){
	                p_slope = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr6, slope_data_type)) total_cells--;
				}
				else{
				    	switch (slope_data_type){
						
						case CELL_TYPE:
	                        p_slope = *(CELL *)ptr6;
	                        break;
	                    case FCELL_TYPE:
	                        p_slope = *(FCELL *)ptr6;
	                        break;
	                    case DCELL_TYPE:
	                        p_slope = *(DCELL *)ptr6;
	                        break;
	                    }
	                }
	                parms.tanslope = tan(p_slope);

	}
	
	if(method==1||method==3){
	
				if (Rast_is_null_value(ptr8, speed_sf_data_type)){
	                p_speed_sf = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr8, speed_sf_data_type)) total_cells--;
				}
				else{
				    	switch (speed_sf_data_type){
						
						case CELL_TYPE:
	                        p_speed_sf = *(CELL *)ptr8;
							p_speed_sf*= SECOND_TO_HOUR;
	                        break;
	                    case FCELL_TYPE:
	                        p_speed_sf = *(FCELL *)ptr8;
							p_speed_sf*= SECOND_TO_HOUR;							
	                        break;
	                    case DCELL_TYPE:
	                        p_speed_sf = *(DCELL *)ptr8;
							p_speed_sf*= SECOND_TO_HOUR;							
	                        break;
	                    }
	                }
	                parms.flow_speeds[1] = p_speed_sf;
					
				if (Rast_is_null_value(ptr9, disp_sf_data_type)){
	                p_disp_sf = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr9, disp_sf_data_type)) total_cells--;
				}
				else{
				    	switch (disp_sf_data_type){
						
						case CELL_TYPE:
	                        p_disp_sf = *(CELL *)ptr9;
							p_disp_sf*= SECOND_TO_HOUR;							
	                        break;
	                    case FCELL_TYPE:
	                        p_disp_sf = *(FCELL *)ptr9;
							p_disp_sf*= SECOND_TO_HOUR;							
	                        break;
	                    case DCELL_TYPE:
	                        p_disp_sf = *(DCELL *)ptr9;
							p_disp_sf*= SECOND_TO_HOUR;							
	                        break;
	                    }
	                }
	                parms.flow_disps[1] = p_disp_sf;
	}
	
	if(method>1){
	
				if (Rast_is_null_value(ptr10, speed_ssf_data_type)){
	                p_speed_ssf = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr10, speed_ssf_data_type)) total_cells--;
				}
				else{
				    	switch (speed_ssf_data_type){
						
						case CELL_TYPE:
	                        p_speed_ssf = *(CELL *)ptr10;
							p_speed_ssf*= DAY_TO_HOUR;
	                        break;
	                    case FCELL_TYPE:
	                        p_speed_ssf = *(FCELL *)ptr10;
							p_speed_ssf*= DAY_TO_HOUR;							
	                        break;
	                    case DCELL_TYPE:
	                        p_speed_ssf = *(DCELL *)ptr10;
							p_speed_ssf*= DAY_TO_HOUR;							
	                        break;
	                    }
	                }
	                parms.flow_speeds[2] = p_speed_ssf;
					
				if (Rast_is_null_value(ptr11, disp_ssf_data_type)){
	                p_disp_ssf = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr11, disp_ssf_data_type)) total_cells--;
				}
				else{
				    	switch (disp_ssf_data_type){
						
						case CELL_TYPE:
	                        p_disp_ssf = *(CELL *)ptr11;
							p_disp_ssf*= DAY_TO_HOUR;							
	                        break;
	                    case FCELL_TYPE:
	                        p_disp_ssf = *(FCELL *)ptr11;
							p_disp_ssf*= DAY_TO_HOUR;							
	                        break;
	                    case DCELL_TYPE:
	                        p_disp_ssf = *(DCELL *)ptr11;
							p_disp_ssf*= DAY_TO_HOUR;							
	                        break;
	                    }
	                }
	                parms.flow_disps[2] = p_disp_ssf;					
					
				if (Rast_is_null_value(ptr12, ksat_data_type)){
	                p_ksat = null_val;
	                if (skip_nulls && !Rast_is_null_value(ptr12, ksat_data_type)) total_cells--;
				}
				else{
				    	switch (ksat_data_type){
						
						case CELL_TYPE:
	                        p_ksat = *(CELL *)ptr12;
	                        break;
	                    case FCELL_TYPE:
	                        p_ksat = *(FCELL *)ptr12;
	                        break;
	                    case DCELL_TYPE:
	                        p_ksat = *(DCELL *)ptr12;
	                        break;
	                    }
	                }
	                parms.ksat = p_ksat;
	}
	
					//Convertit les teneurs en eau en hauteur d'eau mm (m3.m-3 -> mm).
					parms.sat*= parms.depth;
					parms.fc *= parms.depth;
					
					if(method>0)
						parms.pwp*= parms.depth;
					
					//parms.flow_speeds[1] = (5./3.)*p_speed_sf;					
					//parms.flow_disps[1] = (p_speed_sf * p_rh) / (2. * p_slope);	
	
					//parms.flow_speeds[2] = ( (sin(p_slope) * RAD_TO_DEG) * p_ksat ) / (p_sat-p_fc);					
					//parms.flow_disps[2] = ( (sin(p_slope) * RAD_TO_DEG) * p_ksat * p_depth * RES ) / ( (p_sat-p_fc) * (sin(p_slope) * RAD_TO_DEG) * RES );
					
					//Assigne les valeurs 
	                Segment_put(&parms_seg, &parms, row, col);
					
					//Incrémente les pointeurs
	                ptr2 = G_incr_void_ptr(ptr2, sat_dsize);
					ptr3 = G_incr_void_ptr(ptr3, fc_dsize);
					ptr4 = G_incr_void_ptr(ptr4, rum_dsize);
					ptr7 = G_incr_void_ptr(ptr7, depth_dsize);
					
					if(method>0){
						ptr  = G_incr_void_ptr(ptr,  alt_dsize);					
						ptr5 = G_incr_void_ptr(ptr5, pwp_dsize);
						ptr6 = G_incr_void_ptr(ptr6, slope_dsize);
					}
					if(method==1||method==3){
						ptr8 = G_incr_void_ptr(ptr8, speed_sf_dsize);
						ptr9 = G_incr_void_ptr(ptr9, disp_sf_dsize);
					}
					if(method>1){
						ptr10= G_incr_void_ptr(ptr10, speed_ssf_dsize);
						ptr11= G_incr_void_ptr(ptr11, disp_ssf_dsize);					
						ptr12= G_incr_void_ptr(ptr12, ksat_dsize);
					}
	            }
	        }
	
		G_free(sat_cell);
		G_free(fc_cell);
		G_free(rum_cell);
		G_free(depth_cell);
		
		if(method>0){	
			G_free(alt_cell);		
			G_free(pwp_cell);
			G_free(slope_cell);
		}
		if(method==1||method==3){
			G_free(speed_sf_cell);
			G_free(disp_sf_cell);
		}
		if(method>1){
			G_free(speed_ssf_cell);
			G_free(disp_ssf_cell);		
			G_free(ksat_cell);
		}
		G_percent(1, 1, 1);
	}

	/* *********************************************************** */
	/* Construit la liste des méthodes de calcul du bilan hydrique */
	/* *********************************************************** */
	
	static char *build_method_list(void){
		char *buf = G_malloc(1024);
		char *p = buf;
		int i;

		for (i = 0; menu[i].name; i++) {
			char *q;

			if (i)
				*p++ = ',';
			for (q = menu[i].name; *q; p++, q++)
				*p = *q;
		}
		*p = '\0';

		return buf;
	}
	
	/* ********************************************************** */
	/* Construit la liste des sorties du calcul du bilan hydrique */
	/* ********************************************************** */
	
	static char *build_outputs_list(void){
		char *buf = G_malloc(1024);
		char *p = buf;
		int i;

		for (i = 0; menu_outputs[i].name; i++) {
			char *q;

			if (i)
				*p++ = ',';
			for (q = menu_outputs[i].name; *q; p++, q++)
				*p = *q;
		}
		*p = '\0';

		return buf;
	}

	/* *********************************************************************** */
	/* Construit la liste des algorithme de calcul de l'aire de drainage amont */
	/* *********************************************************************** */
	
	static char *build_algorithm_list(void){
		char *buf = G_malloc(1024);
		char *p = buf;
		int i;

		for (i = 0; menu_algorithm[i].name; i++) {
			char *q;

			if (i)
				*p++ = ',';
			for (q = menu_algorithm[i].name; *q; p++, q++)
				*p = *q;
		}
		*p = '\0';

		return buf;
	}
	
	/* *********************************** */
	/* Donne un nom à la couche de sortie  */
	/* *********************************** */
	
	char *make_output_name(const char *output_name){
	static int indice;
	char *buf;
		indice = find_output_name(output_name);
		if(asprintf(&buf, "%s%d%s",menu_outputs[indice].name, n+1, menu[method].suffix)<0)
			fprintf(stderr, "Allocation de memoire lors de la creation du nom de raster de sortie a echoue\n");
		return buf;
	}
	
	/* *************************************************************************** */
	/* Ouvre une couche raster pour la lecture et retourne le fichier descripteur. */
	/* *************************************************************************** */
	
	int openLayer(char *layer){
		int fd;
		const char *mapset;		//Nom du jeu de carte ou se trouve le raster
		
		mapset = G_find_raster(layer, "");	
		if (mapset == NULL){
			G_fatal_error(_("Carte raster <%s> non trouvee"), layer);
		}
		
		fd = Rast_open_old(layer, mapset);
		if (fd < 0){
			G_fatal_error(_("Impossible d'ouvrir la carte raster <%s>"), layer);
		}
		
		return fd;
	}

	/* ********************************************** */
	/* Détecte la méthode de calcul du bilan hydrique */
	/* ********************************************** */
	
	static int find_method(const char *method_name){
		int indice;

			for (indice = 0; menu[indice].name; indice++)
				if (strcmp(menu[indice].name, method_name) == 0)
					return indice;
		
			G_fatal_error(_("Methode <%s> inconnue"), method_name);
		
			return -1;
		}
	
	/* ********************************************** */
	/* Détecte le nom de la couche de sortie désirée  */
	/* ********************************************** */
	
	static int find_output_name(const char *output_name){
		int indice;

			for (indice = 0; menu_outputs[indice].name; indice++)
				if (strcmp(menu_outputs[indice].name,output_name) == 0)
					return indice;
		
			G_fatal_error(_("Nom <%s> inconnu"), output_name);
		
			return -1;
		}

	/* ******************************************************* */
	/* Détecte la méthode de calcul de l'abstraction initiale  */
	/* ******************************************************* */
	
	static int find_ia_method(const char *method_name){
		int indice;

			for (indice = 0; menu_ia[indice].name; indice++)
				if (strcmp(menu_ia[indice].name, method_name) == 0)
					return indice;
		
			G_fatal_error(_("Methode <%s> inconnue"), method_name);
		
			return -1;
		}	

	/* ******************************************************************************************** */
	/* Détecte l'algorithme qui calcule la zone amont contribuant au ruissellement dans une cellule */
	/* ******************************************************************************************** */
	
	static int find_algorithm_method(const char *algorithm_name){
	int indice;
		for (indice = 0; menu_algorithm[indice].name; indice++)
			if (strcmp(menu_algorithm[indice].name, algorithm_name) == 0)
				return indice;
		G_fatal_error(_("Algorithme <%s> inconnu"),algorithm_name);
	
		return -1;
	}	

	/* ********************************************************************************* */
	/* Vérifie que les coordonnées (row,col) de la cellule se situe dans la zone d'étude */
	/* ********************************************************************************* */
	
	int is_OnGrid(int rown, int coln){
		return (rown >= 0 && rown < nrows && coln >=0 && coln < ncols);
	}

	/* ********************************************* */
	/* Calcule l'exposition de la cellule en radians */ 
	/* ********************************************* */
	
	double aspect_on_fly(int row, int col){

	double z[8];
	double dzdx, dzdy, aspect;

	/* calcul de l'exposition de la cellule */
	for(k=0; k<8; k++)
   {
		rown   = row + dy[k];
		coln   = col + dx[k];
		if( !is_OnGrid(rown, coln) )
           return UNDEF;
		   
		Segment_get(&parms_seg, &parms, rown, coln);
		z[k] = parms.altitude;
	}
       dzdx	= ((z[1] + z[0]+z[0] + z[7])-(z[3] + z[4]+z[4] + z[5])) / 8.*RES;
       dzdy 	= ((z[3] + z[2]+z[2] + z[1])-(z[5] + z[6]+z[6] + z[7])) / 8.*RES;

		if(dzdx == 0. && dzdy == 0.)
		   aspect = UNDEF;
		if(dzdx >0. && dzdy >=0.)
		   aspect  = atan(dzdy / dzdx) * M_RAD_TO_DEG + 180.;
		if(dzdx < 0.)
		   aspect = atan(dzdy / dzdx) * M_RAD_TO_DEG + 360.;
		if(dzdx == 0.)
		   aspect = (dzdy>0.) ? 270. : 90.;
		if(dzdx >0. && dzdy <0.)
		   aspect = atan(dzdy / dzdx) * M_RAD_TO_DEG;

	return aspect;
}

	/* *********************************************** */
	/* Calcule la distance en fonction de la direction */
	/* *********************************************** */
	
	double DIST(short dir){
	return (dir%2)==0 ? RES:(RES*M_SQRT2);
	}

	/* ********************************************** */
	/* Methode de calcul de la zone de drainage amont */
	/* ********************************************** */
	
	/*"Deterministic 8:\n"
	"- O'Callaghan, J.F. / Mark, D.M. (1984):\n"
	"    'The extraction of drainage networks from digital elevation data',\n"
	"    Computer Vision, Graphics and Image Processing, 28:323-344\n\n"*/
				  
	void D_8(int row, int col){
	
	Segment_get(&parms_seg, &parms, row, col);

	double z 		= parms.altitude, dz;
	double dzMax	= 0.0;
	short direction = -1;

		/* définit le gradient de pente maximum */
		for(k=0; k<8; k++){
			rown   = row + dy[k];
			coln   = col + dx[k];			
			Segment_get(&parms_seg, &parms, rown, coln);
			if( is_OnGrid(rown, coln) && (dz = ((z - parms.altitude)/DIST(k)) > 0.0) ){
				if( dz > dzMax){
					dzMax 		= dz;
					direction 	= k;
				}
			}
		}
	if(direction>=0){
		landscape[row+dy[direction]][col+dx[direction]].portion[(direction+4)%8]  = 1.;
		landscape[row+dy[direction]][col+dx[direction]].neighbors[(direction+4)%8]= &landscape[row][col];		
	}
	return;
	}

	/*"Deterministic Infinity:\n"
	"- Tarboton, D.G. (1997):\n"
	"    'A new method for the determination of flow directions and upslope areas in grid digital elevation models',\n"
	"    Water Ressources Research, Vol.33, No.2, p.309-319\n\n"*/

	void D_Inf(int row, int col){

	static int k, kk, rowk, colk, rowkk, colkk;
	double portion, Aspect;
	
		Aspect = aspect_on_fly(row, col);
		
		if(Aspect>0.){

			/* calcule la  portion d'aire drainée */
			k          = (int)(Aspect / 45.);
			kk         = k+1;
			portion    = fmod(Aspect, 45.) / 45.;

			/* attribue la portion d'aire drainée à k et k+1 */
			rowk   = row + dy[k];
			colk   = col + dx[k];

			rowkk  = row + dy[kk];
			colkk  = col + dx[kk];

			landscape[rowk][colk].portion[(k+4)%8] 	= 1. - portion;
			landscape[rowkk][colkk].portion[(kk+4)%8] 	= portion;
			
			landscape[rowk][colk].neighbors[(k+4)%8] 	= &landscape[row][col];
			landscape[rowkk][colkk].neighbors[(kk+4)%8]	= &landscape[row][col];
		}
		else
			D_8(row, col);
	return;
	}
 
	/*"Multiple Flow Direction:\n"
	"- Freeman, G.T. (1991):\n"
	"    'Calculating catchment area with divergent flow based on a regular grid',\n"
	"    Computers and Geosciences, 17:413-22\n\n"*/

	void MFD_8(int row, int col){
	
	Segment_get(&parms_seg, &parms, row, col);

	double tanBeta[8];
	double dz, z = parms.altitude;
	double dzSum = 0.0;
	
	/* définit la longueur de contour effective */
	const double l[8] = {0.5,0.354,0.5,0.354,0.5,0.354,0.5,0.354};
			
		/* Calcule et assigne la proportion de flux drainant vers chaque cellule voisine située à l'aval */	
			for(k=0; k<8; k++)
		{
			tanBeta[k] = 0.0;
			rown  = row + dy[k];
			coln  = col + dx[k];

			Segment_get(&parms_seg, &parms, rown, coln);
				if( is_OnGrid(rown, coln) && (dz = (z - parms.altitude)/DIST(k)) > 0.0 )
			{
			    tanBeta[k] = pow(dz, mfd_converge) * l[k];
				dzSum     += tanBeta[k];
			}
			
		}
			if( dzSum > 0.0 )
		{
				for(k=0; k<8; k++)
			{
				if(tanBeta[k]){
					landscape[row+dy[k]][col+dx[k]].portion[(k+4)%8]   = tanBeta[k] / dzSum;
					landscape[row+dy[k]][col+dx[k]].neighbors[(k+4)%8] = &landscape[row][col];
				}
			}
		}
	return;
	}

	/*"Multiple Flow Direction based on maximum downslope gradient:\n"
	"- Qin C, Zhu AX, Pei T, Li B, Zhou C, Yang L (2007):\n"
	" 'An adaptive approach to selecting a flow-partition exponent for a multiple-flow-direction algorithm', \n"
	"International Journal of Geographical Information Science, 21:443-458\n\n"*/
	
	void MFD_md(int row, int col){
	
	Segment_get(&parms_seg, &parms, row, col);

	double tanBeta[8];
	double z = parms.altitude, dz;
	double dzMax, dzSum;
    dzMax = dzSum = 0.;
	
	/* définit la longueur de contour effective */
	const double l[8] = {0.5,0.354,0.5,0.354,0.5,0.354,0.5,0.354};
	
		/* définit le gradient de pente maximum */
			for(k=0; k<8; k++)
	    {
			rown   = row + dy[k];
			coln   = col + dx[k];
			
			Segment_get(&parms_seg, &parms, rown, coln);
				if( is_OnGrid(rown, coln) && (dz = (z - parms.altitude)/DIST(k)) > 0.0 )
			{
				if(dz > dzMax)
					dzMax = dz; 
			}
		}
		
		dzMax	= (dzMax < 1.) ? 8.9 * dzMax + 1.1 : 10.;
		/* Calcule et assigne la proportion de flux drainant vers chaque cellule voisine située à l'aval */	
			for(k=0; k<8; k++)
		{
			tanBeta[k] = 0.0;
			rown  = row + dy[k];
			coln  = col + dx[k];

			Segment_get(&parms_seg, &parms, rown, coln);
				if( is_OnGrid(rown, coln) && (dz = (z - parms.altitude)/DIST(k)) > 0.0 )
			{
			    tanBeta[k] = pow(dz, dzMax) * l[k];
				dzSum     += tanBeta[k];
			}
		}
			if( dzSum > 0.0 )
		{
				for(k=0; k<8; k++)
			{
				if(tanBeta[k]){
					landscape[row+dy[k]][col+dx[k]].portion[(k+4)%8]  = tanBeta[k] / dzSum;
					landscape[row+dy[k]][col+dx[k]].neighbors[(k+4)%8]= &landscape[row][col];
				}
			}
		}
	return;
	}

	/*"Triangular Multiple Flow Direction:\n"
	"- Seibert, J. / McGlynn, B. (2007):\n"
	'A new triangular multiple flow direction algorithm for computing upslope areas from gridded digital elevation models',
	"Water Resources Research, Vol. 43, W04501"*/
	
	void MFD_Inf(int row, int col){

	Segment_get(&parms_seg, &parms, row, col);

	double 	dzSum, dz1, dz2, e[8], e0 = parms.altitude;
	double 	d, s, s_facet[8], d_facet[8];
	double 	cellarea = RES * RES;
	double 	valley[8], portion[8];
	double 	nx, ny, nz, n_norm;
	int  	kk, Dir_inGrid[8];

	/* On commence par calculer la différence d'altitude
	entre la cellule centrale et chacune de ses cellules voisines */

	 for(k=0; k<8; k++)
	{
		rown          = row + dy[k];
		coln          = col + dx[k];
		s_facet[k]    = d_facet[k] = -999.0; // initialise la pente et la direction des facettes à -999.
		if( (Dir_inGrid[k] = is_OnGrid(rown, coln)) ){
			Segment_get(&parms_seg, &parms, rown, coln);
            e[k]  = parms.altitude;
		}
	}
      for(k=0; k<8; k++)
    {
        d = s = -999.0;

		if( Dir_inGrid[k])
		{
		    kk = ( k < 7)? k+1 : 0;

            if(Dir_inGrid[kk] ){

                if(e0>e[k] || e0>e[kk]){

                    nx = ( (e[k]-e0) * dy[kk] - (e[kk]-e0) * dy[k]) * RES;      // nx = z1y2 - z2y1
                    ny = ( (e[k]-e0) * dx[kk] - (e[kk]-e0) * dx[k]) * RES;      // ny = z1x2 - z2x1
                    nz = ( dy[k] * dx[kk] - dy[kk] * dx[k]) * cellarea; 		// nz = y1x2 - y2x1

                    n_norm = sqrt( nx*nx + ny*ny + nz*nz );

                    if( nx == 0.0 )
                    {
                      d = (ny >= 0.0)? 0.0 : M_PI;
                    }
                    else if( nx < 0.0 )
                    {
                      d = M_PI_270 - atan(ny / nx);
                    }
                    else
                    {
                      d = M_PI_090 - atan(ny / nx);
                    }
                    s = tan( acos( nz/n_norm ) );

                    if( d < k * M_PI_045 || d > (k+1) * M_PI_045 ){

                        if( (dz1 = (e0 - e[k])/ DIST(k)) > (dz2 = (e0 - e[kk])/ DIST(kk)) ){
                            d = k * M_PI_045;
                            s = dz1 / DIST(k);
                        }
                        else{
                            d = (k+1) * M_PI_045;
                            s = dz2 / DIST(kk);
                        }
                    }
                }
            }
        s_facet[k] = s;
        d_facet[k] = d;
		}
    }

	/* On continue en calculant la pente et la direction
	 des 8 facettes triangulaires centrées sur la cellule centrale
	 et deux de ses cellules voisines formant un angle de 45° */

    /* On  poursuit en calculant la fraction d'aire drainée par chaque facette (pondération par les pentes) */
	dzSum  = 0.0;
      for(k=0; k<8; k++)
    {
		valley[k]   = 0.0;
		kk = (k < 7)? k+1 : 0;

		if( s_facet[k] > 0.0 )
		{
		    /* Si la direction calculée pointe entre deux cellules voisines */
			if( d_facet[k] >= k * M_PI_045 && d_facet[k] <= (k+1) * M_PI_045 ){
				valley[k] = s_facet[k];// donne la mm pente à la facette
            /* Sinon si la direction calculée est identique entre deux cellules adjacentes */
			}else if( d_facet[k] == d_facet[kk] ){
				valley[k] = s_facet[k]; // donne la mm pente à la facette
            /* Sinon si la pente de la facette suivante n'a pas pu être calculée et que la direction
            de la facette courante pointe vers le voisin suivant */
			}else if( s_facet[kk] == -999.0 && d_facet[k] == (k+1) * M_PI_045 ){
				valley[k] = s_facet[k]; // donne la mm pente à la facette
			}else{
				kk = (k > 0)? k-1 : 7;
			/* Sinon si la pente de la facette précédente n'a pas pu être calculée et que la direction
            de la facette courante pointe vers le voisin courant */
				if( s_facet[kk] == -999.0 && d_facet[k] == k * M_PI_045 ){
					valley[k] = s_facet[k]; // donne la mm pente à la facette
				}
			}
			/* Comme dans l'algorithme MFD8 de Quinn et al. 1991 il est possible de pondérer l'aire de drainage
			et donner plus de poids aux pente les + raides grace à un exposant MFD_Converge */
        valley[k] = pow(valley[k], mfd_converge);
        dzSum    += valley[k];
		}
    /* Si la pente est négative alors la fraction d'aire drainée est nulle */
    portion[k] = 0.0;
    }

    /* On finit par calculer la portion d'aire drainée vers chaque voisin (pondération par les directions des facettes */
      if( dzSum )
    {
		for(k=0; k<8; k++)
		{
			if (k < 7){
				kk = k+1;
			}else{
				kk = 0;
				if( d_facet[k] == 0.0)
					d_facet[k] = M_PI_360;
			}
			/* si la pente de la facette k n'est pas nulle */
			if( valley[k] ){
                // calcule la fraction d'aire drainée par la facette
				valley[k] /= dzSum;
				/* calcule la portion d'aire drainée vers chaque voisin.
				   calcule la différence relative entre la direction de pente de la facette
				   et les directions délimitant cette facette.
				   Attribue une portion de la fraction drainée aux deux cellules constituant la facette */
				portion[k] += valley[k] * ((k+1) * M_PI_045 - d_facet[k]) / M_PI_045;
				portion[kk]+= valley[k] * (d_facet[k] - k * M_PI_045) / M_PI_045;
			}
		}
		for(k=0; k<8; k++)
		{
			landscape[row+dy[k]][col+dx[k]].portion[(k+4)%8]	= portion[k];
			landscape[row+dy[k]][col+dx[k]].neighbors[(k+4)%8]	= &landscape[row][col];
		}
    }
	return;
}

	/* ****************************************************************************** */	
	/* Alloue une nouvelle couche de sol avec des valeurs d'initialisation par défaut */
	/* ****************************************************************************** */
	
	layer *NewLayer(){
		layer *newlayer = (layer *)G_malloc(ncols*sizeof(layer));
		if(!newlayer){	
			fprintf(stderr, "Allocation de memoire pour la couche de sol a echoue\n");
			exit(1);
		}		
		for(col=0;col<ncols;col++){
			/* Alloue un bloc de mémoire pour stocker les valeurs de teneur en eau */
				newlayer[col].swc	= (double *)G_malloc(num_inputs*sizeof(double));
			/* Initialise le pointeur vers l'eau disponible au drainage	et les cellules
			amont contribuant au ruissellement dans la cellule	à leur valeur par défaut
			(i.e. NULL) */
				newlayer[col].raw				= NULL;
				newlayer[col].braw				= NULL;				
				newlayer[col].contribCells		= NULL;
				newlayer[col].UHTsf				= NULL;
				newlayer[col].UHTssf			= NULL;				
			/* Initialise à zéro le nombre de cellules drainant vers la cellule */			
				for(k=0;k<2;k++)				
					newlayer[col].nbContribCells[k]= 0;
			/* Initialise chaque pointeur de cellules et la portion d'aire qu'elles drainent 
			à leur valeur par défaut (i.e. NULL et zéro respectivement) */
				for(k=0;k<8;k++){
					newlayer[col].neighbors[k] 	= NULL;
					newlayer[col].portion[k]	= 0.0;
				}
			/* Initialise les paramètres servant au calcul du ruissellement de surface
			à zéro par défaut */
				newlayer[col].smax 				= 0.0;
				newlayer[col].w1 				= 0.0;
				newlayer[col].w2 				= 0.0;
			/* Initialise la catégorie à laquelle appartient la cellule:
			"surface en eau" et "zone humides" */
				newlayer[col].waterbodies 		= 0;
				newlayer[col].riparian 			= 0;
		}
		
		return newlayer;
	
	}	

	/* ************************************************* */	
	/* Libère la mémoire utilisée par les couches de sol */
	/* ************************************************* */

	void FreeLandscape(){
	
	layer *ptr;	
	
		for(row=0;row<nrows;row++){
		
			ptr = landscape[row];
		 
				for(col=0;col<ncols;col++)
			{
				if(ptr[col].swc)
					G_free(ptr[col].swc);
				if(ptr[col].raw)
					G_free(ptr[col].raw);
				if(ptr[col].braw)
					G_free(ptr[col].braw);
				if(ptr[col].contribCells)
					G_free(ptr[col].contribCells);					
				if(ptr[col].UHTsf)
					free_dvector(ptr[col].UHTsf, 1, num_inputs);			
				if(ptr[col].UHTssf)
					free_dvector(ptr[col].UHTssf, 1, num_inputs);		
			}
			
		}
		G_free(landscape);
		
		return;
		
	}
	
	/* ***************************************************** */
	/* Alloue de la mémoire pour les couches raster d'entrée */
	/* ***************************************************** */
	
	void AllocateMemory(){

		G_message(_("Alloue de la memoire pour les calculs..."));
			
		/* allocation des variables climatiques */
		P 	= G_malloc(num_inputs * sizeof(struct input));
			if(!P)
				G_fatal_error(_("Impossible d allouer de la memoire pour les rasters d entree p="));
		ETP = G_malloc(num_inputs * sizeof(struct input));
			if(!ETP)
				G_fatal_error(_("Impossible d allouer de la memoire pour les rasters d entree etp="));

		/* allocation des rasters de sortie */
		Outputs = G_malloc(num_outputs * sizeof(struct output));
			if(Outputs==NULL)
				G_fatal_error(_("Impossible d allouer de la memoire pour les rasters de sortie"));
		
		/* allocation de la carte d'étude */		
		landscape = (layer **) G_malloc(nrows * sizeof(layer *));
		
			if(flag6){
				waterbodies = (CELL **) G_malloc(nrows * sizeof(CELL *));
				riparian   	= (CELL **) G_malloc(nrows * sizeof(CELL *));
			}
			if(method==1){
				smax 	= (DCELL **) G_malloc(nrows * sizeof(DCELL *));
				w1   	= (DCELL **) G_malloc(nrows * sizeof(DCELL *));
				w2   	= (DCELL **) G_malloc(nrows * sizeof(DCELL *));
			}

		for	(row = 0; row < nrows; row++){
			G_percent(0, nrows, 5);
			
			landscape[row] = NewLayer();
			
				if(flag6){
				waterbodies[row]= (CELL *) G_malloc(ncols * sizeof(CELL));
				riparian[row]   = (CELL *) G_malloc(ncols * sizeof(CELL));		
				}
				if(method==1){
				smax[row] 		= (DCELL *) G_malloc(ncols * sizeof(DCELL));
				w1[row]   		= (DCELL *) G_malloc(ncols * sizeof(DCELL));
				w2[row] 		= (DCELL *) G_malloc(ncols * sizeof(DCELL));
				}
			G_percent(row, nrows - 1, 5);
		}
		
		return;
		
	}

	/* *********************************************** */
	/* Lit les couches raster d'entrée dans la memoire */
	/* *********************************************** */
	
	void ReadInputLayer(){
	
	/* DECLARE */
	int fd;
	int data_size;
	void *cell;
	void *pstructMembers, *ptr2;
	RASTER_MAP_TYPE data_type;
	//struct Options *structMembers[5]   = {parm.waterbodies->answer)?parm.waterbodies:NULL,(parm.riparian->answer)?parm.riparian:NULL,(parm.smax->answer)?parm.smax:NULL,(parm.w->answers)?parm.w[0]:NULL,(parm.w->answers)?parm.w[1]:NULL};
	char *structMembers[5]   			= {parm.waterbodies->answer,parm.riparian->answer,parm.smax->answer,parm.w->answers[0],parm.w->answers[1]};
	void *arraysofPointersToLayers[5] 	= {waterbodies,riparian,smax,w1,w2};
	
	/* INITIALISE */
		for(i=0;i<NUM_ELEM(structMembers);i++){
		
		if(!structMembers[i])
			continue;
		
		// Fichiers descripteurs des couches raster	
		fd 			= openLayer(structMembers[i]);
		// Types de données
		data_type 	= Rast_get_map_type(fd);
		//Taille des données	
		data_size 	= Rast_cell_size(data_type);
		// Alloue de la mémoire et initialise l'adresse des pointeurs
		cell		= Rast_allocate_buf(data_type);

		G_verbose_message(_("Lecture de la carte raster <%s>..."), structMembers[i]);

			for (row = 0; row < nrows; row++) {					
			G_percent(row, nrows, 2);
			
				Rast_get_row(fd, cell, row, data_type);
				
			//if( Rast_get_row(fd, cell, row, data_type) < 0 )
				//G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), structMembers[i]->answer, row);
				ptr2 				= arraysofPointersToLayers[i];
				pstructMembers  	= cell;
				
				for (col = 0; col < ncols; col++){
					// Inscrit les données en mémoire
						switch(data_type){
					
						case CELL_TYPE:					
						*((CELL *)ptr2) = *((CELL *)pstructMembers);
						break;
						
						case FCELL_TYPE:
						*((FCELL *)ptr2) = *((FCELL *)pstructMembers);
						break;
						
						case DCELL_TYPE:
						*((DCELL *)ptr2) = *((DCELL *)pstructMembers);
						break;					
					}				
					// Incrémente les pointeurs
					pstructMembers  = G_incr_void_ptr(pstructMembers,  data_size);
					ptr2  			= G_incr_void_ptr(ptr2,  data_size);
				}		
			}
		G_free(cell);
		Rast_close(fd);		
		G_percent(1, 1, 1);				
		}
	return;
	}

	/* ************************ */
	/* Vide/Nettoye la memoire  */
	/* ************************ */
	
	void Cleanup(){
		for(row = 0; row < nrows; row++)
		{
			if(flag6){
				G_free(waterbodies[row]);
				G_free(riparian[row]);
			}
			if(method==1||method==3){
				G_free(smax[row]);
				G_free(w1[row]);
				G_free(w2[row]);
			}
		}
		if(flag6){
			G_free(waterbodies);
			G_free(riparian);
		}
		if(method==1||method==3){
			G_free(smax);
			G_free(w1);
			G_free(w2);
		}
		return;
	}

	/* ************************************************ */
	/* Retourne un vecteur d'indice des termes non nuls */
	/* ************************************************ */
	
	int *FindNonZeroTermIndices(double *p, int size){
	int *ptr=NULL;
	static int count = 0;	
		for(i=0;i<size;i++){	
			if(p[i]>0.0){
			count++;
			ptr = (int *)G_realloc(ptr, count*sizeof(int));
			ptr[count-1] = i;
			}
		}
	return ptr;
	}

	/* ******************************************************** */
	/* Fonction de réponse à l'échelle d'un chemin d'écoulement */ 
	/* ******************************************************** */

	double FlowPathUnitResponse(node *p, int time_index, int id){
	
		double avg_travel_time 		= p->avg_travel_time[id];
		double var_of_flow_time 	= p->var_of_flow_time[id];
		double sqrt_sigma 			= sqrt(var_of_flow_time);
		double U_t;

		U_t = ( 1.0 / (sqrt_sigma * sqrt( 2.0 * M_PI * pow(time_index,3.0)/pow(avg_travel_time,3.0) ) ) ) * exp( - pow(time_index - avg_travel_time, 2.0)/(2.0 * var_of_flow_time * time_index/avg_travel_time) );
		return U_t;
	}
	
	/* ********************************************* */
	/* Fonction de réponse à l'échelle d'une cellule */ 
	/* ********************************************* */

	double CellOutletResponse(node *p, int time_index, int id){

		double c 	= p->avg_travel_time[id];
		double d 	= p->var_of_flow_time[id];
		double l 	= RES;
		double u_t;

		u_t = ( l / (2.0 * sqrt( M_PI * d * pow(time_index,3.0) ) ) ) * exp( - pow(c * time_index - l, 2.0)/(4.0 * d * time_index) );
		return u_t;
	}	
	
	/* ********************************************** */
	/* Vérifie si un objet est dans la file d'attente */
	/* ********************************************** */
	
	int IsEnqueued(Queue *queue, int rown, int coln){
	
		static int i;
		int check = 0;
		int size = queue->size;
		node *CursorNode;

		for(i=0; i<size; i++)
		{
			CursorNode = Front(queue);
			if(CursorNode->row==rown && CursorNode->col==coln)
			check++;
			DeQueue(queue);
			EnQueue(queue,CursorNode);
		}
		return check;
	}

	/* ******************************************************************** */	
	/* Alloue un nouveau noeud avec des valeurs d'initialisation par défaut */
	/* ******************************************************************** */

	node *NewNode(){
		node *newnode = (node *)G_malloc(sizeof(node));

		if( newnode==NULL ){	
			fprintf(stderr, "Allocation de memoire pour le noeud a echoue\n");
			exit(1);
		}
		/* Initialise les données de ruissellement à zéro */
		for(i=0;i<2;i++){
			newnode->is_visited[i] = 0;
			newnode->travel_time[i] = 0.;
			newnode->avg_travel_time[i] = 0.;
			newnode->var_of_flow_time[i] = 0.;
			newnode->portion[i] = 0.;
		}
		/* Initialise les pointeurs vers les couches adjacentes */
		for(k=0;k<8;k++)
			newnode->neighbors[k] = NULL;
		
		return newnode;
	}

	/* **************************************** */
	/* Récupère un objet dans la file d'attente */
	/* **************************************** */
	
	node *GetNode(Queue *queue, int rown, int coln){

		static int i;
		int size = queue->size;
		node *CursorNode, *TargetNode;

		for(i=0; i<size; i++)
		{
			CursorNode = Front(queue);
			if(CursorNode->row==rown && CursorNode->col==coln)
				TargetNode = CursorNode;

			DeQueue(queue);
			EnQueue(queue,CursorNode);
		}
    return TargetNode;
	}

	/* **************************************************************************** */
	/* Identifie les cellules appartenant au bassin de drainage d'une cellule donné */
	/* **************************************************************************** */
	
	void FindBasin(layer *p){
	
	static int t, kept_processes;
	static double Travel_Time[2];
	static double UpslopeArea[2];

	/* Crée une file d'attente qui va contenir temporairement des pointeurs vers les cellules du réseau */
    Queue *queue = CreateQueue();
	
	/* Crée le premier noeud */
	node *root 	= NewNode();
	root->row 	= row;
	root->col 	= col;

	/* On met le noeud dans la file d'attente */
	EnQueue(queue, root);
		
	/* Pointeur vers la cellule en cours d'analyse. */
	node *CurrentNode = NULL;

		/* Algorithme BFS (Breadth-First-Search) */
		while (!QueueIsEmpty(queue))
		{	
			/* Récupère le premier élément de la file d'attente */
			CurrentNode = Front(queue); 
			DeQueue(queue);
	
			/* Récupère les données sur la cellule */
			Segment_get(&parms_seg, &parms, CurrentNode->row, CurrentNode->col);

			for(k=0;k<8;k++)
			{
				if( landscape[CurrentNode->row][CurrentNode->col].neighbors[k]==NULL)
				continue;
				
				/* Calcule les coordonnées (row, col) de la cellule voisine k */
				rown = CurrentNode->row + dy[k];
				coln = CurrentNode->col + dx[k];

				if(method==3){
				
					Travel_Time[id] 	= (CurrentNode->travel_time[id] + (1.0 /parms.flow_speeds[id])) * DIST(k);
					Travel_Time[id+1] 	= (CurrentNode->travel_time[id+1] + (1.0 /parms.flow_speeds[id+1])) * DIST(k);
					
					if( (Travel_Time[id] < drainage_times[id]) && (Travel_Time[id+1]> drainage_times[id+1]) ){
						kept_processes = 1;				
					}else if( (Travel_Time[id] > drainage_times[id]) && (Travel_Time[id+1] < drainage_times[id+1]) ){
						kept_processes = 2;				
					}else if( (Travel_Time[id] > drainage_times[id]) && (Travel_Time[id+1] > drainage_times[id+1]) ){
						continue;
					}else{					
						kept_processes = 3;				
					}
					
				}else if(method==1){
					Travel_Time[id] 	= (CurrentNode->travel_time[id] + (1.0 /parms.flow_speeds[id])) * DIST(k);				
					kept_processes = 1;				
				}else{
					Travel_Time[id+1] 	= (CurrentNode->travel_time[id+1] + (1.0 /parms.flow_speeds[id+1])) * DIST(k);				
					kept_processes = 2;				
				}
		
				switch(kept_processes){
				
				case 1:
					if(!IsEnqueued(queue, rown, coln)){
						/* Crée un nouveau noeud */
						CurrentNode->neighbors[k] 			= NewNode();
						/* Attribut le numéro de ligne et de colonne */
						CurrentNode->neighbors[k]->row 		= rown;
						CurrentNode->neighbors[k]->col 		= coln;
						/* Marque le noeud comme visité */
						CurrentNode->neighbors[k]->is_visited[id] = 1;						
						/* Calcul le temps de trajet */
						CurrentNode->neighbors[k]->travel_time[id] 	 = Travel_Time[id];
						/* Calcul le temps de trajet moyen */
						CurrentNode->neighbors[k]->avg_travel_time[id] 	 = (CurrentNode->avg_travel_time[id]  + (1.0 /(5./3. * parms.flow_speeds[id]))) * DIST(k);
						/* Calcul de la variance du temps de trajet moyen */
						CurrentNode->neighbors[k]->var_of_flow_time[id]  = (CurrentNode->var_of_flow_time[id]  + 2.0*parms.flow_disps[id]/pow((5./3. * parms.flow_speeds[id]),3.0)) * DIST(k);
						/* Récupèration de la portion d'aire drainant vers la cellule */
						CurrentNode->neighbors[k]->portion[id] = landscape[CurrentNode->row][CurrentNode->col].portion[k];
						/* Ajoute à la file d'attente */
						EnQueue(queue,CurrentNode->neighbors[k]);					
						}else{
						/* Connecte à un noeud existant */
						CurrentNode->neighbors[k] = GetNode(queue, rown, coln);	
							if( Travel_Time[id] < CurrentNode->neighbors[k]->travel_time[id] ){
								/* Recalcul du temps de trajet */							
								CurrentNode->neighbors[k]->travel_time[id] 	= Travel_Time[id];							
								/* Recalcul du temps de trajet moyen */
								CurrentNode->neighbors[k]->avg_travel_time[id]   = (CurrentNode->avg_travel_time[id] + (1.0 /(5./3. * parms.flow_speeds[id]))) * DIST(k);
								/* Recalcul de la variance du temps de trajet moyent */
								CurrentNode->neighbors[k]->var_of_flow_time[id]  = (CurrentNode->var_of_flow_time[id] + 2.0*parms.flow_disps[id]/pow((5./3. * parms.flow_speeds[id]),3.0)) * DIST(k);
								/* Récupération de la portion d'aire drainant vers la cellule */
								CurrentNode->neighbors[k]->portion[id] = landscape[CurrentNode->row][CurrentNode->col].portion[k];
							}
						}
					break;

				case 2:
					if(!IsEnqueued(queue, rown, coln)){
						/* Crée un nouveau noeud */
						CurrentNode->neighbors[k] 			= NewNode();
						/* Attribut le numéro de ligne et de colonne */
						CurrentNode->neighbors[k]->row 		= rown;
						CurrentNode->neighbors[k]->col 		= coln;
						/* Marque le noeud comme visité */
						CurrentNode->neighbors[k]->is_visited[id+1] = 1;						
						/* Calcul le temps de trajet */
						CurrentNode->neighbors[k]->travel_time[id+1] 		= Travel_Time[id+1];						
						/* Calcul le temps de trajet moyen */
						CurrentNode->neighbors[k]->avg_travel_time[id+1]   	= (CurrentNode->avg_travel_time[id+1]  + (1.0 /parms.flow_speeds[id+1])) * DIST(k);
						/* Calcul de la variance du temps de trajet moyen */
						CurrentNode->neighbors[k]->var_of_flow_time[id+1]  	= (CurrentNode->var_of_flow_time[id+1]  + 2.0*parms.flow_disps[id+1]/pow(parms.flow_speeds[id],3.0)) * DIST(k);
						/* Récupèration de la portion d'aire drainant vers la cellule */
						CurrentNode->neighbors[k]->portion[id+1] = landscape[CurrentNode->row][CurrentNode->col].portion[k];
						/* Ajoute à la file d'attente */
						EnQueue(queue,CurrentNode->neighbors[k]);					
						}else{
							/* Connecte à un noeud existant */
							CurrentNode->neighbors[k] = GetNode(queue, rown, coln);				
							if( Travel_Time[id+1] < CurrentNode->neighbors[k]->travel_time[id+1] ){
								/* Recalcul le temps de trajet */							
								CurrentNode->neighbors[k]->travel_time[id+1] 	= Travel_Time[id+1];							
								/* Recalcul du temps de trajet moyen */
								CurrentNode->neighbors[k]->avg_travel_time[id+1]  = (CurrentNode->avg_travel_time[id+1]  + (1.0 /parms.flow_speeds[id+1])) * DIST(k);
								/* Recalcul de la variance du temps de trajet moyen */
								CurrentNode->neighbors[k]->var_of_flow_time[id+1]  = (CurrentNode->var_of_flow_time[id+1]  + 2.0*parms.flow_disps[id+1]/pow(parms.flow_speeds[id+1],3.0)) * DIST(k);
								/* Récupèration de la portion d'aire drainant vers la cellule */
								CurrentNode->neighbors[k]->portion[id+1] = landscape[CurrentNode->row][CurrentNode->col].portion[k];	
							}
						}
					break;
					
				case 3:
					if(!IsEnqueued(queue, rown, coln)){
						/* Crée un nouveau noeud */
						CurrentNode->neighbors[k] 			= NewNode();
						/* Attribut le numéro de ligne et de colonne */
						CurrentNode->neighbors[k]->row 		= rown;
						CurrentNode->neighbors[k]->col 		= coln;
						/* Marque le noeud comme visité */
						CurrentNode->neighbors[k]->is_visited[id] = 1;
						CurrentNode->neighbors[k]->is_visited[id+1] = 1;
						/* Calcul le temps de trajet */
						CurrentNode->neighbors[k]->travel_time[id] 	 	= Travel_Time[id];
						CurrentNode->neighbors[k]->travel_time[id+1] 	= Travel_Time[id+1];						
						/* Calcul le temps de trajet moyen */
						CurrentNode->neighbors[k]->avg_travel_time[id]    = (CurrentNode->avg_travel_time[id] + (1.0 /parms.flow_speeds[id])) * DIST(k);
						CurrentNode->neighbors[k]->avg_travel_time[id+1]  = (CurrentNode->avg_travel_time[id+1] + (1.0 /parms.flow_speeds[id+1])) * DIST(k);						
						/* Calcul de la variance du temps de trajet moyen */
						CurrentNode->neighbors[k]->var_of_flow_time[id]   = (CurrentNode->var_of_flow_time[id] + 2.0*parms.flow_disps[id]/pow(parms.flow_speeds[id],3.0)) * DIST(k);
						CurrentNode->neighbors[k]->var_of_flow_time[id+1] = (CurrentNode->var_of_flow_time[id+1] + 2.0*parms.flow_disps[id+1]/pow(parms.flow_speeds[id+1],3.0)) * DIST(k);
						/* Récupèration de la portion d'aire drainant vers la cellule */
						CurrentNode->neighbors[k]->portion[id] 			  = landscape[CurrentNode->row][CurrentNode->col].portion[k];
						CurrentNode->neighbors[k]->portion[id+1] 		  = landscape[CurrentNode->row][CurrentNode->col].portion[k];
						/* Ajoute à la file d'attente */
						EnQueue(queue,CurrentNode->neighbors[k]);					
						}else{
							/* Connecte à un noeud existant */
							CurrentNode->neighbors[k] = GetNode(queue, rown, coln);				
							if( Travel_Time[id] < CurrentNode->neighbors[k]->travel_time[id] ){
								/* Recalcul le temps de trajet */							
								CurrentNode->neighbors[k]->travel_time[id] 	= Travel_Time[id];							
								/* Recalcul du temps de trajet moyen */
								CurrentNode->neighbors[k]->avg_travel_time[id]   = (CurrentNode->avg_travel_time[id] + (1.0 /(5./3. * parms.flow_speeds[id]))) * DIST(k);
								/* Recalcul de la variance du temps de trajet moyent */
								CurrentNode->neighbors[k]->var_of_flow_time[id]  = (CurrentNode->var_of_flow_time[id] + 2.0*parms.flow_disps[id]/pow((5./3. * parms.flow_speeds[id]),3.0)) * DIST(k);
								/* Récupèration de la portion d'aire drainant vers la cellule */
								CurrentNode->neighbors[k]->portion[id] 			 = landscape[CurrentNode->row][CurrentNode->col].portion[k];
							}
							if( Travel_Time[id+1] < CurrentNode->neighbors[k]->travel_time[id+1]  ){
								/* Recalcul le temps de trajet */							
								CurrentNode->neighbors[k]->travel_time[id+1] 	 = Travel_Time[id+1];							
								/* Recalcul du temps de trajet moyen */
								CurrentNode->neighbors[k]->avg_travel_time[id+1] = (CurrentNode->avg_travel_time[id+1] + (1.0 /parms.flow_speeds[id+1])) * DIST(k);
								/* Recalcul de la variance du temps de trajet moyen */
								CurrentNode->neighbors[k]->var_of_flow_time[id+1]= (CurrentNode->var_of_flow_time[id+1] + 2.0*parms.flow_disps[id+1]/pow(parms.flow_speeds[id+1],3.0)) * DIST(k);
								/* Récupèration de la portion d'aire drainant vers la cellule */
								CurrentNode->neighbors[k]->portion[id+1] 		 = landscape[CurrentNode->row][CurrentNode->col].portion[k];
							}							
						}
					break;
				}
			}
		/* Comptabilise la cellule et mets à jour l'aire de drainage */
		if(CurrentNode->is_visited[id]){
			p->nbContribCells[id]++;
			UpslopeArea[id] += CurrentNode->portion[id];
		}
		if(CurrentNode->is_visited[id+1]){				
			p->nbContribCells[id+1]++;
			UpslopeArea[id] += CurrentNode->portion[id+1];
		}
		/* Réajuste (augmente) la taille du bloc de mémoire contenant les noeuds */
		p->contribCells = (node *)G_realloc(p->contribCells, MAX(p->nbContribCells[id],p->nbContribCells[id+1])*sizeof(node));
		
		/* Ajoute le noeud dans le bloc de mémoire */
		memmove(p->contribCells + (MAX(p->nbContribCells[id],p->nbContribCells[id+1])-1), CurrentNode, sizeof(node));
		}
	/* Calcule la fonction de réponse UHT du bassin de drainage */
	 for(t=1;t<=num_inputs;t++)
	{
		 for(iter=0;iter<MAX(p->nbContribCells[id],p->nbContribCells[id+1]);iter++)
		{
			if(iter<p->nbContribCells[id])
				p->UHTsf[t] += p->contribCells[iter].portion[id]*FlowPathUnitResponse(&p->contribCells[iter],id,t);
			if(iter<p->nbContribCells[id+1])
				p->UHTssf[t] += p->contribCells[iter].portion[id+1]*FlowPathUnitResponse(&p->contribCells[iter],id+1,t);
		}
		if(p->nbContribCells[id])	
			p->UHTsf[t] /= UpslopeArea[id];
		if(p->nbContribCells[id+1])	
			p->UHTssf[t] /= UpslopeArea[id+1];
	}	
		
	/* Détruit la file d'attente */
	DestroyQueue(queue);

	return;
	}
	
	/* ********************************************** */
	/* Initialise la carte avec les options de calcul */
	/* ********************************************** */
	
	void Init(){

		/*****************
		* INITIALISATION *
		******************/
		AllocateMemory();
		ReadInputLayer();

		layer **ptr = landscape;

						
		G_verbose_message(_("Initialisation de la carte en cours..."));

		for (row = 0; row < nrows; row++)
		{
			G_percent(row, nrows, 2);

			for (col = 0; col < ncols; col++)
			{
				Segment_get(&parms_seg, &parms, row, col);
				
				/* Conditions initiales */
				ptr[row][col].swc[0]= parms.sat; /* Initialise la teneur en eau à la saturation */
				ptr[row][col].paw	= parms.rum; /* et la réserve utile à la réserve utile maximale */

				/* Identifie la couche comme une "zone humide" ou une "surface en eau" */
				if(flag6){
					ptr[row][col].waterbodies	= (short)waterbodies[row][col];
					ptr[row][col].riparian 		= (short)riparian[row][col];
				}
				/* Connecte la cellules à ses voisines à travers l'algorithme de calcul de l'aire de drainage amont */
				if(method>0){
					switch(algorithm){
						case 0:
							D_8(row, col);
						break;
						case 1:
							D_Inf(row, col);
						break;
						case 2:
							MFD_8(row, col);
						break;
						case 3:
							MFD_md(row, col);
						break;
						case 4:
							MFD_Inf(row, col);
						break;					
					}
					if(method==1||method==3){
						ptr[row][col].smax 	= smax[row][col];
						ptr[row][col].w1 	= w1[row][col];
						ptr[row][col].w2 	= w2[row][col];
						ptr[row][col].UHTsf = NULL;//dvector(1,num_inputs);
					}
					if(method>1){
						ptr[row][col].raw 	 = (double *)G_calloc(num_inputs, sizeof(double));
						ptr[row][col].braw 	 = NULL;//(double *)G_calloc(num_inputs*sizeof(double));						
						ptr[row][col].UHTssf = NULL;//dvector(1,num_inputs);
					}		
				}
			}
		}
	G_percent(1, 1, 1);
	Cleanup();
	
		if(method>0){
		G_verbose_message(_("Preparation de la carte pour le calcul du ruissellement..."));
				for (row = 0; row < nrows; row++)
			{	
				G_percent(row, nrows, 2);
					for (col = 0; col < ncols; col++)
				{
					Segment_get(&parms_seg, &parms, row, col);				
					FindBasin(&ptr[row][col]);
				}
			}
			G_percent(1, 1, 1);
		}
	
	}	
	
	/* *********************************** */
	/* Exécute le calcul du bilan hydrique */
	/* *********************************** */
	
	void Process(){	

		/**********
		* PROCESS *
		***********/
		
		double aet;
		int init 			= 0;
		double PE, SW, S, Qinsf, Qoutsf, Qinssf, Qoutssf;
		static layer *tmp 	=	NULL;
		static int *ptr		=	NULL;
		static int m, incr;
		
		struct input *prc	= NULL;
		struct input *etp	= NULL;
		struct output *out 	= NULL;
		layer **p 			= landscape;
		
		if( (flag4 && !flag3 && !flag5) || (flag4 && flag5) ){
			month 		= atoi(parm.start->answer)-1;
			sum_days 	= num_days[month];
		}
		
		clock_t start 		= clock();
		G_verbose_message(_("Calcul du bilan hydrique en cours..."));

		/* DEBUT BOUCLE TEMPORELLE (CARTES D ENTREE) */
		for (n = 0; n < num_inputs; n++){
		
		if(num_inputs>1)
			G_percent(n, num_inputs, 2);

		/* Ouvre les cartes d'entrée pour la lecture */
		prc 				= &P[n]; 	
		prc->name 			= parm.prec->answers[n];
		prc->buf 			= Rast_allocate_d_buf();
		prc->fd 			= Rast_open_old(prc->name, "");
		
		etp 				= &ETP[n];
		etp->name 			= parm.etp->answers[n];
		etp->buf 			= Rast_allocate_d_buf();
		etp->fd 			= Rast_open_old(etp->name, "");
		
		/* Ouvre les cartes de sortie à l'écriture */		
		if(options==1 || (options==2 && (n+1)==sum_days) || (options==3 && (n+1)%outiter==0) ){
			char *output_name;		
			for (i = 0; i < num_outputs_names; i++){
				out  			= &Outputs[init+i];
				if(num_outputs_names){
					 output_name = make_output_name(parm.outputs->answers[i]);
				}else{ 
					if(asprintf(&output_name, "%s%d%s", "PAW", n+1, menu[method].suffix)<0)
						fprintf(stderr, "Allocation de memoire lors de la creation du nom de raster de sortie a echoue\n");
				}
				if (G_legal_filename(output_name) < 0)
					G_fatal_error(_("<%s> est un nom de fichier illegal"), output_name);		
				out->name 		= output_name;
				out->buf 		= Rast_allocate_d_buf();
				out->fd 		= Rast_open_new(output_name, DCELL_TYPE);
			}		
		}
	
			/* DEBUT BOUCLE SPATIALE (LIGNES) */ 
			for (row = 0; row < nrows; row++) {

				if(num_inputs==1)
					G_percent(row, nrows, 2);
				
				Rast_get_d_row(P[n].fd, P[n].buf, row);
				Rast_get_d_row(ETP[n].fd, ETP[n].buf, row);

				/* DEBUT BOUCLE SPATIALE (COLONNES) */	
				for (col = 0; col < ncols; col++){
					
					int null = 0;

					double rain 	= (double)P[n].buf[col];
					double etp		= (double)ETP[n].buf[col];
					
					if (Rast_is_d_null_value(&rain) || Rast_is_d_null_value(&etp)){
						null = 1;
					}
					else{

							/* Récupère les données sur la cellule depuis le fichier segmenté */
							Segment_get(&parms_seg, &parms, row, col);

							/*****************************************
							 * Calcul de l'évapotranspiration réelle *
							 *****************************************/	
							 
							if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) ){
						
								p[row][col].p  += rain;
								p[row][col].pet+= etp;
										
								if(flag6){								
									if(p[row][col].waterbodies || rain-etp >= 0.0){
										p[row][col].aet += aet = etp;								
									}else if (parms.rum==0.0){
										p[row][col].aet += aet = EPS;
									}else{
										p[row][col].aet += aet = p[row][col].paw + rain - MAX(p[row][col].paw*exp((rain-etp)/parms.rum), 0.0);
									}
								}
								else if(rain-etp >= 0.0){
									p[row][col].aet += aet = etp;								
								}
								else if (parms.rum==0.0){
									p[row][col].aet += aet = EPS;
								}
								else{
									p[row][col].aet += aet = p[row][col].paw + rain - MAX(p[row][col].paw*exp((rain-etp)/parms.rum), 0.0);
								}
					
							}
							else{		
								p[row][col].p  = rain;
								p[row][col].pet= etp;
								
								if(flag6){								
									if(p[row][col].waterbodies || rain-etp >= 0.0){
										p[row][col].aet = aet = etp;								
									}
									else if (parms.rum==0.0){
										p[row][col].aet = aet = EPS;
									}
									else{
										p[row][col].aet = aet = p[row][col].paw + rain - MAX(p[row][col].paw*exp((rain-etp)/parms.rum), 0.0);
									}
								}else if(rain-etp >= 0.0){
									p[row][col].aet = aet = etp;								
								}
								else if (parms.rum==0.0){
									p[row][col].aet = aet = EPS;
								}else{
									p[row][col].aet = aet = p[row][col].paw + rain - MAX(p[row][col].paw*exp((rain-etp)/parms.rum), 0.0);
								}
							}				
				
								switch(method){

									case 0:
										/*************************************
										 * Calcul de la teneur en eau du sol *
										 *************************************/
										if(flag6){								
											if(p[row][col].waterbodies){
												p[row][col].swc[n]=parms.sat;								
											}else if (p[row][col].riparian){
												p[row][col].swc[n] = MAX( MIN(p[row][col].swc[n] + rain - aet, parms.sat), parms.fc);
											}else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet, parms.sat);
											}
										}else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet, parms.sat);
										}
										/*************************************
										 * Calcul de la réserve utile du sol *
										 *************************************/
										if(flag6){
											if(p[row][col].waterbodies || p[row][col].swc[n] >= parms.fc){
												p[row][col].paw = parms.rum;
											}else{
												p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]), 0.0);
											}
										}else if(p[row][col].swc[n] >= parms.fc){
											p[row][col].paw = parms.rum;
										}else{
											p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]), 0.0);
										}					
									break;
									
									
									case 1:									
										/**************************************
										 * Calcul du ruissellement de surface *
										 **************************************/
										Qinsf = 0.0, Qoutsf = 0.0;
										//static int m, incr;
										/*static layer *tmp 	=	NULL;
										static int *ptr		=	NULL;*/

										/* calcule l'eau disponible en surface */
										SW 	= MAX(p[row][col].swc[n] - parms.pwp, 0.0);
										S 	= p[row][col].smax * (1.0 - SW/(SW+exp(p[row][col].w1 - p[row][col].w2*SW)));

										if(method_ia){						
											if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) )
												p[row][col].sraw += PE = rain>0.05*S ? pow(rain - 0.05*(1.33*pow(S,1.15)),2.0)/(rain + 0.95*(1.33*pow(S,1.15))) : 0.0;
											else 
												p[row][col].sraw = PE = rain>0.05*S ? pow(rain - 0.05*(1.33*pow(S,1.15)),2.0)/(rain + 0.95*(1.33*pow(S,1.15))) : 0.0;
										}else{	
											if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) )
												p[row][col].sraw += PE = rain>0.2*S ? pow(rain - 0.2*S,2.0)/(rain + 0.8*S) : 0.0;
											else
												p[row][col].sraw = PE = rain>0.2*S ? pow(rain - 0.2*S,2.0)/(rain + 0.8*S) : 0.0;
										}
										
										/* calcule le ruissellement entrant et sortant */
										m = 1;
										for(iter=0;iter<p[row][col].nbContribCells[id];iter++){						
											tmp = &landscape[p[row][col].contribCells[iter].row][p[row][col].contribCells[iter].col];
											//ptr = FindNonZeroTermIndices(tmp->sraw, num_inputs);
												if(ptr==NULL) continue;
													//for(n=ptr[0];n<=num_inputs;n++)
													//{
														//for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
														//{
															if(iter==0)
																Qoutsf +=  tmp->sraw * CellOutletResponse(&p[row][col].contribCells[iter], n-m+1, id);											
															else Qinsf +=  tmp->sraw * FlowPathUnitResponse(&p[row][col].contribCells[iter], n-m+1, id);
														//}
													//}
										}
										//for(iter=0;iter<p[row][col].nbContribCells[id];iter++)						
										//	p[row][col].braw[n] += landscape[p[row][col].contribCells[iter]->row][p[row][col].contribCells[iter]->col].raw[n];
										
										//	ptr = FindNonZeroTermIndices(p[row][col].braw[n], num_inputs);
										//		if(ptr){
										//			for(n=ptr[0];n<=num_inputs;n++)
										//			{
										//				for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
										//				{
										//					if(iter==0)
										//						Qoutssf +=  tmp->sraw * CellOutletResponse(p[row][col].contribCells[iter], n-m+1, id);											
										//					else Qinssf +=  p[row][col].braw[m] * UHTsf[n-m+1];
										//				}
										//			} 
										//		}
																				
										if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) ){
											p[row][col].qinsf  += Qinsf;
											p[row][col].qoutsf += Qoutsf;
										}
										else{
											p[row][col].qinsf  = Qinsf;
											p[row][col].qoutsf = Qoutsf;
										}
										//for(iter=1;iter<p[row][col].nbContribCells;iter++)
										//	Qinsf  +=  tmp.sraw * qromb(&FlowPathUnitResponse, p[row][col].contribCells[iter], 0.0, draining_time);
										//  Qoutsf +=  Qinsf * qromb(&CellOutletResponse, p[row][col].contribCells[0], 0.0, draining_time);
										/*************************************
										 * Calcul de la teneur en eau du sol *
										 *************************************/
										if(flag6){								
											if(p[row][col].waterbodies){
												p[row][col].swc[n]=parms.sat;								
											}else if (p[row][col].riparian){
												p[row][col].swc[n] = MAX( MIN(p[row][col].swc[n] + rain - aet + Qinsf - Qoutsf, parms.sat),parms.fc);
											}else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet + Qinsf - Qoutsf, parms.sat);
											}
										}else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet + Qinsf - Qoutsf, parms.sat);
										}
										/*************************************
										 * Calcul de la réserve utile du sol *
										 *************************************/
										if(flag6){
											if(p[row][col].waterbodies || p[row][col].swc[n] >= parms.fc){
												p[row][col].paw = parms.rum;
											}else{
												p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]),0.0);
											}
										}else if(p[row][col].swc[n] >= parms.fc){
											p[row][col].paw = parms.rum; 
										}else{
											p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]),0.0);
										}
									break;
									
									
									case 2:									
										/*****************************************
										 * Calcul du ruissellement de subsurface *
										 *****************************************/
										Qinssf = 0.0, Qoutssf = 0.0;
										//static int m, incr;
										/*static layer *tmp 	=	NULL;
										static int *ptr		=	NULL;*/

										/* calcule le ruissellement entrant et sortant */
										for(iter=0;iter<p[row][col].nbContribCells[id+1];iter++){						
											tmp = &landscape[p[row][col].contribCells[iter].row][p[row][col].contribCells[iter].col];
											ptr = FindNonZeroTermIndices(tmp->raw, num_inputs);
												if(!ptr) continue;
													/*for(n=ptr[0];n<=num_inputs;n++)
													{*/
														for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
														{
															if(iter==0)
																Qoutssf +=  tmp->raw[m] * CellOutletResponse(&p[row][col].contribCells[iter], n-m+1, id+1);											
															else Qinssf +=  tmp->raw[m] * FlowPathUnitResponse(&p[row][col].contribCells[iter], n-m+1, id+1);
														}
													/*}*/
										}
										
										// for(iter=0;iter<p[row][col].nbContribCells[id+1];iter++)						
										//	p[row][col].braw[n] += landscape[p[row][col].contribCells[iter]->row][p[row][col].contribCells[iter]->col].raw[n];
										
										//	ptr = FindNonZeroTermIndices(p[row][col].braw[n], num_inputs);
										//		if(ptr){
										//			for(n=ptr[0];n<=num_inputs;n++)
										//			{
										//				for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
										//				{
										//					if(iter==0)
										//						Qoutssf +=  tmp->raw[m] * CellOutletResponse(p[row][col].contribCells[iter], n-m+1, id+1);											
										//					else Qinssf +=  p[row][col].braw[m] * UHTssf[n-m+1];
										//				}
										//			} 
										//		}


										if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) ){
											p[row][col].qinssf  += Qinssf;
											p[row][col].qoutssf += Qoutssf;
										}
										else{
											p[row][col].qinssf  = Qinssf;
											p[row][col].qoutssf = Qoutssf;
										}

										/*************************************
										 * Calcul de la teneur en eau du sol *
										 *************************************/
										if(flag6){								
											if(p[row][col].waterbodies){
												p[row][col].swc[n] = parms.sat;								
											}
											else if (p[row][col].riparian){
												p[row][col].swc[n] = MAX( MIN(p[row][col].swc[n] + rain - aet + Qinssf - Qoutssf, parms.sat),parms.fc);
											}
											else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet + Qinssf - Qoutssf, parms.sat);
											}
										}
										else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet + Qinssf - Qoutssf, parms.sat);
										}
										/*************************************
										 * Calcul de la réserve utile du sol *
										 *************************************/
										if(flag6){
											if(p[row][col].waterbodies || p[row][col].swc[n] >= parms.fc){
												p[row][col].paw = parms.rum;
											}
											else{
												p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]),0.0);
											}
										}
										else if(p[row][col].swc[n] >= parms.fc){
											p[row][col].paw = parms.rum; 
										}
										else{
											p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]),0.0);
										}				
									break;
									
									
									case 3:									
										/**************************************
										 * Calcul du ruissellement de surface *
										 **************************************/
										Qinsf = 0.0, Qoutsf = 0.0;
										/*static int m, incr;*/
										/*static layer *tmp 	=	NULL;
										static int *ptr		=	NULL;*/

										/* calcule l'eau disponible en surface */
										SW 	= MAX(p[row][col].swc[n] - parms.pwp, 0.0);
										S 	= p[row][col].smax * (1.0 - SW/(SW+exp(p[row][col].w1 - p[row][col].w2*SW)));
										
										if(method_ia){						
											if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) )
												p[row][col].sraw += PE = rain>0.05*S ? pow(rain - 0.05*(1.33*pow(S,1.15)),2.0)/(rain + 0.95*(1.33*pow(S,1.15))) : 0.0;
											else 
												p[row][col].sraw = PE = rain>0.05*S ? pow(rain - 0.05*(1.33*pow(S,1.15)),2.0)/(rain + 0.95*(1.33*pow(S,1.15))) : 0.0;
										}else{	
											if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) )
												p[row][col].sraw += PE = rain>0.2*S ? pow(rain - 0.2*S,2.0)/(rain + 0.8*S) : 0.0;
											else
												p[row][col].sraw = PE = rain>0.2*S ? pow(rain - 0.2*S,2.0)/(rain + 0.8*S) : 0.0;
										}
										
										/* calcule le ruissellement entrant et sortant */
										m = 1;
										for(iter=0;iter<p[row][col].nbContribCells[id];iter++){						
											tmp = &landscape[p[row][col].contribCells[iter].row][p[row][col].contribCells[iter].col];
											//ptr = FindNonZeroTermIndices(tmp->sraw, num_inputs);
												//if(!ptr) continue;
													/*for(n=ptr[0];n<=num_inputs;n++)
													{*/
														for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
														{
															if(iter==0)
																Qoutsf +=  tmp->sraw * CellOutletResponse(&p[row][col].contribCells[iter], n-m+1, id);											
															else Qinsf +=  tmp->sraw * FlowPathUnitResponse(&p[row][col].contribCells[iter], n-m+1, id);
														}
													/*}*/
										}
										if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) ){
											p[row][col].qinsf  += Qinsf;
											p[row][col].qoutsf += Qoutsf;
										}
										else{
											p[row][col].qinsf  = Qinsf;
											p[row][col].qoutsf = Qoutsf;
										}
										//for(iter=1;iter<p[row][col]->nbContribCells;iter++)
										//	Qinsf  +=  tmp->sraw * qromb(&FlowPathUnitResponse, p[row][col]->contribCells[iter], 0.0, draining_time);
										//  Qoutsf +=  Qinsf * qromb(&CellOutletResponse, p[row][col]->contribCells[0], 0.0, draining_time);						
										/*****************************************
										 * Calcul du ruissellement de subsurface *
										 *****************************************/				 
										Qinssf = 0.0, Qoutssf = 0.0;
										//static int m, incr;
										/*static layer *tmp 	=	NULL;
										static int *ptr		=	NULL;*/

										/* calcule le ruissellement entrant et sortant */
										for(iter=0;iter<p[row][col].nbContribCells[id+1];iter++){						
											tmp = &landscape[p[row][col].contribCells[iter].row][p[row][col].contribCells[iter].col];
											ptr = FindNonZeroTermIndices(tmp->raw, num_inputs);
												if(!ptr) continue;
													/*for(n=ptr[0];n<=num_inputs;n++)
													{*/
														for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
														{
															if(iter==0)
																Qoutssf +=  tmp->raw[m] * CellOutletResponse(&p[row][col].contribCells[iter], n-m+1, id+1);
															else Qinssf +=  tmp->raw[m] * FlowPathUnitResponse(&p[row][col].contribCells[iter], n-m+1, id+1);
														}
													/*}*/
										}
										/*for(iter=0;iter<p[row][col].nbContribCells[id+1];iter++)						
											p[row][col].braw[n] += landscape[p[row][col].contribCells[iter]->row][p[row][col].contribCells[iter]->col].raw[n];
										
											ptr = FindNonZeroTermIndices(p[row][col].braw[n], num_inputs);
												if(ptr){
													for(n=ptr[0];n<=num_inputs;n++)
													{
														for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
														{
															if(iter==0)
																Qoutssf +=  tmp->raw[m] * CellOutletResponse(p[row][col].contribCells[iter], n-m+1, id+1);											
															else Qinssf +=  p[row][col].braw[m] * UHTssf[n-m+1];
														}
													} 
												}
										*/						
										
										if( (options==2 && (n+1)<=sum_days) || (options==3 && (n+1)%outiter!=0) ){
											p[row][col].qinssf  += Qinssf;
											p[row][col].qoutssf += Qoutssf;
										}else{
											p[row][col].qinssf  = Qinssf;
											p[row][col].qoutssf = Qoutssf;
										}

										/*************************************
										 * Calcul de la teneur en eau du sol *
										 *************************************/
										if(flag6){								
											if(p[row][col].waterbodies){
												p[row][col].swc[n]=parms.sat;								
											}
											else if (p[row][col].riparian){
												p[row][col].swc[n] = MAX( MIN(p[row][col].swc[n] + rain - aet + Qinsf - Qoutsf + Qinssf - Qoutssf, parms.sat),parms.fc);
											}
											else{
												p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet + Qinsf - Qoutsf + Qinssf - Qoutssf, parms.sat);
											}
										}
										else{
											p[row][col].swc[n] = MIN(p[row][col].swc[n] + rain - aet + Qinsf - Qoutsf + Qinssf - Qoutssf, parms.sat);
										}
										/*************************************
										 * Calcul de la réserve utile du sol *
										 *************************************/
										if(flag6){
											if(p[row][col].waterbodies || p[row][col].swc[n] >= parms.fc){
												p[row][col].paw = parms.rum;
											}
											else{
												p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]),0.0);
											}
										}
										else if(p[row][col].swc[n] >= parms.fc){
											p[row][col].paw = parms.rum; 
										}
										else{
											p[row][col].paw = MIN(parms.rum - (parms.fc - p[row][col].swc[n]),0.0);
										}					
										break;
								}
					}
					// fin du else (null ou pas)					
					/* Inscrit le calcul dans la carte de sortie */
					if(options==1 || (options==2 && (n+1)==sum_days) || (options==3 && (n+1)%outiter==0) ){
						
						double swc; 
						int origin = (options==2) ? n - num_days[(month>12)?month%12:month] : (options==3) ? (n+1) - outiter : 0;
					
						for (i = 0; i < num_outputs_names; i++){
							out  = &Outputs[init+i];
							if(null){
								Rast_set_d_null_value(&out->buf[col], 1);
							}
							else{
								int output_options = find_output_name(parm.outputs->answers[i]);
								switch(output_options){
									case 0:
										out->buf[col] = (options==1) ? (DCELL)(etp-aet):(DCELL)(p[row][col].pet-p[row][col].aet);
									break;
									
									case 1:
										if(options==1 || !origin){
											out->buf[col] = (DCELL)(p[row][col].paw);
										}
										else{ 
											if(flag6){								
												if(p[row][col].waterbodies){
													swc 	= parms.sat;								
												}else if (p[row][col].riparian){
													swc 	= MAX(MIN(p[row][col].swc[origin] + p[row][col].p - p[row][col].aet + p[row][col].qinsf - p[row][col].qoutsf + p[row][col].qinssf - p[row][col].qoutssf, parms.sat), parms.fc);
												}else{
													swc 	= MIN(p[row][col].swc[origin] + p[row][col].p - p[row][col].aet + p[row][col].qinsf - p[row][col].qoutsf + p[row][col].qinssf - p[row][col].qoutssf, parms.sat);
												}
											}else{
													swc 	= MIN(p[row][col].swc[origin] + p[row][col].p - p[row][col].aet + p[row][col].qinsf - p[row][col].qoutsf + p[row][col].qinssf - p[row][col].qoutssf, parms.sat);
												}
											if(flag6){
												if(p[row][col].waterbodies || swc >= parms.fc){
													out->buf[col] = (DCELL)parms.rum;
												}else{
													out->buf[col] = (DCELL) MIN(parms.rum - (parms.fc - swc),0.0);
												}
											}else if(swc >= parms.fc){
												out->buf[col] = (DCELL)parms.rum; 
											}else{
												out->buf[col] = (DCELL)MIN(parms.rum - (parms.fc - swc),0.0);
											}
										}
									break;
									
									case 2:
										if(options==1){
											out->buf[col] 			= (DCELL)p[row][col].swc[origin];
										}
										else if(flag6){								
											if(p[row][col].waterbodies){
												out->buf[col] 	= (DCELL)parms.sat;								
											}
											else if (p[row][col].riparian){
												out->buf[col] 	= (DCELL) (MAX(MIN(p[row][col].swc[origin] + p[row][col].p - p[row][col].aet + p[row][col].qinsf - p[row][col].qoutsf + p[row][col].qinssf - p[row][col].qoutssf, parms.sat), parms.fc));
											}
											else{
												out->buf[col] 	= (DCELL) (MIN(p[row][col].swc[origin] + p[row][col].p - p[row][col].aet + p[row][col].qinsf - p[row][col].qoutsf + p[row][col].qinssf - p[row][col].qoutssf, parms.sat));
											}
										}
										else{
											out->buf[col] 	= (DCELL) (MIN(p[row][col].swc[origin] + p[row][col].p - p[row][col].aet + p[row][col].qinsf - p[row][col].qoutsf + p[row][col].qinssf - p[row][col].qoutssf, parms.sat));
										}
									break;
									
									case 3:
										out->buf[col] = (options==1) ? (DCELL)Qinssf:(DCELL)(p[row][col].qinssf);
									break;
									
									case 4:
										out->buf[col] = (options==1) ? (DCELL)Qoutssf:(DCELL)(p[row][col].qoutssf);
									break;
									
									case 5:
										out->buf[col] = (options==1) ? (DCELL)Qinsf:(DCELL)(p[row][col].qinsf);
									break;
									
									case 6:
										out->buf[col] = (options==1) ? (DCELL)Qoutsf:(DCELL)(p[row][col].qoutsf);
									break;
									
									case 7:
										out->buf[col] = (options==1) ? (DCELL)PE:(DCELL)(p[row][col].sraw);
									break;							
								}
							}
						}
					}
					// fin calcul sortie
				}				
				/* FIN BOUCLE SPATIALE (COLONNES) */
				
				/* Inscrit la ligne dans la carte de sortie */
				if(options==1 || (options==2 && (n+1)==sum_days) || (options==3 && (n+1)%outiter==0) ){
					for (i = 0; i < num_outputs_names; i++)
						Rast_put_d_row(Outputs[init+i].fd, Outputs[init+i].buf);
				}
				if(num_inputs==1)
					G_percent(1, 1, 1);
			}
			/* FIN BOUCLE SPATIALE (LIGNES) */
				
			/* Ferme les cartes d'entrée */	
			Rast_close(P[n].fd);
			Rast_close(ETP[n].fd);
			
			/* Ferme les cartes de sortie */
			if(options==1 || (options==2 && (n+1)==sum_days) || (options==3 && (n+1)%outiter==0) ){
				for (i = 0; i < num_outputs_names; i++)
				{
					out = &Outputs[init+i];
					Rast_close(out->fd);
					Rast_short_history(out->name, "raster", &history);
					Rast_command_history(&history);
					Rast_write_history(out->name, &history);
					G_verbose_message(_("La carte raster <%s> a ete cree"), out->name);
					G_free(out->name);
				}
			}
			
			/* (Re)Définit l'indice mémoire du prochain point d'écriture */
			if(options==1 || (options==2 && (n+1)==sum_days) || (options==3 && (n+1)%outiter==0) )		
				init = init + num_outputs_names;
				
			/* Détermine le prochain point d'écriture */
			if(options==2 && (n+1)>sum_days)
				sum_days += num_days[(month<12)?++month:month%12];

			if(num_inputs>1)
				G_percent(1,1,1);
		}
		/* FIN BOUCLE TEMPORELLE (CARTES D ENTREE) */

		/* Libère la mémoire */
		Segment_close(&parms_seg);
		FreeLandscape();

		clock_t end 		= clock();
		double time_elapsed = (double)(end - start);
		G_verbose_message(_("Temps ecoule pour le calcul: %fsec soit %fmin\n"), time_elapsed/1000.0, time_elapsed/60000.0);
		G_done_msg(_("Le calcul du bilan hydrique est a present termine."));
			
	return;
 
}

	
	/* ____________________________
	   |03      |02      |01      |
	   |        |        |        |
	   |        |  north |        |        
	   |        |        |        |
	   |________|________|________|          
	   |04      |        |00      |
	   |        |        |        |
	   |  east  |    x   |  west  |
	   |        |        |        |
	   |________|________|________|
	   |05      |06      |07      |
	   |        |        |        |
	   |        |  south |        |
	   |        |        |        |
	   |________|________|________|
	 */
		 
		// /* write colors for aspect file */
	// Rast_init_colors(&colors);
	// Rast_read_fp_range(aspect_name, G_mapset(), &range);
	// Rast_get_fp_range_min_max(&range, &min, &max);
	// Rast_make_aspect_fp_colors(&colors, min, max);
	// Rast_write_colors(aspect_name, G_mapset(), &colors);
	
	// Rast_append_format_history(&hist, "aspect map elev = %s", elev_name);
	// Rast_append_format_history(&hist, "zfactor = %.2f", zfactor);
	// Rast_append_format_history(&hist, "min_slope = %f", min_slope);	
	
	
		// /* colortable for slopes */
	//	    struct Colors colors;
	// CELL val1, val2;
	// Rast_init_colors(&colors);
	// val1 = 0;
	// val2 = 2;
	// Rast_add_c_color_rule(&val1, 255, 255, 255, &val2, 255, 255, 0, &colors);
	// val1 = 2;
	// val2 = 5;
	// Rast_add_c_color_rule(&val1, 255, 255, 0, &val2, 0, 255, 0, &colors);
	// val1 = 5;
	// val2 = 10;
	// Rast_add_c_color_rule(&val1, 0, 255, 0, &val2, 0, 255, 255, &colors);
	// val1 = 10;
	// val2 = 15;
	// Rast_add_c_color_rule(&val1, 0, 255, 255, &val2, 0, 0, 255, &colors);
	// val1 = 15;
	// val2 = 30;
	// Rast_add_c_color_rule(&val1, 0, 0, 255, &val2, 255, 0, 255, &colors);
	// val1 = 30;
	// val2 = 50;
	// Rast_add_c_color_rule(&val1, 255, 0, 255, &val2, 255, 0, 0, &colors);
	// val1 = 50;
	// val2 = 90;
	// Rast_add_c_color_rule(&val1, 255, 0, 0, &val2, 0, 0, 0, &colors);
	
	//    Rast_write_colors(parm.dsout, mapset, &colors);
    //Rast_free_colors(&colors);
	
	/* Crée un vecteur de proportion de flux */
	//static double d[8];
	//for(k=0;k<8;k++) double d[k] = MFD_md(k);
				
		//if(CurrentNode->row=row && CurrentNode->col=col) CurrentNode->prop = d[k];

			/* Attribut la même proportion de flux que le noeud courant */
		//CurrentNode->neighbors[k]->prop = CurrentNode->prop;
			/* Ré-Attribut la même proportion de flux que le noeud courant */
			//CurrentNode->neighbors[k]->prop = CurrentNode->prop;
//flowAccum += a[index]->available_runoff * qromb(&FlowPathUnitResponse, a[index]->average_travel_time, a[index]->variance_of_flow_time, 0.0, draining_time);

	/* ********************** */
	/* Ecrit la carte raster  */
	/* ********************** */	
		
	// void writeRaster(){

		// Déclaration
		// double swb_out;
		// out_data_type 	= DCELL_TYPE;
		// int out_dsize 	= G_raster_size(out_data_type);
		// void *out_cell 	= G_allocate_buf(out_data_type);	pointeurs vers un bloc de mémoire
		// void *ptr;										pointeurs pour l'incrementation des colonnes du raster
		
		// Ouvre la couche de bilan hydrique à l'écriture
		// int swb_fd = G_open_raster_new(parm.outputs->answers[n], out_data_type);
		// G_message(_("Ecrit la carte raster <%s>..."), parm.outputs->answers[n]);

			// for (row = 0; row < nrows; row++) {
			// G_percent(row, nrows, 2);

			// INPUT NULL VALUES: ???
			// ptr = out_cell;

				// for (col = 0; col < ncols; col++){
					// if(keep_nulls){
						// if (G_is_null_value(ptr, out_data_type)) {
						// G_set_d_null_value(ptr, 1, out_data_type);
						// ptr = G_incr_void_ptr(ptr, out_dsize);
						// continue;
						// }
					// }				
					// Segment_get(&parms_seg, &parms, row, col);
					// swb_out = parms.swbulated_runoff;
					// if (G_is_d_null_value(&swb_out)) {
						// G_set_d_null_value(ptr, 1, out_data_type);
					// }
					// else{
							// switch (out_data_type){
						  
							// case CELL_TYPE:
								// *(CELL *)ptr = (CELL)(swb_out);
								// break;
							// case FCELL_TYPE:
								// *(FCELL *)ptr = (FCELL)(swb_out);
								// break;
							// case DCELL_TYPE:
								// *(DCELL *)ptr = (DCELL)(swb_out);
								// break;
							// }
						// }
					// ptr = G_incr_void_ptr(ptr, out_dsize);
				// }
				// G_put_row(swb_fd, out_cell, out_data_type)
			// }
		// G_percent(1, 1, 1);
		// G_free(out_cell);
		
		// Segment_close(&parms_seg);   release memory 
		// G_close_cell(swb_fd);
						
		// Ajoute les lignes de commande invoquees dans le fichier historique
		// G_short_history(parm.outputs->answers[n], "raster", &history);
		// G_command_history(&history);
		// G_write_history(parm.outputs->answers[n], &history);
		
		// G_done_msg(_("La carte raster <%s> a ete creee."), parm.outputs->answers[n]);
	
// return;
// }

// /* 		case 2:		
	
		// /* définit le gradient de pente maximum */
			// for(k=0; k<8; k++)
	    // {
			// rown   = row + 2*dy[k];
			// coln   = col + 2*dx[k];
			
			// Segment_get(&parms_seg, &parms, rown, coln);
				// if( is_OnGrid(rown, coln) && (dz = z - parms.altitude) > 0.0 )
			// {
				// if( (double) (dz / DIST(k)) > dzMax) dzMax = (double) (dz / DIST(k)); 
			// }
		
		// }
	
		// /* Calcule et assigne la proportion de flux drainant vers chaque cellule voisine située à l'aval */	
			// for(k=0; k<8; k++)
		// {
			// tanBeta[k] = 0.0;
			// rown  = row + 2*dy[k];
			// coln  = col + 2*dx[k];

			// Segment_get(&parms_seg, &parms, rown, coln);
				// if( is_OnGrid(rown, coln) && (dz = z - parms.altitude) > 0.0 )
			// {
				// dzSum += tanBeta[k] = ( pow( (double) (dz / DIST(k)), MFD_md_pow(dzMax)) ) * l[k];
			// }
				
		// }
			// if( dzSum > 0.0 )
		// {
			// rown  = row + 2*dy[(dir+8/2)%8];
			// coln  = col + 2*dx[(dir+8/2)%8];
			
			// if( is_OnGrid(rown, coln) && tanBeta[(dir+8/2)%8] > 0.0)
			// {
				// p = tanBeta[(dir+8/2)%8] / dzSum;
			// } else p = 0.0;
		// } else p = 0.0;	 */				
	/*}*/
	
	
	
						// /* calcule le ruissellement sortant */
						// p[col]->raw[n] 	= MAX(p[col]->swc[n] - parms.fc, 0.0);
						// ptr				= FindNonZeroTermIndices(p[col]->raw, num_inputs);	
						// if(ptr!=NULL){
							// for(n=ptr[0];n<=num_inputs;n++)
								// for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
									// Qoutssf +=  p[col]->raw[m] * CellOutlaetesponse(p[col]->contribCells[0], n+m-1);
							// ;	
						// }
						// p[col]->qoutssf += Qoutssf;
						
						/*double Ds 		= (p[col]->swc<=parms.fc)? 0.0 : parms.depth * ((p[col]->swc-parms.fc)/(parms.sat-parms.fc));
						p[col]->qoutssf = MIN((parms.ksat*(Ds * RES)*parms.tanslope*1000.)/(RES*RES), p[col]->raw);*/
						

						// /* calcule le ruissellement sortant */
						// p[col]->raw[n] 	= MAX(p[col]->swc[n] - parms.fc, 0.0);
						// ptr				= FindNonZeroTermIndices(p[col]->raw, num_inputs);	
						// if(ptr!=NULL){
							// for(n=ptr[0];n<=num_inputs;n++)
								// for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
									// Qoutssf +=  p[col]->raw[m] * CellOutlaetesponse(p[col]->contribCells[0], n+m-1);
							// ;	
						// }
						// p[col]->qoutssf += Qoutssf;
						
						/*double Ds 		= (p[col]->swc<=parms.fc)? 0.0 : parms.depth * ((p[col]->swc-parms.fc)/(parms.sat-parms.fc));
						p[col]->qoutssf = MIN((parms.ksat*(Ds * RES)*parms.tanslope*1000.)/(RES*RES), p[col]->raw);*/						
	

						// /* calcule le ruissellement sortant */
						// ptr				= FindNonZeroTermIndices(p[col]->sraw, num_inputs);	
						// if(ptr!=NULL){
							// for(n=ptr[0];n<=num_inputs;n++)
								// for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
									// Qoutsf +=  p[col]->sraw[m] * CellOutlaetesponse(p[col]->contribCells[0], n+m-1);
							// ;	
						// }
						// p[col]->qoutsf += Qoutsf;	
						
						// /* calcule le ruissellement sortant */
						// ptr				= FindNonZeroTermIndices(p[col]->sraw, num_inputs);	
						// if(ptr!=NULL){
							// for(n=ptr[0];n<=num_inputs;n++)
								// for(m=ptr[0],incr=0;m<=n;m=ptr[incr++])
									// Qoutsf +=  p[col]->sraw[m] * CellOutlaetesponse(p[col]->contribCells[0], n+m-1);
							// ;	
						// }
						// p[col]->qoutsf += Qoutsf;

	
	

	
	// void createSEGMENT(){

    // G_verbose_message(_("Cree un fichier temporaire..."));
	
    // if (Segment_open(&parms_seg, G_tempfile(), nrows, ncols, srows, scols, sizeof(layer), segments_in_memory) != 1)
		// G_fatal_error(_("Ne peux pas creer le fichier temporaire"));
	
	// /* DECLARE */
	// int fd;
	// int data_size;
	// int skip_nulls;
	// double dval;
	// void *cell;
	// void *ptr;
	// RASTER_MAP_TYPE data_type;
	// struct Options *structMembers[10] 			= {parm.altitude,parm.slope,parm.depth,parm.sat,parm.fc,parm.pwp,parm.rum,parm.ksat,parm.flow_speeds,parm.flow_disps};
	// double *segParms[NUM_ELEM(structMembers)+2] = {&parms.altitude,&parms.tanslope,&parms.depth,&parms.sat,&parms.fc,&parms.pwp,&parms.rum,&parms.ksat,&parms.flow_speeds[0],&parms.flow_speeds[1],&parms.flow_disps[0],&parms.flow_disps[1]};

	// /* INITIALISE */
		// for(i=0;i<NUM_ELEM(structMembers);i++){

		//Fichiers descripteurs des couches raster	
	// FD: fd 			= if(i>8) ? openLayer(structMembers[i]->answers[id]) : openLayer(structMembers[i]->answer);
		//Types de données
		// data_type 	= G_get_raster_map_type(fd);
		//Taille des données	
		// data_size 	= G_cell_size(data_type);
		//Alloue de la mémoire et initialise l'adresse des pointeurs
		// cell		= G_allocate_buf(data_type);
		//Fixe les valeurs nulles
		// G_set_null_value(&nullval, 1, data_type);
		//Esquive les valeurs nulles ??
		// skip_nulls 	= G_is_null_value(&null_val,data_type);

		// dval = 0.0;
		// total_cells = nrows * ncols;

			// for (row = 0; row < nrows; row++) {
					
			// G_percent(row, nrows, 2);

			// if( G_get_row(fd, cell, row, data_type) < 0 )
				// G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), structMembers[i]->answer, row);

			// ptr  = cell;
			
				// for (col = 0; col < ncols; col++){		
					// if (G_is_null_value(ptr, data_type)){
						// dval = null_val;
						//if (skip_nulls && !G_is_null_value(ptr, data_type)) total_cells--;
					// }
					// else{
						// switch (data_type){
						
						// case CELL_TYPE:
							// dval = *(CELL *)ptr;
							// dval = (i==2) ? tan(dval*DEG_TO_RAD) : (i==8) ? dval*DAY_TO_HOUR : dval;							
							// break;
						// case FCELL_TYPE:
							// dval = *(FCELL *)ptr;
							// dval = (i==2) ? tan(dval*DEG_TO_RAD) : (i==8) ? dval*DAY_TO_HOUR : dval;							
							// break;
						// case DCELL_TYPE:
							// dval = *(DCELL *)ptr;
							// dval = (i==2) ? tan(dval*DEG_TO_RAD) : (i==8) ? dval*DAY_TO_HOUR : dval;
							// break;
						// }
					// }
				// *segParms[(id>0)?(i+1):i] =  dval; 
				//Assigne les valeurs 
				// Segment_put(&parms_seg, segParms[(id>0)?(i+1):i], row, col);
				//Incrémente les pointeurs
				// ptr  = G_incr_void_ptr(ptr,  data_size);					
				// }
			// }
		// if(i>8){
			// if(id>0){
				// id--;
				// goto END;
			// }else{
				// id++;
				// goto FD;
			// }
		// }			
	// END:G_free(cell);
		// G_close_cell(fd);		
		// G_percent(1, 1, 1);				
		// }
	// return;
	// }
	
	
	
	// void ReadInputLayer(){

	// if(flag6){
	 
		// DECLARATION
		// int wetlands_fd, riparian_fd;			Fichiers descripteurs des couches rasters
		// int wetlands_dsize, riparian_dsize;		Taille du type de donnees raster
		// void *wetlands_cell, *riparian_cell;	Pointeur vers les blocs de (mémoire) lignes du raster
		// RASTER_MAP_TYPE wetlands_data_type, riparian_data_type;
		
		// INITIALISATION
		// wetlands_fd = openLayer(parm.wetlands->answer);
		// riparian_fd = openLayer(parm.riparian->answer);

		// wetlands_data_type 	= G_get_raster_map_type(wetlands_fd);
		// riparian_data_type 	= G_get_raster_map_type(riparian_fd);
		
		// wetlands_dsize 		= G_raster_size(wetlands_data_type);
		// riparian_dsize 		= G_raster_size(riparian_data_type);
		
		// wetlands_cell 		= G_allocate_raster_buf(wetlands_data_type);
		// riparian_cell 		= G_allocate_raster_buf(riparian_data_type);
		
		// void *ptr_wetlands, *ptr_riparian;					//pointeur pour l'incrementation des colonnes du raster
		// G_important_message(_("Lecture des cartes rasters <%s> et <%s>..."), parm.wetlands->answer,parm.riparian->answer);
		// for (row = 0; row < nrows; row++){
			// G_percent(row, nrows, 2);	
			// if (G_get_raster_row(wetlands_fd, wetlands_cell, row, wetlands_data_type) < 0)
				// G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.wetlands->answer, row);
			// if (G_get_raster_row(riparian_fd, riparian_cell, row, riparian_data_type) < 0)
				// G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.riparian->answer, row);
				
			// ptr_wetlands = wetlands_cell;
			// ptr_riparian = riparian_cell;
			// for (col = 0; col < ncols; col++){
				// wetlands[row][col] = G_get_raster_value_d(ptr_wetlands, wetlands_data_type);
				// riparian[row][col] = G_get_raster_value_d(ptr_riparian, riparian_data_type);
				
				// ptr_wetlands = G_incr_void_ptr(ptr_wetlands, wetlands_dsize);	//avance le pointeur de cellule
				// ptr_riparian = G_incr_void_ptr(ptr_riparian, riparian_dsize);
			// }
		// }
		// G_percent(1, 1, 1);
		
		// cleanup
		// G_close_cell(wetlands_fd);
		// G_close_cell(riparian_fd);
		// G_free(wetlands_cell);
		// G_free(riparian_cell);
	// }
	
		// if(method==2||method==4){
		
		// DECLARATION
		// int smax_dsize, w1_dsize, w2_dsize;
		// int smax_fd, w1_fd, w2_fd;
		// void *smax_cell, *w1_cell, *w2_cell;
		// RASTER_MAP_TYPE smax_data_type, w1_data_type, w2_data_type;
		
		// INITIALISATION
		// smax_fd 		= openLayer(parm.smax->answer);
		// w1_fd 			= openLayer(parm.w->answers[1]);
		// w2_fd 			= openLayer(parm.w->answers[2]);
		
		// smax_data_type 	= G_get_raster_map_type(smax_fd);
		// w1_data_type 	= G_get_raster_map_type(w1_fd);
		// w2_data_type 	= G_get_raster_map_type(w2_fd);
		
		// smax_dsize 		= G_raster_size(smax_data_type);
		// w1_dsize 		= G_raster_size(w1_data_type);
		// w2_dsize 		= G_raster_size(w2_data_type);
		
		// smax_cell 		= G_allocate_raster_buf(smax_data_type);
		// w1_cell 		= G_allocate_raster_buf(w1_data_type);
		// w2_cell 		= G_allocate_raster_buf(w2_data_type);
		
		// void *ptr_smax, *ptr_w1, *ptr_w2;					//pointeur pour l'incrementation des colonnes du raster
		// G_important_message(_("Lecture des cartes rasters <%s>, <%s> et <%s>..."), parm.smax->answer,parm.w1->answer,parm.w2->answer);
		
		// for (row = 0; row < nrows; row++){
			// G_percent(row, nrows, 2);	
			// if (G_get_raster_row(smax_fd, smax_cell, row, smax_data_type) < 0)
				// G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.smax->answer, row);
			// if (G_get_raster_row(w1_fd, w1_cell, row, w1_data_type) < 0)
				// G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.w1->answer, row);
			// if (G_get_raster_row(w2_fd, w2_cell, row, w2_data_type) < 0)
				// G_fatal_error(_("Impossible de lire la carte raster <%s> ligne %d"), parm.w2->answer, row);
				
			// ptr_smax 	= smax_cell;
			// ptr_w1 		= w1_cell;
			// ptr_w2 		= w2_cell;
			// for (col = 0; col < ncols; col++){
				// smax[row][col] 	= G_get_raster_value_d(ptr_smax, smax_data_type);
				// w1[row][col] 	= G_get_raster_value_d(ptr_w1, w1_data_type);
				// w2[row][col] 	= G_get_raster_value_d(ptr_w2, w2_data_type);
				
				// ptr_smax 	= G_incr_void_ptr(ptr_smax, smax_dsize);	//avance le pointeur de cellule
				// ptr_w1 		= G_incr_void_ptr(ptr_w1, w1_dsize);
				// ptr_w2 		= G_incr_void_ptr(ptr_w2, w2_dsize);
			// }
		// }
		// G_percent(1, 1, 1);
	   // }
	// return;
	// }

	
/*	void MFD_Inf(int row, int col){

	Segment_get(&parms_seg, &parms, row, col);
	
	int   	kk;
	double 	z = parms.altitude, dzSum, dz[8];
	double 	hs, hd, s_facet[8], d_facet[8];
	double 	valley[8], portion[8];
	double 	cellarea = RES * RES;
	double 	nx, ny, nz, n_norm;
	int  	Dir_inGrid[8];

	 for(k=0; k<8; k++)
	{
		rown          = row + dy[k];
		coln          = col + dx[k];
		Dir_inGrid[k] = 1;
		s_facet[k] = d_facet[k] = -999.0;
		
		Segment_get(&parms_seg, &parms, rown, coln);

		if( is_OnGrid(rown, coln) )
			dz[k] = z - parms.altitude;
		else{
			dz[k]   		= 0.0;
			Dir_inGrid[k] 	= 0;
		}
	}
	  
      for(k=0; k<8; k++)
    {
		hd = hs = -999.0;
		
		if( Dir_inGrid[k] )
		{
			kk = ( k < 7)? k+1 : 0;
			  
			if( Dir_inGrid[kk] )
			{
				nx = ( dz[kk] * dy[k] - dz[k] * dy[kk]) * RES;
				ny = ( dz[k] * dx[kk] - dz[kk] * dx[k]) * RES;
				nz = ( dx[k] * dy[kk] - dx[kk] * dy[k]) * cellarea;
				
				n_norm = sqrt( nx*nx + ny*ny + nz*nz );

				if( nx == 0.0 )
				{
				  hd = (ny >= 0.0)? 0.0 : M_PI;
				} 
				else if( nx < 0.0 )
				{
				  hd = M_PI_270 - atan(ny / nx);
				} 
				else
				{
				  hd = M_PI_090 - atan(ny / nx);
				}
				hs = -tan( acos( nz/n_norm ) );
				//SHOULD IT BE LIKE THIS: (( hr <= i * M_PI_045 || hr >= ii * M_PI_045 )  OR AS BELOW???
				if( hd < k * M_PI_045 || hd > (k+1) * M_PI_045 )
				{
					if( dz[k] > dz[kk] ){
						hd = k * M_PI_045;
						hs = dz[k] / DIST(k);
					}else{
						hd = kk * M_PI_045;
						hs = dz[kk] / DIST(kk);                               
					}
				}
			}
			else if( dz[k] > 0.0 ){
				hd = k * M_PI_045;
				hs = dz[k] / DIST(k);
			}
		 s_facet[k] = hs;
		 d_facet[k] = hd;
		}
    }
      
	dzSum  = 0.0;
      for(k=0; k<8; k++)
    {           
		valley[k]   = 0.0;
		kk = (k < 7)? k+1 : 0;
            
		if( s_facet[k] > 0.0 )
		{
			if( d_facet[k] > k * M_PI_045 && d_facet[k] < (k+1) * M_PI_045 ){
				valley[k] = s_facet[k];
			}else if( d_facet[k] == d_facet[kk] ){
				valley[k] = s_facet[k];
			}else if( s_facet[kk] == -999.0 && d_facet[k] == (k+1) * M_PI_045){
				valley[k] = s_facet[k];
			}else{
				kk = (k > 0)? k-1 : 7;
				if( s_facet[kk] == -999.0 && d_facet[k] == k * M_PI_045 ){
					valley[k] = s_facet[k];
				}
			}
		dzSum += valley[k] = pow(valley[k], MFD_Converge);
		} 
    portion[k] = 0.0;
    }

      if( dzSum )
    {
		for(k=0; k<8; k++)
		{
			if (k < 7){
				kk = k+1;
			}else{
				kk = 0;
				if( d_facet[k] == 0.0) 
					d_facet[k] = M_PI_360;
			}
			if( valley[k] ){
				valley[k] /= dzSum;
				portion[k] += valley[k] * ((k+1) * M_PI_045 - d_facet[k]) / M_PI_045;
				portion[kk]+= valley[k] * (d_facet[k] - k * M_PI_045) / M_PI_045;
			}
		}
		for(k=0; k<8; k++)
		{
			landscape[row][col].prop[k]							= portion[k];
			landscape[row+dy[k]][col+dx[k]].neighbors[(k+4)%8]	= &landscape[row][col];
		}
    }
}
*/	
/*
	DCELL aspect_fly(DCELL * n, DCELL * c, DCELL * s, double d)
{
    double xslope = ((n[-1] + c[-1] + c[-1] + s[-1]) -
		     (n[1] + c[1] + c[1] + s[1])) / (8 * d);
    double yslope = ((s[-1] + s[0] + s[0] + s[1]) -
		     (n[-1] + n[0] + n[0] + n[1])) / (8 * region.ns_res);
    double asp;

    if (!yslope)
	if (!xslope)
	    asp = UNDEF;
	else if (xslope > 0)
	    asp = parm.up ? 270. : 90.;
	else
	    asp = parm.up ? 90. : 270.;
    else if ((asp = atan2(xslope, yslope) / DEG2RAD) < 0.)
	asp += 360.;

    return asp;
}
*/
	
	
// SinSlopedeg * KhSAT * Profondeur * 50 / ((ThetaS - ThetaFC) * SinSlopedeg * 50)
// SinSlopedeg * KhSAT / (ThetaS - ThetaFC)
						