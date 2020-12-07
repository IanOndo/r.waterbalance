/***************************************************************************************************************************************************************************************************************************
 *
 * MODULE:       r.flowrouting
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

/***********************************************************************************************
 *
 *				FlowRouting.h
 *				Ce fichier d'en-tête déclare les fonctions, variables et structures des données
 *				utilisées par la fonction principale du programme du module r.flowrouting
 *
 ***********************************************************************************************/
 
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <grass/gis.h>
#include "Queue.h"

#ifndef _HEAD_H
#define _HEAD_H

/*******************************************************
 * DECLARATIONS DES TYPES, VARIABLES, STRUCTURES...ETC *
 *******************************************************/
 
// Prototype de structure
typedef struct SoilLayer layer;		   				/* définit le type layer qui a la structure SoilLayer */
typedef struct Node  node;  							/* définit le type node qui a la structure Node */

// Définit la structure d'une couche de sol
struct SoilLayer
{
	// Quantité d'eau contenue dans la couche
	double *swc;
	
	// Quantité d'eau disponible pour les plantes ou le ruissellement dans la couche et dans le bassin versant
	double paw, *raw, *braw, sraw;

	// Quantité d'eau précipitée ou evapotranspirée
	double p, pet, aet;

	// Quantité d'eau drainée vers/depuis la couche par ruissellement de surface/subsurface
	double qinsf, qinssf, qoutsf, qoutssf;
	
	// Paramètres pour le calcul du ruissellement de surface
	double smax, w1, w2;
	
	// Portion de l'aire de drainage amont de la cellule
	double portion[8];

	// Pointeur vers les cellules contribuant au ruissellement dans la cellule en amont
	layer *neighbors[8];
	
	// Nombre de cellules contribuant au ruissellement dans la cellule
	int nbContribCells[2];
	
	// Cellules contribuant au ruissellement dans la cellule 
	node *contribCells;
	
	// Pointeurs vers les ordonnées à l'origine de l'hydrographe du bassin versant drainé en amont par la cellule 
	double *UHTsf, *UHTssf;
	
 	// Status
	short waterbodies, riparian;
};

struct Parm
{
	double altitude;
	double tanslope, depth;
	double sat, fc, pwp, rum;
	double ksat, flow_speeds[2], flow_disps[2];
}parms;

struct Node
{
	int row, col;
	int is_visited[2];
    double avg_travel_time[2];									/* Temps de trajet moyen de la couche jusqu'à un exutoire situé en aval */
    double var_of_flow_time[2]; 								/* Variance du temps d'écoulement le long d'un trajet */
	double travel_time[2];
	double portion[2];
	node *neighbors[8];
};

struct input
{
    const char *name;
    int fd;
    DCELL *buf;
};

struct output
{
    char *name;
    int fd;
    DCELL *buf;
};

struct Cell_head window;									/* Stocke les informations sur la région et les informations d'en-tête des couches rasters */
extern struct Cell_head window;
struct GModule *module;										/* Module GRASS pour les arguments d'analyse */
struct Flag *flag, *flag2, *flag3, *flag4, *flag5, *flag6;	/* Drapeau GRASS pour spécifier des options supplémentaires */
struct History history;     								/* Contient les méta-données (titres, commentaires,...) */
struct
{	
	struct Option *prec, *etp;
	struct Option *altitude, *slope, *depth;								
	struct Option *sat, *fc, *pwp, *rum, *ksat;
	struct Option *flow_speeds, *flow_disps;
	struct Option *smax, *w;
	struct Option *waterbodies, *riparian;	
	struct Option *outputs;									
	struct Option *drainage_times;//, *upper, *lower;
	struct Option *convergence;	
	struct Option *start;
	struct Option *method, *algorithm;
	struct Option *init_abs;
	struct Option *outiter;
	struct Option *mem;
} parm;	

struct menu
{	
	char	*suffix;					/* suffix pour les couches de sortie */
    char 	*name;                  /* nom de la méthode  */
    char 	*text;                  /* Affichage du menu - description complète */
} menu[] = {
    {"_CLM", 	"climat",    			"bilan hydrique climatique i.e P + ETP"},
    {"_SF",		"surface_account",   	"bilan hydrique climatique + ruissellement de surface"},
    {"_SSF",	"subsurface_account",   "bilan hydrique climatique + ruissellement de subsurface"},
    {"_FULL",	"full_account",       	"bilan hydrique complet"},
    {NULL,		NULL,      				NULL}
};

struct menu_outputs
{
    char 	*name;                  /* nom de la couche de sortie */
    char 	*text;                  /* Affichage du menu - description complète */
} menu_outputs[] = {
    {"DE", 		"deficit of evapotranspiration"},
    {"PAW",		"plant available water"},
    {"SWC",		"soil water content"},
    {"QINSSF",	"lateral subsurface inflow"},
	{"QOUTSSF",	"lateral subsurface outflow"},
	{"QINSF",	"lateral surface inflow"},
	{"QOUTSF",	"lateral surface outflow"},
	{"PE",		"precipitation in excess"},
    {NULL,		NULL}
};

struct menu_algorithm
{
    char 	*name;                  /* nom de la couche de sortie */
    char 	*text;                  /* Affichage du menu - description complète */
} menu_algorithm[] = {
    {"D8", 		"Deterministic 8:\n - O'Callaghan, J.F. / Mark, D.M. (1984):\n 'The extraction of drainage networks from digital elevation data',\nComputer Vision, Graphics and Image Processing, 28:323-344\n\n"},
    {"DInf",	"Deterministic Infinity:\n - Tarboton, D.G. (1997):\n 'A new method for the determination of flow directions and upslope areas in grid digital elevation models',\n Water Ressources Research, Vol.33, No.2, p.309-319\n\n"},
    {"MFD8",	"Multiple Flow Direction:\n - Freeman, G.T. (1991):\n 'Calculating catchment area with divergent flow based on a regular grid',\n Computers and Geosciences, 17:413-22\n\n"},
    {"MFDmd",	"Multiple Flow Direction based on maximum downslope gradient:\n - Qin C, Zhu AX, Pei T, Li B, Zhou C, Yang L (2007):\n 'An adaptive approach to selecting a flow-partition exponent for a multiple-flow-direction algorithm', \n International Journal of Geographical Information Science, 21:443-458\n\n"},
	{"MFDInf",	"Triangular Multiple Flow Direction:\n - Seibert, J. / McGlynn, B. (2007):\n 'A new triangular multiple flow direction algorithm for computing upslope areas from gridded digital elevation models', Water Resources Research, Vol. 43, W04501"},
    {NULL,		NULL}
};

struct menu_ia
{	
    char 	*name;                  /* nom de la méthode */
    char 	*text;                  /* Affichage du menu - description complète */
} menu_ia[] = {
    {"traditional",    	"Ia = 0.2 * S -{USDA SCS 1972}-"},
    {"alternative",   	"Ia = 0.05 * S -{Woodward et al. 2003}- "},
    {NULL,      		NULL}
};

int i, n;
short k;
long iter;
int id = 0;
int row;                    								/* Numéro de la ligne */
int col;                    								/* Numéro de la colonne */
double RES;													/* Résolution de la carte */
int rown;													/* Ligne voisine */
int coln;													/* Colonne voisine */
int nrows;                  								/* Nombre total de lignes */
int ncols;                  								/* Nombre total de colonnes */
int srows, scols;
const int num_days[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
const short dy[8] = {0,-1,-1,-1,0,1,1,1};
const short dx[8] = {1,1,0,-1,-1,-1,0,1};
layer **landscape = NULL;
struct input *P = NULL;
struct input *ETP = NULL;
struct output *Outputs = NULL;
DCELL *Prec =NULL, *Etp = NULL;
CELL **waterbodies = NULL, **riparian = NULL;
DCELL **smax = NULL, **w1 = NULL, **w2 = NULL;
void *ptr, *ptr2, *ptr3, *ptr4, *ptr5, *ptr6, *ptr7, *ptr7, *ptr8, *ptr9, *ptr10, *ptr11, *ptr12;
double p_sat, p_fc, p_rum, p_depth, p_alt, p_pwp, p_slope, p_speed_sf, p_disp_sf, p_speed_ssf, p_disp_ssf, p_ksat;
int method, method_ia;
int algorithm;
int outiter;
int month, sum_days;
double mfd_converge, drainage_times[2];
int num_inputs;
int num_outputs, num_outputs_names;
int options;
int keep_nulls = 1;
double null_val, dnullval;
double disk_mb, mem_mb, pq_mb;
int nseg;
int maxmem;
int segments_in_memory;
int total_cells;

SEGMENT parms_seg;

/* Type de données d'entrée (CELL/FCELL/DCELL [entier,décimale,double décimale]) */	
RASTER_MAP_TYPE alt_data_type, speed_sf_data_type, disp_sf_data_type, speed_ssf_data_type, disp_ssf_data_type, sat_data_type, fc_data_type, rum_data_type, pwp_data_type, slope_data_type, depth_data_type, ksat_data_type; 



/******************************
 * DECLARATIONS DES FONCTIONS *
 ******************************/
 
void parseOptions(int argc, char *argv[]);
void createSEGMENT();
static char *build_method_list(void);
static char *build_outputs_list(void);
static char *build_algorithm_list(void);
char *make_output_name(const char *output_name);
int  openLayer(char *layer);
int is_OnGrid(int rown, int coln);
static int find_method(const char *method_name);
static int find_output_name(const char *output_name);
static int find_ia_method(const char *method_name);
static int find_algorithm_method(const char *algorithm_name);
double aspect_on_fly(int row, int col);
void D_8(int row, int col);
void D_Inf(int row, int col);
void MFD_8(int row, int col);
void MFD_Inf(int row, int col);
void MFD_md(int row, int col);
void AllocateMemory();
void ReadInputLayer();
void Cleanup();
layer *NewLayer();
void FreeLandscape();
int *FindNonZeroTermIndices(double *p, int size);
double FlowPathUnitResponse(node *p, int time_index, int id);
double CellOutletResponse(node *p, int time_index, int id);
double DIST(short dir);
int IsEnqueued(Queue *queue, int rown, int coln);
node *NewNode();
node *GetNode(Queue *queue, int rown, int coln);
void FindBasin(layer *a);
void Init();
void Process();

 #endif