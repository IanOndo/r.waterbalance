#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"

#define END 1
#define FREE_ARG char*


double dstdDev(double data[], int n)
{
    double sum = 0.0, stdDev = 0.0, mean;

    for(int i=0; i<n; ++i)
        sum += data[i];

    mean = sum/n;

    for(i=0; i<n; ++i)
        stdDev += pow(data[i] - mean, 2);

    return sqrt(stdDev/n);
}




float *vector(long nl, long nh)
/* alloue un vecteur de décimaux allant de [nl..à..nh]*/
{
	float *ptr;
	ptr = (float *)malloc((size_t) ((nh-nl+1+END)*sizeof(float)));
	if(!ptr){
	fprintf(stderr, "Allocation de memoire pour le vecteur a echoue\n");
	exit(1);
	}
	return ptr-nl+END;
}

double *dvector(long nl, long nh)
/* alloue un vecteur de décimaux allant de [nl..à..nh]*/
{
	double *ptr;
	ptr = (double *)malloc((size_t) ((nh-nl+1+END)*sizeof(double)));
	if(!ptr){
	fprintf(stderr, "Allocation de memoire pour le vecteur a echoue\n");
	exit(1);
	}
	return ptr-nl+END;
}

int *ivector(long nl, long nh)
/* alloue un vecteur d'entier allant de [nl..à..nh]*/
{
	int *ptr;
	ptr = (int *)malloc((size_t) ((nh-nl+1+END)*sizeof(int)));
	if(!ptr){
	fprintf(stderr, "Allocation de memoire pour le vecteur a echoue\n");
	exit(1);
	}
	return ptr-nl+END;
}

unsigned char *cvector(long nl, long nh)
/* alloue un vecteur de caractère allant de [nl..à..nh]*/
{
	unsigned char *ptr;
	ptr = (unsigned char *)malloc((size_t) ((nh-nl+1+END)*sizeof(unsigned char)));
	if(!ptr){
	fprintf(stderr, "Allocation de memoire pour le vecteur a echoue\n");
	exit(1);
	}
	return ptr-nl+END;
}

unsigned long *lvector(long nl, long nh)
/* alloue un vecteur d'entier long allant de [nl..à..nh]*/
{
	unsigned long *ptr;
	ptr = (unsigned long *)malloc((size_t) ((nh-nl+1+END)*sizeof(unsigned long)));
	if(!ptr){
	fprintf(stderr, "Allocation de memoire pour le vecteur a echoue\n");
	exit(1);
	}
	return ptr-nl+END;
}

void free_vector(float *v, long nl, long nh)
/* Déalloue la mémoire alloué avec la fonction vector() */
{
	free((FREE_ARG) (v+nl-END));
	return;
}
	
void free_ivector(int *v, long nl, long nh)
/* Déalloue la mémoire alloué avec la fonction ivector() */
{
	free((FREE_ARG) (v+nl-END));
	return;
}
				
void free_cvector(unsigned char *v, long nl, long nh)
/* Déalloue la mémoire alloué avec la fonction cvector() */
{
	free((FREE_ARG) (v+nl-END));
	return;
}
	
void free_lvector(unsigned long *v, long nl, long nh)
/* Déalloue la mémoire alloué avec la fonction lvector() */
{
	free((FREE_ARG) (v+nl-END));
	return;
}
	
void free_dvector(double *v, long nl, long nh)
/* Déalloue la mémoire alloué avec la fonction dvector() */
{
	free((FREE_ARG) (v+nl-END));
	return;			
}

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
/* Connaissant les vecteurs xa[1..n] et ya[1..n], et une valeur x, la routine retourne une valeur y, et une erreur d'estimation dy.
Si P(x) est le polynôme de degré N-1 tel que P(xai) = yai, i = 1, . . . , n, alors la valeur retournée y est y = P(x).*/
{
	int i, m, ns=1;			
	double den, dif, dift, ho, hp, w;
	double *c, *d;
	dif = fabs(x-xa[1]);
	c = dvector(1,n);
	d = dvector(1,n);
	
	for (i=1;i<=n;i++) { 									/* Ici on trouve l'index ns de la table d'entrée la plus proche, et itnitialise le tableau des valeurs de c et de d.*/
		if ( (dift=fabs(x-xa[i])) < dif) {
		ns = i;
		dif = dift;
		}
		c[i] = ya[i]; 
		d[i] = ya[i];
	}
	*y = ya[ns--]; 											/* Approximation initiale de y. */
	for (m=1;m<n;m++) { 									/* Pour chaque colonne du tableau, on boucle sur les c et d actuels et on les met à jour */
		for (i=1;i<=n-m;i++) { 
			ho = xa[i]-x; 
			hp = xa[i+m]-x;
			w = c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {fprintf(stderr,"Erreur dans la routine polint\n");exit(1);}
			/* L'erreur se produit uniquement si les deux valeurs xa d'entrée sont identiques. */
			den = w/den;
			d[i] = hp*den; 									/* Ici les c et les d sont mis à jour. */
			c[i] = ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
		/* After each column in the tableau is completed, we decide which correction, c or d,
		we want to add to our accumulating value of y, i.e., which path to take through the
		tableau—forking up or down. We do this in such a way as to take the most “straight
		line” route through the tableau to its apex, updating ns accordingly to keep track of
		where we are. This route keeps the partial approximations centered (insofar as possible)
		on the target x. The last dy added is thus the error indication. */
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}

#define FUNC(a, b, x) (*func)(a, b, x)
#define EPS 1.0e-5
#define JMAX 20

/**************************
* Méthode du point milieu *
***************************/

 double midpnt (double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n)
 /* La routine calcule le nième niveau de reprécision de la règle étendue du point médian. func est un pointeur d'entrée vers la fonction à intégrer entre
les limites a et b qui sont aussi des valeurs d'entrée. Pour n=1, la routine retourne l'estimation la plus grossière de l'intégrale.
Les appels aux valeurs n=2,3,...(par ordre séquentiel) améliorera la précision en ajoutant (2/3)*3^n-1 points supplémentaires pour le calcul. s ne doit pas être modifiée
entre les appels successifs */
{
	double x, tnm, sum, del, ddel;
	static double s;
	int it, j;
	
	if(n == 1){
	return (s=(b-a)*FUNC(time, sigma, 0.5*(a+b)));
	} else {
	for(it=1,j=1;j<n-1;j++) it *= 3;
	tnm = it;
	del = (b-a)/(3.0*tnm);					/* les espaces entre les points ajoutés alternent entre del et ddel */
	ddel = del+del;
	x = a+0.5*del;
	sum = 0.0;
	for(j=1;j<=it;j++){
		sum += FUNC(time, sigma, x);
		x +=ddel;
		sum += FUNC(time, sigma, x);
		x += del;
	}
	s = (s+(b-a)*sum/tnm)/3.0;				/* la nouvelle somme est combinée avec l'ancienne intégrale pour donner une intégrale plus réajustée */
	return s;
	}
}

/*********************
* Méthode du trapèze *
**********************/

double trapzd(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n)
/* La routine calcule le nieme niveau de reprécision de la règle trapézoidale étendue. func est un pointeur d'entrée vers la fonction à intégrer
entre les limites a et b, qui sont aussi des valeurs d'entrée. Pour n=1, la routine retourne l'estimation la plus grossière de l'intégrale.
Les appels aux valeurs n=2,3,...(par ordre séquentiel) améliorera la précision en ajoutant 2^n-2 points supplémentaires pour le calcul. */
{
	double x,tnm,sum,del;
	static double s;
	int it,j;
	
	if(n == 1) {
		return (s=0.5*(b-a)*(FUNC(time, sigma, a)+FUNC(time, sigma, b)));
	} else {
	for(it=1,j=1;j<n-1;j++) it <<= 1;
	tnm = it;
	del = (b-a)/tnm;		/* Espacement entre les points à ajouter */
	x = a+0.5*del;
	for(sum=0.0,j=1;j<=it;j++,x+=del) sum+= FUNC(time, sigma, x);
	s = 0.5*(s+(b-a)*sum/tnm); /* Remplace s par sa valeur réajustée */
	return s;
	}
}

double qtrap(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b)
/* Retourne l'intégrale de la fonction func de a à b. Le paramètre EPS peut être fixé à la précision désirée et JMAX tel que 2 à la puissance JMAX-1 est le 
nombre de pas maximum autorisé. L'intégration est effectuée par la règle du trapèze. */
{
	double trapzd(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n);
	int j;
	double s, olds;
	
	olds = -1.0e-30;			/* N'importe quel nombre improbable d'etre la moyenne de la fonction à ces points terminaux */
	for(j=1;j<=JMAX;j++){
	s = trapzd(func, time, sigma, a, b, j);
	if(j>5)						/* Evite une convergence ennuyeuse trop précoce */
		if(fabs(s-olds) < EPS*fabs(olds) || (s == 0.0 && olds == 0.0)) return s;
	olds = s;
	}
	fprintf(stderr,"Trop d iteration dans la routine qtrap\n");exit(1);
	return 0.0;
}
#undef JMAX

#define JMAX 10
double qtrap_modif(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b)
/* Retourne l'intégrale de la fonction func de a à b. Le paramètre EPS peut être fixé à la précision désirée et JMAX tel que 2 à la puissance JMAX-1 est le 
nombre de pas maximum autorisé. L'intégration est effectuée par la règle du trapèze. */
{
	double midpnt (double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n);
	int j;
	double s, olds;
	
	olds = -1.0e-30;			/* N'importe quel nombre improbable d'etre la moyenne de la fonction à ces points terminaux */
	for(j=1;j<=JMAX;j++){
	s = trapzd(func, time, sigma, a, b, j);
	if(j>5)						/* Evite une convergence ennuyeuse trop précoce */
		if(fabs(s-olds) < EPS*fabs(olds) || (s == 0.0 && olds == 0.0)) return s;
	olds = s;
	}
	fprintf(stderr,"Trop d iteration dans la routine qtrap modifiee\n");exit(1);
	return 0.0;
}
#undef JMAX
#undef EPS

#define JMAX 20
#define EPS 1.0e-5

/***********
 * SIMPSON *
 ***********/
	
double qsimp(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b)
/* Retourne l'intégrale de la fonction func de a à b. Le paramètre EPS peut être fixé à la précision désirée et JMAX tel que 2 à la puissance JMAX-1 est le 
nombre de pas maximum autorisé. L'intégration est effectuée par la règle de Simpson. */
{
	double trapzd(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n);
	int j;
	double s, st, ost=0.0, os=0.0;
	
	for(j=1;j<=JMAX;j++){
	st = trapzd(func, time, sigma, a, b, j);
	s = (4.0*st-ost)/3.0;
	if(j>5)						/* Evite une convergence ennuyeuse trop précoce */
		if(fabs(s-os) < EPS*fabs(os) || (s == 0.0 && os == 0.0)) return s;
	os = s;
	ost = st;
	}
	fprintf(stderr,"Trop d iteration dans la routine qsimp\n");exit(1);
	return 0.0;
}
#undef JMAX

#define JMAX 15
double qsimp_modif(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b)
/* Retourne l'intégrale de la fonction func de a à b. Le paramètre EPS peut être fixé à la précision désirée et JMAX tel que 2 à la puissance JMAX-1 est le 
nombre de pas maximum autorisé. L'intégration est effectuée par la règle de Simpson avec une la fonction d'intégration du point milieu qui remplace
celle de la méthode du trapèze dans la version classique. L'étape d'extrapolation est également changé. */
{
	double midpnt(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n);
	int j;
	double s, st, ost=0.0, os=0.0;
	
	for(j=1;j<=JMAX;j++){
	st = trapzd(func, time, sigma, a, b, j);
	s = (9.0*st-ost)/8.0;		/* remplace (4.0*st-ost)/3.0 */
	if(j>5)						/* Evite une convergence ennuyeuse trop précoce */
		if(fabs(s-os) < EPS*fabs(os) || (s == 0.0 && os == 0.0)) return s;
	os = s;
	ost = st;
	}
	fprintf(stderr,"Trop d iteration dans la routine qsimp modifiee\n");exit(1);
	return 0.0;
}
#undef JMAX

#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

/* Ici EPS est la précision désirée, comme déterminée par l'estimation de l'erreur d'extrapolation ;
JMAX limite le nombre total de pas; K est le nombre de points utilisés dans l'extrapolation. */

double qromb(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b)
/* Retourne l'intégrale de la fonction func de a à b. L'intégration est effectué par la méthode de Romberg d'ordre 2K,
 où, par exemple K=2 est la règle de Simpson.*/
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, int n);
	double ss, dss, s[JMAXP], h[JMAXP+1]; 			/* Stocke les approximations trapezoidales successives et leur longueurs de pas successives.*/
	int j; 
	
	h[1] = 1.0;
	for (j=1;j<=JMAX;j++) {
		s[j] = trapzd(func, time, sigma, a, b, j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1] = 0.25*h[j];
		/* C'est l'étape clé: le facteur est 0.25 même si la longueur du pas est diminuée par 0.5 seulement. 
		Cela permet l'extrapolation d'un polynôme en h2, pas juste un polynôme en h.*/
	}
	fprintf(stderr,"Trop d iterations dans la routine qromb\n");exit(1);
	return 0.0; 
}
#undef JMAX

#define JMAX 14

double qromo(double (*func)(double time, double sigma, double t), double time, double sigma, double a, double b, double (*choose)(double(*)(double, double, double), double, double, double, double, int))
/* Intégration par la méthode de Romberg sur un interval ouvert. Retourne l'intégrale de la fonction func de a à b, en utilisant une fonction d'intégration
spécifiée choose et la méthode de Romberg. Normalement choose sera une formule ouverte, n'évaluant pas la fonction à ses points terminaux.
On considère que choose triple le nombre de pas à chaque appel. Les routines midpnt, midinf, midsql, midsqu, midexp, sont les choix possibles pour choose.
Les paramètres ont la même signification que dans la routine qromb. */
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	int j;
	double ss, dss, h[JMAXP+1], s[JMAXP];
	
	h[1] = 1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func, time, sigma, a, b, j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1] = h[j]/9.0; /* C'est ici que l'hypothèse du triplage du pas et de l'érreur de série est utilisée. */
	} 
	fprintf(stderr, "Trop d iteration dans la routine qromo\n");exit(1);
	return 0.0; 
}
#undef EPS
#undef JMAX
