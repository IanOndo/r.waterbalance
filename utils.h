#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef _UTILS_H
#define _UTILS_H

float *vector(long nl, long nh);						/* Alloue un vecteur de décimaux avec une précision simple allant de nl à nh [nl..nh].*/
int *ivector(long nl, long nh); 						/* Alloue un vecteur d'entier allant de nl à nh [nl..nh].*/
unsigned char *cvector(long nl, long nh); 				/* Alloue un vecteur de caractère non signé allant de nl à nh [nl..nh].*/
unsigned long *lvector(long nl, long nh); 				/* Alloue un vecteur d'entier long non signé allant de nl à nh [nl..nh].*/
double *dvector(long nl, long nh); 						/* Alloue un vecteur de décimaux avec une précision double allant de nl à nh [nl..nh].*/
void free_vector(float *v, long nl, long nh);			/* Déalloue un vecteur de décimaux avec une précision simple allant de nl à nh [nl..nh].*/
void free_ivector(int *v, long nl, long nh);			/* Déalloue un vecteur d'entier allant de nl à nh [nl..nh].*/
void free_cvector(unsigned char *v, long nl, long nh);	/* Déalloue un vecteur de caractère non signé allant de nl à nh [nl..nh].*/
void free_lvector(unsigned long *v, long nl, long nh);	/* Déalloue un vecteur d'entier long non signé allant de nl à nh [nl..nh].*/
void free_dvector(double *v, long nl, long nh);			/* Déalloue un vecteur de décimaux avec une précision double allant de nl à nh [nl..nh].*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

double trapzd(double (*func)(double, double, double), double, double, double a, double b, int n);
double midpnt(double (*func)(double, double, double), double, double, double a, double b, int n);
double qtrap(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b);
double qtrap_modif(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b);
double qsimp(double (*func) (double, double, double), double, double, double a, double b);
double qsimp_modif(double (*func) (double time, double sigma, double t), double time, double sigma, double a, double b);
double qromb(double (*func)(double, double, double), double, double, double a, double b);
double qromo(double (*func)(double, double, double), double, double, double a, double b, double (*choose)(double(*)(double, double, double), double, double, double, double, int));
 #endif