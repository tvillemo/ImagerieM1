/*
 * Function.h
 *
 *  Created on: 11 févr. 2014
 *      Author: Tanguy
 */

#ifndef FUNCTION_H_
#define FUNCTION_H_


#include "opencv/cv.h"
#include "opencv/highgui.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>

#define NBGRIS 256

typedef struct{
	float Energie;
	float Uniformite;
	float Entropie;
	float homoLocal;
	float corelation;
}Haralick;

unsigned int** cooccurrence(double **imgIn,int nbLg,int nbCol,int direction);
void img2Op(double** imgB,double** imgG,double** imgR,int nbLg, int nbCol);
double** img2Hue(double** imgB,double** imgG,double** imgR,int nbLg, int nbCol);
void saveImage(string fileName,IplImage* imgBGR);
float* seekMax(double** imgB,double** imgG,double** imgR,int nbLg,int nbCol);
float* seekMin(double** imgB,double** imgG,double** imgR,int nbLg,int nbCol);
double** planCouleur(IplImage* img, int numPlan);
IplImage* reconstruireImage8Bits(double** img,int nbLg, int nbCol);
IplImage* reconstruireImageCouleur(double** plan1, double** plan2, double** plan3,int nbLg, int nbCol);
float Max1Plan(double** img,int nbLg,int nbCol);
float Min1Plan(double** img,int nbLg,int nbCol);
float** normalize(unsigned int** mat);
void writeCoo(unsigned int** mat,string filename);
Haralick primitiveCoo(float** mat);
void createTabHaralick(double **imgHue,int nbLg, int nbCol, std::string name);
void ecrirePrimitive(Haralick primitive0,Haralick primitive45,Haralick primitive90,Haralick primitive135, std::string name);
void img2NOp(double** NOp1,double** NOp2,double** imgR,double** imgG,double** imgB,int nbLg, int nbCol);
double*** calcLBP(double** imgR,double** imgG,double** imgB,int nbLg, int nbCol);
void retVoisins(int i,int j,double** voisins,double** img);

#endif /* FUNCTION_H_ */
