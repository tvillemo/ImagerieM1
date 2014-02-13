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

#endif /* FUNCTION_H_ */
