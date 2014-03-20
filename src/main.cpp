/*
 * main.cpp
 *
 *  Created on: 11 f�vr. 2014
 *      Author: Tanguy
 */

#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include "Function.h"

using namespace cv;

int main(){
	string filename="Images/Leaves_3.jpg";
	IplImage* img=cvLoadImage(const_cast<char*>(filename.c_str()));
	if (!img)
		cout<<"erreur"<<endl;

	double** imgB = planCouleur(img,1);
	double** imgG = planCouleur(img,2);
	double** imgR = planCouleur(img,3);
	int nbCol=cvGetSize(img).width; //largeur (nombre de colonnes)
	int nbLg=cvGetSize(img).height;// hauteur (nombre de lignes)
	double** imgHue = img2Hue(imgB, imgG, imgR,nbLg,nbCol);
	//IplImage* imgColHue = reconstruireImage8Bits(imgHue,nbLg,nbCol);
	String out_img_name="imgLBPHue.jpg";
	/*saveImage(out_img_name,imgColHue);
	createTabHaralick(imgHue,nbLg,nbCol,filename);*/
	calcLBP(imgB,imgG,imgR,nbLg,nbCol,"Opp");
	system("pause");
	return 0;
}

