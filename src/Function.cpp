/*
 * Function.cpp
 *
 *  Created on: 11 f�vr. 2014
 *      Author: Tanguy
 */

#include "Function.h"
using namespace std;
using namespace cv;


unsigned int** cooccurrence(double **imgIn,int nbLg,int nbCol,int direction){

	unsigned int** cooccurrence=new unsigned int* [NBGRIS];
	for (int i=0;i<NBGRIS;i++){
		cooccurrence[i]=new unsigned int[NBGRIS];
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			cooccurrence[i][j]=0;
		}
	}
	int ** img=new int*[nbLg];
	for (int i=0;i<nbLg;i++){
		img[i]=new int[nbCol];
	}
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			img[j][i]=(int)round(imgIn[j][i]);
		}
	}
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			if (direction==0){
				if(i!=nbCol-1){
					cooccurrence[img[j][i]][img[j][i+1]]++;
				}
				if(i!=0){
					cooccurrence[img[j][i]][img[j][i-1]]++;
				}
			}
			else if(direction==45){
				if(i!=nbCol-1 && j!=nbLg-1){
					cooccurrence[img[j][i]][img[j+1][i+1]]++;
				}
				if(i!=0&&j!=0){
					cooccurrence[img[j][i]][img[j-1][i-1]]++;
				}
			}
			else if(direction==90){
				if(j!=nbLg-1){
					cooccurrence[img[j][i]][img[j+1][i]]++;
				}
				if(j!=0){
					cooccurrence[img[j][i]][img[j-1][i]]++;
				}
			}
			else if(direction==135){
				if(i!=nbCol-1&&j!=0){
					cooccurrence[img[j][i]][img[j-1][i+1]]++;
				}
				if(i!=0&&j!=nbLg-1){
					cooccurrence[img[j][i]][img[j+1][i-1]]++;
				}
			}
		}
	}
	return cooccurrence;
}

void img2Op(double** imgB,double** imgG,double** imgR,int nbLg,int nbCol){

	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			imgB[j][i]=(imgR[j][i]-imgG[j][i])/sqrt(2.0); //R-G/sqrt(2)
			imgG[j][i]=(imgR[j][i]+imgG[j][i]-2.0*imgB[j][i])/sqrt(6.0);
			imgR[j][i]=(imgR[j][i]+imgG[j][i]+imgB[j][i])/sqrt(3.0);
		}
	}
	float* min=seekMin(imgB,imgG,imgR,nbLg,nbCol);
	float* max=seekMax(imgB,imgG,imgR,nbLg,nbCol);
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			imgB[j][i]=(((imgB[j][i]-min[0]))/(max[0]-min[0]))*255.0;
			imgG[j][i]=(((imgG[j][i]-min[1]))/(max[1]-min[1]))*255.0;
			imgR[j][i]=(((imgR[j][i]-min[2]))/(max[2]-min[2]))*255.0;
		}
	}
}

//fonction d'enregistrement d'une image en format JPEG
void saveImage(string fileName,IplImage* imgBGR)
{
	if(!cvSaveImage(const_cast<char*>(fileName.c_str()), imgBGR)){
		cout<<"Could not save: "<<fileName<<endl;
	}
}

float* seekMax(double** imgB,double** imgG,double** imgR,int nbLg,int nbCol){
	float* max= new float[3];
	for (int i=0;i<3;i++){
		max[i]=0.0;
	}
	for (int i=0;i<nbLg;i++){
		for (int j=0;j<nbCol;j++){

			if(max[0]<imgB[i][j]){
				max[0]=imgB[i][j];
			}
			if(max[1]<imgG[i][j]){
				max[1]=imgG[i][j];
			}
			if(max[2]<imgR[i][j]){
				max[2]=imgR[i][j];
			}
		}
	}
	return max;
}

float* seekMin(double** imgB,double** imgG,double** imgR,int nbLg,int nbCol){
	float* min= new float[3];
	for (int i=0;i<3;i++){
		min[i]=numeric_limits<float>::max();
	}
	for (int i=0;i<nbLg;i++){
		for (int j=0;j<nbCol;j++){

			if(min[0]>imgB[i][j]){
				min[0]=imgB[i][j];
			}
			if(min[1]>imgG[i][j]){
				min[1]=imgG[i][j];
			}
			if(min[2]>imgR[i][j]){
				min[2]=imgR[i][j];
			}
		}
	}
	return min;
}

double** planCouleur(IplImage* img, int numPlan)
{
	int nbCol=cvGetSize(img).width; //largeur (nombre de colonnes)
	int nbLg=cvGetSize(img).height;// hauteur (nombre de lignes)

	double** imgCoul=new double*[nbLg];
	for(int i=0;i<nbLg;i++)
		imgCoul[i] = new double[nbCol];


	for(int i=0;i<nbLg;i++)
	{
		for(int j=0;j<nbCol;j++)
		{
			CvScalar color=cvGet2D(img,i,j);
			imgCoul[i][j]=(int)color.val[numPlan-1]; //je fais un cast pour passer du double en int
		}
	}

	return imgCoul;
}

IplImage* reconstruireImageCouleur(double** plan1, double** plan2, double** plan3,int nbLg, int nbCol)
{
	IplImage* img=cvCreateImage(cvSize(nbCol,nbLg),IPL_DEPTH_8U,3);

	for(int i=0;i<nbLg;i++)
	{
		for(int j=0;j<nbCol;j++)
		{
			CvScalar color;
			color.val[0]=plan1[i][j];
			color.val[1]=plan2[i][j];
			color.val[2]=plan3[i][j];
			color.val[3]=0;
			cvSet2D(img,i,j,color);
		}
	}

	return img;
}

double** img2Hue(double** imgB,double** imgG,double** imgR,int nbLg, int nbCol){
	img2Op(imgB,imgG,imgR,nbLg,nbCol);
	double** imgC=new double*[nbLg];
	for (int i=0;i<nbLg;i++){
		imgC[i]=new double[nbCol];
	}
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			if (imgG[j][i]==0.0){
				imgG[j][i]=M_PI/2.0;
			}
			imgC[j][i]=atan(imgB[j][i]/imgG[j][i]);

		}
	}
	float max= Max1Plan(imgC, nbLg, nbCol);
	float min= Min1Plan(imgC, nbLg, nbCol);
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			imgC[j][i]=(((imgC[j][i]-min))/(max-min))*255.0;
		}
	}
	return imgC;
}

IplImage* reconstruireImage8Bits(double** img,int nbLg, int nbCol)
{
	IplImage* img8Bits=cvCreateImage(cvSize(nbCol,nbLg),IPL_DEPTH_8U,1);

	for(int i=0;i<nbLg;i++)
	{
		for(int j=0;j<nbCol;j++)
		{
			CvScalar color;
			color.val[0]=img[i][j];
			color.val[1]=0;
			color.val[2]=0;
			color.val[3]=0;
			cvSet2D(img8Bits,i,j,color);
		}
	}

	return img8Bits;
}


float Max1Plan(double** img,int nbLg,int nbCol){
	float max=0.0;
	for (int i=0;i<nbLg;i++){
		for (int j=0;j<nbCol;j++){

			if(max<img[i][j]){
				max=img[i][j];
			}
		}
	}
	return max;
}

float Min1Plan(double** img,int nbLg,int nbCol){
	float min;
	min=numeric_limits<float>::max();
	for (int i=0;i<nbLg;i++){
		for (int j=0;j<nbCol;j++){

			if(min>img[i][j]){
				min=img[i][j];
			}
		}
	}
	return min;
}

void writeCoo(unsigned int** mat,string filename){
	fstream fileHandle;
	fileHandle.open(filename.c_str(), fstream::out | ios::binary);
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			fileHandle<<mat[j][i]<<" ";
		}
		fileHandle<<"\n";
	}
}

float** normalize(unsigned int** mat){
	float n=0;
	float** matOut=new float*[NBGRIS];
	for (int i=0;i<NBGRIS;i++){
		matOut[i]=new float[NBGRIS];
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			n+=(float)mat[j][i];
		}
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			matOut[j][i]=(float)mat[j][i]/n;
		}
	}
	return matOut;
}

Haralick primitiveCoo(float** mat){
	Haralick primitive;
	primitive.Energie=0.0;
	primitive.Entropie=0.0;
	primitive.Uniformite=0.0;
	primitive.corelation=0.0;
	primitive.homoLocal=0.0;
	double* mux=new double[NBGRIS];
	double* muy=new double[NBGRIS];
	double* sigmax=new double[NBGRIS];
	double* sigmay=new double[NBGRIS];
	double* Mx=new double [NBGRIS];
	double* My=new double [NBGRIS];
	for (int i=0; i<NBGRIS;i++){
		mux[i]=0.0;
		muy[i]=0.0;
		sigmax[i]=0.0;
		sigmay[i]=0.0;
		Mx[i]=0.0;
		My[i]=0.0;
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			primitive.Energie+=pow(mat[j][i],2);
			if (mat[j][i]!=0.0)
				primitive.Entropie-=mat[j][i]*log10(mat[j][i]);
			primitive.homoLocal+=mat[j][i]/(1+(pow((double)(j-i),2)));
			Mx[j]+=mat[j][i];
			My[i]+=mat[j][i];
		}
		primitive.Uniformite+=pow(mat[i][i],2);

	}
	for (int i=0;i<NBGRIS;i++){
		mux[i]+=i*Mx[i];
		muy[i]+=i*My[i];
		sigmax[i]+=pow(i-mux[i],2)*Mx[i];
		sigmay[i]+=pow(i-muy[i],2)*My[i];
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			if(sigmax[j]!=0.0&&sigmay[i]!=0.0){
				primitive.corelation+=((j-mux[j])*(i-muy[i])*(mat[j][i]))/(sqrt(sigmax[j])*sqrt(sigmay[i]));
			}
		}
	}
	return primitive;
}
