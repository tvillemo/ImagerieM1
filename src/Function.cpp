/*
 * Function.cpp
 *
 *  Created on: 11 févr. 2014
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
	double **OB=new double *[nbLg];
	double **OG=new double *[nbLg];
	double **OR=new double *[nbLg];
	for (int i=0;i<nbLg;i++){
		OB[i]=new double[nbCol];
		OG[i]=new double[nbCol];
		OR[i]=new double[nbCol];
	}

	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			OB[j][i]=(imgR[j][i]-imgG[j][i])/sqrt(2.0); //R-G/sqrt(2)
			OG[j][i]=(imgR[j][i]+imgG[j][i]-2.0*imgB[j][i])/sqrt(6.0);
			OR[j][i]=(imgR[j][i]+imgG[j][i]+imgB[j][i])/sqrt(3.0);
		}
	}
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			imgB[j][i]=OB[j][i];
			imgG[j][i]=OG[j][i];
			imgR[j][i]=OR[j][i];
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
	double mux=0.0;
	double muy=0.0;
	double sigmax=0.0;
	double sigmay=0.0;
	double* Mx=new double [NBGRIS];
	double* My=new double [NBGRIS];
	for (int i=0; i<NBGRIS;i++){
		Mx[i]=0.0;
		My[i]=0.0;
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			primitive.Energie+=pow(mat[j][i],2);
			if (mat[j][i]!=0.0)
				primitive.Entropie-=mat[j][i]*log10(mat[j][i]);
			primitive.homoLocal+=mat[j][i]/(1+(pow((double)(j-i),2.0)));
			Mx[j]+=mat[j][i];
			My[i]+=mat[j][i];
		}
		primitive.Uniformite+=pow(mat[i][i],2.0f);
	}
	for (int i=0;i<NBGRIS;i++){
		mux+=(float)i*Mx[i];
		muy+=(float)i*My[i];
	}
	for (int i=0;i<NBGRIS;i++){
		sigmax+=pow((float)i-mux,2.0)*Mx[i];
		sigmay+=pow((float)i-muy,2.0)*My[i];
	}
	for (int i=0;i<NBGRIS;i++){
		for (int j=0;j<NBGRIS;j++){
			primitive.corelation+=((i-mux)*(j-muy)*(mat[j][i]));
		}
	}
	primitive.corelation/=(sqrt(sigmax*sigmay));
	return primitive;
}

void createTabHaralick(double **imgHue,int nbLg, int nbCol, string name){
	unsigned int** cooc0 = cooccurrence(imgHue,nbLg,nbCol,0);
	unsigned int** cooc45 = cooccurrence(imgHue,nbLg,nbCol,45);
	unsigned int** cooc90 = cooccurrence(imgHue,nbLg,nbCol,90);
	unsigned int** cooc135 = cooccurrence(imgHue,nbLg,nbCol,135);
	float **NormMat0=normalize(cooc0);
	float **NormMat45=normalize(cooc45);
	float **NormMat90=normalize(cooc90);
	float **NormMat135=normalize(cooc135);
	Haralick primitive0 = primitiveCoo(NormMat0);
	Haralick primitive45 = primitiveCoo(NormMat45);
	Haralick primitive90 = primitiveCoo(NormMat90);
	Haralick primitive135 = primitiveCoo(NormMat135);
	writeCoo(cooc0,"coo0.txt");
	writeCoo(cooc45,"coo45.txt");
	writeCoo(cooc90,"coo90.txt");
	writeCoo(cooc135,"coo135.txt");
	ecrirePrimitive(primitive0,primitive45,primitive90,primitive135,name);
}

void ecrirePrimitive(Haralick primitive0,Haralick primitive45,Haralick primitive90,Haralick primitive135,string name){
	String filename="Primitive.csv";
	fstream fileHandle;
	fileHandle.open(filename.c_str(), fstream::out|fstream::app | ios::binary);
	fileHandle<<name<<"\n";
	fileHandle<<";"<<"Entropie"<<";"<<"Energie"<<";"<<"Uniformite"<<";"<<"HomoLocal"<<";"<<"correlation"<<"\n";
	fileHandle<<"0;"<<primitive0.Entropie<<";"<<primitive0.Energie<<";"<<primitive0.Uniformite<<";"<<primitive0.homoLocal<<";"<<primitive0.corelation<<"\n";
	fileHandle<<"45;"<<primitive45.Entropie<<";"<<primitive45.Energie<<";"<<primitive45.Uniformite<<";"<<primitive45.homoLocal<<";"<<primitive45.corelation<<"\n";
	fileHandle<<"90;"<<primitive90.Entropie<<";"<<primitive90.Energie<<";"<<primitive90.Uniformite<<";"<<primitive90.homoLocal<<";"<<primitive90.corelation<<"\n";
	fileHandle<<"135;"<<primitive135.Entropie<<";"<<primitive135.Energie<<";"<<primitive135.Uniformite<<";"<<primitive135.homoLocal<<";"<<primitive135.corelation<<"\n";
}

void img2NOp(double** NOp1,double** NOp2,double** imgB,double** imgG,double** imgR,int nbLg, int nbCol){
	img2Op(imgR,imgG,imgB,nbLg,nbCol);
	for (int i=0;i<nbCol;i++){
		for (int j=0;j<nbLg;j++){
			NOp1[j][i] = imgB[j][i]==0 ? 0 : imgR[j][i]/imgB[j][i];
			NOp1[j][i] = imgB[j][i]==0 ? 0 : imgG[j][i]/imgB[j][i];
		}
	}
}

double*** calcLBP(double** imgB,double** imgG,double** imgR,int nbLg, int nbCol,string type){
	double*** imgOut;
	double** voisins=new double*[3];
	for (int i =0;i<3;i++){
		voisins[i]=new double[3];
	}
	if (type=="nOpp"){
		imgOut=new double**[2];
		for (int i=0;i<2;i++){
			imgOut[i]=new double*[nbLg];
			for (int j=0;j<nbLg;j++){
				imgOut[i][j]=new double [nbCol];
			}
		}
		double** NOp1=new double*[nbLg];
		double** NOp2=new double*[nbLg];
		for (int i=0;i<nbLg;i++){
			NOp1[i]=new double[nbCol];
			NOp2[i]=new double[nbCol];
		}
		img2NOp(NOp1,NOp2,imgB,imgG,imgR,nbLg,nbCol);
		double S;
		for (int i=1;i<nbCol-1;i++){
			for (int j=1;j<nbLg-1;j++){
				retVoisins(i,j,voisins,NOp1);
				double resInt=0.0;
				for (int x=0;x<2;x++){
					for (int y=0;y<2;y++){
						S=(voisins[x][y]-voisins[1][1])==0 ? 0 : 1;
						resInt+=S*pow(2,7.0);
					}
				}
				imgOut[0][j][i]=resInt;

				retVoisins(i,j,voisins,NOp2);
				double resInt=0.0;
				for (int x=0;x<2;x++){
					for (int y=0;y<2;y++){
						S=voisins[x][y]-voisins[1][1]==0 ? 0 : 1;
						resInt+=S*pow(2,7.0);
					}
				}
				imgOut[1][j][i]=resInt;
			}
		}
	}
	else if(type=="Opp"){
		imgOut=new double**[3];
		for (int i=0; i<3; i++){
			imgOut[i]=new double*[nbLg];
			for(int j=0; j<nbLg; j++){
				imgOut[i][j]=new double[nbCol];
			}
		}
		double** Op1=new double*[nbLg];
		double** Op2=new double*[nbLg];
		double** Op3=new double*[nbLg];
		for (int i=0; i<nbLg; i++){
			Op1[i]= new double[nbCol];
			Op2[i]= new double[nbCol];
			Op3[i]= new double[nbCol];
		}
		double S;
		img2Op(imgB,imgG,imgR,nbLg,nbCol);
		for (int i=1; i<nbCol-1; i++){
			for(int j=1; j<nbLg-1; j++){
				retVoisins(i,j,voisins,Op1);
				double resInt=0.0;
				for (int x=0;x<2;x++){
					for (int y=0;y<2;y++){
						S=voisins[x][y]-voisins[1][1]==0 ? 0 : 1;
						resInt+=S*pow(2,7.0);
					}
				}
				imgOut[0][j][i]=resInt;

				retVoisins(i,j,voisins,Op2);
				double resInt=0.0;
				for (int x=0;x<2;x++){
					for (int y=0;y<2;y++){
						S=voisins[x][y]-voisins[1][1]==0 ? 0 : 1;
						resInt+=S*pow(2,7.0);
					}
				}
				imgOut[1][j][i]=resInt;

				retVoisins(i,j,voisins,Op3);
				double resInt=0.0;
				for (int x=0;x<2;x++){
					for (int y=0;y<2;y++){
						S=voisins[x][y]-voisins[1][1]==0 ? 0 : 1;
						resInt+=S*pow(2,7.0);
					}
				}
				imgOut[2][j][i]=resInt;
			}
		}


	}
	else{
		imgOut=new double**[1];
		double** imgHue=new double*[nbLg];
		imgOut[0]=new double*[nbLg];
		for (int j=0;j<nbLg;j++){
			imgOut[0][j]=new double [nbCol];
			imgHue[j]=new double[nbCol];
		}
		double S;
		imgHue=img2Hue(imgB,imgG,imgR,nbLg,nbCol);
		for (int i=1;i<nbCol-1;i++){
			for (int j=1;j<nbLg-1;j++){
				retVoisins(i,j,voisins,imgHue);
				double resInt=0.0;
				for (int x=0;x<2;x++){
					for (int y=0;y<2;y++){
						S=voisins[x][y]-voisins[1][1]==0 ? 0 : 1;
						resInt+=S*pow(2,7.0);
					}
				}
				imgOut[0][j][i]=resInt;
			}
		}
	}
	return imgOut;
}

void retVoisins(int i,int j,double** voisins,double** img){
	voisins[0][0]=img[j-1][i-1];
	voisins[1][0]=img[j][i-1];
	voisins[2][0]=img[j+1][i-1];
	voisins[0][1]=img[j-1][i];
	voisins[1][1]=img[j][i];
	voisins[2][1]=img[j+1][i];
	voisins[0][2]=img[j-1][i+1];
	voisins[1][2]=img[j][i+1];
	voisins[2][2]=img[j+1][i+1];

}

