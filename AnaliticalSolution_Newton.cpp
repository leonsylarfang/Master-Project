#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
//#include <GL/glut.h>
#include <Eigen/Dense>  

FILE *in0,*in1,*out;
using namespace std; 

int ciggj(double **,int,double *);
void init();

const double pi=3.1415926;

//const int Nt=5;
const int Nt = 12;

#define Nunknow 1+2*Nt
#define Nunknow1 2+4*Nt
#define Nframe 30;

void init(void) 
{
#if 1
	/*in0 = fopen("D:/BU Lab/Publication Plan/third paper 2014/Project1/Scan_curves/Scan_04_60.dat", "rt");
	in1 = fopen("D:/BU Lab/Publication Plan/third paper 2014/Project1/Scan_curves/Scan_04_60.dat", "rt");*/
//	double frame = 0.5;
	double frame = 0.2;
//	in0 = fopen("in_laugh_cheek&mouth_without_nose.dat", "rt");
//	in1 = fopen("in_reference_cheek&mouth_without_nose.dat", "rt");
//	in0 = fopen("in_laugh_outer_face_line&inner_face_line.dat", "rt");
//	in1 = fopen("in_reference_outer_face_line&inner_face_line.dat", "rt");
//	in0 = fopen("in_laugh.txt", "rt");
//	in1 = fopen("in_reference.txt", "rt");
	in0 = fopen("in_laugh_outer_face_line&inner_face_line_1.dat", "rt");
	in1 = fopen("in_reference_outer_face_line&inner_face_line_1.dat", "rt");
    out = fopen("out.dat","wt");
	fstream OriginalC("origCurve.mel", ios::out);
	fstream CalculatC("keyframes_nt_outer_face_line&inner_face_line_1.mel", ios::out);
	fstream  UnifyC("UnifyC.mel", ios::out);
	fstream RotateC("RotateC.mel", ios::out);
	fstream FC("FC.mel", ios::out);
#endif
#if 0
	in0 = fopen("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/AtoB&BtoA/Curve01.dat", "rt");
	in1 = fopen("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/AtoB&BtoA/Curve06.dat", "rt");
	out = fopen("out.dat","wt");
	fstream OriginalC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Define_parameters/origCurve.mel", ios::out);
	fstream CalculatC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Define_parameters/CalcCurve.mel", ios::out);
	fstream UnifyC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Define_parameters/UnifyC.mel", ios::out);
	fstream RotateC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Define_parameters/RotateC.mel", ios::out);
	fstream FC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Define_parameters/FC.mel", ios::out);
#endif
#if 0 //elephant
	in0 = fopen("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/leg.txt", "rt");
	in1 = fopen("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/leg.txt", "rt");
	out = fopen("out.dat", "wt");
	fstream OriginalC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/origCurve.mel", ios::out);
	fstream CalculatC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/CalcCurve.mel", ios::out);
	fstream UnifyC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/UnifyC.mel", ios::out);
	fstream RotateC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/RotateC.mel", ios::out);
	fstream FC("D:/BU Lab/Publication Plan/SIGGRAPH_Latex/Review_comments_04-07-2016/hybridskindeformation/acm-tog/Comparison_figure materials/Elephant/data/FC.mel", ios::out);
#endif
	int i,ii,j,k,kk,iq,n,curve,curve0,icurve,i0,k0;
	int nn = Nframe; 
	fscanf(in0, "%d %d\n", &curve, &k0);
	printf("curve,k0= %d %d\n", curve, k0);
//	exit(1);
	icurve = curve;
	const int It = k0;

	double dv,rio,An,a1q,a2q,vi,t,c0[3],c0Bar[3],dt,f;
	double **cnBar1 = new double *[Nt];
	double **cnBar2 = new double *[Nt];
	double **dnBar1 = new double *[Nt];
	double **dnBar2 = new double *[Nt];
	for(i=0;i<Nt;i++)
	{
		cnBar1[i]=new double[3];
		cnBar2[i]=new double[3];
		dnBar1[i]=new double[3];
		dnBar2[i]=new double[3];
	}
	double *bb = new double [Nunknow];
	double *bbx = new double [Nunknow];
	double **aa=new double *[Nunknow];
	double **aax=new double *[Nunknow];
	for(i=0;i<Nunknow;i++)
	{
		aa[i]=new double[Nunknow];
		aax[i]=new double[Nunknow];
	}
	double *bbn = new double [Nunknow1];
	double *bbnx = new double [Nunknow1];
	double **aan=new double *[Nunknow1];
	double **aanx=new double *[Nunknow1];
	for(i=0;i<Nunknow1;i++)
	{
		aan[i]=new double[Nunknow1];
		aanx[i]=new double[Nunknow1];
	}
	double **qi0=new double *[It];
	double **qi1=new double *[It];
	double **qt=new double *[It];
	for(i=0;i<It;i++)
	{
		qi0[i]=new double[3];
		qi1[i]=new double[3];
		qt[i]=new double[3];
	}

  
   double*** inData0 = new double**[curve]; // 
   double*** inData0L = new double** [curve]; // 
   for(i=0;i<curve;i++)
   {
      inData0[i] = new double*[It];
	  inData0L[i] = new double* [It];
	  for (short j = 0; j < It; j++)
	  {
		  inData0[i][j] = new double[3]; // 
		  inData0L[i][j] = new double[3]; // 
	  }
   }
   double*** inData1 = new double**[curve]; // 
   double*** inData1L = new double** [curve]; // 
   for(i=0;i<curve;i++)
   {
      inData1[i] = new double*[It];
	  inData1L[i] = new double* [It];
	  for (short j = 0; j < It; j++)
	  {
		  inData1[i][j] = new double[3]; // 
		  inData1L[i][j] = new double[3]; // 
	  }
   }
   double*** inData8 = new double**[curve]; // 
   for(i=0;i<curve;i++)
   {
      inData8[i] = new double*[It];
      for(short j=0;j<It;j++)
         inData8[i][j] = new double[3]; // 
   }

   double*** iniuni = new double**[curve]; // 
   for(i=0;i<curve;i++)
   {
      iniuni[i] = new double*[It];
      for(short j=0;j<It;j++)
         iniuni[i][j] = new double[3]; // 
   }
   double*** inruni = new double**[curve]; //  
   for(i=0;i<curve;i++)
   {
      inruni[i] = new double*[It];
      for(short j=0;j<It;j++)
         inruni[i][j] = new double[3]; // 
   }
   double*** induni = new double**[curve]; // 
   for(i=0;i<curve;i++)
   {
      induni[i] = new double*[It];
      for(short j=0;j<It;j++)
         induni[i][j] = new double[3]; // 
   }

   double*** incuni = new double**[curve]; // 
   for (i = 0; i<curve; i++)
   {
	   incuni[i] = new double*[It];
	   for (short j = 0; j<It; j++)
		   incuni[i][j] = new double[3]; // 
   }
   //double*** incuni = new double**[It]; // 
   //for(i=0;i<It;i++)
   //{
   //   incuni[i] = new double*[curve];
   //   for(short j=0;j<curve;j++)
   //      incuni[i][j] = new double[3]; // 
   //}
     
   double ****AnimData = new double ***[nn]; // 
   for(i=0;i<nn;i++)
   {
      AnimData[i] = new double **[curve];
      for(short j=0;j<curve;j++)
	  {
         AnimData[i][j] = new double *[It]; 
		 for(short k=0;k<It;k++)
	         AnimData[i][j][k] = new double[3]; 
	  }
   }

   double*** qn = new double**[curve]; // 
   for(i=0;i<curve;i++)
   {
      qn[i] = new double*[Nunknow];
      for(short j=0;j<Nunknow;j++)
         qn[i][j] = new double[3]; // 
   }
   double*** rn = new double**[curve]; // 
   for(i=0;i<curve;i++)
   {
      rn[i] = new double*[Nunknow];
      for(short j=0;j<Nunknow;j++)
         rn[i][j] = new double[3]; // 
   }
   double*** dn = new double**[curve]; // 
   for(i=0;i<curve;i++)
   {
      dn[i] = new double*[Nunknow];
      for(short j=0;j<Nunknow;j++)
         dn[i][j] = new double[3]; // 
   }

//inpput the data
	for(i=0;i<curve;i++)
   {
		for(j=0;j<k0;j++)
		{
			fscanf(in0,"%le %le %le\n",&inData0[i][j][0],&inData0[i][j][1],&inData0[i][j][2]);
//			printf("j,inData0= %d %le %le %le\n", j, inData0[i][j][0], inData0[i][j][1], inData0[i][j][2]);
		}
   }
//	exit(1);
    fscanf(in1, "%d %d\n", &curve0, &k0);
	for(i=0;i<curve;i++)
   {
		for(int j=0;j<k0;j++)
		{
			fscanf(in1,"%le %le %le\n",&inData1[i][j][0],&inData1[i][j][1],&inData1[i][j][2]);
//			fscanf(in1, "%le %le %le\n", &inData1[i][k0-j-1][0], &inData1[i][k0 - j - 1][1], &inData1[i][k0 - j - 1][2]);
		}
   }
//	exit(1);
//Generate Mel file of input curves
	 for(i=0;i<curve;i++)
		 //for (i = 0; i<1; i++)
	 {
         OriginalC << "curve -d 3\n";         
         for(j=0;j<k0;j++)
         {
             OriginalC << "-p " << inData0[i][j][0] << " " << inData0[i][j][1] << " " << inData0[i][j][2] << "\n";
			 inData0L[i][j][0] = inData0[i][j][0];
			 inData0L[i][j][1] = inData0[i][j][1];
			 inData0L[i][j][2] = inData0[i][j][2];
		 }
		 OriginalC << ";\n"; 
	 }
/*
	 for (i = 0; i < curve; i++)
	 {
		 for (j = 0; j < k0; j++)
		 {
			 fprintf(out, "i,j,in0= %d %d %15.12f %15.12f %15.12f\n", 
				 i,j, inData0[i][j][0], inData0[i][j][1], inData0[i][j][2]);
			 fprintf(out, "i,j,in0= %d %d %15.12f %15.12f %15.12f\n", 
				 i, j, inData0L[i][j][0], inData0L[i][j][1], inData0L[i][j][2]);
		 }
	 }
*/
//printf("first curve,It= %d %d\n",curve,It);
	 for(i=0;i<curve;i++)
		// //for (i = 0; i<1; i++)
	 {
         OriginalC << "curve -d 3\n";         
         for(j=0;j<k0;j++)
         {
             OriginalC << "-p " << inData1[i][j][0] << " " << inData1[i][j][1] << " " << inData1[i][j][2] << "\n";
			 inData1L[i][j][0] = inData1[i][j][0];
			 inData1L[i][j][1] = inData1[i][j][1];
			 inData1L[i][j][2] = inData1[i][j][2];
		 }
		 OriginalC << ";\n"; 
	 }

	 for (i = 0; i < curve; i++)
	 {
		 for (j = 0; j < k0; j++)
		 {
			 fprintf(out, "i,j,in1= %d %d %15.12f %15.12f %15.12f\n",
				 i, j, inData1[i][j][0], inData1[i][j][1], inData1[i][j][2]);
			 fprintf(out, "i,j,in1= %d %d %15.12f %15.12f %15.12f\n",
				 i, j, inData1L[i][j][0], inData1L[i][j][1], inData1L[i][j][2]);
		 }
	 }

	 double ti, dti, x, y, z;
//	dv=1.0/100;						
	dti = 1.0 / nn;
	for (k = 0; k <= nn; k++)
	{
		ti = k * dti;
//  ******************  linear interpolation  ************************//
/*
		printf("ti= %le\n", ti);
		for (i = 0; i < curve; i++)
		{
			CalculatC << "curve -d 3\n";
			for (j = 0; j < k0; j++)
			{
				x = inData1L[i][j][0] + ti * (inData0L[i][j][0] - inData1L[i][j][0]);
				y = inData1L[i][j][1] + ti * (inData0L[i][j][1] - inData1L[i][j][1]);
				z = inData1L[i][j][2] + ti * (inData0L[i][j][2] - inData1L[i][j][2]);
				CalculatC << "-p " << x << " " << y << " " << z << "\n";
				fprintf(out,"j,ti,in1-0-x= %d %4.2f %15.12f %15.12f %15.12f\n", j, ti, inData1L[i][j][0], inData0L[i][j][0],x);
			}
			CalculatC << ";\n";
		}
	}
*/
//  ******************  nonlinear interpolation from Newton's second law ************************//

		printf("ti= %le\n", ti);
		for (i = 0; i < curve; i++)
		{
			CalculatC << "curve -d 3\n";
			for (j = 0; j < k0; j++)
			{
				x = inData1L[i][j][0] + ti * ti * (inData0L[i][j][0] - inData1L[i][j][0]);
				y = inData1L[i][j][1] + ti * ti * (inData0L[i][j][1] - inData1L[i][j][1]);
				z = inData1L[i][j][2] + ti * ti * (inData0L[i][j][2] - inData1L[i][j][2]);
				CalculatC << "-p " << x << " " << y << " " << z << "\n";
				fprintf(out,"j,ti,in1-0-x= %d %4.2f %15.12f %15.12f %15.12f\n", j, ti, inData1L[i][j][0], inData0L[i][j][0],x);
			}
			CalculatC << ";\n";
		}
	}	
	
	//fstream C_P4("C_P4.txt", ios::out);
	//for (j = 0; j <= It - 1; j++)
	//{
	//	for (ii = 0; ii<curve; ii++)
	//	{
	//		C_P4 << incuni[ii][j][0] << " " << incuni[ii][j][1] << " " << incuni[ii][j][2] << "\n";
	//	}
	//}
//getchar();
#if 0
	FILE *C_p42 = fopen("C_p42.txt", "rt");
	for (j = 0; j<k0; j++)
	{
		for (i = 0; i<icurve; i++)
		{
			fscanf(C_p42, "%le %le %le\n", &inruni[i][j][0], &inruni[i][j][1], &inruni[i][j][2]);
		}
	}
	fstream C_p42_mel("C_p42.mel", ios::out);
	for (i = 0; i<icurve; i++)
	{
		C_p42_mel << "curve -d 3\n";
		for (j = 1; j <= It - 1; j++)
		{
			C_p42_mel << "-p " << inruni[i][j][0] << " " << inruni[i][j][1] << " " << inruni[i][j][2] << "\n";
		}
		C_p42_mel << ";\n";
	}
#endif
#if 0
	FILE *C_p43 = fopen("C_p43.txt", "rt");
	for (j = 0; j<k0; j++)
	{
		for (i = 0; i<icurve; i++)
		{
			fscanf(C_p43, "%le %le %le\n", &inruni[i][j][0], &inruni[i][j][1], &inruni[i][j][2]);
		}
	}
	fstream C_p43_mel("C_p43.mel", ios::out);
	for (i = 0; i<icurve; i++)
	{
		C_p43_mel << "curve -d 3\n";
		for (j = 1; j <= It - 1; j++)
		{
			C_p43_mel << "-p " << inruni[i][j][0] << " " << inruni[i][j][1] << " " << inruni[i][j][2] << "\n";
		}
		C_p43_mel << ";\n";
	}
#endif
exit(1);
	delete [] bb;
	delete [] bbx;
	delete [] bbn;
	delete [] bbnx;
	for(i=0; i<Nunknow; i++)
    {
		delete [] aa[i];
		delete [] aax[i];
	}	
	delete [] aa;
	delete [] aax;
	for(i=0; i<Nunknow1; i++)
    {
		delete [] aan[i];
		delete [] aanx[i];
	}	
	delete [] aan;
	delete [] aanx;
	for(i=0; i<Nt; i++)
    {
		delete [] cnBar1[i];
		delete [] cnBar2[i];
		delete [] dnBar1[i];
		delete [] dnBar2[i];
		delete [] qi0[i];
		delete [] qi1[i];
		delete [] qt[i];
	}	
	delete [] cnBar1;
	delete [] cnBar2;
	delete [] dnBar1;
	delete [] dnBar2;
	delete [] qi0[i];
	delete [] qi1[i];
	delete [] qt[i];

	for(i=0;i<icurve;i++)
	{
		for(j=0;j<It;j++)
		{
			delete [] inData0[i][j];
			delete [] inData1[i][j];
			delete[] inData0L[i][j];
			delete[] inData1L[i][j];
			delete [] inData8[i][j];
			delete [] iniuni[i][j];
			delete [] inruni[i][j];
			delete [] induni[i][j];
		}
		delete [] inData0[i];
		delete [] inData1[i];
		delete[] inData0L[i];
		delete[] inData1L[i];
		delete [] inData8[i];
		delete [] iniuni[i];
		delete [] inruni[i];
		delete [] induni[i];
	}
	delete [] inData0;
	delete [] inData1;
	delete[] inData0L;
	delete[] inData1L;
	delete [] inData8;
	delete [] iniuni;
	delete [] inruni;
	delete [] induni;

 	for(i=0;i<nn;i++)
	{
		for(j=0;j<curve;j++)
		{
			for(k=0;k<It;k++)
			{
				delete [] AnimData[i][j][k];
			}
			delete [] AnimData[i][j];
		}
		delete [] AnimData[i];
	}
	delete [] AnimData;

 	for(i=0; i<curve; i++)
    {
		for(j=0;j<Nunknow;j++)
		{
			delete [] qn[i][j];
			delete [] rn[i][j];
			delete [] dn[i][j];

		}
		delete [] qn[i];
		delete [] rn[i];
		delete [] dn[i];

	}	
	delete [] qn;
	delete [] rn;
	delete [] rn;


   for(i=0;i<It;i++)
   {
      for(j=0;j<curve;j++)
	  {
		  delete [] incuni[i][j];
	  }
	  delete [] incuni[i];
   }
   delete [] incuni;

//   glEndList();
   fclose(out);
   fclose(in0);
   fclose(in1);
   printf("end\n");
}
/*
void display(void)
{

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   glPushMatrix();

   glEnable(GL_COLOR_MATERIAL);
   glEnable (GL_LIGHTING);
   glShadeModel (GL_SMOOTH);
   glColor3f(0.5, 0.0, 1.0);
   glTranslatef(-0.0, -1.0, 0.0);
   glPushMatrix();
   glScalef(1.0, 1.0, 1.0);
   glRotatef(-70.0, 1.0, 0.0, 0.0);
   glCallList(startList);
   glPopMatrix();

   glFlush();
}

void reshape (int w, int h)
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   if (w <= h)
      glOrtho(-2.5, 2.5, -2.5*(GLfloat)h/(GLfloat)w,
         2.5*(GLfloat)h/(GLfloat)w, -10.0, 10.0);
   else
      glOrtho(-2.5*(GLfloat)w/(GLfloat)h,
         2.5*(GLfloat)w/(GLfloat)h, -2.5, 2.5, -10.0, 10.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
      case 27:
         exit(0);
         break;
   }
}
*/
int main(int argc, char** argv)
{
/*   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize(500, 500); 
   glutInitWindowPosition(100, 100);
   glutCreateWindow(argv[0]);
*/   init();
/*   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutKeyboardFunc(keyboard);
   glutMainLoop();
*/   return 0;
}
int ciggj(double **a,int n,double *b)      
{
	int i, j, k, is, u, v;
	int *js = new int [n];
	double d, t;
	for (k=0;k<=n-1;k++)
	{
		d=0.0;
		for (i=k;i<=n-1;i++)
		for (j=k;j<=n-1;j++)
			{t=fabs(a[i][j]);			
			if(t>d) {d=t;js[k]=j;is=i;}
			}
			if(d+1.0==1.0)
				{delete js;printf("fail\n");return(0);}
			if (is!=k)
				{for(j=k;j<=n-1;j++)
					{u=k*n+j;v=is*n+j;
					 t=a[k][j];a[k][j]=a[is][j];a[is][j]=t;
					 }
				  t=b[k];b[k]=b[is];b[is]=t;
				 }
			if (js[k]!=k)
			for (i=0;i<=n-1;i++)
				{u=i*n+k;v=i*n+js[k];
				 t=a[i][k];a[i][k]=a[i][js[k]];a[i][js[k]]=t;
				}
			t=a[k][k];
			for (j=k+1;j<=n-1;j++)
				{u=k*n+j;
				 if (a[k][j]!=0.0)a[k][j]=a[k][j]/t;
				}
			b[k]=b[k]/t;
			for (j=k+1;j<=n-1;j++)
				{u=k*n+j;
				 if(a[k][j]!=0.0)
				 	{for (i=0;i<=n-1;i++)
				 		{v=i*n+k;
				 		 if((i!=k) && (a[i][k]!=0.0))
				 		 	{is=i*n+j;
				 		 	 a[i][j]=a[i][j]-a[i][k]*a[k][j];
				 		 	}
				 		 }
				 	}
				 }
				 for(i=0;i<=n-1;i++)
				 	{u=i*n+k;
				 	 if((i!=k) && (a[i][k]!=0.0))
				 	 	b[i]=b[i]-a[i][k]*b[k];
				 	 }
			}
			for (k=n-1;k>=0;k--)
				if(k!=js[k])
					{t=b[k];b[k]=b[js[k]];b[js[k]]=t;}
				delete js;
				return(1);
			}

