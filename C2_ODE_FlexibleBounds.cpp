
#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio>
#include <cmath>
#include <ctime>
#include <vector>

FILE *in, *out, *OriginalC, *CalculatC;
const double pi = 3.1415926;

int ciggj(double **, int, double *);
//int de(double, double, double, double &, double &, double *, double *);
int de(double, double, double, double &, double &, std::vector <double> &, std::vector <double> &);
int flip1 = 0;
int flip2 = 0;
int addCurves = 11;
double Param_cal = 0.001;

void main()
{
	in = fopen("frame31_outer_face_line&inner_face_line.dat", "rt");
	out = fopen("out.dat", "wt");
	OriginalC = fopen("OrigCurve_frame31_outer_face_line&inner_face_line.mel", "wt");
	CalculatC = fopen("out_frame31_outer_face_line&inner_face_line.mel", "wt");
	int i, j, k, n, curve, icurve, nu, nv, points;
	double ksi1, ksi2, a1, a2, a3, du, dv, ui, vi, x, y, z, a13[3][3], ksi[2][3], x1, y1, z1, L01[200], L12[200], L23[200], L34[200], L45[200], u2[200], u3[200], u5[200], u6[200], Li, L, vj[6][200], surfaceVbound[6][102][3];

	//************** Define the parameters_array d10v and d10 for d1 to d10 *****************//
	/*double* d10 = new double[10];
	double** d10v = new double*[10];
	for (i = 0; i<10; i++)
	{
		d10v[i] = new double[3];
	}*/
	std::vector <double> d10(10);
	std::vector <std::vector <double> > d10v(10, std::vector <double>(3));
	//************** Define the parameters_array e13v and e13 for e1 to e13 *****************//
	/*double* e13 = new double[13];
	double** e13v = new double*[13];
	for (i = 0; i<13; i++)
	{
		e13v[i] = new double[3];
	}*/
	std::vector <double> e13(13);
	std::vector <std::vector <double> > e13v(13, std::vector <double>(3));
	//************** Specifying the values of vector-valued shape parameters *****************//
	for (i = 0; i<3; i++)
	{
		for (j = 0; j<3; j++)
		{
			//	SP1
			/*
			if(i==0)a13[i][j]=1.0;
			if(i==1)a13[i][j]=3.0;
			if(i==2)a13[i][j]=1.0;
			*/
			//	SP2
			if (i == 0)a13[i][j] = 0.0001;
			if (i == 1)a13[i][j] = 3.0;
			if (i == 2)a13[i][j] = 0.0001;
		}
	}

	//************** Input the number of boundary constraints and the number of points on each boundary curve *****************//
	fscanf(in, "%d %d\n", &curve, &points);

	//************** Define the array to keep c1(v) to c6(v) *****************//
	//double*** c16 = new double**[curve]; // 
	//for (i = 0; i<curve; i++)
	//{
	//	c16[i] = new double*[points];
	//	for (j = 0; j<points; j++)
	//		c16[i][j] = new double[3]; // 
	//}
	//double*** temp = new double**[curve - 2]; // 
	//for (i = 0; i<curve - 2; i++)
	//{
	//	temp[i] = new double*[points];
	//	for (j = 0; j<points; j++)
	//		temp[i][j] = new double[3]; // 
	//}
	// Define 3D vector
	std::vector <std::vector <std::vector <double> > > c16;
	// Initialization 
	c16.resize(curve);
	for (int i = 0; i < curve; ++i) {
		c16[i].resize(points);

		for (int j = 0; j < points; ++j)
			c16[i][j].resize(3);
	}

	std::vector <std::vector <std::vector <double> > > temp;
	temp.resize(curve-2);
	for (int i = 0; i < curve-2; ++i) {
		temp[i].resize(points);

		for (int j = 0; j < points; ++j)
			temp[i][j].resize(3);
	}
	//************** Input c1(v) to c6(v) *****************//
	for (i = 0; i<curve; i++) //different from the first approach
	{
		for (int j = 0; j<points; j++)
		{
			fscanf(in, "%le %le %le\n", &c16[i][j][0], &c16[i][j][1], &c16[i][j][2]);
		}
	}

	//************** Generate Mel file of input curves *****************//
	for (i = 0; i<curve; i++)
	{
		fprintf(OriginalC, "curve -d 3\n");
		for (j = 0; j<points; j++)
		{
			fprintf(OriginalC, "-p %le %le %le \n", c16[i][j][0], c16[i][j][1], c16[i][j][2]);
		}
		fprintf(OriginalC, ";\n");
	}

	//************** Re-parametrization *****************//
	for (i = 0; i<curve; i++)
	{
		L = 0.0;
		for (int j = 0; j<points - 1; j++)
		{
			L = L + sqrt((c16[i][j + 1][0] - c16[i][j][0])*(c16[i][j + 1][0] - c16[i][j][0]) + (c16[i][j + 1][1] - c16[i][j][1])*(c16[i][j + 1][1] - c16[i][j][1]) + (c16[i][j + 1][2] - c16[i][j][2])*(c16[i][j + 1][2] - c16[i][j][2]));
		}
		vj[i][0] = 0.0;
		Li = 0.0;
		for (int j = 0; j<points - 1; j++)
		{
			Li = Li + sqrt((c16[i][j + 1][0] - c16[i][j][0])*(c16[i][j + 1][0] - c16[i][j][0]) + (c16[i][j + 1][1] - c16[i][j][1])*(c16[i][j + 1][1] - c16[i][j][1]) + (c16[i][j + 1][2] - c16[i][j][2])*(c16[i][j + 1][2] - c16[i][j][2]));
			vj[i][j + 1] = Li / L;
		}
	}

	//************** Interpolate 6 curves *****************//

	//************** Calculate the arc length between the 6 curves  (here some questions) *****************// 
	for (int j = 0; j<points; j++)
	{
		L01[j] = sqrt((c16[1][j][0] - c16[0][j][0])*(c16[1][j][0] - c16[0][j][0]) + (c16[1][j][1] - c16[0][j][1])*(c16[1][j][1] - c16[0][j][1]) + (c16[1][j][2] - c16[0][j][2])*(c16[1][j][2] - c16[0][j][2]));
		L12[j] = sqrt((c16[2][j][0] - c16[1][j][0])*(c16[2][j][0] - c16[1][j][0]) + (c16[2][j][1] - c16[1][j][1])*(c16[2][j][1] - c16[1][j][1]) + (c16[2][j][2] - c16[1][j][2])*(c16[2][j][2] - c16[1][j][2]));
		L23[j] = sqrt((c16[3][j][0] - c16[2][j][0])*(c16[3][j][0] - c16[2][j][0]) + (c16[3][j][1] - c16[2][j][1])*(c16[3][j][1] - c16[2][j][1]) + (c16[3][j][2] - c16[2][j][2])*(c16[3][j][2] - c16[2][j][2]));
		L34[j] = sqrt((c16[4][j][0] - c16[3][j][0])*(c16[4][j][0] - c16[3][j][0]) + (c16[4][j][1] - c16[3][j][1])*(c16[4][j][1] - c16[3][j][1]) + (c16[4][j][2] - c16[3][j][2])*(c16[4][j][2] - c16[3][j][2]));
		L45[j] = sqrt((c16[5][j][0] - c16[4][j][0])*(c16[5][j][0] - c16[4][j][0]) + (c16[5][j][1] - c16[4][j][1])*(c16[5][j][1] - c16[4][j][1]) + (c16[5][j][2] - c16[4][j][2])*(c16[5][j][2] - c16[4][j][2]));
		u2[j] = L01[j] / (L01[j] + L12[j] + L23[j] + L34[j] + L45[j]);
		u3[j] = (L01[j] + L12[j]) / (L01[j] + L12[j] + L23[j] + L34[j] + L45[j]);
		/*
		u5[j]=(L01[j]+L12[j]+L23[j])/(L01[j]+L12[j]+L23[j]+L34[j]+L45[j]);
		u6[j]=(L01[j]+L12[j]+L23[j]+L34[j])/(L01[j]+L12[j]+L23[j]+L34[j]+L45[j]);
		*/
		u5[j] = (L34[j] + L45[j]) / (L01[j] + L12[j] + L23[j] + L34[j] + L45[j]);
		u6[j] = L45[j] / (L01[j] + L12[j] + L23[j] + L34[j] + L45[j]);
	}

	//************** Calculate boundary tangent and boundary curves from 6 curves *****************//
	for (int j = 0; j<points; j++)
	{
		for (i = 0; i<3; i++)
		{
//			temp[0][j][i] = (c16[1][j][i] - c16[0][j][i]) / u2[j];
//			temp[1][j][i] = (u2[j] * c16[2][j][i] - u3[j] * c16[1][j][i] + (u3[j] - u2[j])*c16[0][j][i]) / u2[j] / u2[j] / (u3[j] - u2[j]);
			temp[0][j][i] = Param_cal *(c16[1][j][i] - c16[0][j][i]) / u2[j];
			temp[1][j][i] = Param_cal *(u2[j] * c16[2][j][i] - u3[j] * c16[1][j][i] + (u3[j] - u2[j])*c16[0][j][i]) / u2[j] / u2[j] / (u3[j] - u2[j]);
			/*
			temp[2][j][i]=(c16[5][j][i]-c16[4][j][i])/(1.0-u6[j]);
			temp[3][j][i]=((u6[j]-u5[j])*c16[5][j][i]+(1.0-u6[j])*c16[3][j][i]-(1.0-u5[j])*c16[4][j][i])/(1.0-u6[j])/(1.0-u6[j])/(u6[j]-u5[j]);
			*/
//			temp[2][j][i] = -(c16[4][j][i] - c16[5][j][i]) / u6[j];
//			temp[3][j][i] = (-u5[j] * c16[4][j][i] + u6[j] * c16[3][j][i] + (u5[j] - u6[j])*c16[5][j][i]) / u6[j] / u6[j] / (u5[j] - u6[j]);
			temp[2][j][i] = -Param_cal *(c16[4][j][i] - c16[5][j][i]) / u6[j];
			temp[3][j][i] = Param_cal *(-u5[j] * c16[4][j][i] + u6[j] * c16[3][j][i] + (u5[j] - u6[j])*c16[5][j][i]) / u6[j] / u6[j] / (u5[j] - u6[j]);

			if (flip1 == 1)
			{
				temp[0][j][i] = -temp[0][j][i];
				temp[1][j][i] = -temp[1][j][i];
			}
			if (flip2 == 1)
			{
				temp[2][j][i] = -temp[2][j][i];
				temp[3][j][i] = -temp[3][j][i];
			}
			

		}
	}
	//*******************FlexibleControlBounds*********************************************************//
	std::vector <std::vector <double> > Boundary1(points, std::vector <double>(3));
	std::vector <std::vector <double> > Tangent1(points, std::vector <double>(3));
	std::vector <std::vector <double> > Curvature1(points, std::vector <double>(3));
	std::vector <std::vector <double> > Boundary2(points, std::vector <double>(3));	
	std::vector <std::vector <double> > Tangent2(points, std::vector <double>(3));	
	std::vector <std::vector <double> > Curvature2(points, std::vector <double>(3));
	for (int j = 0; j<points; j++)
	{
		for (i = 0; i<3; i++)
		{
			Boundary1[j][i] = c16[0][j][i];
			Tangent1[j][i] = temp[0][j][i];
			Curvature1[j][i] = temp[1][j][i];

			Boundary2[j][i] = c16[5][j][i];
			Tangent2[j][i] = temp[2][j][i];
			Curvature2[j][i] = temp[3][j][i];
		}
	}

	int tan_controler1= 1.0,     // the parameters used to define tangent and curvature should be paid attention to 
		tan_controler2 = 1.0,
		cur_controler1 = 3.0,
		cur_controler2= 1.0;
	for (int j = 0; j<points; j++)
	{
		for (i = 0; i<3; i++)
		{
			c16[0][j][i] = Boundary1[j][i];
			c16[1][j][i] = tan_controler1 * Tangent1[j][i];
			c16[2][j][i] = cur_controler1 * Curvature1[j][i];

			c16[3][j][i] = Boundary2[j][i];
			c16[4][j][i] = tan_controler2 * Tangent2[j][i];
			c16[5][j][i] = cur_controler2 * Curvature2[j][i];
		}
	}
	//*************************************************************************************************//
	//for (int j = 0; j<points; j++)       //original fix Bounds 
	//{
	//	for (i = 0; i<3; i++)
	//	{
	//		c16[3][j][i] = c16[5][j][i];
	//		c16[1][j][i] = temp[0][j][i];
	//		c16[2][j][i] = temp[1][j][i];
	//		c16[4][j][i] = temp[2][j][i];
	//		c16[5][j][i] = temp[3][j][i];
	//	}
	//}

	//************** Print the boundary curves, boundary tangents and boundary curvature for checking purpose *****************//
	for (i = 0; i<curve; i++)
	{
		for (int j = 0; j<points; j++)
		{
			fprintf(out, "i,j,x,y,z= %d %d %le %le %le\n", i, j, c16[i][j][0], c16[i][j][1], c16[i][j][2]);
		}
	}

	//************** X component *****************//
	a1 = a13[0][0];
	a2 = a13[1][0];
	a3 = a13[2][0];
	de(a1, a2, a3, ksi1, ksi2, d10, e13);
	ksi[0][0] = ksi1;
	ksi[1][0] = ksi2;
	for (i = 0; i<10; i++)
	{
		d10v[i][0] = d10[i];
	}
	for (i = 0; i<13; i++)
	{
		e13v[i][0] = e13[i];
	}

	//************** Y component *****************//
	a1 = a13[0][1];
	a2 = a13[1][1];
	a3 = a13[2][1];
	de(a1, a2, a3, ksi1, ksi2, d10, e13);
	ksi[0][1] = ksi1;
	ksi[1][1] = ksi2;
	for (i = 0; i<10; i++)
	{
		d10v[i][1] = d10[i];
	}
	for (i = 0; i<13; i++)
	{
		e13v[i][1] = e13[i];
	}

	//************** Z component *****************//
	a1 = a13[0][2];
	a2 = a13[1][2];
	a3 = a13[2][2];
	de(a1, a2, a3, ksi1, ksi2, d10, e13);
	ksi[0][2] = ksi1;
	ksi[1][2] = ksi2;
	for (i = 0; i<10; i++)
	{
		d10v[i][2] = d10[i];
	}
	for (i = 0; i<13; i++)
	{
		e13v[i][2] = e13[i];
	}

	//************** Create surfaces *****************//
	double* bb = new double[6]; // 
	double** aa = new double*[6]; // 
	for (i = 0; i<6; i++)
	{
		aa[i] = new double[6];
	}

	//nu=100;
	//nu = 15;
	nu = addCurves;
	nv = points;
	du = 1.0 / nu;
	dv = 1.0 / nv;

	for (i = 0; i <= nu; i++)
	{
		ui = i * du;
		fprintf(CalculatC, "curve -d 3\n");
		int h = 0;
		for (j = 0; j<nv; j++)
		{
			// vi=j*dv;
			vi = vj[0][j];
			int m, n;
			double a11, a12, a21, a22, b1v, b2v, a3, a4, a5, a6, c1, c2, c3, c4, c5, c6, bc1, bc2, bc3, bc4, bc5, bc6;
			//	x component
			ksi1 = ksi[0][0];
			ksi2 = ksi[1][0];
			for (m = 0; m<6; m++)
			{
				bb[m] = 0.0;
				for (n = 0; n<6; n++)
				{
					aa[m][n] = 0.0;
				}
			}
			aa[0][0] = 1.0;
			aa[0][2] = 1.0;
			aa[0][5] = 1.0;
			aa[1][1] = ksi1;
			aa[1][3] = ksi2;
			aa[1][4] = 1.0;
			aa[2][0] = -ksi1 * ksi1;
			aa[2][2] = -ksi2 * ksi2;
			aa[3][0] = cos(ksi1);
			aa[3][1] = sin(ksi1);
			aa[3][2] = cos(ksi2);
			aa[3][3] = sin(ksi2);
			aa[3][4] = 1.0;
			aa[3][5] = 1.0;
			aa[4][0] = -ksi1 * sin(ksi1);
			aa[4][1] = ksi1 * cos(ksi1);
			aa[4][2] = -ksi2 * sin(ksi2);
			aa[4][3] = ksi2 * cos(ksi2);
			aa[4][4] = 1.0;
			aa[5][0] = -ksi1 * ksi1*cos(ksi1);
			aa[5][1] = -ksi1 * ksi1*sin(ksi1);
			aa[5][2] = -ksi2 * ksi2*cos(ksi2);
			aa[5][3] = -ksi2 * ksi2*sin(ksi2);
			bb[0] = c16[0][j][0];
			bb[1] = c16[1][j][0];
			bb[2] = c16[2][j][0];
			bb[3] = c16[3][j][0];
			bb[4] = c16[4][j][0];
			bb[5] = c16[5][j][0];
			ciggj(aa, 6, bb);
			c1 = bb[0];
			c2 = bb[1];
			c3 = bb[2];
			c4 = bb[3];
			c5 = bb[4];
			c6 = bb[5];
			x1 = c1 * cos(ksi1*ui) + c2 * sin(ksi1*ui) + c3 * cos(ksi2*ui) + c4 * sin(ksi2*ui) + c5 * ui + c6;
			//	y component
			ksi1 = ksi[0][1];
			ksi2 = ksi[1][1];
			for (m = 0; m<6; m++)
			{
				bb[m] = 0.0;
				for (n = 0; n<6; n++)
				{
					aa[m][n] = 0.0;
				}
			}
			aa[0][0] = 1.0;
			aa[0][2] = 1.0;
			aa[0][5] = 1.0;
			aa[1][1] = ksi1;
			aa[1][3] = ksi2;
			aa[1][4] = 1.0;
			aa[2][0] = -ksi1 * ksi1;
			aa[2][2] = -ksi2 * ksi2;
			aa[3][0] = cos(ksi1);
			aa[3][1] = sin(ksi1);
			aa[3][2] = cos(ksi2);
			aa[3][3] = sin(ksi2);
			aa[3][4] = 1.0;
			aa[3][5] = 1.0;
			aa[4][0] = -ksi1 * sin(ksi1);
			aa[4][1] = ksi1 * cos(ksi1);
			aa[4][2] = -ksi2 * sin(ksi2);
			aa[4][3] = ksi2 * cos(ksi2);
			aa[4][4] = 1.0;
			aa[5][0] = -ksi1 * ksi1*cos(ksi1);
			aa[5][1] = -ksi1 * ksi1*sin(ksi1);
			aa[5][2] = -ksi2 * ksi2*cos(ksi2);
			aa[5][3] = -ksi2 * ksi2*sin(ksi2);
			bb[0] = c16[0][j][1];
			bb[1] = c16[1][j][1];
			bb[2] = c16[2][j][1];
			bb[3] = c16[3][j][1];
			bb[4] = c16[4][j][1];
			bb[5] = c16[5][j][1];
			ciggj(aa, 6, bb);
			c1 = bb[0];
			c2 = bb[1];
			c3 = bb[2];
			c4 = bb[3];
			c5 = bb[4];
			c6 = bb[5];
			y1 = c1 * cos(ksi1*ui) + c2 * sin(ksi1*ui) + c3 * cos(ksi2*ui) + c4 * sin(ksi2*ui) + c5 * ui + c6;
			//	z component
			ksi1 = ksi[0][2];
			ksi2 = ksi[1][2];
			for (m = 0; m<6; m++)
			{
				bb[m] = 0.0;
				for (n = 0; n<6; n++)
				{
					aa[m][n] = 0.0;
				}
			}
			aa[0][0] = 1.0;
			aa[0][2] = 1.0;
			aa[0][5] = 1.0;
			aa[1][1] = ksi1;
			aa[1][3] = ksi2;
			aa[1][4] = 1.0;
			aa[2][0] = -ksi1 * ksi1;
			aa[2][2] = -ksi2 * ksi2;
			aa[3][0] = cos(ksi1);
			aa[3][1] = sin(ksi1);
			aa[3][2] = cos(ksi2);
			aa[3][3] = sin(ksi2);
			aa[3][4] = 1.0;
			aa[3][5] = 1.0;
			aa[4][0] = -ksi1 * sin(ksi1);
			aa[4][1] = ksi1 * cos(ksi1);
			aa[4][2] = -ksi2 * sin(ksi2);
			aa[4][3] = ksi2 * cos(ksi2);
			aa[4][4] = 1.0;
			aa[5][0] = -ksi1 * ksi1*cos(ksi1);
			aa[5][1] = -ksi1 * ksi1*sin(ksi1);
			aa[5][2] = -ksi2 * ksi2*cos(ksi2);
			aa[5][3] = -ksi2 * ksi2*sin(ksi2);
			bb[0] = c16[0][j][2];
			bb[1] = c16[1][j][2];
			bb[2] = c16[2][j][2];
			bb[3] = c16[3][j][2];
			bb[4] = c16[4][j][2];
			bb[5] = c16[5][j][2];
			ciggj(aa, 6, bb);
			c1 = bb[0];
			c2 = bb[1];
			c3 = bb[2];
			c4 = bb[3];
			c5 = bb[4];
			c6 = bb[5];
			z1 = c1 * cos(ksi1*ui) + c2 * sin(ksi1*ui) + c3 * cos(ksi2*ui) + c4 * sin(ksi2*ui) + c5 * ui + c6;
			fprintf(CalculatC, "-p %le %le %le\n", x1, y1, z1);

			/*//  Save V direction boundaries: v1 to v6
			if (j==0||j==1||j==2||j==(nv-3)||j==(nv-2)||j==(nv-1))
			{
			surfaceVbound[h][i][0]=x1;
			surfaceVbound[h][i][1]=y1;
			surfaceVbound[h][i][2]=z1;
			h++;
			}*/
		}
		fprintf(CalculatC, ";\n");
	}
	/*//  draw V direction boundaries: v1 to v6
	for(int h=0;h<=5;h++)
	{
	fprintf(outToMayaout,"curve -d 3\n");
	for(i=0;i<=nu;i++)
	{
	fprintf(outToMayaout,"-p %le %le %le\n",surfaceVbound[h][i][0],surfaceVbound[h][i][1],surfaceVbound[h][i][2]);
	}
	fprintf(outToMayaout,";\n");
	}*/
}

int de(double a1, double a2, double a3, double &ksi1, double &ksi2, std::vector <double> &d10, std::vector <double> &e13)
{
	int i;
	ksi1 = sqrt(a2*(1.0 + sqrt(1.0 - 4.0*a1*a3 / a2 / a2)) / 2.0 / a1);
	ksi2 = sqrt(a2*(1.0 - sqrt(1.0 - 4.0*a1*a3 / a2 / a2)) / 2.0 / a1);
	fprintf(out, "inside ksi1,ksi2= %le %le\n", ksi1, ksi2);
	e13[0] = ksi1 * (cos(ksi1) - 1.0) + ksi1 * ksi1*sin(ksi1)*(1.0 - cos(ksi2)) / ksi2 / sin(ksi2);
	e13[1] = sin(ksi1) - ksi1 + ksi1 * ksi1*sin(ksi1)*(1.0 / sin(ksi2) - 1.0 / ksi2) / ksi2;
	e13[2] = cos(ksi1) - 1.0 + ksi1 * ksi1*(1.0 - cos(ksi2)) / ksi2 / ksi2 + ksi1 * ksi1*(cos(ksi2) - cos(ksi1))
		*(1.0 / ksi2 - 1.0 / sin(ksi2)) / ksi2;
	e13[3] = ksi1 * (-sin(ksi1) + ksi1 * sin(ksi2) / ksi2) + ksi1 * ksi1*(cos(ksi2) - cos(ksi1))*(cos(ksi2)
		- 1.0) / ksi2 / sin(ksi2);
	e13[4] = (1.0 / ksi2 - cos(ksi2) / sin(ksi2)) / ksi2;
	e13[5] = (sin(ksi2) + cos(ksi2)*cos(ksi2) / sin(ksi2) - cos(ksi2) / sin(ksi2)) / ksi2;
	e13[6] = (1.0 / sin(ksi2) - 1.0 / ksi2) / ksi2;
	e13[7] = (1.0 - cos(ksi2)) / ksi2 / sin(ksi2);
	e13[10] = ksi1 * ksi1 / ksi2 / ksi2 - 1.0;
	e13[11] = cos(ksi2) / ksi2 / sin(ksi2);
	e13[12] = 1.0 / ksi2 / sin(ksi2);
	d10[0] = e13[0] / (e13[0] * e13[2] - e13[1] * e13[3]);
	d10[1] = -e13[1] / (e13[0] * e13[2] - e13[1] * e13[3]);
	d10[2] = e13[2] / (e13[0] * e13[2] - e13[1] * e13[3]);
	d10[3] = -e13[3] / (e13[0] * e13[2] - e13[1] * e13[3]);
	e13[8] = ksi1 * ksi1*(d10[3] * sin(ksi1) - d10[0] * (cos(ksi2) - cos(ksi1))) / ksi2 / ksi2 / sin(ksi2);
	e13[9] = ksi1 * ksi1*(d10[2] * sin(ksi1) - d10[1] * (cos(ksi2) - cos(ksi1))) / ksi2 / ksi2 / sin(ksi2);
	d10[4] = d10[0] * e13[4] + d10[1] * e13[5];
	d10[5] = d10[0] * e13[6] + d10[1] * e13[7];
	d10[6] = d10[3] * e13[4] + d10[2] * e13[5];
	d10[7] = d10[3] * e13[6] + d10[2] * e13[7];
	d10[8] = (d10[0] + d10[1])*(e13[10] + 1.0);
	d10[9] = -1.0 / ksi2 / ksi2 + (e13[10] + 1.0)*d10[4];
	fprintf(out, "ksi2,d10[4],d10[9],e13[10]= %le %le %le %le\n", ksi2, d10[4], d10[9], e13[10]);
	return(1);
}

// ciggj solver
int ciggj(double **a, int n, double *b)
{
	int i, j, k, is, u, v;
	int *js = new int[n];
	double d, t;
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				t = fabs(a[i][j]);
				if (t>d) { d = t; js[k] = j; is = i; }
			}
		if (d + 1.0 == 1.0)
		{
			delete js; printf("fail\n"); return(0);
		}
		if (is != k)
		{
			for (j = k; j <= n - 1; j++)
			{
				u = k * n + j; v = is * n + j;
				t = a[k][j]; a[k][j] = a[is][j]; a[is][j] = t;
			}
			t = b[k]; b[k] = b[is]; b[is] = t;
		}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k; v = i * n + js[k];
				t = a[i][k]; a[i][k] = a[i][js[k]]; a[i][js[k]] = t;
			}
		t = a[k][k];
		for (j = k + 1; j <= n - 1; j++)
		{
			u = k * n + j;
			if (a[k][j] != 0.0)a[k][j] = a[k][j] / t;
		}
		b[k] = b[k] / t;
		for (j = k + 1; j <= n - 1; j++)
		{
			u = k * n + j;
			if (a[k][j] != 0.0)
			{
				for (i = 0; i <= n - 1; i++)
				{
					v = i * n + k;
					if ((i != k) && (a[i][k] != 0.0))
					{
						is = i * n + j;
						a[i][j] = a[i][j] - a[i][k] * a[k][j];
					}
				}
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			u = i * n + k;
			if ((i != k) && (a[i][k] != 0.0))
				b[i] = b[i] - a[i][k] * b[k];
		}
	}
	for (k = n - 1; k >= 0; k--)
		if (k != js[k])
		{
			t = b[k]; b[k] = b[js[k]]; b[js[k]] = t;
		}
	delete js;
	return(1);
}

