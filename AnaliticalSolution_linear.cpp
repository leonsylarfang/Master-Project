#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio>
#include <cmath>
#include <ctime>

FILE* in0, * in1, * outToMaya;

int main()
{
	int i, j, curve, points;
	double dv = 5;
	in0 = fopen("in_laugh_outer_face_line&inner_face_line_6.dat", "rt");
	in1 = fopen("in_reference_outer_face_line&inner_face_line_6.dat", "rt");
	outToMaya = fopen("keyframes_outer_face_line&inner_face_line_6.mel", "wt");
	fscanf(in0, "%d %d\n", &curve, &points);
	fscanf(in1, "%d %d\n", &curve, &points);
	double*** c16_1 = new double** [curve]; // 
	double*** c16_2 = new double** [curve];
	for (i = 0; i < curve; i++)
	{
		c16_1[i] = new double* [points];
		for (j = 0; j < points; j++)
			c16_1[i][j] = new double[3]; // 
	}

	for (i = 0; i < curve; i++)
	{
		c16_2[i] = new double* [points];
		for (j = 0; j < points; j++)
			c16_2[i][j] = new double[3]; // 
	}



	for (i = 0; i < curve; i++)
	{
		for (int j = 0; j < points; j++)
		{
			fscanf(in0, "%le %le %le\n", &c16_1[i][j][0], &c16_1[i][j][1], &c16_1[i][j][2]);
		}
	}

	for (i = 0; i < curve; i++)
	{
		for (int j = 0; j < points; j++)
		{
			fscanf(in1, "%le %le %le\n", &c16_2[i][j][0], &c16_2[i][j][1], &c16_2[i][j][2]);
		}
	}


	double nn = 1 / dv;
	double ti = 0;

	for (i = 0; i <= dv; i++)
	{
		ti = i * nn;
		fprintf(outToMaya, "curve -d 3\n");
		for (int j = 0; j < points; j++)
		{
			double xn = 0;
			double yn = 0;
			double zn = 0;
			xn = c16_1[0][j][0] + (c16_2[0][j][0] - c16_1[0][j][0]) * ti;
			yn = c16_1[0][j][1] + (c16_2[0][j][1] - c16_1[0][j][1]) * ti;
			zn = c16_1[0][j][2] + (c16_2[0][j][2] - c16_1[0][j][2]) * ti;
			fprintf(outToMaya, "-p %le %le %le\n", xn, yn, zn);
		}
		fprintf(outToMaya, ";\n");
	}
	return 0;
}