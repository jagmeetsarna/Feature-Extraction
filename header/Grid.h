#pragma once

#include <vector>
#include "VectorAttributes.h"
#include "ScalarAttributes.h"
#include <stdio.h>
#include <string>
#include <Eigen/Dense>

using namespace std;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

class Grid
{
public:
	Grid(int P, int C) : scalars(P), pointNormals(P), faceNormals(C), pointsZ(P)															//Resize the array for the points and cells for the grid
	{
		pointsX.resize(P);
		pointsY.resize(P);
		pointsZ.resize(P);
		cells.resize(3 * C);

		for (int i = 0; i < 10; i++) {
			A1hist.push_back(0);
			D1hist.push_back(0);
			D2hist.push_back(0);
			D3hist.push_back(0);
			D4hist.push_back(0);
			A3hist.push_back(0);
		}
	}

	float calculateSurfaceArea();
	float calculateEccentricity();

	float calculateBoundingBoxVol();
	float calculateDiameter();

	float calculateVolume();

	float calculateSphericity();

	double calculateAngleBetweenPoints();



	void setExtremes(float mix, float max, float miy, float may, float miz, float maz) {
		minX = mix;
		maxX = max;
		minY = miy;
		maxY = may;
		minZ = miz;
		maxZ = maz;
	}

	int numPoints()
	{
		return pointsX.size();
	}

	int numCells()
	{
		return cells.size() / 3;
	}

	void getPoint(int i, float* p);

	void setPoint(int i, vector<float> p);

	void setCell(int cell, vector<int> vertices);

	void setClass(string newcls);

	string getClass();

	int	 getCell(int cell, int* vertices);

	int findCell(float* p);

	void normalize();

	void computeFaceNormals();

	void computeVertexNormals();

	void computeCovarianceMatrix();

	void computeEigenvectors();

	VectorAttributes& getFaceNormals()
	{
		return faceNormals;
	}

	VectorAttributes& getPointNormals()
	{
		return pointNormals;
	}

	Eigen::Matrix3f& getCovarianceMatrix() {
		return covarianceMatrix;
	}

	vector<vector<float>> getEigenvectors() {
		vector<vector<float>> vectors;
		vectors.push_back(eigenVec1);
		vectors.push_back(eigenVec2);
		vectors.push_back(eigenVec3);

		return vectors;
	}

	vector<Point3d> getCellCentroids();

	void momentTest();

	int sgn(float x) {
		if (x > 0) return 1;
		if (x < 0) return -1;
		return 1;
	}
	void PCARotation();

	void calculateD1();
	vector<float> getD1hist() {
		return D1hist;
	}

	void calculateD2(int n);
	vector<float> getD2hist() {
		return D2hist;
	}

	void calculateD3(int n);
	vector<float> getD3hist() {
		return D3hist;
	}

	void calculateD4(int n);
	vector<float> getD4hist() {
		return D4hist;
	}

	void calculateAngleBetweenPoints(int n);
	vector<float> getA3hist() {
		return A3hist;
	}


protected:

	ScalarAttributes	scalars;

	Eigen::Matrix3f		covarianceMatrix;

	vector<float>		eigenVec1, eigenVec2, eigenVec3;

	vector<Point3d>		cellCentroids;

	vector<float>		pointsX, pointsY, pointsZ;
	vector<int>			cells;
	VectorAttributes    pointNormals;
	VectorAttributes    faceNormals;
	std::string	cls;

	float				surfaceArea;
	float				volume;
	float				compactness;
	float				sphericity;
	float				minX, maxX, minY, maxY, minZ, maxZ;
	float				boundingBoxVolume;
	float				diameter;
	float				eccentricity;
	vector<float>		A1hist, D1hist, D2hist, D3hist, D4hist,A3hist;

};



