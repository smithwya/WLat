#pragma once
#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;

class lattice
{

public:
	struct position {
		position();
		position(int a, int b, int c, int d);
		int t;
		int x;
		int y;
		int z;
	};

	position shift(position n, int dir, int amount);
	lattice();
	lattice(int nColors, int spatialExtent, int timeExtent, double coupling, double anisotropy);
	double getAnisotropy();
	double getBeta();
	int getNc();
	double getBareAnisotropy(double xiR);
	void ColdStart();
	void readFile(string filename);
	void writeFile(string filename);
	MatrixXcd getLink(position n, int dir, bool fwd);
	void update(int numUpdates);
	MatrixXcd matrixPow(MatrixXcd G, double omega, int cutoff);
	void fixCoulombGauge(double tolerance);
	void fixCoulombGaugeSubgroups(double tolerance);
	double checkCoulombGauge();
	vector<double> measureObservable(function<vector<double>(position)> func);
	double getAvgPlaq();
	double getCorrelator(int mu, int nu, int R, int T);
	double getLoop(int mu, int nu, int R, int T);
	MatrixXcd getCoulombLoop(position n, int mu, int r, int t);
	MatrixXcd getLine(position n, int dir, int length);
	MatrixXcd getLoop(position n, int dir1, int dir2, int len1, int len2);
	vector<double> getGreensiteCorrelator(int rmax);
	vector<double> getImprovedCorrelator(int rmax, int t);
	vector<double> getAverageSSLoop(int rmax, int t);
	vector<double> getAverageTSLoop(int rmax, int t);
	vector<double> getPolyakovLoop(int rmax);
	vector<double> getSquareLoop(int rmax);
	vector<double> getAverageWilsonLoop(int rmax, int t);

private:
	int Nc;
	int N;
	int T;
	double beta;
	double xi_R;
	double xi_0;
	vector<vector<vector<vector<vector<MatrixXcd>>>>> configuration;
	void transformLat(function<MatrixXcd(position,int)>);
	MatrixXcd genStaple(position n, int mu);
	MatrixXcd genLinkSU2(MatrixXcd A);
	MatrixXcd genLinkSU3(position n, int mu);
	MatrixXcd getSubgroup(MatrixXcd U, int i);
	MatrixXcd SU3mat(MatrixXcd U,int i);
	MatrixXcd reunitarizeSU2(MatrixXcd U);
	MatrixXcd reunitarizeSU3(MatrixXcd U);

	double randNum();

};

