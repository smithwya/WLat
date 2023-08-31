#include "lattice.h"
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <ctime>

using namespace std;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
typedef std::complex<double> comp;
//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//mt19937 generator(seed);
// seed
random_device rd;
seed_seq sd{rd(), rd(), rd(), rd()};
auto generator = mt19937(sd);


//returns random number between -1/2 and 1/2

double lattice::randNum() {
	return (double)generator() / generator.max() - 0.5;
}

//position constructor defaults to (0,0,0,0)
lattice::position::position() {
	t = 0;
	x = 0;
	y = 0;
	z = 0;
}

//Position constructor for (a,b,c,d)
lattice::position::position(int a, int b, int c, int d) {
	t = a;
	x = b;
	y = c;
	z = d;
}

//shifts a position n in direction dir by amount
lattice::position lattice::shift(position n, int dir, int amount)
{
	switch (dir) {
	case 0: return position((n.t + amount + T) % T, n.x, n.y, n.z);

	case 1: return position(n.t, (n.x + amount + N) % N, n.y, n.z);

	case 2: return position(n.t, n.x, (n.y + amount + N) % N, n.z);

	case 3: return position(n.t, n.x, n.y, (n.z + amount + N) % N);
	}
	return n;
}

//default constructor for lattice sets all fields to 0
lattice::lattice()
{
	configuration = {};
	Nc = 0;
	N = 0;
	T = 0;
	beta = 0;
	xi_R = 0;
	xi_0 = 0;
}

//constructor for lattice
lattice::lattice(int nColors, int spatialExtent, int timeExtent, double coupling, double anisotropy)
{
	Nc = nColors;
	N = spatialExtent;
	T = timeExtent;
	beta = coupling;
	xi_R = anisotropy;
	xi_0 = setAnisotropy(xi_R);
	ColdStart();
}

//returns renormalized anisotropy
double lattice::getAnisotropy()
{
	return xi_R;
}

//returns beta
double lattice::getBeta()
{
	return beta;
}


//returns number of colors
int lattice::getNc()
{
	return Nc;
}

//sets the bare anisotropy from renorm anisotropy using JPAC parameterization
double lattice::setAnisotropy(double xiR)
{

	if (xiR == 1.0) return 1.0;

	double c0 = 0.3398;
	double c1 = -0.3068;
	double c2 = -0.0449;

	double eta1 = c0 + c1 / xiR + c2 / (pow(xiR, 3));

	double a1 = -1.0238;
	double a2 = -1.8921;

	return xiR / (1 + (4.0 / 6.0) * (1 + a1 / beta) / (1 + a2 / beta) * (eta1 / beta));
}

//initializes lattice to identity matrices
void lattice::ColdStart()
{
	MatrixXcd U = MatrixXcd::Identity(Nc, Nc);

	vector<MatrixXcd> directions(4, U);
	vector<vector<MatrixXcd>> Z(N, directions);
	vector<vector<vector<MatrixXcd>>> Y(N, Z);
	vector<vector<vector<vector<MatrixXcd>>>> X(N, Y);
	vector<vector<vector<vector<vector<MatrixXcd>>>>> Ts(T, X);
	configuration = Ts;
	return;
}

//reads in a lattice configuration from a text file
void lattice::readFile(string filename) {

	std::fstream infile(filename);

	for (int i = 0; i < T; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				for (int l = 0; l < N; l++) {

					for (int mu = 0; mu < 4; mu++) {

						comp a, b;
						infile >> a;
						infile >> b;

						configuration[i][j][k][l][mu] << a, b, -conj(b), conj(a);

					}
				}

			}

		}
	}
	infile.close();
}

//writes a text file for the current config
void lattice::writeFile(string filename)
{
	std::ofstream outfile(filename);

	for (int i = 0; i < T; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				for (int l = 0; l < N; l++) {

					position n = position(i, j, k, l);

					for (int mu = 0; mu < 4; mu++) {

						MatrixXcd U = getLink(n, mu, true);

						outfile << U(0, 0) << " " << U(0, 1) << endl;

					}
				}
			}
		}
	}

	outfile.close();
}

//returns the link at position n in direction (fwd*dir)
MatrixXcd lattice::getLink(position n, int dir, bool fwd)
{
	if (fwd) {
		return configuration[n.t][n.x][n.y][n.z][dir];
	}
	else {
		position m = shift(n, dir, -1);
		return getLink(m, dir, true).adjoint();
	}
}

//averages a function f(n)=vector<double> over the entire lattice
vector<double> lattice::measureObservable(function<vector<double>(lattice::position)> func) {

	int len = func(lattice::position()).size();
	double norm = T * pow(N, 3);
	vector<double> results = vector<double>(len, 0);

	for (int i = 0; i < T; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				for (int l = 0; l < N; l++) {
					position n = position(i, j, k, l);
					vector<double> tempres = func(n);
					for (int r = 0; r < len; r++) {
						results[r] += tempres[r];


					}
				}
			}
		}
	}

	std::transform(results.begin(), results.end(), results.begin(), [norm](double& c) { return c / norm; });

	return results;
}

//replaces U_mu(n) with f(n,mu) on entire lattice
void lattice::transformLat(function<MatrixXcd(position, int)> func) {

	for (int i = 0; i < T; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				for (int l = 0; l < N; l++) {
					position n = position(i, j, k, l);
					for (int mu = 0; mu < 4; mu++) {
						configuration[i][j][k][l][mu] = func(n, mu);
					}
				}
			}
		}
	}

}

//generates the sum of staples
MatrixXcd lattice::genStaple(position n, int mu)
{
	MatrixXcd A = MatrixXcd::Zero(Nc, Nc);
	MatrixXcd B, C, D;

	for (int nu = 0; nu <= 3; nu++) {
		if (nu != mu) {
			B = getLink(shift(n, mu, 1), nu, true) * getLink(shift(shift(n, mu, 1), nu, 1), mu, false) * getLink(shift(n, nu, 1), nu, false);
			C = getLink(shift(n, mu, 1), nu, false) * getLink(shift(shift(n, mu, 1), nu, -1), mu, false) * getLink(shift(n, nu, -1), nu, true);
			D = B + C;

			if (mu == 0 || nu == 0) D = D * xi_0;
			else D = D / xi_0;

			A = A + D;


		}
	}

	return A;
}

//generates a new SU(2) U_mu(n) using the heatbath algorithm
MatrixXcd lattice::genLink(position n, int mu)
{

	MatrixXcd A = genStaple(n, mu);

	double det = sqrt(A.determinant().real());

	double r1, r2, r3;
	double x0, x1, x2, x3;
	double lambdasq = 0;

	MatrixXcd V = A / det;

	while (true) {
		r1 = 0.5 + randNum();
		r2 = 0.5 + randNum();
		r3 = 0.5 + randNum();
		lambdasq = (-1 / (2 * det * beta)) * (log(r1) + pow(cos(2 * 3.14159265 * r2), 2) * log(r3));
		double r = 0.5 + randNum();
		if (pow(r, 2) <= 1 - lambdasq) break;
	}
	x0 = 1 - 2 * lambdasq;

	while (true) {
		r1 = 2 * randNum();
		r2 = 2 * randNum();
		r3 = 2 * randNum();
		if (r1 * r1 + r2 * r2 + r3 * r3 <= 1) break;
	}

	double length = sqrt(pow(r1, 2) + pow(r2, 2) + pow(r3, 2));
	double norm = sqrt(1 - pow(x0, 2));

	x1 = r1 * norm / length;
	x2 = r2 * norm / length;
	x3 = r3 * norm / length;

	MatrixXcd X = MatrixXcd::Zero(Nc, Nc);
	X << comp(x0, x3), comp(x2, x1), comp(-x2, x1), comp(x0, -x3);

	return reunitarize(X * V.adjoint());
}

//reunitarizes SU(2) matrix
MatrixXcd lattice::reunitarize(MatrixXcd U)
{
	MatrixXcd newU = MatrixXcd::Zero(Nc, Nc);
	comp x = sqrt(U(0,0)*conj(U(0,0))+U(0,1)*conj(U(0,1)));
	
	newU << U(0, 0), U(0, 1), -conj(U(0, 1)), conj(U(0, 0));
	return newU/x;
}

//thermalizes the lattice by performing numUpdates sweeps of the heatbath algorithm
void lattice::update(int numUpdates)
{

	auto func = [&](position n, int mu) {
		return genLink(n, mu);
	};


	for (int i = 0; i < numUpdates; i++) {
		transformLat(func);
	}

	return;
}

//returns the omega'th power of matrix G, up to cutoff terms in a series
MatrixXcd lattice::matrixPow(MatrixXcd G, double omega, int cutoff) {

	MatrixXcd Gw = MatrixXcd::Zero(Nc,Nc);
	MatrixXcd GminusI = G-MatrixXcd::Identity(Nc,Nc);

	for (int i = 0; i <= cutoff; i++) {

		MatrixXcd temp = MatrixXcd::Identity(Nc, Nc);

		for (int j = 1; j <= i; j++) {
			temp = GminusI*temp;
		}

		double c = tgamma(omega + 1) / (tgamma(omega + 1 - i) * tgamma(i + 1));

		Gw = Gw+c*temp;
	}

	return reunitarize(Gw);
}


//returns F
double lattice::checkCoulombGauge() {

	auto F = [&](position n) {
		double sum = 0;
		
		for (int mu = 1; mu < 4; mu++) {
			sum+=getLink(n,mu,true).trace().real();
		}
		vector<double> dat = {sum/3.0};
		return dat;
	};
	double Fval= measureObservable(F)[0]/4.0;
	
	return Fval;
}

//fixes the coulomb gauge up to tolerance dF<=tolerance
void lattice::fixCoulombGauge(double tolerance) {
	double omega = 1.75;

	auto gsweep = [&](position n) {
		MatrixXcd K = MatrixXcd::Zero(Nc, Nc);

		//1 to 3 gives Coulomb gauge, 0 to 3 gives Landau gauge
		for (int dir = 1; dir < 4; dir++) {
			K += getLink(n, dir, true);
			K += getLink(shift(n, dir, -1), dir, true).adjoint();
		}

		MatrixXcd G = K.adjoint();
		comp a = G.determinant();
		G = G / sqrt(a);

		//overrelaxation by replacing gauge transf g with g^omega
		G = matrixPow(G, omega, 2);

		for (int mu = 0; mu < 4; mu++) {
			configuration[n.t][n.x][n.y][n.z][mu] = G * getLink(n, mu, true);

			position m = shift(n, mu, -1);
			configuration[m.t][m.x][m.y][m.z][mu] = getLink(m, mu, true) * (G.adjoint());
		}
		vector<double> count = {1.0};
		return count;
	};
	
	int nsweeps = 0;
	double last = 0;
	double current = checkCoulombGauge();
	
	while(true){
		last = current;
		nsweeps += measureObservable(gsweep)[0];
		current = checkCoulombGauge();

		if(abs(current-last)<=tolerance) break;
	}
	return;
}



//returns the product of matrices in one direction
MatrixXcd lattice::getLine(position n, int dir, int length) {
	MatrixXcd L = MatrixXcd::Identity(Nc,Nc);
	
	for (int i = 0; i < length; i++) {
		L = L*getLink(shift(n, dir, i), dir, true);
	}
	return L;
}

//returns a loop of matrices
MatrixXcd lattice::getLoop(position n, int dir1, int dir2, int length1, int length2) {

	MatrixXcd L1 = getLine(n, dir1, length1);
	MatrixXcd L2 = getLine(shift(n, dir1, length1), dir2, length2);
	MatrixXcd L3 = getLine(shift(n, dir2, length2), dir1, length1);
	MatrixXcd L4 = getLine(n, dir2, length2);
	return L1*L2*(L3.adjoint())*(L4.adjoint());
}

MatrixXcd lattice::getCoulombLoop(position n, int mu, int r, int t){
	
	MatrixXcd L2 = getLine(shift(n, mu, r), 0, t);
	MatrixXcd L4 = getLine(n, 0, t);
	return L2*(L4.adjoint());
	
}

vector<double> lattice::getGreensiteCorrelator(int rmax){
	
	auto greenFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 0; r <= rmax; r++){
			double sum = 0;
			for(int mu = 1; mu <4; mu++){
				sum+=getCoulombLoop(n,mu,r,1).trace().real();
			}
			data.push_back(sum/3);
		}
		
		return data;
	};
	
	return measureObservable(greenFunc);
	
}

vector<double> lattice::getImprovedCorrelator(int rmax, int t){
	
	auto impFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 0; r <= rmax; r++){
			double sum = 0;
			for(int mu = 1; mu <4; mu++){
				sum+=getCoulombLoop(n,mu,r,t).trace().real();
			}
			data.push_back(sum/3);
		}
		
		return data;
	};
	
	return measureObservable(impFunc);
	
}


vector<double> lattice::getAverageSSLoop(int rmax, int t){
	
	auto impFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 0; r <= rmax; r++){
			double sum = 0;
			
			sum+=getLoop(n,1,2,r,t).trace().real();
			sum+=getLoop(n,1,3,r,t).trace().real();
			sum+=getLoop(n,2,3,r,t).trace().real();
			
			data.push_back(sum/3);
		}
		
		return data;
	};
	
	return measureObservable(impFunc);
	
}

vector<double> lattice::getAverageTSLoop(int rmax, int t){
	
	auto impFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 0; r <= rmax; r++){
			double sum = 0;
			
			for(int mu = 1; mu <4; mu++){
				sum+=getLoop(n,mu,0,r,t).trace().real();
			}
			data.push_back(sum/3);
		}
		
		return data;
	};
	
	return measureObservable(impFunc);
	
}

//returns the average plaquette
double lattice::getAvgPlaq() {

	auto plaq = [&](position n) {

		MatrixXcd P = MatrixXcd::Zero(Nc, Nc);

		for (int mu = 0; mu < 4; mu++) {
			for (int nu = 0; nu < mu; nu++) {
				P+=getLoop(n,mu,nu,1,1);			}
		}
		
		vector<double> v = {};
		v.push_back(P.trace().real() / 6.0);
		return v;
	};
	
	return measureObservable(plaq)[0];

}

vector<double> lattice::getPolyakovLoop(int rmax){
	
	auto impFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 0; r <= rmax; r++){
			double sum = 0;
			for(int mu = 1; mu <4; mu++){
				position m = shift(n,mu,r);
				double p1 = getLine(n,0,T).trace().real();
				double p2 = getLine(m,0,T).adjoint().trace().real();
				sum+=p1*p2;
			}
			data.push_back(sum/3);
		}
		
		return data;
	};
	
	return measureObservable(impFunc);
	
}

vector<double> lattice::getSquareLoop(int rmax){
	
	auto impFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 0; r <= rmax; r++){
			double sum = 0;
			for(int mu = 0; mu <4; mu++){
				for(int nu =0; nu < mu; nu++){
					sum+=getLoop(n,mu,nu,r,r).trace().real();
				}
			}
			data.push_back(sum/6);
		}
		
		return data;
	};
	
	return measureObservable(impFunc);
	
}



