#include "lattice.h"
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <ctime>
#include <iomanip>

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
	xi_R = 1;
	xi_0 = 1;
}

//constructor for lattice
lattice::lattice(int nColors, int spatialExtent, int timeExtent, double coupling, double anisotropy)
{
	Nc = nColors;
	N = spatialExtent;
	T = timeExtent;
	beta = coupling;
	xi_R = anisotropy;
	xi_0 = getBareAnisotropy(xi_R);
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
double lattice::getBareAnisotropy(double xiR)
{

	if (xiR == 1.0) return 1.0;

	if(Nc==2){
		double c0 = 0.3398;
		double c1 = -0.3068;
		double c2 = -0.0449;

		double eta1 = c0 + c1 / xiR + c2 / (pow(xiR, 3));

		double a1 = -1.0238;
		double a2 = -1.8921;

		return xiR / (1 + (4.0 / 6.0) * (1 + a1 / beta) / (1 + a2 / beta) * (eta1 / beta));
	}
	if(Nc==3){
		//Klassen parameterization
		double eta1 = (1.002503 + 0.391*pow(xiR,-1) + 1.47130*pow(xiR,-2) - 0.19231*pow(xiR,-3))/(1 + 0.26287*pow(xiR,-1) + 1.59008*pow(xiR,-2) - 0.18224*pow(xiR,-3));
		double eta = 1+(1-1/xiR)*(eta1/beta)*((beta-6.0*0.55055)/(beta-0.77810*6.0));
		return xiR/eta;
	}
	
	return 1.0;
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
		
	
	auto func2 = [&](position n, int mu) {
		comp a, b;
		infile >> a;
		infile >> b;
		MatrixXcd U {{a,b},{-conj(b),conj(a)}};
		return U;
	};
	
	
		
	auto func3 = [&](position n, int mu) {
		comp a, b, c, d, e, f;
		infile >> a;
		infile >> b;
		infile >> c;
		infile >> d;
		infile >> e;
		infile >> f;
		MatrixXcd U {{a, b, c}, {d, e, f}, {comp(0, 0), comp(0, 0), comp(0, 0)}};
		return reunitarizeSU3(U);
	};
		
	
	if(Nc==2){
		transformLat(func2);
	}
	if(Nc==3){
		transformLat(func3);
	}
	
	infile.close();
	return;
}

//writes a text file for the current config
void lattice::writeFile(string filename)
{
	std::ofstream outfile(filename);
	outfile<<std::fixed<<setprecision(std::numeric_limits<double>::digits10 + 1)<<endl;


	for (int i = 0; i < T; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				for (int l = 0; l < N; l++) {

					position n = position(i, j, k, l);

					for (int mu = 0; mu < 4; mu++) {

						MatrixXcd U = getLink(n, mu, true);
						if(Nc==2){
							outfile << U(0, 0) << " " << U(0, 1) << endl;
						}
						if(Nc==3){
							
							outfile<< U(0,0) << " " << U(0,1) <<" "<<U(0,2) <<" "<< U(1,0)<<" "<<U(1,1)<<" "<<U(1,2)<<endl;
						}

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
MatrixXcd lattice::genLinkSU2(MatrixXcd A)
{

	MatrixXcd Ap = A-A.adjoint() + (A.adjoint().trace())*MatrixXcd::Identity(2,2);
	double q = 0.5*sqrt(Ap.determinant()).real();
	Ap = Ap/(2.0*q);
	double betaprime = beta*q/Nc;
	double alpha = 2*betaprime;


	double R, Rp, Rpp, Rppp;
	double X, Xp, B, C;
	double deltabar;

	while (true) {
		R = 0.5 + randNum();
		Rp = 0.5 + randNum();
		
		X = -(log(R)/alpha);
		Xp =-(log(Rp)/alpha); 
		
		Rpp = 0.5 + randNum();
		C=pow(cos(2*3.14159265*Rpp),2);
		
		B = X*C;
		
		deltabar = Xp + B;
		
		Rppp = 0.5+randNum();
		
		if(pow(Rppp,2)<= 1-0.5*deltabar) break;
	}
	double x0 = 1-deltabar;
	
	double r1, r2, r3;
	while (true) {
		r1 = 2 * randNum();
		r2 = 2 * randNum();
		r3 = 2 * randNum();
		if (r1 * r1 + r2 * r2 + r3 * r3 <= 1) break;
	}

	double length = sqrt(pow(r1, 2) + pow(r2, 2) + pow(r3, 2));
	double norm = sqrt(1 - pow(x0, 2));

	double x1 = r1 * norm / length;
	double x2 = r2 * norm / length;
	double x3 = r3 * norm / length;
	

	MatrixXcd U{{comp(x0, x3), comp(x2, x1)},{comp(-x2, x1), comp(x0, -x3)}};

	return reunitarizeSU2(U * (Ap.adjoint()));
}

MatrixXcd lattice::getSubgroup(MatrixXcd U, int i){
	switch (i){
		case 1: return MatrixXcd{{U(0,0),U(0,1)},{U(1,0),U(1,1)}};
		case 2: return MatrixXcd{{U(0,0),U(0,2)},{U(2,0),U(2,2)}};
		case 3: return MatrixXcd{{U(1,1),U(1,2)},{U(2,1),U(2,2)}};
		
	}
	
	return MatrixXcd::Identity(2,2);
}

MatrixXcd lattice::SU3mat(MatrixXcd U, int i){
	switch(i){
		case 1: return MatrixXcd{{U(0,0),U(0,1),0},{U(1,0),U(1,1),0},{0,0,1}};
		case 2: return MatrixXcd{{U(0,0),0,U(0,1)},{0,1,0},{U(1,0),0,U(1,1)}};
		case 3: return MatrixXcd{{1,0,0},{0,U(0,0),U(0,1)},{0,U(1,0),U(1,1)}};
	}
	
	return MatrixXcd::Identity(3,3);
}


//generates a new SU(3) link using SU(2) subgroups
MatrixXcd lattice::genLinkSU3(position n, int mu){
	
	MatrixXcd U = getLink(n,mu,true);

	MatrixXcd W = U*genStaple(n,mu);
	
	MatrixXcd R2 = genLinkSU2(getSubgroup(W,1));
	
	MatrixXcd R = SU3mat(R2,1);
	
	W = R*W;
	MatrixXcd S2 = genLinkSU2(getSubgroup(W,2));
	MatrixXcd S = SU3mat(S2,2);
	
	W = S*W;
	MatrixXcd T2 = genLinkSU2(getSubgroup(W,3));
	MatrixXcd T = SU3mat(T2,3);
	
	return reunitarizeSU3(T*S*R*U);
}

//reunitarizes SU(2) matrix
MatrixXcd lattice::reunitarizeSU2(MatrixXcd U)
{
	MatrixXcd newU = MatrixXcd::Zero(2, 2);
	comp x = sqrt(U(0,0)*conj(U(0,0))+U(0,1)*conj(U(0,1)));
	
	newU << U(0, 0), U(0, 1), -conj(U(0, 1)), conj(U(0, 0));
	return newU/x;
}
//reunitarizes SU(3) matrix
MatrixXcd lattice::reunitarizeSU3(MatrixXcd U)
{
	Eigen::Vector3cd u(U(0,0),U(0,1),U(0,2));

	u=u/sqrt(u.dot(u));

	
	Eigen::Vector3cd v(U(1,0),U(1,1),U(1,2));
	v=v-(u.adjoint()*v)*u;
	v=v/sqrt(v.dot(v));

	Eigen::Vector3cd w = (u.conjugate()).cross(v.conjugate());
	MatrixXcd Up = MatrixXcd::Identity(3,3);
	Up.row(0) = u;
	Up.row(1) = v;
	Up.row(2) = w.conjugate();
	return Up;
}

//thermalizes the lattice by performing numUpdates sweeps of the heatbath algorithm
void lattice::update(int numUpdates)
{


	auto func2 = [&](position n, int mu) {
			MatrixXcd A = genStaple(n, mu);
			return genLinkSU2(A);
		};
		
	auto func3 = [&](position n, int mu) {
			return genLinkSU3(n,mu);
		};
		

	for (int i = 0; i < numUpdates; i++) {
		if(Nc==2){
			transformLat(func2);

			
		}
		if(Nc==3){
			transformLat(func3);

		}
	}
	
	return;
}

//returns the omega'th power of matrix G, up to cutoff terms in a series
MatrixXcd lattice::matrixPow(MatrixXcd G, double omega, int cutoff) {
	int dim = G.rows();

	MatrixXcd Gw = MatrixXcd::Zero(dim,dim);
	MatrixXcd GminusI = G-MatrixXcd::Identity(dim,dim);

	for (int i = 0; i <= cutoff; i++) {

		MatrixXcd temp = MatrixXcd::Identity(dim, dim);

		for (int j = 1; j <= i; j++) {
			temp = GminusI*temp;
		}

		double c = tgamma(omega + 1) / (tgamma(omega + 1 - i) * tgamma(i + 1));

		Gw = Gw+c*temp;
	}
	if(dim==2) return reunitarizeSU2(Gw);
	return reunitarizeSU3(Gw);
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

//fixes the coulomb gauge up to tolerance dF<=tolerance
void lattice::fixCoulombGaugeSubgroups(double tolerance) {
	double omega = 1.75;
	
	auto gsweep = [&](position n) {

		//for each of the 3 su2 subgroups
		for(int s =1; s<=3; s++){
			MatrixXcd K = MatrixXcd::Zero(2, 2);
			//construct K for the subgrp (1-3 = Coulombg, 0-3 = Landau)
			for (int dir = 1; dir < 4; dir++) {
				K += getSubgroup(getLink(n,dir,true),s);
				K += getSubgroup(getLink(shift(n, dir, -1), dir, true).adjoint(),s);
			}
			MatrixXcd G = K.adjoint();
			comp a = G.determinant();
			G = G / sqrt(a);
			//overrelaxation by replacing gauge transf g with g^omega
			G = matrixPow(G, omega, 2);
			
			MatrixXcd G3 = SU3mat(G,s);
			
			
			for (int mu = 0; mu < 4; mu++) {
				configuration[n.t][n.x][n.y][n.z][mu] = reunitarizeSU3(G3 * getLink(n, mu, true));

				position m = shift(n, mu, -1);
				configuration[m.t][m.x][m.y][m.z][mu] = reunitarizeSU3(getLink(m, mu, true) * (G3.adjoint()));
			}
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
				sum+=getCoulombLoop(n,mu,r,1).trace().real()/Nc;
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
		v.push_back(P.trace().real() / (6.0*Nc));
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

vector<double> lattice::getAverageWilsonLoop(int rmax, int t){
	
	auto impFunc = [&](position n){
		vector<double> data = {};
		
		for(int r = 1; r <= rmax; r++){
			double sum = 0;
			
			for(int mu = 0; mu <4; mu++){
				for(int nu=0; nu<mu; nu++){
					sum+=getLoop(n,mu,nu,r,t).trace().real();
				}
			}
			data.push_back(sum/(Nc*6));
		}
		
		return data;
	};
	
	return measureObservable(impFunc);
	
}

