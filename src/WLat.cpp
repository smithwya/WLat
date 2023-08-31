
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "lattice.h"
#include <chrono>
using namespace std;



void saveData(string dataname, vector<double> v) {
    std::ofstream datfile(dataname,ios::app);
    for (double x : v) {
        datfile << x << " ";
    }

    datfile << endl;
    datfile.close();
    return;
}



int main(int argc, char ** argv)
{
	//code only works for SU(2) currently
    int Nc = 2;
    
    //takes in arguments from command line
    int jobnum = atoi(argv[1]);
    double beta = atof(argv[2]);
    string betaname = argv[2];
    int N = atoi(argv[3]);
    int T = atoi(argv[4]);
    double xi_R = atof(argv[5]);
    int hotStart = atoi(argv[6]);
    int nSweeps = atoi(argv[7]);
    double gTol = atof(argv[8]);
    string configname = (string)argv[9]+"/run"+to_string(jobnum)+".txt";
    string dataname = (string)argv[10]+"/run"+to_string(jobnum);
    string infoname = (string)argv[10]+"/config.txt";

	//start the clock
	auto begin = std::chrono::steady_clock::now();

	//make lattice object
	lattice lat = lattice(Nc,N,T,beta,xi_R);
	
	//read in previously saved config
	if(hotStart==1){
		lat.readFile(configname);
	}
	
	
	lat.update(nSweeps);
	//get time elapsed for thermalization
	auto mid = std::chrono::steady_clock::now();
	
	for(int i = 0; i <20; i++){
		lat.update(100);
		//gauge fix
		if(gTol!=0){
			lat.fixCoulombGauge(gTol);
		}
		
		//measure data and save to file
		//renormalization
		saveData(dataname+"greenbake",lat.getGreensiteCorrelator(8));
		//saveData(dataname+"WRR",lat.getSquareLoop(8));
		
	}
	//save the file
	lat.writeFile(configname);
	
	auto end = std::chrono::steady_clock::now();
	
	if(jobnum>1) return 0;
	
	//print some stuff to config file
	std::ofstream infofile;
	
	if(hotStart==1){
		infofile.open(infoname, std::ios::app);
		} else{
		infofile.open(infoname,std::ios::out);
		}
	infofile <<"Thermalized lattice with "<<nSweeps<<" sweeps in"<<std::chrono::duration_cast<std::chrono::minutes>(mid - begin).count() << " min" <<endl;
	infofile << "Finished 10 measurements separated by 100 sweeps and gauge fixing to dF<"<<gTol<<" in "<<std::chrono::duration_cast<std::chrono::minutes>(end - mid).count() << " min" <<endl;
	infofile.close();
    return 0;
}
