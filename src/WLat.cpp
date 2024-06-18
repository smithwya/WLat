
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
    string infoname = (string)argv[9]+"/config.txt";
    string suffix = argv[11];
    int nMeasurements = atoi(argv[12]);
    int sweeps_per_meas = atoi(argv[13]);
	int Nc = atoi(argv[14]);


	auto begin = std::chrono::steady_clock::now();
	lattice lat = lattice(Nc,N,T,beta,xi_R);
	if(hotStart==1){
		lat.readFile(configname);
	}
	lat.update(nSweeps);
	auto mid = std::chrono::steady_clock::now();
	
	

	for(int i = 1; i <nMeasurements; i++){
		if(gTol!=0){
			if(Nc==3){
				lat.fixCoulombGaugeSubgroups(gTol);
			}
			if(Nc==2){
				lat.fixCoulombGauge(gTol);
			}
		}

		configname = (string)argv[9]+"/run"+to_string(i+1)+".txt;
		lat.update(sweeps_per_meas);
		lat.writeFile(configname);
	}

	auto end = std::chrono::steady_clock::now();
	
	if(jobnum>1) return 0;

	std::ofstream infofile;
	
	if(hotStart==1){
		infofile.open(infoname, std::ios::app);
		} else{
		infofile.open(infoname,std::ios::out);
		}/*
	infofile <<"Thermalized lattice with "<<nSweeps<<" sweeps in"<<std::chrono::duration_cast<std::chrono::minutes>(mid - begin).count() << " min" <<endl;
	infofile << "Measurement name:"<<suffix<<". Finished "<<nMeasurements<<" measurements separated by "<<sweeps_per_meas<<" sweeps and gauge fixing to dF<"<<gTol<<" in "<<std::chrono::duration_cast<std::chrono::minutes>(end - mid).count() << " min" <<endl;*/
	infofile <<"Split into "<<nMeasurements<<" streams separated by "<<sweeps_per_meas<<" sweeps"<<endl;
	infofile.close(); 
	 
    return 0;
    
}
