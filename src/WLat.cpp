#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <fstream>
#include "lattice.h"
#include <chrono>
#include <ctime>
using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::seconds;


void saveData(string dataname, vector<double> v) {
    std::ofstream datfile(dataname,ios::app);
	datfile<<std::fixed<<setprecision(std::numeric_limits<double>::digits10 + 1);
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

	cout<<"[INIT] Starting job "<<jobnum<<endl;
	cout<<"Beta = "<<beta<<endl;
	cout<<"L = "<<N<<endl;
	cout<<"T = "<<T<<endl;
	cout<<"xi_R = "<<xi_R<<endl;

	//start the clock
	auto t_s = high_resolution_clock::now();

	//make lattice object
	lattice lat = lattice(Nc,N,T,beta,xi_R);
	auto t_e = high_resolution_clock::now();
	cout<<"[INIT] Lattice initialized (" << duration_cast<seconds>(t_e-t_s).count()<< "s)\n";
	//read in previously saved config
	if(hotStart==1){
		cout<<"[IO] Loading config: "<<configname<<endl;
		t_s=high_resolution_clock::now();
		lat.readFile(configname);
		t_e=high_resolution_clock::now();
		cout<<"...done ("<< duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
	}
	if(nSweeps>0){
		cout<<"[UPDATE] Thermalizing "<<nSweeps<<" times...";
		t_s = high_resolution_clock::now();
		lat.update(nSweeps);
		t_e = high_resolution_clock::now();
		cout<<"done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
		t_s = high_resolution_clock::now();

		cout<<"[IO] Saving thermalized  lattice to "<<configname<<endl;
		lat.writeFile(configname);
		t_e = high_resolution_clock::now();
		cout<<"done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
	}
	cout<<"[MEASURE] Performing "<<nMeasurements<<" measurements of G(R,T)\n";
	//measure data and save to file

	for(int i = 1; i <= nMeasurements; i++){
		cout<<"[MEASURE] Beginning measurement #"<<i<<endl;
		vector<vector<double>> dat_pt_list = {};
		t_s = high_resolution_clock::now();
		lat.update(sweeps_per_meas);
		t_e = high_resolution_clock::now();
		cout<<"Performed "<< sweeps_per_meas<<" updates (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";

		if(gTol!=0){
			t_s = high_resolution_clock::now();
			cout<<"Gauge fixing to dF<"<<gTol<<"...";
			if(Nc==3){
				lat.fixCoulombGaugeSubgroups(gTol);
			}
			if(Nc==2){
				lat.fixCoulombGauge(gTol);
			}
			t_e = high_resolution_clock::now();
			cout<<"done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
		}
		t_s = high_resolution_clock::now();
		cout<<"Measuring G(R,T) for R,T<="<<N/2<<"...";
		for(int j = 1; j <=N/2; j++){
			//saveData(dataname+suffix,lat.getAverageTSLoop(N/2,xi_R*j));
			dat_pt_list.push_back(lat.getImprovedCorrelator(N/2,j));
		}
		t_e = high_resolution_clock::now();
		cout<<"done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";

		for(vector<double> dl : dat_pt_list){
		saveData(dataname+suffix,dl);
		}

		configname=(string)argv[9]+"/run"+to_string(i+jobnum)+".txt";
		t_s = high_resolution_clock::now();
		cout<<"Saving configuration to "<<configname<<endl;
		lat.writeFile(configname);
		cout<<"done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";

	}

    return 0;

}
