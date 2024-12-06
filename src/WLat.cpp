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
    string cnfg = argv[1];
    string dat = argv[2];
    int cnfgStart = atoi(argv[3]);
    int cnfgEnd = atoi(argv[4]);
    int cnfgInc = atoi(argv[5]);
    
    int Nc = atoi(argv[6]);
    double beta = atof(argv[7]);
    int L = atoi(argv[8]);
    int T = atoi(argv[9]);
    double xi_R = atof(argv[10]);
    
    bool multiGen = atoi(argv[11]);
    int nTherm = atoi(argv[12]);
    double gTol = atof(argv[13]);
	bool measure = atoi(argv[14]);
	int Rmax = atoi(argv[15]);
	int Tmax = atoi(argv[16]);

	cout<<"[WLAT] Initializing"<<endl;
	cout<<"Nc = "<<Nc<<endl;
	cout<<"Beta = "<<beta<<endl;
	cout<<"L = "<<L<<endl;
	cout<<"T = "<<T<<endl;
	cout<<"xi_R = "<<xi_R<<endl;
	
	auto t_s = high_resolution_clock::now();
	lattice lat = lattice(Nc,L,T,beta,xi_R);
	auto t_e = high_resolution_clock::now();

	cout<<"[WLAT] Lattice object initialized ("<< duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
	
	if(multiGen){
		cout<<"[WLAT] Generating initial configurations ONLY"<<endl;
		 
		cout<<"[UPDATE] Thermalizing initial configuration from Cold Start with "<<nTherm<<" heatbath sweeps";
		t_s = high_resolution_clock::now();
		lat.ColdStart();
		lat.update(nTherm);
		t_e = high_resolution_clock::now();
		cout<<" ...done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
		
		cout<<"[UPDATE] Generating configurations in range ["<<cnfgStart<<", "<<cnfgEnd<<", "<<cnfgInc<<"] ";
		
		for(int i = cnfgStart; i <=cnfgEnd; i=i+cnfgInc){
		
			string fname=cnfg+"."+std::to_string(i);
			lat.update(cnfgInc);
			lat.writeFile(fname);
			cout<<"     Saving configuration "<<fname<<endl;
			
		}
		
		t_e = high_resolution_clock::now();
		cout<<" ...done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
		cout<<"[WLAT] Exiting"<<endl;
		return 0;
	}
	
	
	
	for(int i = cnfgStart; i <=cnfgEnd; i=i+cnfgInc){
	
		string fname=cnfg+"."+std::to_string(i);
		
		
		t_s=high_resolution_clock::now();
		lat.readFile(fname);
		t_e=high_resolution_clock::now();
		cout<<" ...done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
		if(gTol != 0){
			cout<<"[GFIX] Coulomb Gauge fixing to df<"<<gTol;
			
			t_s = high_resolution_clock::now();
		
			if(Nc==3){
				lat.fixCoulombGaugeSubgroups(gTol);
			}
			if(Nc==2){
				lat.fixCoulombGauge(gTol);
			}

			lat.writeFile(fname);
			t_e = high_resolution_clock::now();
			cout<<" ...done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";
		}
		
		
		if(measure){
		cout<<"[MEASURE] Measuring G(R,T) with R_max = "<<Rmax<<" T_max="<<Tmax;
		
			vector<vector<double> > dat_pt_list = {};
			
			t_s = high_resolution_clock::now();
			for(int j = 1; j <=Tmax; j++){
				dat_pt_list.push_back(lat.getImprovedCorrelator(Rmax,j));
			}
			t_e = high_resolution_clock::now();
			
			cout<<" ...done (" << duration_cast<seconds>(t_e-t_s).count()<<"s)\n";

			for(vector<double> dl : dat_pt_list){
				saveData(dat,dl);
			}
			cout<<"[MEASURE] Data saved to "<<dat<<endl;

		}
	
	}
		

    return 0;

}
