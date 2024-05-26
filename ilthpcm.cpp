#include "ilthpcm.h"

int main(int argc,char* argv[])
{
	/*Declarations*/
	int nt, noOfElements;
	Surface surface;
	vector<Layer> layers;
	vector<Weather> weatherData;
	vector<double> x,dx,alpha,T,Tnew;
	vector<double> a,b,c,d;
	double solarrad, qirr, qconv, qrad, albedo;
	double dt, xi = -0.1; //dt in seconds, xi in C/m
	ofstream fout;
	double fl0; //liquid fraction initial update
	double dfl; //liquid fraction derivative
	vector<double> fl, deltaH; //liquid fraction vector
	double Sc,Sp, underrelax_factor; // Source term in RHS and Sorce Term in LHS
	double k_eff; //Composite thermal conductivity
	double c_eff; //Composite heat capacity
	double rho_eff; //Composite density

	//Variables for MEPDG format solution
	ofstream fMEPDG;
	vector<double> xPCC, TPCC;
	int timestamp;

	/*Read data*/
	//Read input file
	readInputFile(nt, underrelax_factor, surface, layers, argv[1]);
	//Read weather data
	readWeatherData(weatherData, argv[2]);

	//Define mesh (x and dx or each element, also diffusivity of each element)
	noOfElements = 0;
	for(int i=0;i<layers.size();i++){
		noOfElements += layers[i].numLayerElements;
	}
	
	defineMesh(x, dx, layers, noOfElements);

	a.assign(noOfElements,0.0);
	b.assign(noOfElements,0.0);
	c.assign(noOfElements,0.0);
	d.assign(noOfElements,0.0);
	alpha.assign(noOfElements,0.0);
	deltaH.assign(noOfElements,0.0);
	

	dt = 3600.0/(double)nt;

	//Define stiffness matrix
	update_liquid_fraction(layers, T);
	update_thermal_properties(layers);
	assign_layer_to_element(alpha, deltaH, layers, noOfElements);
	stiffnessmat(a, b, c, x, dx, alpha, dt, noOfElements);

	//Create output file
	fout.open(argv[3],ios::trunc);
	fMEPDG.open("ThermalPCC_ILTH.dat",ios::trunc);
	fMEPDG << fixed << setprecision(1) << showpoint;
	
	//Write headers in output file
	fout << "Year" << "," << "Month" << "," << "Day" << "," << "Hour" << ",";
	for(int j=0;j<noOfElements;j++){
		fout << x[j]*1000.0 << ",";
	}
	fout << endl;

	//Initialize temperature field
	T.assign(noOfElements,weatherData[0].AirTemp);
	Tnew.assign(noOfElements,0.0);

	//Begin loop for each weather case 
	for(int i=0;i<weatherData.size();i++){
		cout << "Analyzing hour " << i+1 << " of " << weatherData.size() << endl;
		//Calculate solar energy ( assumed constant over the hour)
		solarrad = solar(weatherData[i]);

		//Iterate nt times to cover the hour
		for(int t=0;t<nt;t++){
			//Calculate surface energy balance
			/*Update fl*/
			/*Calculate alpha from fl*/
			/*Calculate a,b,c,d*/
			
			qirr = longwave(weatherData[i], T[0], surface.isothermalEmissivity);
			qconv = convection(weatherData[i], T[0]);

			if(surface.surfaceCode == 0)
			{
				albedo = surface.isothermalAlbedo;
			}
			else
			{
				albedo = thermochromic_albedo(surface, T[0]);
			}

			qrad = (solarrad*(1.0-albedo) + qirr + qconv)/(layers[0].solidMatrixDensity*layers[0].solidMatrixHeatCapacity);

			//Assume new liquid fraction to be as old liquid fraction - 0th iteration
                        fl = flnewfinal

			// Iterative process
			while (true) {	
				//Define RHS [d] - liquid fraction should be an input here
				rhsvector(d, T, x, dx, alpha, dt, qrad, xi, noOfElements,fl);
				
				//Solve the system of equations and store result in Tnew
				solve(Tnew, a, b, c, d, noOfElements);
				
				//Swap solution for 1st iteration of temperature
				T1 = Tnew;
	
				//Calculate the first iteration of liquid fraction 
				flnew = fl + under-relaxation*(update liquid fraction(layers,T))
				if (flnew > 1) {
	   			flnew = 1;
				}
				if (flnew < 0) {
				flnew = 0;
				}
			
				//Define RHS [d]  - liquid fraction should be an input here
				rhsvector(d, T, x, dx, alpha, dt, qrad, xi, noOfElements,flnew);
				
				//Solve the system of equations and store result in Tnew
				solve(Tnew, a, b, c, d, noOfElements);
				
				//Swap solution for 2nd iteration of temperature
				T2 = Tnew;

				//Update
				if(!maxDifferenceGreaterThan(T1, T2, 1e-4)){
				Tnewfinal = T2;
				flnewfinal = flnew;
				break;
				}

				else{
				T = T2;
				fl = flnew;
				}
			}
				
			//update_liquid_fraction(layers, T); This may not be needed here.
			update_thermal_properties(layers);
			assign_layer_to_element(alpha, deltaH, layers, noOfElements);
			stiffnessmat(a, b, c, x, dx, alpha, dt, noOfElements);
		}

		//Print current case to file
		fout << weatherData[i].Year << "," << weatherData[i].Month << "," << weatherData[i].Day << "," << weatherData[i].Hour << ","; 
		for(int j=0;j<noOfElements;j++){
			fout << T[j] << ",";
		}
		fout << endl;
		
		//Write to MEPDG format (top stabilizedLayer only)
		xPCC = vector<double>(x.begin(),x.begin()+layers[0].numLayerElements);
		TPCC = vector<double>(T.begin(),T.begin()+layers[0].numLayerElements);
		timestamp = weatherData[i].Year*1E6 + weatherData[i].Month*1E4 + weatherData[i].Day*1E2 + weatherData[i].Hour;
		
		WriteMEPDG(xPCC,TPCC,layers[0].numLayerElements,layers[0].thickness,11,2,fMEPDG,timestamp);
	}

	fout.close();
	fMEPDG.close();

	return 0;
}
