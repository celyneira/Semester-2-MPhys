#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1D.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/MediumSiliconCarbide.hh"

#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Plotting.hh"

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int read_txt(std::string filename, std::vector<double>& col){
  std::ifstream inputFile(filename); // Replace "data.csv" with your actual file name
  
  std::string line, cell;
  if (!inputFile) {

    std::cerr << "Unable to open file: "<< filename<<std::endl;

    return 1;

  }

  double value1;

  while (inputFile >> value1) {
      col.push_back(value1);

  }
  inputFile.close();

  return 0; 
}

int read_csv(std::string filename, std::vector<double>& column1, std::vector<double>& column2){
  std::ifstream inputFile(filename); // Replace "data.csv" with your actual file name
  
  std::string line, cell;
  if (!inputFile) {

    std::cerr << "Unable to open file: "<< filename << std::endl;

    return 1;

  }

  double value1, value2;

  while (inputFile >> value1 >> value2) {
      column1.push_back(value1);

      column2.push_back(value2);

  }
  inputFile.close();

  return 0; 
}



std::vector<std::vector<double>>  read_table(std::string filename)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

  std::string line;
  std::vector<std::vector<double>> data;

  while(getline(file, line)){
    std::stringstream ss(line);
    std::vector<double> row;
    double value;
    while (ss >> value){
        
      row.push_back(value);
    }
    data.push_back(row);
  }
  file.close();
  return data;
}

double look_up(int z_indice, int r_indice, std::vector<std::vector<double>> data)
{
  if (z_indice <= 0 || z_indice > data.size() || r_indice <= 0 || r_indice > data[0].size()) {
        std::cout<<"z size "<<data.size()<<std::endl; 
        std::cout<<"z "<<z_indice<<" r "<< r_indice<<std::endl;
        std::cerr << "Invalid row or column specified." << std::endl;
        return 1;
    }
  
  double number = data[z_indice - 1][r_indice - 1];
  
  return number;
}

int main(int argc, char * argv[]) {

  constexpr unsigned int nEvents = 1;
  
  std::vector<double> voxel_position_vals;
  
  read_txt("/Users/celyn1/Documents/Sem2MPhys/Missing_charges/linspace.txt", voxel_position_vals);
  std::vector<double> time_step, transI;
  read_csv("/Users/celyn1/Documents/Sem2MPhys/GarfieldSim_SiC/impulse_response.txt", time_step, transI);
  //std::cout<<time_step.size();

  
  TApplication app("app", &argc, argv);
  constexpr bool plotVelocity = false;
  constexpr bool plotSignal = true;
  constexpr bool plotDrift = false;
  constexpr bool plotField = false;
  constexpr bool plotWeightingField = false;
  constexpr bool writeDeleted = false;

  constexpr bool includeDiffusion = true;
  constexpr bool addNoise = false;
  constexpr bool transfer_func = true;
  constexpr bool custom_file_name = false;
  constexpr bool voltage_signal = false;
  // Particle type
  const std::string particle{"ELECTRON"};

  //file name
  std::string filename = "build/temp_data/test.txt" ;
  ;
  //voxel information
  double voxel_z_max = 50e-4; //cm
  double voxel_r_max = 25e-4;
  
  int voxel_z_refinement = 500;
  int voxel_r_refinement = 251;

  // Sensor thickness [cm]
  const double thickness = 50e-4;
  // Sensor temperature [K]
  const double temperature{300};
  //Bias voltage [V]
  const double v_bias{-200};

  const double sigma {0};
  //extracting normalisation constants
  std::vector<double> normalisations;
  read_txt("/Users/celyn1/Documents/Sem2MPhys/GarfieldSim_SiC/voxels/normalisation/normalisation_constants_0.000000e+00.txt", normalisations);

  
  //number of voxel positions
  std::vector<double> voxel_position_values;
  read_txt("/Users/celyn1/Documents/Sem2MPhys/GarfieldSim_SiC/voxels/positions/voxel_positions_0.000000e+00.txt", voxel_position_values);
  std::cout << "Voxel position values:" << std::endl;
  for (const auto& value : voxel_position_values) {
    std::cout << value << std::endl;
}
  const int n_positions = 1;


  //Voxel position
  const std::vector<double> voxel_locs =voxel_position_values;


  // Flag to save the signal to a file.
  constexpr bool writeSignal = true;
  for (unsigned int i = 0; i < nEvents; ++i) {
  for (unsigned int locs = 0; locs < n_positions; ++locs) {
    double momentum = 1e6;
    
    //getting the correct voxel mesh
    char voxel_data[50];
    std::cout<<"Opening position: "<< voxel_position_values[locs]<<std::endl;
    sprintf(voxel_data, "/Users/celyn1/Documents/Sem2MPhys/GarfieldSim_SiC/voxels/charge_carrier_mesh_%e_%e.txt", voxel_position_values[locs],sigma);

    std::cout<<voxel_data<<std::endl;
    std::vector<std::vector<double>> data = read_table(voxel_data);

    double Normalisation = normalisations[locs];
    double NO_GENERATED =   10000; 
      double d = thickness;
      
        // Define the medium.
        MediumSiC si;
        si.SetTemperature(temperature);

        //std::cout << "The temperature is: " << variable_temperature[temp] << " with index " << temp << "\n";
        char doping_type;
        double doping_val;
        
          double total_number_of_carriers = 0;
          double total_drifts = 0;
          // double e0 = photon_energy[0];
          double stopping_power;

          // Make a plot of the drift velocities.
          plottingEngine.SetDefaultStyle();
          if (plotVelocity) {
            si.PlotVelocity("eh", new TCanvas("cM", "", 600, 600));
          }
        
          // Make a component with constant drift field and weighting field.
          // Bias voltage [V]
          double vbias = v_bias;
          ComponentConstant uniformField;
          uniformField.SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
          uniformField.SetMedium(&si);
          uniformField.SetElectricField(0, vbias / d, 0);
          uniformField.SetWeightingField(0, -1. / d, 0, "");

          // Depletion voltage [V]
          constexpr double vdep = -0.;
          // Make a component with linear drift field.
          auto eLinear = [d,vbias,vdep](const double /*x*/, const double y, 
                                        const double /*z*/,
                                        double& ex, double& ey, double& ez) {
            ex = ez = 0.;
            ey = (vbias - vdep) / d; 
          };

          ComponentUser linearField;
          linearField.SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
          linearField.SetMedium(&si);
          linearField.SetElectricField(eLinear);

          // Make a component with analytic weighting field for a strip or pixel.
          ComponentAnalyticField wField;
          wField.SetMedium(&si);
          wField.AddPlaneY(0, vbias, "back");
          wField.AddPlaneY(d, 0, "front");
          wField.AddStripOnPlaneY('z', d, -2*d,  2*d, "strip");

          // Create a sensor. 
          Sensor sensor;
          sensor.AddComponent(&linearField); 
          const std::string label = "strip";
          sensor.AddElectrode(&wField, label);
          if (transfer_func)
            sensor.SetTransferFunction(time_step, transI);
          
          // Plot the drift field if requested.
          if (plotField) {
            ViewField* fieldView = new ViewField(&sensor); 
            fieldView->SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
            fieldView->PlotContour("ey");
          }
          // Plot the weighting potential if requested.
          if (plotWeightingField) {
            ViewField* wfieldView = new ViewField(&wField); 
            wfieldView->SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
            wfieldView->PlotContourWeightingField("strip", "v");
          }

          // Set the time bins.
          const unsigned int nTimeBins = 10000;
          const double tmin =  -16.5;
          const double tmax = 16.5;
          const double tstep = (tmax - tmin) / nTimeBins;
          sensor.SetTimeWindow(tmin, tstep, nTimeBins);
          si.GetDoping(doping_type, doping_val);
          //std::cout<<doping_type<<" = type \n" << doping_val << "= val\n";
          // Set up Heed with a photon
          TrackHeed track;
          track.SetSensor(&sensor);

          track.SetParticle(particle);
          track.SetMomentum(momentum);
          
          // Simulate electron/hole drift lines using MC integration.
          AvalancheMC drift(&sensor);
          
          // Use steps of 1 micron.
          if (!includeDiffusion)
            drift.DisableDiffusion();
          
          

          drift.SetDistanceSteps(1.e-4);
        
          // Plot the signal if requested.
          TCanvas* cSignal = nullptr;
          if (plotSignal) { 
            cSignal = new TCanvas("cSignal", "", 600, 600);
          }

          ViewDrift* driftView = nullptr;
          TCanvas* cDrift = nullptr;
          if (plotDrift) {
            cDrift = new TCanvas("cDrift", "", 600, 600);
            driftView = new ViewDrift();
            driftView->SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
            driftView->SetCanvas(cDrift);
            drift.EnablePlotting(driftView);
          }

          if (plotDrift) driftView->Clear();
            // Reset the signal.
            sensor.ClearSignal();
            if (i % 100 == 0) std::cout << i << "/" << nEvents << "\n"; 
            // Simulate a charged-particle track.
            double xt = 0.;
            // if (smearx) xt = -0.5 * pitch + RndmUniform() * pitch;
            track.NewTrack(xt, 0, 0, 0, 0, 1, 0);

          double voxel_location = voxel_locs[locs]*d;

          for (int z_ind = 1 ; z_ind < voxel_z_refinement; z_ind ++){
          for (int r_ind = 1 ; r_ind < voxel_r_refinement; r_ind ++){
            
            double number_charge_carriers = look_up(z_ind, r_ind, data); 
            double z_position = (z_ind - 0.5) * voxel_z_max  / voxel_z_refinement;
            
            if (z_position < 0) continue;
            if (z_position > d) continue;
            
            double r_position = (r_ind - 0.5) * voxel_r_max  / voxel_r_refinement - voxel_r_max/2;
            total_number_of_carriers += number_charge_carriers;
            number_charge_carriers = std::round(number_charge_carriers);
            for (int carrier_indice = 0; carrier_indice < number_charge_carriers; ++ carrier_indice)
            {
              total_drifts += 1;
              double the = RndmUniform()*2*Pi;
              drift.DriftElectron(r_position*cos(the),z_position,r_position*sin(the),0.);
              drift.DriftHole(r_position*cos(the),z_position,r_position*sin(the),0.);
              //if (carrier_indice % 100 == 0) {
              //std::cout << "Current index: " << carrier_indice << std::endl;}
            }
          }}
          if (addNoise)
            sensor.AddWhiteNoise(label, 0.009067864135419285/Normalisation*1e6*10000, true, 1);
         
          sensor.ConvoluteSignals();
            if (plotSignal) {
              sensor.PlotSignal(label, cSignal);
              cSignal->Update();
              gSystem->ProcessEvents();
            }
            if (plotDrift) {
              constexpr bool twod = true;
              driftView->Plot(twod);
              cDrift->Update();
              gSystem->ProcessEvents();
            }
          
            // Save the induced current signal to a file.
            if (writeSignal) {
            //if(!custom_file_name){
                char filename[50];
                sprintf(filename, "data/signal_SiC_%04d_%.0f_%.4f_%.0fV_%.4fm.txt", i, 
                  temperature, thickness, v_bias, voxel_locs[locs]*100);
              //  }
                std::ofstream outfile;
              outfile.open(filename, std::ios::out);
              if (voltage_signal){
                for (unsigned int j = 0; j < nTimeBins; ++j) {
                  const double t = (j + 0.5) * tstep;
                  const double f = sensor.GetSignal(label, j)/NO_GENERATED*50;
                  const double fe = sensor.GetElectronSignal(label, j)/NO_GENERATED*50;
                  const double fh = sensor.GetIonSignal(label, j)/NO_GENERATED*50;
                  outfile << t << "  " << f << "  " << fe << "  " << fh << "\n";
                }}
              if (!voltage_signal){
                for (unsigned int j = 0; j < nTimeBins; ++j) {
                    const double t = (j + 0.5) * tstep;
                    const double f = sensor.GetSignal(label, j)/NO_GENERATED;
                    const double fe = sensor.GetElectronSignal(label, j)/NO_GENERATED;
                    const double fh = sensor.GetIonSignal(label, j)/NO_GENERATED;
                    outfile << t << "  " << f << "  " << fe << "  " << fh << "\n";
                  }
              }
              outfile.close();
            }

            if (writeDeleted) {
              char s_filename[50];
              sprintf(s_filename, "data/drifted_deleted_%02d_%.4f.txt", i, voxel_locs[locs]*100);
              std::ofstream outfile_2;
              outfile_2.open(s_filename, std::ios::out);
              for (unsigned int j = 0; j < nTimeBins; ++j) {
                outfile_2 << total_number_of_carriers << "  "<<total_drifts<<"\n";
              }
              outfile_2.close();
            }
            std::cout<<"The total number of drifted e-h pairs is: "<<total_number_of_carriers<<std::endl;
            std::cout<<total_drifts << " pairs were drifted\n";
        
      }
    
  std::cout<<"Done!";
  if (plotVelocity || plotSignal || plotDrift || 
            plotField || plotWeightingField) {
          app.Run();
            } 
}}
