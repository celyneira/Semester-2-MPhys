#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1D.h>

#include "Garfield/MediumSiliconCarbide.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"

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

    std::cerr << "Unable to open file";

    return 1;

  }

  double value1;

  while (inputFile >> value1) {
      col.push_back(value1);

  }
  inputFile.close();

  return 0; 
}
int main(int argc, char * argv[]) {

  constexpr unsigned int nEvents = 1;
    
  std::vector<double> variable;
  read_txt("/Users/celyn1/Documents/MPhys/sil_sens/linspace.txt", variable);
  TApplication app("app", &argc, argv);
  constexpr bool plotVelocity = true;
  constexpr bool plotSignal = true;
  constexpr bool plotDrift = true;
  constexpr bool plotField = true;
  constexpr bool plotWeightingField = true;
  constexpr bool writeStoppingPower = false;

  constexpr bool includeDiffusion = true;

  // Particle type
  const std::string particle{"ELECTRON"};

  // Particle momentum [eV]
  const std::vector<double> particle_momentum{0.5e8, 0.5e8, 0.5e9, 0.5e10};
  constexpr int nMomt = 1;

  // Sensor thickness [cm]
  const std::vector<double> thickness{150.e-4, 200e-4, 250e-4,300e-4};
  constexpr int nThicc = 1;

  // Sensor temperature [K]
  //const std::vector<double> variable_temperature{100, 200, 300, 400, 500};
  const std::vector<double> variable_temperature = {300};
  constexpr int nTemps = 1;

  //Bias voltage [V]
  const std::vector<double> variable_v_bias{-50, -100, -75};
  constexpr int nVoltages = 1;

  //Strip thickness [cm]
  const std::vector<double> variable_strip{600e-4, 150e-4};
  //, 100e-4, 150e-4, 200e-4, 250e-4, 300e-4};
  constexpr int nStrip_T = 1;

  //Strip centre-point
  const std::vector<double> strip_centre{0, -100e-4, -50e-4, 0, 150e-4, 100e-4, 50e-4};
  const int nCentres = 1;

  // Flag to save the signal to a file.
  constexpr bool writeSignal = true;

  for (unsigned int centres = 0; centres < nCentres; ++centres) {
  for (unsigned int volt = 0; volt < nVoltages; ++volt) {
  for (unsigned int momt = 0; momt < nMomt; ++momt) {
    std::cout<<momt<<std::endl<<nMomt<<std::endl;

    double momentum = particle_momentum[momt];

    std::cout << "The momentum is " << momentum << " with index " << momt << "\n";

    for (unsigned int strip_thicc = 0; strip_thicc < nStrip_T; ++ strip_thicc){
    for (unsigned int thicc = 0; thicc < nThicc; ++thicc) {
      
      double d = thickness[thicc];

      std::cout << "The thickness is " << thickness[thicc] << " with index " << thicc << "\n";


      for (unsigned int temp = 0; temp < nTemps; ++temp) {
        std::cout<<variable_temperature[temp];
        // Define the medium.
        MediumSiC si;
        si.SetTemperature(variable_temperature[temp]);

        std::cout << "The temperature is: " << variable_temperature[temp] << " with index " << temp << "\n";
        char doping_type;
        double doping_val;
        
        for (unsigned int i = 0; i < nEvents; ++i) {

          // double e0 = photon_energy[0];
          double stopping_power;

          // Make a plot of the drift velocities.
          plottingEngine.SetDefaultStyle();
          if (plotVelocity) {
            si.PlotVelocity("eh", new TCanvas("cM", "", 600, 600));
          }
        
          // Make a component with constant drift field and weighting field.
          // Bias voltage [V]
          double vbias = variable_v_bias[volt];
          ComponentConstant uniformField;
          uniformField.SetArea(-2 * d, 0., -2 * d, 2 * d, d, 2 * d);
          uniformField.SetMedium(&si);
          uniformField.SetElectricField(0, vbias / d, 0);
          uniformField.SetWeightingField(0, -1. / d, 0, "");

          // Depletion voltage [V]
          constexpr double vdep = 0.;
          // Make a component with linear drift field.
          auto eLinear = [d,vbias,vdep](const double /*x*/, const double y, 
                                        const double /*z*/,
                                        double& ex, double& ey, double& ez) {
            ex = ez = 0.;
            ey = (vbias - vdep) / d ; 
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
          wField.AddStripOnPlaneY('z', d, strip_centre[centres]-variable_strip[strip_thicc]/2,  strip_centre[centres]+variable_strip[strip_thicc]/2, "strip");

          // Create a sensor. 
          Sensor sensor;
          sensor.AddComponent(&linearField); 
          const std::string label = "strip";
          sensor.AddElectrode(&wField, label);

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
          const double tmin =  0.;
          const double tmax = 10.;
          const double tstep = (tmax - tmin) / nTimeBins;
          sensor.SetTimeWindow(tmin, tstep, nTimeBins);
          si.GetDoping(doping_type, doping_val);
          std::cout<<doping_type<<" = type \n" << doping_val << "= val\n";
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
            if (i % 10 == 0) std::cout << i << "/" << nEvents << "\n"; 
            // Simulate a charged-particle track.
            double xt = 0.;
            // if (smearx) xt = -0.5 * pitch + RndmUniform() * pitch;
            track.NewTrack(xt, 0, 0, 0, 0, 1, 0);
            
            // Drift an electron hole pair
            drift.DriftElectron(0, d/2, 0,0);
            drift.DriftHole(0, d/2, 0,0);

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
              char filename[50];
              sprintf(filename, "results/signal_%04d_%.0f_%.4f_%.0f_%.0fV_%.4f_%.4f.txt", i, 
                variable_temperature[temp], thickness[thicc], particle_momentum[momt], variable_v_bias[volt], variable_strip[strip_thicc], strip_centre[centres]);
              std::ofstream outfile;
              outfile.open(filename, std::ios::out);
              for (unsigned int j = 0; j < nTimeBins; ++j) {
                const double t = (j + 0.5) * tstep;
                const double f = sensor.GetSignal(label, j);
                const double fe = sensor.GetElectronSignal(label, j);
                const double fh = sensor.GetIonSignal(label, j);
                outfile << t << "  " << f << "  " << fe << "  " << fh << "\n";
              
              }
              outfile.close();
            }

            if (writeStoppingPower) {
              char s_filename[50];
              sprintf(s_filename, "stopping_power_%02d.txt", i);
              std::ofstream outfile_2;
              outfile_2.open(s_filename, std::ios::out);
              for (unsigned int j = 0; j < nTimeBins; ++j) {
                outfile_2 << stopping_power << std::endl;
              }
              outfile_2.close();

            }
        }
        
      }
    }
  }
  }
  }}
  if (plotVelocity || plotSignal || plotDrift || 
            plotField || plotWeightingField) {
          app.Run();
            }
  
}
