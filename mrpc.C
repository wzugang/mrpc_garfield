#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TApplication.h>

#include "ViewSignal.hh"
#include "ComponentAnalyticField.hh"
#include "MediumMagboltz.hh"
#include "SolidTube.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "TrackHeed.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewField.hh"
#include "ViewGeometry.hh"

using namespace std;
using namespace Garfield;

// Finds ion creation time based on distance from the electron cluster 
// to the wire (found by fitting the ion creation time of many avalanches)
double IonTiming(double dist) {

  double p0 =  1.49880e-13;   
  double p1 =  2.09250e+02;
  double p2 =  2.61998e+02;
  double p3 = -1.24766e+02;
  
  return p0 + p1 * dist + p2 * dist * dist + p3 * dist * dist * dist;

}

int main(int argc, char *argv[]) {

  TApplication app("app",&argc, argv);

  plottingEngine.SetDefaultStyle();
 
//  // Variables used to describe the geometry
//  const double dWire = 50.e-4;
//  const double rAnode = 1.46;

/*
  // Variables describing the signal histogram
  const double tmin = 0.;
  const double tstep = 1.;
  const int nTimeBins = 1000; 

  // Average ion creation point
  const double xIon = 0.00253;

  // Location of track
  double x0 = 1.2;
  double y0 = -(sqrt(rAnode * rAnode - x0 * x0) - 0.002);
  double z0 = 0.;
  double t0 = 0.;

  // Variables used to describe the avalanche
  double gain = 1190;
  double theta = 0.2654;
  */
  
  // Make a gas medium
  MediumMagboltz* gas = new MediumMagboltz();
  // Set the temperature [K] and pressure [Torr]
  const double pressure = 1 * AtmosphericPressure;
  const double temperature = 293.15;
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("C2F4H2", 95., "iC4H10", 4.5, "SF6", 0.5);
  gas->LoadGasFile("r134a_95_iso_4.5_sf6_0.5.gas"); //for electrons
  //gas->LoadIonMobility("/Users/chiu/software/garfield/Data/IonMobility_CO2+_CO2.txt"); // for ions
 
  // Build the geometry.
  GeometrySimple* geo = new GeometrySimple();

  // Dimensions of 1 gas gap [cm]
  const Double_t lx = 10/2.;   // half lengths
  const Double_t ly = 0.02/2.;   
  const Double_t lz = 10/2.;  

  SolidBox *tube = new SolidBox(0.,0.,0.,lx,ly,lz);

  // Add  the solid to the geometry, together with the medium inside
  geo->AddSolid(tube, gas);

  // Make a component with analytic electric field.
  ComponentAnalyticField* comp = new ComponentAnalyticField();
  comp->SetGeometry(geo);

  // Wire radius [cm] (100 um gold plate CuBe wire) and voltage
  const Double_t vPlane = 2000.; // volts
  const Double_t vGround =    0.;

  // Add the ends.
  comp->AddPlaneY(ly, vPlane, "t");
  comp->AddPlaneY(-ly, vGround, "g");

  comp->AddReadout("g");

  // Make a sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(comp);
  //sensor->AddElectrode(comp, "s");
  sensor->SetArea(-lx, -ly, -lz, lx, ly, lz);
  //sensor->SetTimeWindow(tmin, tstep, nTimeBins);
  //sensor->ClearSignal();
  std::cout << "Num Electrodes = " << sensor->GetNumberOfElectrodes() << std::endl;
 
  // Plot isopotential contours
  ViewField* fView = new ViewField();
  fView->SetSensor(sensor);
  fView->SetArea(-lx,-ly,lx,ly);
  fView->SetVoltageRange(-10., 2200.);
  fView->PlotContour();
  gPad->SaveAs("voltage.png");


  /*
  string junk;
  cin >> junk;

  return 0;
  */

  /*
  ViewGeometry* view = new ViewGeometry();
  view->SetGeometry(geo);
  view->Plot();
  */

  // MC integration
  /*
  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);
  drift->EnableSignalCalculation();
  */
 
  // Setup HEED
  TrackHeed* track = new TrackHeed();
  track->SetParticle("muon");
  track->SetMomentum(3.e9); // in eV?
  double x0 = 0;
  double y0 = ly;
  double z0 = 0.;
  double t0 = 0.;   // Initial position
  double dx0 = 0., dy0 = -1., dz0 = 0.;   // Direction of travel
  track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
  track->SetSensor(sensor);

  // Start simulation 
  sensor->ClearSignal();
  sensor->NewSignal();
  track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);

  // Cluster 
  double xc, yc, zc, tc, ec, extra;
  int nc;

  double time;
  double location;
  int count = 0;

  // Loop over all clusters created by muon
  while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    count++;
    // Find radial location of cluster
    location = sqrt(xc * xc + zc * zc);

    // Find the creation time of the ions created by the electron cluster
    //time = IonTiming(location);
    std::cout << "Cluster: " << count << std::endl;
    std::cout << "Size of Cluster: " << nc << std::endl;
    std::cout << "Location of Cluster: " << yc << "\tr\t" << location << std::endl;
    std::cout << xc << "\t" << yc << "\t" << zc << "\t" << tc << "\t" << nc << "\t" << ec << "\t" << extra << std::endl;

      // For each electron in the cluster:
      /*
      for (int i = 0; i < nc; i++) {
        // Get the number of electron tracks.
        //double np = gain * RndmPolya(theta);
        // Scale the effect of the ion induced signal by the number of electron tracks.
        //drift->SetIonSignalScalingFactor(np);
        // Drift 1 highly charged ion to estimate the signal created by the avalanche from 1 electron
        //drift->DriftIon(xIon, 0., 0., time);
      }
      */
  }

/*
  std::cout << "Plotting Signal..." << std::endl;
  ViewSignal* signalView = new ViewSignal();
  signalView->SetSensor(sensor);
  TCanvas* c1 = new TCanvas("c1", "Signal", 21, 28, 500, 527);
  signalView->SetCanvas(c1);
  signalView->PlotSignal("s");
  // Save the plotted driftlines to a pdf
  c1->SaveAs("ApproxSignal.pdf");
  */
  
  app.Run(kTRUE);
}

  
