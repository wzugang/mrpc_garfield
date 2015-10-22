#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TApplication.h>
#include <TFile.h>
#include <TH1.h>
#include <TGraph.h>
#include <TString.h>

#include "ViewSignal.hh"
#include "ComponentAnalyticField.hh"
#include "ComponentConstant.hh"
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

const int NUM_GASES = 2;
TH1 *h_primaries[NUM_GASES];
TGraph* g_townsend[NUM_GASES];    // townsend coefficient
TGraph* g_attachment[NUM_GASES];  // attachment coefficient
TGraph* g_effective[NUM_GASES];   // effective gain coefficient
TGraph* g_drift[NUM_GASES];       // drift velocity

int main(int argc, char *argv[]) {

  TApplication app("app",&argc, argv);
  TFile *savefile = new TFile("primaries.root","RECREATE");

  plottingEngine.SetDefaultStyle();
 
  TString name;
  TString title;

  // Make a gas medium
  MediumMagboltz* gas = new MediumMagboltz();
  // Set the temperature [K] and pressure [Torr]
  const double pressure = 1 * AtmosphericPressure;  // 760 Torr
  const double temperature = ZeroCelsius + 20.0;    // 293.15 C
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);


  for (int igas=0; igas<NUM_GASES; igas++)
  {
  
    //== Set Gas Composition
    if ( igas==0 ) {
      gas->SetComposition("C2H2F4", 95., "iC4H10", 4.5, "SF6", 0.5);
      gas->LoadGasFile("r134a_95_iso_4.5_sf6_0.5.gas"); //for electrons
      title = "r134a_95_iso_4.5_sf6_0.5";
    } else if (igas==1 ) {
      gas->SetComposition("C2H2F4", 97., "SF6", 3.0);
      gas->LoadGasFile("r134a_97_sf6_3.0.gas");
      title = "r134a_97_sf6_3.0";
    }

    //gas->LoadIonMobility("/Users/chiu/software/garfield/Data/IonMobility_CO2+_CO2.txt"); // for ions
    /*
    gas->SetComposition("CO2", 98., "iC4H10", 8.0);
    gas->LoadGasFile("co2_92_iso_8.gas"); //for electrons
    */

    // Make histograms
    name = "h_primaries"; name += igas;
    h_primaries[igas] = new TH1F(name,title,31,-0.5,30.5);
    h_primaries[igas]->SetLineColor(igas+2);
    h_primaries[igas]->SetXTitle("num primaries");

    // Retrieve the points at which the Townsend coefficient was computed
    std::vector<double> efields,  bfields, angles;
    gas->GetFieldGrid(efields, bfields, angles);
    cout << "num grid points\t" << efields.size() << "\t" << bfields.size() << "\t" << angles.size() << endl;

    double e[100];
    double alpha[100];
    double attach[100];
    double effective[100];
    double drift[100];
    // Retrieve the Townsend coefficient and prepare plot vectors
    for (unsigned long i=0; i < efields.size(); i++) {
      e[i] = efields[i];

      double logalpha;
      gas->GetElectronTownsend(i, 0, 0, logalpha);
      alpha[i] = exp(logalpha);

      double logattach;
      gas->GetElectronAttachment(i, 0, 0, logattach);
      attach[i] = exp(logattach);

      effective[i] = alpha[i] - attach[i];

      double tempdrift;
      gas->GetElectronVelocityE(i, 0, 0, tempdrift);
      drift[i] = tempdrift;

      cout << "townsend/attach " << i << "\t" << e[i] << "\t"
        << alpha[i] << "\t" << attach[i] << "\t" << drift[i] << endl;

    }

    // Fill the graphs
    name = "g_townsend"; name += igas;
    g_townsend[igas] = new TGraph(efields.size(), e, alpha);
    g_townsend[igas]->SetName(name);
    g_townsend[igas]->SetLineColor(igas+1);
    g_townsend[igas]->SetMarkerColor(igas+1);
    g_townsend[igas]->SetLineStyle(2);

    name = "g_attachment"; name += igas;
    g_attachment[igas] = new TGraph(efields.size(), e, attach);
    g_attachment[igas]->SetName(name);
    g_attachment[igas]->SetLineColor(igas+1);
    g_attachment[igas]->SetMarkerColor(igas+1);
    g_attachment[igas]->SetLineStyle(9);

    name = "g_effective"; name += igas;
    g_effective[igas] = new TGraph(efields.size(), e, effective);
    g_effective[igas]->SetName(name);
    g_effective[igas]->SetLineColor(igas+1);
    g_effective[igas]->SetMarkerColor(igas+1);

    name = "g_drift"; name += igas;
    g_drift[igas] = new TGraph(efields.size(), e, drift);
    g_drift[igas]->SetName(name);
    g_drift[igas]->SetLineColor(igas+1);
    g_drift[igas]->SetMarkerColor(igas+1);

    // Build the geometry.
    // Build the geometry.
    GeometrySimple* geo = new GeometrySimple();

    // Dimensions of 1 gas gap [cm]
    const Double_t lx = 10/2.;   // half lengths
    const Double_t ly = 0.0105/2.;   
    const Double_t lz = 10/2.;  

    SolidBox *tube = new SolidBox(0.,0.,0.,lx,ly,lz);

    // Add  the solid to the geometry, together with the medium inside
    geo->AddSolid(tube, gas);

    ComponentConstant* comp = new ComponentConstant();
    comp->SetGeometry(geo);
    comp->SetElectricField(0., 152000., 0.);

    // Make a sensor.
    Sensor* sensor = new Sensor();
    sensor->AddComponent(comp);
    //sensor->AddElectrode(comp, "s");
    sensor->SetArea(-lx, -ly, -lz, lx, ly, lz);
    //sensor->SetTimeWindow(tmin, tstep, nTimeBins);
    //sensor->ClearSignal();
    std::cout << "Num Electrodes = " << sensor->GetNumberOfElectrodes() << std::endl;

    // Setup HEED
    TrackHeed* track = new TrackHeed();
    track->SetSensor(sensor);
    track->SetParticle("pi");
    track->SetMomentum(3.e9); // in eV
    // Switch on debugging to print out some information (stopping power, W value, ...)
    track->EnableDebugging();

    for (int ievt=0; ievt<10000; ievt++)
    {
      double x0 = 0, y0 = ly, z0 = 0., t0 = 0.;   // Initial position
      double dx0 = 0., dy0 = -1., dz0 = 0.;       // Direction of travel
      track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);

      // Cluster coordinates
      double xc = 0., yc = 0., zc = 0., tc = 0.;
      // Number of electrons produced in a collision
      int nc = 0;
      // Energy loss in a collision
      double ec = 0.;
      // Dummy variable (not used at present)
      double extra = 0.;
      // Total energy loss along the track
      double esum = 0.;
      // Total number of electrons produced along the track
      int nsum = 0;
      // Loop over the clusters. 
      // GetCluster returns false if the list of clusters is exhausted 
      // (or if the calculation of the track failed). 
      int nloop = 0;
      while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
        esum += ec;
        nsum += nc;
        cout << nloop << "\t" << xc << "\t" << yc << "\t" << zc << "\t" << tc
          << "\t" << nc << "\t" << ec << endl;
        for (int i = 0; i < nc; ++i) {
          // Coordinates of the conduction electron
          double xe, ye, ze, te;
          // Energy and direction of the electron (these are not provided by Heed, 
          // but the parameters are included for compatibility with other Track classes). 
          double ee, dxe, dye, dze;
          track->GetElectron(i, xe, ye, ze, te, ee, dxe, dye, dze);
          // Calculate electron and ion drift lines starting from xe, ye, ze
          cout << ":" << i << "\t" << xe << "\t" << ye << "\t" << ze << "\t" << te
            << "\t" << ee << "\t" << dxe << "\t" << dye << "\t" << dze << endl;
        }
      }

      h_primaries[igas]->Fill(nsum);
      //cout << "Found all tracks " << nsum << endl;
    }

    delete tube;
    delete comp;
    delete sensor;
    delete track;
  }

  TCanvas *c_properties = new TCanvas("c_properties","townsend/attachment coeff",550,425);
  g_townsend[0]->Draw("alp");
  g_townsend[0]->GetHistogram()->SetYTitle("Townsend/Attachment Coeffs (1/cm)");
  g_townsend[0]->GetHistogram()->SetXTitle("E-field (V/cm)");
  g_attachment[0]->Draw("lp");
  g_effective[0]->Draw("lp");
  for (int igas=1; igas<NUM_GASES; igas++)
  {
    g_townsend[igas]->Draw("lp");
    g_attachment[igas]->Draw("lp");
    g_effective[igas]->Draw("lp");
  }
  // Get Maximum of all plots
  /*
  Double_t gain_max = -1e9;
  Double_t gain_min = 1e9;
  cout << "maxmin " << gain_max << "\t" << gain_min << endl;
  for (int igas=0; igas<NUM_GASES; igas++)
  {
    if ( g_townsend[igas]->GetMaximum() > gain_max ) gain_max = g_townsend[igas]->GetMaximum();
    if ( g_attachment[igas]->GetMaximum() > gain_max ) gain_max = g_attachment[igas]->GetMaximum();
    if ( g_effective[igas]->GetMaximum() > gain_max ) gain_max = g_effective[igas]->GetMaximum();
    if ( g_townsend[igas]->GetMinimum() < gain_min ) gain_min = g_townsend[igas]->GetMinimum();
    if ( g_attachment[igas]->GetMinimum() < gain_min ) gain_min = g_attachment[igas]->GetMinimum();
    if ( g_effective[igas]->GetMinimum() < gain_min ) gain_min = g_effective[igas]->GetMinimum();
  cout << "maxmin " << gain_max << "\t" << gain_min << endl;
  }
  g_townsend[0]->SetMaximum( gain_max + fabs(gain_max*0.1) );
  g_townsend[0]->SetMinimum( gain_min - fabs(gain_min*0.1) );
  cout << "maxmin " << gain_max << "\t" << gain_min << endl;
  */
  g_townsend[0]->SetMaximum( 10000. );
  g_townsend[0]->SetMinimum( -10000. );
  gPad->Modified();
  gPad->Update();
  c_properties->SaveAs("c_properties.png");

  TCanvas *c_drift = new TCanvas("c_drift","drift",550,425);
  g_drift[0]->Draw("alp");
  g_drift[0]->GetHistogram()->SetXTitle("E-field (V/cm)");
  for (int igas=1; igas<NUM_GASES; igas++)
  {
    g_drift[igas]->Draw("lp");
  }
  c_drift->SaveAs("c_drift.png");

  TCanvas *c_primaries = new TCanvas("c_primaries","primary electrons",550,425);
  for (int igas=0; igas<NUM_GASES; igas++)
  {
    Double_t nevents = h_primaries[igas]->GetEntries();
    h_primaries[igas]->Sumw2();
    h_primaries[igas]->Scale(1.0/nevents);
  }
  h_primaries[0]->GetXaxis()->SetRangeUser(0,10);
  h_primaries[0]->Draw("ehist");
  for (int igas=1; igas<NUM_GASES; igas++) {
    h_primaries[igas]->Draw("ehistsame");
  }
  c_primaries->SaveAs("c_primaries.png");

  savefile->Write();
  //savefile->Close();

  app.Run(kTRUE);

  // Wait for input to prevent crashing at end
  /*
  cout << "End of gas_properties" << endl;
  string junk;
  cin >> junk;
  */

  return 0;
}


