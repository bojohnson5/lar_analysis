#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

double qf(double en);
double nr_f90(double en, TRandom3 &rng);

void cohar750_neutrons(string sim_files, int nbins, int ly_low, int ly_high,
                       string output_name, double sim_area, double meas_flux) {
  auto CENNS = new TChain("CENNS");
  CENNS->Add(sim_files.c_str());

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  auto f = TFile::Open("../data/cohar750_cevns.root");
  auto fit_func = (TF1 *)f->Get("fit_func");

  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  string title = "Neutron Sim;Detected PEs;Counts";

  const double alpha = -1.25;
  const double A =
      (alpha + 1.0) / (pow(40.0, alpha + 1) - pow(2.5, alpha + 1)) * meas_flux;
  const double delta_e = (300.0 - 0.0) / 300.0;

  auto init_en = TH1F("init_en", "Simulated neutron energies", 300, 0, 300);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    init_en.Fill(ev->fPrimaryEnergy);
  }

  auto neutrons_f90 =
      TH1F(output_name.c_str(), title.c_str(), nbins, ly_low, ly_high);
  auto name = output_name + "_nof90";
  auto neutrons_nof90 =
      TH1F(name.c_str(), title.c_str(), nbins, ly_low, ly_high);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    double top_pes = 0.0;
    double bot_pes = 0.0;
    for (int j = 0; j < ev->fNumChans; j++) {
      // smear the PEs from a channel
      for (const auto pe : ev->fvChannels[j]) {
        double r = rng.Gaus(1.0, sigma);
        if (r < 0.0)
          r = 0.0;
        if (pe.pos[2] > 0)
          top_pes += r;
        else
          bot_pes += r;
      }
    }
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      int b = init_en.FindBin(ev->fPrimaryEnergy);
      double w = sim_area * A * pow(ev->fPrimaryEnergy, alpha) /
                 init_en.GetBinContent(b);
      double reco_en = fit_func->Eval(top_pes + bot_pes);
      double en_quenched = reco_en * qf(reco_en);
      neutrons_nof90.Fill(en_quenched, w);
      double reco_f90 = nr_f90(en_quenched, rng);
      if (reco_f90 > 0.5 && reco_f90 < 0.9) {
        neutrons_f90.Fill(en_quenched, w);
      }
    }
  }
  int b1 = neutrons_f90.FindBin(5);
  int b2 = neutrons_f90.FindBin(40);
  cout << "Integral[f90]: " << neutrons_f90.Integral(b1, b2) << endl;
  cout << "Integral[nof90]: " << neutrons_nof90.Integral(b1, b2) << endl;

  stringstream ss;
  ss << "../data/" << output_name << ".root";
  auto fout = new TFile(ss.str().c_str(), "recreate");
  neutrons_f90.Write();
  neutrons_nof90.Write();
  fout->Close();
  delete fout;
  delete CENNS;
}

double qf(double en) {
  if (en < 100 && en > 0) {
    return 0.251 + 0.000752 * en;
  } else {
    return 0.3262;
  }
}

double nr_f90(double en, TRandom3 &rng) {
  double sigma = 0.0266 + 0.473 / en;
  double mu = 0;
  if (en < 40) {
    mu = 0.765 - 1.014 / en - 0.0000112 / (en * en);
  } else {
    mu = 0.767 - 0.000264 * en;
  }

  return rng.Gaus(mu, sigma);
}
