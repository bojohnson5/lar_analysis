#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

void ly(string filename, int nbins, int ly_low, int ly_high, string plotname,
        double sim_area, double meas_flux) {
  auto CENNS = new TChain("CENNS");
  CENNS->Add(filename.c_str());

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  auto f = TFile::Open("../data/cohar750_ar39.root");
  auto fit_func = (TF1 *)f->Get("fit_func");

  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  string title = plotname + " Sim;Detected PEs;Counts";

  const double alpha = -1.25;
  const double A =
      (alpha + 1.0) / (pow(40.0, alpha + 1) - pow(2.5, alpha + 1)) * meas_flux;
  const double delta_e = (300.0 - 0.0) / 300.0;

  auto init_en = TH1F("init_en", "Simulated neutron energies", 300, 0, 300);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    init_en.Fill(ev->fPrimaryEnergy);
  }

  auto ly = TH1F("ly", title.c_str(), nbins, ly_low, ly_high);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    double top_pes = 0.0;
    double bot_pes = 0.0;
    for (int j = 0; j < ev->fNumChans; j++) {
      // smear the PEs from a channel
      double smeared_pes = 0.0;
      for (const auto pe : ev->fvChannels[j]) {
        if (pe.ischerenkov)
          continue;
        double r = rng.Gaus(1.0, sigma);
        if (r < 0.0)
          r = 0.0;
        smeared_pes += r;
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
      ly.Fill(top_pes + bot_pes, w);
    }
  }

  auto c = new TCanvas();
  ly.Draw("hist");
  // int b1 = ly.FindBin(0);
  int b1 = ly.FindBin(30);
  int b2 = ly.FindBin(240);
  cout << "CEvNS ROI counts: " << ly.Integral(b1, b2) << endl;
  c->SetLogy();
  replace(plotname.begin(), plotname.end(), ' ', '_');
  c->SaveAs(("../plots/" + plotname + ".png").c_str());
  delete c;

  delete CENNS;
}
