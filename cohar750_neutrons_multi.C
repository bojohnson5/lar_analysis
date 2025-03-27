#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TRandom3.h"
#include "progressbar.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

double qf(double en);
double nr_f90(double en, TRandom3 &rng);

void cohar750_neutrons(vector<string> sim_files, int nbins, int ly_low,
                       int ly_high, string output_name, double sim_area,
                       double meas_flux) {
  auto CENNS = new TChain("CENNS");
  for (auto fn : sim_files) {
    CENNS->Add(fn.c_str());
  }
  // CENNS->Add(sim_files.c_str());

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
  const double SE = (pow(300.0, alpha + 1) - pow(0.01, alpha + 1)) /
                    (pow(40.0, alpha + 1) - pow(2.5, alpha + 1));
  const double delta_e = (300.0 - 0.0) / 300.0;

  auto init_en = TH1F("init_en", "Simulated neutron energies", 30000, 0, 300);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    init_en.Fill(ev->fPrimaryEnergy);
  }
  // auto g = new TGraph();
  // for (int i = 1; i <= init_en.GetNbinsX(); i++) {
  //   double en = init_en.GetBinCenter(i);
  //   std::cout << en << std::endl;
  //   int n = init_en.GetBinContent(i);
  //   if (n == 0)
  //     n = 1;
  //   std::cout << n << std::endl;
  //   g->AddPoint(en, pow(en, alpha) / n);
  // }
  // auto cg = new TCanvas();
  // cg->SetLogy();
  // g->Draw();
  // cg->SaveAs("../plots/graph.png");
  // auto ch = new TCanvas();
  // init_en.Draw("hist");
  // ch->SaveAs("../plots/hist.png");
  // return;

  auto neutrons_f90 =
      TH1F(output_name.c_str(), title.c_str(), nbins, ly_low, ly_high);
  auto name = output_name + "_nof90";
  auto neutrons_nof90 =
      TH1F(name.c_str(), title.c_str(), nbins, ly_low, ly_high);
  // auto name_ly = output_name + "_ly";
  // auto neutrons_ly_f90 =
  //     TH1F(name_ly.c_str(), title.c_str(), nbins, ly_low, ly_high);
  // name_ly = name_ly + "_nof90";
  // auto neutrons_ly_nof90 =
  //     TH1F(name_ly.c_str(), title.c_str(), nbins, ly_low, ly_high);
  ProgressBar<long long> bar(CENNS->GetEntries());
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    bar.update(i);
    CENNS->GetEntry(i);
    double top_pes = ev->fDetPhotonsTop;
    double bot_pes = ev->fDetPhotonsBottom;
    double top_smear = 0.0;
    double bot_smear = 0.0;
    // bool isNR = ev->fEDepLAr == ev->fEDepNR;
    for (int j = 0; j < top_pes; j++) {
      double r = rng.Gaus(1.0, sigma);
      if (r < 0.0)
        r = 0.0;
      top_smear += r;
    }
    for (int j = 0; j < bot_pes; j++) {
      double r = rng.Gaus(1.0, sigma);
      if (r < 0.0)
        r = 0.0;
      bot_smear += r;
    }
    // if (top_smear >= 2.0 && bot_smear >= 2.0 && isNR) {
    if (top_smear >= 2.0 && bot_smear >= 2.0) {
      int b = init_en.FindBin(ev->fPrimaryEnergy);
      // double w = 1.0;
      double w = sim_area * A * pow(ev->fPrimaryEnergy, alpha) /
                 init_en.GetBinContent(b);
      // neutrons_ly_nof90.Fill(top_smear + bot_smear);
      double reco_en = fit_func->Eval(top_smear + bot_smear);
      double en_quenched = reco_en * qf(reco_en);
      neutrons_nof90.Fill(en_quenched, w);
      double reco_f90 = nr_f90(en_quenched, rng);
      if (reco_f90 > 0.5 && reco_f90 < 0.9) {
        neutrons_f90.Fill(en_quenched, w);
        // neutrons_ly_f90.Fill(top_smear + bot_smear);
      }
    }
  }
  bar.finish();
  int b1 = neutrons_f90.FindBin(5);
  int b2 = neutrons_f90.FindBin(40);
  cout << "Integral[f90]: " << neutrons_f90.Integral(b1, b2) << endl;
  cout << "Integral[nof90]: " << neutrons_nof90.Integral(b1, b2) << endl;

  auto c = new TCanvas();
  neutrons_f90.Draw("hist");
  neutrons_nof90.SetLineColor(kRed);
  neutrons_nof90.Draw("hist;same");
  auto plot_name = "../plots/" + output_name + ".png";
  c->SaveAs(plot_name.c_str());
  delete c;

  stringstream ss;
  ss << "../data/" << output_name << ".root";
  auto fout = new TFile(ss.str().c_str(), "recreate");
  neutrons_f90.Write();
  neutrons_nof90.Write();
  // neutrons_ly_f90.Write();
  // neutrons_ly_nof90.Write();
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
