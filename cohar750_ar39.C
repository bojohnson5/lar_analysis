#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace ::std;

double f90(double en);
double er_f90(double en);
double full_f90(double *x, double *pars);
bool is_sharing_energy(CENNSEvent *ev);

void cohar750_ar39() {
  ifstream f("../data/ar_39spec.csv");
  vector<double> xs, ys;
  double x, y;

  while (!f.eof()) {
    f >> x >> y;
    xs.push_back(x);
    ys.push_back(y);
  }
  f.close();

  auto c1 = new TCanvas("c1", "plot", 1000, 2000);
  c1->Divide(3, 3);

  auto g1 = new TGraph(xs.size(), xs.data(), ys.data());
  double norm = g1->Integral();
  cout << "ar39 spec norm: " << norm << endl;
  g1->Scale(476.0 / norm);
  c1->cd(1);
  g1->Draw("ALP");

  auto h_ar39spec =
      new TH1F("h_ar39spec",
               "COH-Ar-750 Ar39 #beta spectrum; True energy [keV];Counts [/s]",
               100, 0, 200);
  h_ar39spec->Sumw2();
  double bin_width = h_ar39spec->GetBinWidth(1);
  double nbins = h_ar39spec->GetNbinsX();

  for (int i = 1; i < nbins + 1; i++) {
    double center = h_ar39spec->GetBinCenter(i);
    h_ar39spec->SetBinContent(i, g1->Eval(center) * bin_width);
    h_ar39spec->SetBinError(i, TMath::Sqrt(g1->Eval(center) * bin_width));
  }
  c1->cd(2);
  h_ar39spec->Draw("hist");

  auto CENNS = new TChain("CENNS");
  /*CENNS->Add("~/ceem_coherent/CENNS750/sim_out/5227976/test_*.root");*/
  CENNS->Add("~/ceem_coherent/CENNS750/sim_out/5359845/test_*.root");
  Long_t num_ent = CENNS->GetEntries();

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  // auto in_f = TFile::Open("../data/cohar750_ar39.root", "update");
  // TH1F *h_sim = (TH1F *)in_f->Get("h_sim");
  TH1F *h_sim = new TH1F(
      "h_sim", "Simulated event energy;Energy [keVee];Counts", 100, 0, 200);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    if (i % 15000 == 0)
      cout << "[h_sim]on entry " << i << " of " << num_ent << endl;
    CENNS->GetEntry(i);
    h_sim->Fill(ev->fPrimaryEnergy * 1000);
  }
  // h_sim->Write();
  c1->cd(3);
  h_sim->Draw("hist");

  TH1F *h3 = (TH1F *)h_ar39spec->Clone();
  h3->Divide(h_sim);

  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  // TH2F *h_resp = (TH2F *)in_f->Get("h_resp_ar39");
  auto h_resp = new TH2F(
      "h_resp_ar39", "COH-Ar-750 e- response;PEs (smear); True energy[keVee]",
      333, 0, 1000, 100, 0, 200);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    if (i % 15000 == 0)
      cout << "[h_resp]on entry " << i << " of " << num_ent << endl;
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
    // if (top_pes >= 2.0 && bot_pes >= 2.0 && is_sharing_energy(ev)) {
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      h_resp->Fill(top_pes + bot_pes, ev->fPrimaryEnergy * 1000);
    }
  }

  // TF1 *fit_func = (TF1 *)in_f->Get("fit_func");
  auto profx = (TH1F *)h_resp->ProfileX();
  auto fit_func = new TF1("fit_func", "[0] * x", 0, 1000);
  profx->Fit(fit_func, "R0Q");
  c1->cd(4);
  h_resp->Draw("colz");
  fit_func->Draw("same");

  auto h_ar39_no_f90 =
      new TH1F("h_ar39_no_f90",
               "Ar39 Signal Prediction (no F90);Reconstructed energy "
               "[keVee];Counts [/s]",
               40, 0, 40);

  auto f_bkg = TFile::Open("../data/nowater_f90.root");
  auto bkg = (TH2F *)f_bkg->Get("bkg");
  auto h_ar39_adj1_f90 =
      new TH1F("h_ar39_adj1_f90",
               "Ar39 Signal Prediction (Data-driven F90);Reconstructed energy "
               "[keVee];Counts [/s]",
               40, 0, 40);

  auto f_bkg2 = TFile::Open("../data/PDF_noWaterData.root");
  auto PDFTree = (TTree *)f_bkg2->Get("PDFTree");
  PDFTree->Draw("F90:(4.5*energy)>>bkg2(500, 0, 1000, 50, 0, 1)", "!isBeam",
                "goff");
  auto bkg2 = (TH2F *)gDirectory->Get("bkg2");
  auto h_ar39_adj2_f90 =
      new TH1F("h_ar39_adj2_f90",
               "Ar39 Signal Prediction (Data-driven F90);Reconstructed energy "
               "[keVee];Counts [/s]",
               40, 0, 40);

  auto h_ar39_adj3_f90 =
      new TH1F("h_ar39_adj3_f90",
               "Ar39 Signal Prediction (Data-driven F90);Reconstructed energy "
               "[keVee];Counts [/s]",
               40, 0, 40);
  for (int i = 0; i < num_ent; i++) {
    if (i % 15000 == 0)
      cout << "on entry " << i << " of " << num_ent << endl;
    CENNS->GetEntry(i);
    double top_pes = 0.0;
    double bot_pes = 0.0;
    for (int j = 0; j < ev->fNumChans; j++) {
      // smear the PEs from a channel
      for (auto pe : ev->fvChannels[j]) {
        double r = rng.Gaus(1.0, sigma);
        if (r < 0.0)
          r = 0.0;
        if (pe.pos[2] > 0)
          top_pes += r;
        else
          bot_pes += r;
      }
    }
    // if (top_pes >= 2.0 && bot_pes >= 2.0 && is_sharing_energy(ev)) {
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      double reco_en = fit_func->Eval(top_pes + bot_pes);
      int weight_bin = h3->FindBin(ev->fPrimaryEnergy * 1000);
      double weight = h3->GetBinContent(weight_bin);
      h_ar39_no_f90->Fill(reco_en, weight);
      int b1 = bkg->ProjectionX()->FindBin(reco_en);
      double f901 = bkg->ProjectionY("", b1, b1 + 1)->GetRandom();
      if (f901 > 0.5 && f901 < 0.9)
        h_ar39_adj1_f90->Fill(reco_en, weight);
      int b2 = bkg2->ProjectionX()->FindBin(top_pes + bot_pes);
      double f902 = bkg2->ProjectionY("", b2, b2 + 1)->GetRandom();
      if (f902 > 0.5 && f902 < 0.9)
        h_ar39_adj2_f90->Fill(reco_en, weight);
      double reco_f90 = er_f90(reco_en);
      if (reco_f90 > 0.5 && reco_f90 < 0.9)
        h_ar39_adj3_f90->Fill(reco_en, weight);
    }
  }
  c1->cd(5);
  h_ar39_no_f90->Draw("hist");
  c1->cd(6);
  h_ar39_adj1_f90->Draw("hist");
  c1->cd(7);
  h_ar39_adj2_f90->Draw("hist");
  c1->cd(8);
  h_ar39_adj3_f90->Draw("hist");

  c1->Update();
  c1->SaveAs("../plots/cohar750_ar39.png");

  int b1 = h_ar39_no_f90->FindBin(5);
  int b2 = h_ar39_no_f90->FindBin(39);
  // int b2 = h_ar39_no_f90->FindBin(40);
  cout << "bin1: " << b1 << " bin2: " << b2 << endl;

  cout << "Integral nof90: " << h_ar39_no_f90->Integral(b1, b2) << endl;
  cout << "Integral 1: " << h_ar39_adj1_f90->Integral(b1, b2) << endl;
  cout << "Integral 2: " << h_ar39_adj2_f90->Integral(b1, b2) << endl;
  cout << "Integral 3: " << h_ar39_adj3_f90->Integral(b1, b2) << endl;

  // in_f->Close();

  auto fout = new TFile("../data/cohar750_ar39.root", "recreate");
  g1->Write();
  h_ar39spec->Write();
  h3->Write();
  h_sim->Write();
  h_resp->Write();
  fit_func->Write();
  h_ar39_no_f90->Write();
  h_ar39_adj1_f90->Write();
  h_ar39_adj2_f90->Write();
  h_ar39_adj3_f90->Write();
  fout->Close();
  delete fout;
  delete CENNS;
}

double f90(double en) {
  double mu = 0.268 + 1.09 / en - 0.633 / (en * en);
  double sigma = 0.0199 + 0.631 / en;

  TRandom3 rng;
  return rng.Gaus(mu, sigma);
}

double er_f90(double en) {
  double sigma = 0.0199 + 0.631 / en;
  double mu = 0;
  if (en < 40) {
    mu = 0.268 + 1.09 / en - 0.633 / (en * en);
  } else {
    mu = 0.292 - 0.000186 * en;
  }

  TRandom3 rng;
  return rng.Gaus(mu, sigma);
}

double full_f90(double *x, double *pars) {
  double en = x[0];
  double xx = x[1];
  double Aer = pars[0] + pars[1] * en + pars[2] * (en * en);
  double Anr = pars[3] + pars[4] * exp(-pars[5] * en);
  double muer = pars[6] + pars[7] / en + pars[8] / (en * en);
  double munr = pars[9] + pars[10] / en + pars[11] / (en * en);
  double sigmaer = pars[12] + pars[13] / en;
  double sigmanr = pars[14] + pars[15] / en;

  double gaus1 = Aer / (Aer + Anr) * 1.0 / (sigmaer * sqrt(2.0 * TMath::Pi())) *
                 exp(-0.5 * pow((xx - muer) / sigmaer, 2.0));
  double gaus2 = Anr / (Aer + Anr) * 1.0 / (sigmanr * sqrt(2.0 * TMath::Pi())) *
                 exp(-0.5 * pow((xx - munr) / sigmanr, 2.0));

  return gaus1 + gaus2;
}

bool is_sharing_energy(CENNSEvent *ev) {
  const int num_chan = ev->fNumChans;
  double chans[num_chan];
  int q1 = ceil(0.25 * num_chan);
  int q3 = ceil(0.75 * num_chan);
  for (int i = 0; i < num_chan; i++) {
    chans[i] = ev->fvChannels[i].size();
  }

  double iqr = chans[q3] - chans[q1];
  double outlier = chans[q3] + 1.5 * iqr;
  long count = 0;
  for (int i = 0; i < num_chan; i++) {
    if (chans[i] > outlier) {
      count++;
    }
  }

  if (count > 0.1 * num_chan)
    return false;
  else
    return true;
}
