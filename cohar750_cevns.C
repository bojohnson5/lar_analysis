#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <vector>

using namespace ::std;

double qf(double en);
double f90(double en);
double er_f90(double en);
double nr_f90(double en);
double full_f90(double *x, double *pars);
bool is_sharing_energy(vector<vector<double>> *channels);

void cohar750_cevns() {
  /* signal prediction from dukecevns */
  ifstream f("../data/sns_diff_rates-cohar750conv-Ar-klein-Ar40.out");
  vector<double> xs, ys;
  double x, y;

  while (!f.eof()) {
    f >> x >> y;
    xs.push_back(x * 1000.0);
    ys.push_back(y / 1000.0);
  }
  f.close();

  auto g1 = new TGraph(xs.size(), xs.data(), ys.data());

  ifstream fprompt("../data/sns_diff_rates-cohar750prompt-Ar-klein-Ar40.out");
  vector<double> xsprompt, ysprompt;

  while (!fprompt.eof()) {
    fprompt >> x >> y;
    xsprompt.push_back(x * 1000.0);
    ysprompt.push_back(y / 1000.0);
  }
  fprompt.close();

  auto gprompt = new TGraph(xsprompt.size(), xsprompt.data(), ysprompt.data());

  ifstream fdelayed("../data/sns_diff_rates-cohar750delayed-Ar-klein-Ar40.out");
  vector<double> xsdelayed, ysdelayed;

  while (!fdelayed.eof()) {
    fdelayed >> x >> y;
    xsdelayed.push_back(x * 1000.0);
    ysdelayed.push_back(y / 1000.0);
  }
  fdelayed.close();

  auto gdelayed =
      new TGraph(xsdelayed.size(), xsdelayed.data(), ysdelayed.data());

  auto c1 = new TCanvas("c1", "canvas", 3000, 1000);
  c1->Divide(2, 3);
  c1->cd(1);
  g1->Draw("ALP");

  auto h_dukecevns = new TH1F(
      "h_dukecevns", "COHAr750 dukecevns prediction;True energy [keVr];Events",
      100, 0, 200);
  h_dukecevns->Sumw2();
  auto h_dukecevnsprompt =
      new TH1F("h_dukecevnsprompt",
               "COHAr750 dukecevns prompt prediction;True energy [keVr];Events",
               100, 0, 200);
  h_dukecevnsprompt->Sumw2();
  auto h_dukecevnsdelayed = new TH1F(
      "h_dukecevnsdelayed",
      "COHAr750 dukecevns delayed prediction;True energy [keVr];Events", 100, 0,
      200);
  h_dukecevnsdelayed->Sumw2();
  double bin_width = h_dukecevns->GetBinWidth(1);
  int nbins = h_dukecevns->GetNbinsX();
  for (int i = 1; i < nbins + 1; i++) {
    double center = h_dukecevns->GetBinCenter(i);
    h_dukecevns->SetBinContent(i, g1->Eval(center) * bin_width);
    h_dukecevns->GetBinError(i, sqrt(g1->Eval(center) * bin_width));
    h_dukecevnsprompt->SetBinContent(i, gprompt->Eval(center) * bin_width);
    h_dukecevnsprompt->GetBinError(i, sqrt(gprompt->Eval(center) * bin_width));
    h_dukecevnsdelayed->SetBinContent(i, gdelayed->Eval(center) * bin_width);
    h_dukecevnsdelayed->GetBinError(i,
                                    sqrt(gdelayed->Eval(center) * bin_width));
  }
  cout << "Events before cuts: " << h_dukecevns->Integral() << endl;
  cout << "Events before cuts[prompt]: " << h_dukecevnsprompt->Integral()
       << endl;
  cout << "Events before cuts[delayed]: " << h_dukecevnsdelayed->Integral()
       << endl;
  c1->cd(2);
  h_dukecevns->Draw("hist");

  /* simulation results */
  auto fin = TFile::Open("../data/cohar750_cevns_sim_output.root");
  auto CENNS = fin->Get<TTree>("CENNS");

  double en;
  vector<int> *top = 0;
  vector<vector<double>> *chans = 0;
  CENNS->SetBranchAddress("en", &en);
  CENNS->SetBranchAddress("top", &top);
  CENNS->SetBranchAddress("chans", &chans);

  TH1F *h_sim = new TH1F(
      "h_sim", "Simulated event energy;Energy [keVee];Counts", 100, 0, 200);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    h_sim->Fill(en * 1000);
  }
  c1->cd(3);
  h_sim->Draw("hist");

  TH1F *h3 = (TH1F *)h_dukecevns->Clone();
  TH1F *h4 = (TH1F *)h_dukecevnsprompt->Clone();
  TH1F *h5 = (TH1F *)h_dukecevnsdelayed->Clone();
  h3->Divide(h_sim);
  h4->Divide(h_sim);
  h5->Divide(h_sim);

  /* determine detector response function */
  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  auto h_resp = new TH2F(
      "h_resp", "COH-Ar-750 Ar40 response;PEs (smear); True energy [keVr]", 250,
      0, 500, 100, 0, 200);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    double top_pes = 0.0;
    double bot_pes = 0.0;
    double total_pes = 0.0;
    for (int j = 0; j < chans->size(); j++) {
      // smear the PEs from a channel
      double smeared_pes = 0.0;
      for (auto pe : chans->at(j)) {
        double r = rng.Gaus(1.0, sigma);
        if (r < 0.0)
          r = 0.0;
        smeared_pes += r;
      }
      total_pes += smeared_pes;
      if (top->at(j) == 1) {
        top_pes += smeared_pes;
      } else {
        bot_pes += smeared_pes;
      }
    }
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      h_resp->Fill(total_pes, en * 1000);
    }
  }

  auto profx = (TH1F *)h_resp->ProfileX();
  auto fit_func =
      new TF1("fit_func", "[0] * x + sqrt([1] * x^2 + [2] * x)", 0, 500);
  fit_func->SetParameters(0.25, 0.25, 0.25);
  fit_func->SetParLimits(0, 0.1, 5.0);
  fit_func->SetParLimits(1, -1.0, 1.0);
  fit_func->SetParLimits(2, 0.0, 200.0);
  profx->Fit(fit_func, "R0Q");
  c1->cd(4);
  h_resp->Draw("colz");
  fit_func->Draw("same");

  auto h_cevns =
      new TH1F("h_cevns",
               "COHAr750 CEvNS Prediction (no F90);Reconstructed energy "
               "[keVee];Counts",
               40, 0, 40);
  auto h_cevnsadjf901 = new TH1F(
      "h_cevnsadjf901",
      "COHAr750 CEvNS Prediction (Data-driven F90);Reconstructed energy "
      "[keVee];Counts",
      40, 0, 40);
  auto h_cevnsprompt =
      new TH1F("h_cevnsprompt",
               "COHAr750 Prompt CEvNS Prediction (no F90);Reconstructed energy "
               "[keVee];Counts",
               40, 0, 40);
  auto h_cevnsadjf901prompt = new TH1F(
      "h_cevnsadjf901prompt",
      "COHAr750 Prompt CEvNS Prediction (Data-driven F90);Reconstructed energy "
      "[keVee];Counts",
      40, 0, 40);
  auto h_cevnsdelayed = new TH1F(
      "h_cevnsdelayed",
      "COHAr750 Delayed CEvNS Prediction (no F90);Reconstructed energy "
      "[keVee];Counts",
      40, 0, 40);
  auto h_cevnsadjf901delayed =
      new TH1F("h_cevnsadjf901delayed",
               "COHAr750 Delayed CEvNS Prediction (Data-driven "
               "F90);Reconstructed energy "
               "[keVee];Counts",
               40, 0, 40);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    double top_pes = 0.0;
    double bot_pes = 0.0;
    double total_pes = 0.0;
    for (int j = 0; j < chans->size(); j++) {
      // smear the PEs from a channel
      double smeared_pes = 0.0;
      for (auto pe : chans->at(j)) {
        double r = rng.Gaus(1.0, sigma);
        if (r < 0.0)
          r = 0.0;
        smeared_pes += r;
      }
      total_pes += smeared_pes;
      if (top->at(j) == 1) {
        top_pes += smeared_pes;
      } else {
        bot_pes += smeared_pes;
      }
    }
    // if (top_pes >= 2.0 && bot_pes >= 2.0 && is_sharing_energy(chans)) {
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      double reco_en = fit_func->Eval(total_pes);
      int weight_bin = h3->FindBin(en * 1000);
      double weight = h3->GetBinContent(weight_bin);
      int weight_binprompt = h4->FindBin(en * 1000);
      double weightprompt = h4->GetBinContent(weight_bin);
      int weight_bindelayed = h5->FindBin(en * 1000);
      double weightdelayed = h5->GetBinContent(weight_bin);
      double en_quench = reco_en * qf(reco_en);
      h_cevns->Fill(en_quench, weight);
      h_cevnsprompt->Fill(en_quench, weightprompt);
      h_cevnsdelayed->Fill(en_quench, weightdelayed);
      double reco_f90 = nr_f90(en_quench);
      if (reco_f90 > 0.5 && reco_f90 < 0.9) {
        h_cevnsadjf901->Fill(en_quench, weight);
        h_cevnsadjf901prompt->Fill(en_quench, weightprompt);
        h_cevnsadjf901delayed->Fill(en_quench, weightdelayed);
      }
    }
  }
  c1->cd(5);
  h_cevns->Draw("hist");
  c1->cd(6);
  h_cevnsadjf901->Draw("hist");

  c1->Update();
  int b1 = h_cevns->FindBin(5);
  int b2 = h_cevns->FindBin(39.5);
  cout << "Bins: " << b1 << " and " << b2 << endl;

  cout << "Events after cuts[full]: " << h_cevns->Integral(b1, b2) << endl;
  cout << "Events after cuts and F90[full]: "
       << h_cevnsadjf901->Integral(b1, b2) << endl;
  cout << "Events after cuts[prompt]: " << h_cevnsprompt->Integral(b1, b2)
       << endl;
  cout << "Events after cuts and F90[prompt]: "
       << h_cevnsadjf901prompt->Integral(b1, b2) << endl;
  cout << "Events after cuts[delayed]: " << h_cevnsdelayed->Integral(b1, b2)
       << endl;
  cout << "Events after cuts and F90[delayed]: "
       << h_cevnsadjf901delayed->Integral(b1, b2) << endl;
  auto fout = new TFile("../data/cohar750_cevns.root", "recreate");
  g1->Write();
  h_dukecevns->Write();
  h_dukecevnsprompt->Write();
  h_dukecevnsdelayed->Write();
  h3->Write();
  h4->Write();
  h5->Write();
  h_resp->Write();
  fit_func->Write();
  h_cevns->Write();
  h_cevnsadjf901->Write();
  h_cevnsprompt->Write();
  h_cevnsadjf901prompt->Write();
  h_cevnsdelayed->Write();
  h_cevnsadjf901delayed->Write();
  fout->Close();
  delete fout;
}

double qf(double en) {
  if (en < 100 && en > 0) {
    return 0.251 + 0.000752 * en;
  } else {
    return 0.3262;
  }
}

double f90(double en) {
  double mu = 0.765 - 1.014 / en;
  double sigma = 0.0266 + 0.473 / en;

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

double nr_f90(double en) {
  double sigma = 0.0266 + 0.473 / en;
  double mu = 0;
  if (en < 40) {
    mu = 0.765 - 1.014 / en - 0.0000112 / (en * en);
  } else {
    mu = 0.767 - 0.000264 * en;
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

bool is_sharing_energy(vector<vector<double>> *channels) {
  const int num_chan = channels->size();
  double chans[num_chan];
  int q1 = ceil(0.25 * num_chan);
  int q3 = ceil(0.75 * num_chan);
  for (int i = 0; i < num_chan; i++) {
    chans[i] = channels->at(i).size();
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
