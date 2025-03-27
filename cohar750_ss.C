#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

using namespace ::std;

bool is_sharing_energy(CENNSEvent *ev);

void cohar750_ss(string sim_files, string output_name) {
  auto CENNS = new TChain("CENNS");
  /*CENNS->Add("~/ceem_coherent/CENNS750/sim_out/5227976/test_*.root");*/
  CENNS->Add(sim_files.c_str());
  Long_t num_ent = CENNS->GetEntries();

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  auto f = TFile::Open("../data/cohar750_ar39.root");
  auto fit_func = (TF1 *)f->Get("fit_func");

  auto f_bkg = TFile::Open("../data/PDF_noWaterData.root");
  auto PDFTree = (TTree *)f_bkg->Get("PDFTree");
  PDFTree->Draw("F90:(4.5*energy)>>bkg2(500, 0, 1000, 50, 0, 1)", "!isBeam",
                "goff");
  auto bkg = (TH2F *)gDirectory->Get("bkg2");
  auto h_f90 =
      new TH1F(output_name.c_str(),
               "Signal Prediction (Data-driven F90);Reconstructed energy "
               "[keVee];Counts [/s]",
               40, 0, 40);
  auto name = output_name + "nof90";
  auto h_nof90 = new TH1F(name.c_str(),
                          "Signal Prediction (no F90);Reconstructed energy "
                          "[keVee];Counts [/s]",
                          40, 0, 40);
  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  cout << "Total events: " << num_ent << endl;
  for (int i = 0; i < num_ent; i++) {
    CENNS->GetEntry(i);
    double top_pes = ev->fDetPhotonsTop;
    double bot_pes = ev->fDetPhotonsBottom;
    double top_smear = 0.0;
    double bot_smear = 0.0;
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
    if (i % 1000000 == 0) {
      cout << "\ttop_pes: " << top_pes << endl;
      cout << "\tbot_pes: " << bot_pes << endl;
      cout << "\ttop_smear: " << top_smear << endl;
      cout << "\tbot_smear: " << bot_smear << endl;
    }
    if (top_smear >= 2.0 && bot_smear >= 2.0) {
      double reco_en = fit_func->Eval(top_smear + bot_smear);
      h_nof90->Fill(reco_en);
      int b1 = bkg->ProjectionX()->FindBin(top_smear + bot_smear);
      double f901 = bkg->ProjectionY("", b1, b1 + 1)->GetRandom(&rng);
      if (f901 > 0.5 && f901 < 0.9)
        h_f90->Fill(reco_en);
    }
  }
  int b1 = h_f90->FindBin(5);
  int b2 = h_f90->FindBin(39.5);
  cout << "Integral[f90]: " << h_f90->Integral(b1, b2) << endl;
  cout << "Integral[nof90]: " << h_nof90->Integral(b1, b2) << endl;

  stringstream ss;
  ss << "../data/" << output_name << ".root";
  auto fout = new TFile(ss.str().c_str(), "recreate");
  h_f90->Write();
  h_nof90->Write();
  fout->Close();
  delete fout;
  delete h_f90;
  delete h_nof90;
  delete CENNS;
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
