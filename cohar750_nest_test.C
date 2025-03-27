#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include <vector>

using namespace ::std;

void cohar750_nest_test() {
  /*auto fin =*/
  /*    TFile::Open("~/ceem_coherent/CENNS750/sim_out/5215865/output.root");*/
  /*auto fin =*/
  /*    TFile::Open("~/ceem_coherent/CENNS750/sim_out/5012974/output.root");*/
  auto fin =
      TFile::Open("~/ceem_coherent/CENNS750/sim_out/5219907/output.root");
  auto CENNS = fin->Get<TTree>("CENNS");

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  auto h_sim = new TH1F("h_sim", "Energy;True energy [keVee];Counts", 100, 0, 200);
  auto h_resp = new TH2F(
      "h_resp", "COHAr750 e- response (NEST);PEs (smear);True energy [keVee]",
      333, 0, 1000, 100, 0, 200);
  auto f90ven =
      new TH2F("f90ven", "COHAr750 e- response (NEST);True energy [keVee];F90",
               100, 0, 200, 50, 0, 1);
  for (int i = 0; i < CENNS->GetEntries(); i++) {
    CENNS->GetEntry(i);
    double top_pes = 0.0, bot_pes = 0.0, total_pes = 0.0, f90_pes = 0.0;
    h_sim->Fill(ev->fPrimaryEnergy * 1000);
    for (int j = 0; j < ev->fNumChans; j++) {
      double smeared_pes = 0.0;
      for (const auto pe : ev->fvChannels[j]) {
        if (pe.ischerenkov)
          continue;
        double r = rng.Gaus(1.0, sigma);
        if (r < 0.0)
          r = 0.0;
        smeared_pes += r;
        if (pe.time < 90.0)
          f90_pes += r;
      }
      total_pes += smeared_pes;
      if (smeared_pes > 0.0 && ev->fvChannels[j][0].pos[2] > 0) {
        top_pes += smeared_pes;
      } else if (smeared_pes > 0.0) {
        bot_pes += smeared_pes;
      }
    }
    double f90 = 0.0;
    if (total_pes != 0.0)
      f90 = f90_pes / total_pes;
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      h_resp->Fill(total_pes, ev->fPrimaryEnergy * 1000);
      f90ven->Fill(ev->fPrimaryEnergy * 1000, f90);
    }
  }
  auto c1 = new TCanvas();
  h_resp->Draw("colz");
  c1->Draw();
  auto c2 = new TCanvas();
  f90ven->Draw("colz");
  c2->Draw();

  /*auto fout = new TFile("../data/cohar750_nest_test.root", "recreate");*/
  /*h_resp->Write();*/
  /*f90ven->Write();*/
  /*h_sim->Write();*/
  /*fout->Close();*/
}
