#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TRandom3.h"

using namespace std;

void pb_ly(int nbins, int low, int high) {
  auto CENNS = new TChain("CENNS");
  CENNS->Add("~/ceem_coherent/CENNS750/sim_out/5337801/test_*");
  CENNS->Add("~/ceem_coherent/CENNS750/sim_out/5337806/test_*");

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  auto ly = new TH1F("ly", "Light yield", nbins, low, high);
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
      ly->Fill(top_pes + bot_pes);
    }
  }

  auto c = new TCanvas();
  ly->Draw("hist");
  c->SetLogy();
  c->Draw();
  c->SaveAs("../plots/pb_ly.png");
}
