#include "/N/u/bojohn/Quartz/cohar750_sim/cenns/io/CENNSEvent.hh"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TRandom3.h"
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

void ly(string filename, int nbins, int low, int high, string plotname) {
  auto CENNS = new TChain("CENNS");
  CENNS->Add(filename.c_str());

  CENNSEvent *ev = nullptr;
  CENNS->SetBranchAddress("Event", &ev);

  TRandom3 rng;
  rng.SetSeed(1234);
  double sigma = 0.44;
  string title = plotname + " Sim;Detected PEs;Counts";
  auto ly = TH1F("ly", title.c_str(), nbins, low, high);
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
        if (pe.pos[2] > 0) {
          // top_pes += r;
          top_pes += 1.0;
        } else {
          // bot_pes += r;
          bot_pes += 1.0;
        }
      }
    }
    if (top_pes >= 2.0 && bot_pes >= 2.0) {
      ly.Fill(top_pes + bot_pes);
    }
  }
  delete CENNS;

  auto c = new TCanvas();
  ly.Draw("hist");
  // int b1 = ly.FindBin(0);
  // int b1 = ly.FindBin(30);
  // int b2 = ly.FindBin(240);
  // cout << "CEvNS ROI counts: " << ly.Integral(b1, b2) << endl;
  cout << "CEvNS ROI counts: " << ly.Integral(0, 240) << endl;
  // c->SetLogy();
  replace(plotname.begin(), plotname.end(), ' ', '_');
  c->SaveAs(("../plots/" + plotname + ".png").c_str());
  delete c;
}
