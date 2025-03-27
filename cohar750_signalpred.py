#!/usr/bin/env python
# coding: utf-8

# In[14]:


import ROOT


# In[15]:


# scaling factors
sns_live = 6480  # sec/SNS-year
hog_scale = 4.176226
pb_scale = 9718.48336907423
concrete_scale = 151.52
brn_scale = 2 * 5000  # MW * hr/SNS-year


# In[16]:


f_cevns = ROOT.TFile.Open("../data/cohar750_cevns.root")
h_cevns = f_cevns.Get("h_cevns")
h_cevnsf90 = f_cevns.Get("h_cevnsadjf901")


# In[17]:


c = ROOT.TCanvas()
h_cevns.Draw("hist")
h_cevnsf90.SetLineColorAlpha(ROOT.kRed)
h_cevnsf90.Draw("hist;same")
c.Draw()


# In[18]:


cevns_rate = h_cevns.Integral(6, 39)
cevnsf90_rate = h_cevnsf90.Integral(6, 39)
print(f"Yearly CEvNS: {cevns_rate:.2f}")
print(f"Yearly CEvNS w/F90: {cevnsf90_rate:.2f}")


# In[19]:


f_ar39 = ROOT.TFile.Open("../data/cohar750_ar39.root")
h_ar39 = f_ar39.Get("h_ar39_no_f90")
h_ar39f90 = f_ar39.Get("h_ar39_adj2_f90")
c = ROOT.TCanvas()
h_ar39.Scale(sns_live)
h_ar39f90.Scale(sns_live)
h_ar39.GetYaxis().SetRangeUser(0, 3.5 * sns_live)
h_ar39.Draw("hist")
h_ar39f90.SetLineColorAlpha(ROOT.kRed)
h_ar39f90.Draw("hist;same")
c.Draw()


# In[20]:


ar39_rate = h_ar39.Integral(6, 39)
ar39f90_rate = h_ar39f90.Integral(6, 39)
print(f"Yearly Ar39: {ar39_rate:.2f}")
print(f"Yearly Ar39 w/F90: {ar39f90_rate:.2f}")
print(f"Ratio: {ar39_rate / ar39f90_rate}")


# In[21]:


f_hog_phase1 = ROOT.TFile.Open("../data/hog_phase1.root")
hog_phase1 = f_hog_phase1.Get("hog_phase1")
hog_phase1nof90 = f_hog_phase1.Get("hog_phase1nof90")
c = ROOT.TCanvas()
hog_phase1.Scale(1/hog_scale*sns_live)
hog_phase1nof90.Scale(1/hog_scale*sns_live)
hog_phase1nof90.GetYaxis().SetRangeUser(0, 580/hog_scale*sns_live)
hog_phase1nof90.Draw("hist")
hog_phase1.SetLineColorAlpha(ROOT.kRed)
hog_phase1.Draw("hist;same")
c.Draw()


# In[22]:


hog_phase1nof90_rate = hog_phase1nof90.Integral(6, 39)
hog_phase1_rate = hog_phase1.Integral(6, 39)
print(f"Yearly HOG Phase 1: {hog_phase1nof90_rate:.2f}")
print(f"Yearly HOG Phase 1 w/F90: {hog_phase1_rate:.2f}")
print(f"Ratio: {hog_phase1nof90_rate / hog_phase1_rate}")


# In[29]:


f_hog_phase2 = ROOT.TFile.Open("../data/hog_phase2_morestats.root")
hog_phase2 = f_hog_phase2.Get("hog_phase2_morestats")
hog_phase2nof90 = f_hog_phase2.Get("hog_phase2_morestatsnof90")
c = ROOT.TCanvas()
hog_scale = 1858.95576481144
hog_phase2nof90.Scale(1/hog_scale*sns_live)
hog_phase2.Scale(1/hog_scale*sns_live)
hog_phase2nof90.GetYaxis().SetRangeUser(0, 10000)
hog_phase2nof90.Draw("hist")
hog_phase2.SetLineColorAlpha(ROOT.kRed)
hog_phase2.Draw("hist;same")
c.Draw()


# In[30]:


hog_phase2nof90_rate = hog_phase2nof90.Integral(6, 39)
hog_phase2_rate = hog_phase2.Integral(6, 39)
print(f"Yearly HOG Phase 2: {hog_phase2nof90_rate:.2f}")
print(f"Yearly HOG Phase 2 w/F90: {hog_phase2_rate:.2f}")
print(f"Ratio: {hog_phase2nof90_rate / hog_phase2_rate}")


# In[31]:


f_pb210 = ROOT.TFile.Open("../data/pb210_morestats.root")
pb210 = f_pb210.Get("h_f90")
pb210nof90 = f_pb210.Get("h")
c = ROOT.TCanvas()
pb210nof90.Scale(1/pb_scale*sns_live)
pb210.Scale(1/pb_scale*sns_live)
pb210nof90.GetYaxis().SetRangeUser(0, 600)
pb210nof90.Draw("hist")
pb210.SetLineColorAlpha(ROOT.kRed)
pb210.Draw("hist;same")
c.Draw()


# In[32]:


pb210nof90_rate = pb210nof90.Integral(6, 39)
pb210_rate = pb210.Integral(6, 39)
print(f"Yearly Pb210: {pb210nof90_rate:.2f}")
print(f"Yearly Pb210 w/F90: {pb210_rate:.2f}")
print(f"Ratio: {pb210nof90_rate / pb210_rate}")


# In[33]:


f_k40wall = ROOT.TFile.Open("../data/k40_wall_morestats.root")
k40wall = f_k40wall.Get("k40_wall_morestats")
k40wallnof90 = f_k40wall.Get("k40_wall_morestatsnof90")
f_k40floor = ROOT.TFile.Open("../data/k40_floor_morestats.root")
k40floor = f_k40floor.Get("k40_floor_morestats")
k40floornof90 = f_k40floor.Get("k40_floor_morestatsnof90")
c = ROOT.TCanvas()
concrete_scale = 75757.576
k40wallnof90.Scale(1/concrete_scale*sns_live)
k40wall.Scale(1/concrete_scale*sns_live)
concrete_scale = 90909.091
k40floornof90.Scale(1/concrete_scale*sns_live)
k40floor.Scale(1/concrete_scale*sns_live)
k40wallnof90 += k40floornof90
k40wall += k40floor
k40wallnof90.GetYaxis().SetRangeUser(0, 1000)
k40wallnof90.Draw("hist")
k40wall.SetLineColorAlpha(ROOT.kRed)
k40wall.Draw("hist;same")
c.Draw()


# In[34]:


k40wallnof90_rate = k40wallnof90.Integral(6, 39)
k40wall_rate = k40wall.Integral(6, 39)
print(f"Yearly K40: {k40wallnof90_rate:.2f}")
print(f"Yearly K40 w/F90: {k40wall_rate:.2f}")
print(f"Ratio: {k40wallnof90_rate / k40wall_rate}")


# In[35]:


f_u238wall = ROOT.TFile.Open("../data/u238_wall_morestats.root")
u238wall = f_u238wall.Get("u238_wall_morestats")
u238wallnof90 = f_u238wall.Get("u238_wall_morestatsnof90")
f_u238floor = ROOT.TFile.Open("../data/u238_floor_morestats.root")
u238floor = f_u238floor.Get("u238_floor_morestats")
u238floornof90 = f_u238floor.Get("u238_floor_morestatsnof90")
c = ROOT.TCanvas()
concrete_scale = 24242.424
u238wallnof90.Scale(1/concrete_scale*sns_live)
u238wall.Scale(1/concrete_scale*sns_live)
concrete_scale = 9090.909
u238floornof90.Scale(1/concrete_scale*sns_live)
u238floor.Scale(1/concrete_scale*sns_live)
u238wallnof90 += u238floornof90
u238wall += u238floor
u238wallnof90.GetYaxis().SetRangeUser(0, 6500)
u238wallnof90.Draw("hist")
u238wall.SetLineColorAlpha(ROOT.kRed)
u238wall.Draw("hist;same")
c.Draw()


# In[36]:


u238wallnof90_rate = u238wallnof90.Integral(6, 39)
u238wall_rate = u238wall.Integral(6, 39)
print(f"Yearly U238: {u238wallnof90_rate:.2f}")
print(f"Yearly U238 w/F90: {u238wall_rate:.2f}")
print(f"Ratio: {u238wallnof90_rate / u238wall_rate}")


# In[37]:


f_th232wall = ROOT.TFile.Open("../data/th232_wall_morestats.root")
th232wall = f_th232wall.Get("th232_wall_morestats")
th232wallnof90 = f_th232wall.Get("th232_wall_morestatsnof90")
f_th232floor = ROOT.TFile.Open("../data/th232_floor_morestats.root")
th232floor = f_th232floor.Get("th232_floor_morestats")
th232floornof90 = f_th232floor.Get("th232_floor_morestatsnof90")
c = ROOT.TCanvas()
concrete_scale = 18181.818
th232wallnof90.Scale(1/concrete_scale*sns_live)
th232wall.Scale(1/concrete_scale*sns_live)
concrete_scale = 6060.606
th232floornof90.Scale(1/concrete_scale*sns_live)
th232floor.Scale(1/concrete_scale*sns_live)
th232wallnof90 += th232floornof90
th232wall += th232floor
th232wallnof90.GetYaxis().SetRangeUser(0, 8000)
th232wallnof90.Draw("hist")
th232wall.SetLineColorAlpha(ROOT.kRed)
th232wall.Draw("hist;same")
c.Draw()


# In[38]:


th232wallnof90_rate = th232wallnof90.Integral(6, 39)
th232wall_rate = th232wall.Integral(6, 39)
print(f"Yearly Th232: {th232wallnof90_rate:.2f}")
print(f"Yearly Th232 w/F90: {th232wall_rate:.2f}")
print(f"Ratio: {th232wallnof90_rate / th232wall_rate}")


# In[39]:


concretenof90 = k40wallnof90 + u238wallnof90 + th232wallnof90
concrete = k40wall + u238wall + th232wall
c = ROOT.TCanvas()
concretenof90.GetYaxis().SetRangeUser(0, 15000)
concretenof90.Draw("hist")
concrete.SetLineColorAlpha(ROOT.kRed)
concrete.Draw("hist;same")
c.Draw()


# In[40]:


concretenof90_rate = concretenof90.Integral(6, 39)
concrete_rate = concrete.Integral(6, 39)
print(f"Yearly Concrete: {concretenof90_rate:.2f}")
print(f"Yearly Concrete w/F90: {concrete_rate:.2f}")
print(f"Ratio: {concretenof90_rate / concrete_rate}")


# In[41]:


f_neutrons = ROOT.TFile.Open("../data/neutrons_side_morestats.root")
neutrons = f_neutrons.Get("neutrons_side_morestats")
neutronsnof90 = f_neutrons.Get("neutrons_side_morestats_nof90")
c = ROOT.TCanvas()
neutronsnof90.Scale(brn_scale)
neutrons.Scale(brn_scale)
neutronsnof90.GetYaxis().SetRangeUser(0, 2000)
neutronsnof90.Draw("hist")
neutrons.SetLineColorAlpha(ROOT.kRed)
neutrons.Draw("hist;same")
c.Draw()


# In[42]:


neutronsnof90_rate = neutronsnof90.Integral(6, 39)
neutrons_rate = neutrons.Integral(6, 39)
print(f"Yearly Neutrons: {neutronsnof90_rate:.2f}")
print(f"Yearly Neutrons w/F90: {neutrons_rate:.2f}")
print(f"Ratio: {neutronsnof90_rate / neutrons_rate}")


# In[51]:


hs = ROOT.THStack("hs", "COH-Ar-750 Signal Predictions;Energy [keVee];Counts/SNS-year")
h_cevnsf90.SetFillColor(ROOT.kRed+1)
hs.Add(h_cevnsf90)
h_ar39f90.SetFillColor(ROOT.kBlue+1)
hs.Add(h_ar39f90)
hog_phase2.SetFillColor(ROOT.kGreen+1)
hs.Add(hog_phase2)
pb210.SetFillColor(ROOT.kYellow+1)
hs.Add(pb210)
concrete.SetFillColor(ROOT.kOrange+1)
hs.Add(concrete)
neutrons.SetFillColor(ROOT.kSpring+1)
hs.Add(neutrons)
leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
leg.AddEntry(h_cevnsf90, "CEvNS", "f")
leg.AddEntry(h_ar39f90, "Ar39", "f")
leg.AddEntry(hog_phase2, "HOG Phase 2", "f")
leg.AddEntry(pb210, "Pb210", "f")
leg.AddEntry(concrete, "Concrete", "f")
leg.AddEntry(neutrons, "BRNs", "f")
c = ROOT.TCanvas()
c.SetLogy()
hs.Draw("hist")
hs.SetMinimum(1)
leg.Draw()
c.Draw()


# In[ ]:




