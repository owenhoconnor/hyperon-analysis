#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamCCNC()
{
//=========Macro generated from canvas: c3/True CCNC
//=========  (Mon Aug 10 13:57:24 2026) by ROOT version 6.28/12
   TCanvas *c3 = new TCanvas("c3", "True CCNC",0,0,2000,2000);
   gStyle->SetOptStat(0);
   c3->Range(-0.75,-442.1813,1.75,3979.631);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);
   
   TH1F *hTrueCCNC__3 = new TH1F("hTrueCCNC__3","",2,-0.5,1.5);
   hTrueCCNC__3->SetBinContent(1,3369);
   hTrueCCNC__3->SetBinContent(2,931);
   hTrueCCNC__3->SetMinimum(0);
   hTrueCCNC__3->SetEntries(4300);
   hTrueCCNC__3->SetDirectory(nullptr);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 1179;
   color = new TColor(ci, 0.2, 0.6, 1, " ", 0.7);
   hTrueCCNC__3->SetFillColor(ci);
   hTrueCCNC__3->SetLineWidth(2);
   hTrueCCNC__3->GetXaxis()->SetTitle("Interaction Current");
   hTrueCCNC__3->GetXaxis()->SetBinLabel(1,"CC");
   hTrueCCNC__3->GetXaxis()->SetBinLabel(2,"NC");
   hTrueCCNC__3->GetXaxis()->CenterTitle(true);
   hTrueCCNC__3->GetXaxis()->SetLabelFont(42);
   hTrueCCNC__3->GetXaxis()->SetTitleOffset(1);
   hTrueCCNC__3->GetXaxis()->SetTitleFont(42);
   hTrueCCNC__3->GetYaxis()->SetTitle("Interactions");
   hTrueCCNC__3->GetYaxis()->CenterTitle(true);
   hTrueCCNC__3->GetYaxis()->SetLabelFont(42);
   hTrueCCNC__3->GetYaxis()->SetTitleFont(42);
   hTrueCCNC__3->GetZaxis()->SetLabelFont(42);
   hTrueCCNC__3->GetZaxis()->SetTitleOffset(1);
   hTrueCCNC__3->GetZaxis()->SetTitleFont(42);
   hTrueCCNC__3->Draw("HIST TEXT0");
   c3->Modified();
   c3->cd();
   c3->SetSelected(c3);
}
