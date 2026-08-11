#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamCCNC()
{
//=========Macro generated from canvas: c3/True CCNC
//=========  (Tue Aug 11 11:24:54 2026) by ROOT version 6.40.02
   TCanvas *c3 = new TCanvas("c3", "True CCNC", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c3->Range(-0.75,-164.1938,1.75,1477.744);
   c3->SetFillColor(0);
   c3->SetFillStyle(1001);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   
   TH1F *hTrueCCNC__3 = new TH1F("hTrueCCNC", "", 2, -0.5, 1.5);
   hTrueCCNC__3->SetBinContent(1,1251);
   hTrueCCNC__3->SetBinContent(2,257);
   hTrueCCNC__3->SetMinimum(0);
   hTrueCCNC__3->SetEntries(1508);
   hTrueCCNC__3->SetDirectory(nullptr);
   hTrueCCNC__3->SetFillColor(TColor::GetColor("#3399ffb2"));
   hTrueCCNC__3->SetFillStyle(1001);
   hTrueCCNC__3->SetLineWidth(2);
   hTrueCCNC__3->GetXaxis()->SetTitle("Interaction Current");
   hTrueCCNC__3->GetXaxis()->SetBinLabel(1, "CC");
   hTrueCCNC__3->GetXaxis()->SetBinLabel(2, "NC");
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
   c3->SetSelected(c3);
}
