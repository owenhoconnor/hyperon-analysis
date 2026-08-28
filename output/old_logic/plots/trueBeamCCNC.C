#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamCCNC()
{
//=========Macro generated from canvas: c3/True CCNC
//=========  (Wed Aug 12 17:54:13 2026) by ROOT version 6.40.02
   TCanvas *c3 = new TCanvas("c3", "True CCNC", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c3->Range(-0.75,-113.7938,1.75,1024.144);
   c3->SetFillColor(0);
   c3->SetFillStyle(1001);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   
   TH1F *hTrueCCNC__2 = new TH1F("hTrueCCNC", "", 2, -0.5, 1.5);
   hTrueCCNC__2->SetBinContent(1,867);
   hTrueCCNC__2->SetBinContent(2,86);
   hTrueCCNC__2->SetMinimum(0);
   hTrueCCNC__2->SetEntries(953);
   hTrueCCNC__2->SetFillColor(TColor::GetColor("#3399ffb2"));
   hTrueCCNC__2->SetFillStyle(1001);
   hTrueCCNC__2->SetLineWidth(2);
   hTrueCCNC__2->GetXaxis()->SetTitle("Interaction Current");
   hTrueCCNC__2->GetXaxis()->SetBinLabel(1, "CC");
   hTrueCCNC__2->GetXaxis()->SetBinLabel(2, "NC");
   hTrueCCNC__2->GetXaxis()->CenterTitle(true);
   hTrueCCNC__2->GetXaxis()->SetLabelFont(42);
   hTrueCCNC__2->GetXaxis()->SetTitleOffset(1);
   hTrueCCNC__2->GetXaxis()->SetTitleFont(42);
   hTrueCCNC__2->GetYaxis()->SetTitle("Interactions");
   hTrueCCNC__2->GetYaxis()->CenterTitle(true);
   hTrueCCNC__2->GetYaxis()->SetLabelFont(42);
   hTrueCCNC__2->GetYaxis()->SetTitleFont(42);
   hTrueCCNC__2->GetZaxis()->SetLabelFont(42);
   hTrueCCNC__2->GetZaxis()->SetTitleOffset(1);
   hTrueCCNC__2->GetZaxis()->SetTitleFont(42);
   hTrueCCNC__2->Draw("HIST TEXT0");
   c3->Modified();
   c3->SetSelected(c3);
}
