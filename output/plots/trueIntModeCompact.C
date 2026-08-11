#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueIntModeCompact()
{
//=========Macro generated from canvas: c4/True Interaction Mode Compact
//=========  (Mon Aug 10 18:14:50 2026) by ROOT version 6.28/12
   TCanvas *c4 = new TCanvas("c4", "True Interaction Mode Compact",0,0,2000,2000);
   gStyle->SetOptStat(0);
   c4->Range(-0.125,-0.2625,1.125,2.3625);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   c4->SetFrameBorderMode(0);
   c4->SetFrameBorderMode(0);
   
   TH1F *hTrueIntModeCompact__4 = new TH1F("hTrueIntModeCompact__4","",1,0,1);
   hTrueIntModeCompact__4->SetBinContent(1,1);
   hTrueIntModeCompact__4->SetBinError(1,1);
   hTrueIntModeCompact__4->SetEntries(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTrueIntModeCompact__4->SetLineColor(ci);
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(1,"");
   hTrueIntModeCompact__4->GetXaxis()->SetBit(TAxis::kLabelsVert);
   hTrueIntModeCompact__4->GetXaxis()->SetLabelFont(42);
   hTrueIntModeCompact__4->GetXaxis()->SetTitleOffset(1);
   hTrueIntModeCompact__4->GetXaxis()->SetTitleFont(42);
   hTrueIntModeCompact__4->GetYaxis()->SetLabelFont(42);
   hTrueIntModeCompact__4->GetYaxis()->SetTitleFont(42);
   hTrueIntModeCompact__4->GetZaxis()->SetLabelFont(42);
   hTrueIntModeCompact__4->GetZaxis()->SetTitleOffset(1);
   hTrueIntModeCompact__4->GetZaxis()->SetTitleFont(42);
   hTrueIntModeCompact__4->Draw("HIST TEXT0");
   c4->Modified();
   c4->cd();
   c4->SetSelected(c4);
}
