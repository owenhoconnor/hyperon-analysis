#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueIntModeCompact()
{
//=========Macro generated from canvas: c4/True Interaction Mode Compact
//=========  (Mon Aug 10 13:57:24 2026) by ROOT version 6.28/12
   TCanvas *c4 = new TCanvas("c4", "True Interaction Mode Compact",0,0,2000,2000);
   gStyle->SetOptStat(0);
   c4->Range(-0.5,-236.775,4.5,2130.975);
   c4->SetFillColor(0);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   c4->SetFrameBorderMode(0);
   c4->SetFrameBorderMode(0);
   
   TH1F *hTrueIntModeCompact__4 = new TH1F("hTrueIntModeCompact__4","",4,0,4);
   hTrueIntModeCompact__4->SetBinContent(1,1804);
   hTrueIntModeCompact__4->SetBinContent(2,1620);
   hTrueIntModeCompact__4->SetBinContent(3,367);
   hTrueIntModeCompact__4->SetBinContent(4,10);
   hTrueIntModeCompact__4->SetBinError(1,42.47352);
   hTrueIntModeCompact__4->SetBinError(2,40.24922);
   hTrueIntModeCompact__4->SetBinError(3,19.15724);
   hTrueIntModeCompact__4->SetBinError(4,3.162278);
   hTrueIntModeCompact__4->SetEntries(4);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTrueIntModeCompact__4->SetLineColor(ci);
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(1,"");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(2,"");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(3,"");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(4,"");
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
