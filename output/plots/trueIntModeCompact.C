#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueIntModeCompact()
{
//=========Macro generated from canvas: c4/True Interaction Mode Compact
//=========  (Tue Aug 11 11:24:54 2026) by ROOT version 6.40.02
   TCanvas *c4 = new TCanvas("c4", "True Interaction Mode Compact", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c4->Range(-0.625,-97.12501,5.625,874.125);
   c4->SetFillColor(0);
   c4->SetFillStyle(1001);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   
   TH1F *hTrueIntModeCompact__4 = new TH1F("hTrueIntModeCompact", "", 5, 0, 5);
   hTrueIntModeCompact__4->SetBinContent(1,446);
   hTrueIntModeCompact__4->SetBinContent(2,740);
   hTrueIntModeCompact__4->SetBinContent(3,133);
   hTrueIntModeCompact__4->SetBinContent(4,3);
   hTrueIntModeCompact__4->SetBinContent(5,186);
   hTrueIntModeCompact__4->SetBinError(1,21.11871208194287);
   hTrueIntModeCompact__4->SetBinError(2,27.20294101747089);
   hTrueIntModeCompact__4->SetBinError(3,11.5325625946708);
   hTrueIntModeCompact__4->SetBinError(4,1.732050807568877);
   hTrueIntModeCompact__4->SetBinError(5,13.63818169698586);
   hTrueIntModeCompact__4->SetEntries(5);
   hTrueIntModeCompact__4->SetFillColor(0);
   hTrueIntModeCompact__4->SetFillStyle(1001);
   hTrueIntModeCompact__4->SetLineColor(TColor::GetColor("#000099"));
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(1, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(2, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(3, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(4, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(5, "");
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
   c4->SetSelected(c4);
}
