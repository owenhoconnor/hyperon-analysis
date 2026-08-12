#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueIntModeCompact()
{
//=========Macro generated from canvas: c4/True Interaction Mode Compact
//=========  (Wed Aug 12 12:33:09 2026) by ROOT version 6.40.02
   TCanvas *c4 = new TCanvas("c4", "True Interaction Mode Compact", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c4->Range(-0.625,-78.75001,5.625,708.75);
   c4->SetFillColor(0);
   c4->SetFillStyle(1001);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   
   TH1F *hTrueIntModeCompact__4 = new TH1F("hTrueIntModeCompact", "", 5, 0, 5);
   hTrueIntModeCompact__4->SetBinContent(1,157);
   hTrueIntModeCompact__4->SetBinContent(2,600);
   hTrueIntModeCompact__4->SetBinContent(3,85);
   hTrueIntModeCompact__4->SetBinContent(4,1);
   hTrueIntModeCompact__4->SetBinContent(5,110);
   hTrueIntModeCompact__4->SetBinError(1,12.52996408614167);
   hTrueIntModeCompact__4->SetBinError(2,24.49489742783178);
   hTrueIntModeCompact__4->SetBinError(3,9.219544457292887);
   hTrueIntModeCompact__4->SetBinError(4,1);
   hTrueIntModeCompact__4->SetBinError(5,10.48808848170152);
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
