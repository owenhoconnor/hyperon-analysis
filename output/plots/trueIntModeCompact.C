#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueIntModeCompact()
{
//=========Macro generated from canvas: c4/True Interaction Mode Compact
//=========  (Wed Aug 12 17:54:13 2026) by ROOT version 6.40.02
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
   
   TH1F *hTrueIntModeCompact__3 = new TH1F("hTrueIntModeCompact", "", 5, 0, 5);
   hTrueIntModeCompact__3->SetBinContent(1,157);
   hTrueIntModeCompact__3->SetBinContent(2,600);
   hTrueIntModeCompact__3->SetBinContent(3,85);
   hTrueIntModeCompact__3->SetBinContent(4,1);
   hTrueIntModeCompact__3->SetBinContent(5,110);
   hTrueIntModeCompact__3->SetBinError(1,12.52996408614167);
   hTrueIntModeCompact__3->SetBinError(2,24.49489742783178);
   hTrueIntModeCompact__3->SetBinError(3,9.219544457292887);
   hTrueIntModeCompact__3->SetBinError(4,1);
   hTrueIntModeCompact__3->SetBinError(5,10.48808848170152);
   hTrueIntModeCompact__3->SetEntries(5);
   hTrueIntModeCompact__3->SetStats(0);
   hTrueIntModeCompact__3->SetFillColor(0);
   hTrueIntModeCompact__3->SetFillStyle(1001);
   hTrueIntModeCompact__3->SetLineColor(TColor::GetColor("#000099"));
   hTrueIntModeCompact__3->GetXaxis()->SetBinLabel(1, "QE");
   hTrueIntModeCompact__3->GetXaxis()->SetBinLabel(2, "RES");
   hTrueIntModeCompact__3->GetXaxis()->SetBinLabel(3, "DIS");
   hTrueIntModeCompact__3->GetXaxis()->SetBinLabel(4, "Coh");
   hTrueIntModeCompact__3->GetXaxis()->SetBinLabel(5, "MEC");
   hTrueIntModeCompact__3->GetXaxis()->SetBit(TAxis::kLabelsVert);
   hTrueIntModeCompact__3->GetXaxis()->SetLabelFont(42);
   hTrueIntModeCompact__3->GetXaxis()->SetTitleOffset(1);
   hTrueIntModeCompact__3->GetXaxis()->SetTitleFont(42);
   hTrueIntModeCompact__3->GetYaxis()->SetLabelFont(42);
   hTrueIntModeCompact__3->GetYaxis()->SetTitleFont(42);
   hTrueIntModeCompact__3->GetZaxis()->SetLabelFont(42);
   hTrueIntModeCompact__3->GetZaxis()->SetTitleOffset(1);
   hTrueIntModeCompact__3->GetZaxis()->SetTitleFont(42);
   hTrueIntModeCompact__3->Draw("HIST TEXT0");
   c4->Modified();
   c4->SetSelected(c4);
}
