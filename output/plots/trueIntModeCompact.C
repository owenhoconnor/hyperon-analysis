#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueIntModeCompact()
{
//=========Macro generated from canvas: c4/True Interaction Mode Compact
//=========  (Thu Jul 30 10:01:53 2026) by ROOT version 6.40.02
   TCanvas *c4 = new TCanvas("c4", "True Interaction Mode Compact", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(1111);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c4->Range(-0.7500001,-5503.182,6.75,49528.63);
   c4->SetFillColor(0);
   c4->SetFillStyle(1001);
   c4->SetBorderMode(0);
   c4->SetBorderSize(2);
   
   TH1F *hTrueIntModeCompact__4 = new TH1F("hTrueIntModeCompact", "", 6, 0, 6);
   hTrueIntModeCompact__4->SetBinContent(1,41929);
   hTrueIntModeCompact__4->SetBinContent(2,21958);
   hTrueIntModeCompact__4->SetBinContent(3,6679);
   hTrueIntModeCompact__4->SetBinContent(4,238);
   hTrueIntModeCompact__4->SetBinContent(5,3);
   hTrueIntModeCompact__4->SetBinContent(6,9255);
   hTrueIntModeCompact__4->SetBinError(1,204.765719787273);
   hTrueIntModeCompact__4->SetBinError(2,148.1823201330037);
   hTrueIntModeCompact__4->SetBinError(3,81.72514912803769);
   hTrueIntModeCompact__4->SetBinError(4,15.42724862054151);
   hTrueIntModeCompact__4->SetBinError(5,1.732050807568877);
   hTrueIntModeCompact__4->SetBinError(6,96.2029105588807);
   hTrueIntModeCompact__4->SetEntries(6);
   
   TPaveStats *ptstats = new TPaveStats(0.78, 0.775, 0.98, 0.935, "brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetFillStyle(1001);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_text12 = ptstats->AddText("hTrueIntModeCompact");
   ptstats_text12->SetTextSize(0.03680000081658363);
   TText *ptstats_text13 = ptstats->AddText("Entries = 6      ");
   TText *ptstats_text14 = ptstats->AddText("Mean  =      0");
   TText *ptstats_text15 = ptstats->AddText("Std Dev   =      0");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->SetParent(hTrueIntModeCompact__4);
   hTrueIntModeCompact__4->GetListOfFunctions()->Add(ptstats);
   hTrueIntModeCompact__4->SetFillColor(0);
   hTrueIntModeCompact__4->SetFillStyle(1001);
   hTrueIntModeCompact__4->SetLineColor(TColor::GetColor("#000099"));
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(1, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(2, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(3, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(4, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(5, "");
   hTrueIntModeCompact__4->GetXaxis()->SetBinLabel(6, "");
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
