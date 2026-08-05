#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamCCNC()
{
//=========Macro generated from canvas: c3/True CCNC
//=========  (Thu Jul 30 10:01:53 2026) by ROOT version 6.40.02
   TCanvas *c3 = new TCanvas("c3", "True CCNC", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(1111);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c3->Range(-0.1375,-7737.319,1.1375,69635.87);
   c3->SetFillColor(0);
   c3->SetFillStyle(1001);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   
   TH1F *hTrueCCNC__3 = new TH1F("hTrueCCNC", "", 10, -0.01, 1.01);
   hTrueCCNC__3->SetBinContent(1,58951);
   hTrueCCNC__3->SetBinContent(10,21113);
   hTrueCCNC__3->SetEntries(80064);
   hTrueCCNC__3->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78, 0.775, 0.98, 0.935, "brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetFillStyle(1001);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_text8 = ptstats->AddText("hTrueCCNC");
   ptstats_text8->SetTextSize(0.03680000081658363);
   TText *ptstats_text9 = ptstats->AddText("Entries = 80064  ");
   TText *ptstats_text10 = ptstats->AddText("Mean  = 0.2637");
   TText *ptstats_text11 = ptstats->AddText("Std Dev   = 0.4406");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->SetParent(hTrueCCNC__3);
   hTrueCCNC__3->GetListOfFunctions()->Add(ptstats);
   hTrueCCNC__3->SetFillColor(0);
   hTrueCCNC__3->SetFillStyle(1001);
   hTrueCCNC__3->SetLineColor(TColor::GetColor("#000099"));
   hTrueCCNC__3->GetXaxis()->SetRange(1, 10);
   hTrueCCNC__3->GetXaxis()->SetLabelFont(42);
   hTrueCCNC__3->GetXaxis()->SetTitleOffset(1);
   hTrueCCNC__3->GetXaxis()->SetTitleFont(42);
   hTrueCCNC__3->GetYaxis()->SetLabelFont(42);
   hTrueCCNC__3->GetYaxis()->SetTitleFont(42);
   hTrueCCNC__3->GetZaxis()->SetLabelFont(42);
   hTrueCCNC__3->GetZaxis()->SetTitleOffset(1);
   hTrueCCNC__3->GetZaxis()->SetTitleFont(42);
   hTrueCCNC__3->Draw("HIST");
   c3->Modified();
   c3->SetSelected(c3);
}
