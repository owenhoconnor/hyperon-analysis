#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamCCNC()
{
//=========Macro generated from canvas: c3/True CCNC
//=========  (Thu Aug  6 10:24:12 2026) by ROOT version 6.28/12
   TCanvas *c3 = new TCanvas("c3", "True CCNC",0,0,2000,2000);
   c3->Range(-0.125,-442.1813,1.125,3979.631);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);
   
   TH1F *hTrueCCNC__3 = new TH1F("hTrueCCNC__3","",10,0,1);
   hTrueCCNC__3->SetBinContent(1,3369);
   hTrueCCNC__3->SetBinContent(10,931);
   hTrueCCNC__3->SetEntries(4300);
   hTrueCCNC__3->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("hTrueCCNC");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 4300   ");
   ptstats_LaTex = ptstats->AddText("Mean  = 0.2165");
   ptstats_LaTex = ptstats->AddText("Std Dev   = 0.4119");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hTrueCCNC__3->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hTrueCCNC__3);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTrueCCNC__3->SetLineColor(ci);
   hTrueCCNC__3->GetXaxis()->SetRange(1,10);
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
   c3->cd();
   c3->SetSelected(c3);
}
