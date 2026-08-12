#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamIntType()
{
//=========Macro generated from canvas: c1/True Interaction Type
//=========  (Wed Aug 12 12:33:09 2026) by ROOT version 6.40.02
   TCanvas *c1 = new TCanvas("c1", "True Interaction Type", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(1111);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c1->Range(986.4087,-113.7938,1112.621,1024.144);
   c1->SetFillColor(0);
   c1->SetFillStyle(1001);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   
   TH1F *hTrueIntType__1 = new TH1F("hTrueIntType", "", 10, 999.03, 1100);
   hTrueIntType__1->SetBinContent(1,867);
   hTrueIntType__1->SetBinContent(10,86);
   hTrueIntType__1->SetEntries(953);
   hTrueIntType__1->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78, 0.775, 0.98, 0.935, "brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetFillStyle(1001);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_text0 = ptstats->AddText("hTrueIntType");
   ptstats_text0->SetTextSize(0.03680000081658363);
   TText *ptstats_text1 = ptstats->AddText("Entries = 953    ");
   TText *ptstats_text2 = ptstats->AddText("Mean  =   1011");
   TText *ptstats_text3 = ptstats->AddText("Std Dev   =  25.49");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->SetParent(hTrueIntType__1);
   hTrueIntType__1->GetListOfFunctions()->Add(ptstats);
   hTrueIntType__1->SetFillColor(0);
   hTrueIntType__1->SetFillStyle(1001);
   hTrueIntType__1->SetLineColor(TColor::GetColor("#000099"));
   hTrueIntType__1->GetXaxis()->SetRange(1, 10);
   hTrueIntType__1->GetXaxis()->SetLabelFont(42);
   hTrueIntType__1->GetXaxis()->SetTitleOffset(1);
   hTrueIntType__1->GetXaxis()->SetTitleFont(42);
   hTrueIntType__1->GetYaxis()->SetLabelFont(42);
   hTrueIntType__1->GetYaxis()->SetTitleFont(42);
   hTrueIntType__1->GetZaxis()->SetLabelFont(42);
   hTrueIntType__1->GetZaxis()->SetTitleOffset(1);
   hTrueIntType__1->GetZaxis()->SetTitleFont(42);
   hTrueIntType__1->Draw("HIST");
   c1->Modified();
   c1->SetSelected(c1);
}
