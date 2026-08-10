#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamIntType()
{
//=========Macro generated from canvas: c1/True Interaction Type
//=========  (Thu Aug  6 10:24:11 2026) by ROOT version 6.28/12
   TCanvas *c1 = new TCanvas("c1", "True Interaction Type",0,0,2000,2000);
   c1->Range(987.5,-513.8438,1112.5,4624.594);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *hTrueIntType__1 = new TH1F("hTrueIntType__1","",10,1000,1100);
   hTrueIntType__1->SetBinContent(1,3915);
   hTrueIntType__1->SetBinContent(2,8);
   hTrueIntType__1->SetBinContent(10,377);
   hTrueIntType__1->SetEntries(4300);
   hTrueIntType__1->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("hTrueIntType");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 4300   ");
   ptstats_LaTex = ptstats->AddText("Mean  =   1010");
   ptstats_LaTex = ptstats->AddText("Std Dev   =  25.31");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hTrueIntType__1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hTrueIntType__1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTrueIntType__1->SetLineColor(ci);
   hTrueIntType__1->GetXaxis()->SetRange(1,10);
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
   c1->cd();
   c1->SetSelected(c1);
}
