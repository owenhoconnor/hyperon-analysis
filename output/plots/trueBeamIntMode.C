#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamIntMode()
{
//=========Macro generated from canvas: c2/True Interaction Mode
//=========  (Mon Aug 10 13:57:23 2026) by ROOT version 6.28/12
   TCanvas *c2 = new TCanvas("c2", "True Interaction Mode",0,0,2000,2000);
   c2->Range(-1.25,-236.775,11.25,2130.975);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   
   TH1F *hTrueIntMode__2 = new TH1F("hTrueIntMode__2","",10,0,10);
   hTrueIntMode__2->SetBinContent(1,1804);
   hTrueIntMode__2->SetBinContent(2,1620);
   hTrueIntMode__2->SetBinContent(3,367);
   hTrueIntMode__2->SetBinContent(4,10);
   hTrueIntMode__2->SetBinContent(11,499);
   hTrueIntMode__2->SetEntries(4300);
   hTrueIntMode__2->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("hTrueIntMode");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 4300   ");
   ptstats_LaTex = ptstats->AddText("Mean  = 0.6272");
   ptstats_LaTex = ptstats->AddText("Std Dev   = 0.6654");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hTrueIntMode__2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hTrueIntMode__2);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTrueIntMode__2->SetLineColor(ci);
   hTrueIntMode__2->GetXaxis()->SetBinLabel(1,"QE");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(2,"RES");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(3,"DIS");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(4,"Coh");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(5,"Coh Elastic");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(6,"e Scatter");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(7,"IMD Annih");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(8,"Inverse Beta");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(9,"Glashow Res");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(10,"AMNuGamma");
   hTrueIntMode__2->GetXaxis()->SetRange(1,10);
   hTrueIntMode__2->GetXaxis()->SetLabelFont(42);
   hTrueIntMode__2->GetXaxis()->SetTitleOffset(1);
   hTrueIntMode__2->GetXaxis()->SetTitleFont(42);
   hTrueIntMode__2->GetYaxis()->SetLabelFont(42);
   hTrueIntMode__2->GetYaxis()->SetTitleFont(42);
   hTrueIntMode__2->GetZaxis()->SetLabelFont(42);
   hTrueIntMode__2->GetZaxis()->SetTitleOffset(1);
   hTrueIntMode__2->GetZaxis()->SetTitleFont(42);
   hTrueIntMode__2->Draw("HIST");
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);
}
