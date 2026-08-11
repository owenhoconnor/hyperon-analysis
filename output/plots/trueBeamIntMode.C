#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamIntMode()
{
//=========Macro generated from canvas: c2/True Interaction Mode
//=========  (Mon Aug 10 18:14:50 2026) by ROOT version 6.28/12
   TCanvas *c2 = new TCanvas("c2", "True Interaction Mode",0,0,2000,2000);
   c2->Range(-0.275,-0.13125,2.475,1.18125);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   
   TH1F *hTrueIntMode__2 = new TH1F("hTrueIntMode__2","",10,0,2.2);
   hTrueIntMode__2->SetBinContent(5,1);
   hTrueIntMode__2->SetEntries(1);
   hTrueIntMode__2->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_LaTex = ptstats->AddText("hTrueIntMode");
   ptstats_LaTex->SetTextSize(0.0368);
   ptstats_LaTex = ptstats->AddText("Entries = 1      ");
   ptstats_LaTex = ptstats->AddText("Mean  =      1");
   ptstats_LaTex = ptstats->AddText("Std Dev   =      0");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hTrueIntMode__2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hTrueIntMode__2);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 1179;
   color = new TColor(ci, 0.2, 0.6, 1, " ", 0.7);
   hTrueIntMode__2->SetFillColor(ci);
   hTrueIntMode__2->SetLineWidth(2);
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
