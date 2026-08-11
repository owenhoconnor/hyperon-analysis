#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamIntMode()
{
//=========Macro generated from canvas: c2/True Interaction Mode
//=========  (Tue Aug 11 11:24:54 2026) by ROOT version 6.40.02
   TCanvas *c2 = new TCanvas("c2", "True Interaction Mode", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(1111);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c2->Range(-1.375,-97.12501,11.375,874.125);
   c2->SetFillColor(0);
   c2->SetFillStyle(1001);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   
   TH1F *hTrueIntMode__2 = new TH1F("hTrueIntMode", "", 10, -0.1, 10.1);
   hTrueIntMode__2->SetBinContent(1,446);
   hTrueIntMode__2->SetBinContent(2,740);
   hTrueIntMode__2->SetBinContent(3,133);
   hTrueIntMode__2->SetBinContent(4,3);
   hTrueIntMode__2->SetBinContent(10,186);
   hTrueIntMode__2->SetEntries(1508);
   hTrueIntMode__2->SetDirectory(nullptr);
   
   TPaveStats *ptstats = new TPaveStats(0.78, 0.775, 0.98, 0.935, "brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetFillStyle(1001);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *ptstats_text4 = ptstats->AddText("hTrueIntMode");
   ptstats_text4->SetTextSize(0.03680000081658363);
   TText *ptstats_text5 = ptstats->AddText("Entries = 1508   ");
   TText *ptstats_text6 = ptstats->AddText("Mean  =  1.906");
   TText *ptstats_text7 = ptstats->AddText("Std Dev   =  3.092");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->SetParent(hTrueIntMode__2);
   hTrueIntMode__2->GetListOfFunctions()->Add(ptstats);
   hTrueIntMode__2->SetFillColor(TColor::GetColor("#3399ffb2"));
   hTrueIntMode__2->SetFillStyle(1001);
   hTrueIntMode__2->SetLineWidth(2);
   hTrueIntMode__2->GetXaxis()->SetBinLabel(1, "QE");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(2, "RES");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(3, "DIS");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(4, "Coh");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(5, "Coh Elastic");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(6, "e Scatter");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(7, "IMD Annih");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(8, "Inverse Beta");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(9, "Glashow Res");
   hTrueIntMode__2->GetXaxis()->SetBinLabel(10, "AMNuGamma");
   hTrueIntMode__2->GetXaxis()->SetRange(1, 10);
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
   c2->SetSelected(c2);
}
