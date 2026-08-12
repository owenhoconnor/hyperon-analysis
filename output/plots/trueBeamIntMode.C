#ifdef __CLING__
#pragma cling optimize(0)
#endif
void trueBeamIntMode()
{
//=========Macro generated from canvas: c2/True Interaction Mode
//=========  (Wed Aug 12 17:54:12 2026) by ROOT version 6.40.02
   TCanvas *c2 = new TCanvas("c2", "True Interaction Mode", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   c2->Range(-2.25,-78.75001,15.25,708.75);
   c2->SetFillColor(0);
   c2->SetFillStyle(1001);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   
   TH1F *hTrueIntMode__1 = new TH1F("hTrueIntMode", "", 14, -0.5, 13.5);
   hTrueIntMode__1->SetBinContent(1,157);
   hTrueIntMode__1->SetBinContent(2,600);
   hTrueIntMode__1->SetBinContent(3,85);
   hTrueIntMode__1->SetBinContent(4,1);
   hTrueIntMode__1->SetBinContent(11,110);
   hTrueIntMode__1->SetEntries(953);
   hTrueIntMode__1->SetFillColor(TColor::GetColor("#3399ffb2"));
   hTrueIntMode__1->SetFillStyle(1001);
   hTrueIntMode__1->SetLineWidth(2);
   hTrueIntMode__1->GetXaxis()->SetBinLabel(1, "QE");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(2, "RES");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(3, "DIS");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(4, "Coh");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(5, "Coh Elastic");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(6, "e Scatter");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(7, "IMD Annih");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(8, "Inverse Beta");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(9, "Glashow Res");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(10, "AMNuGamma");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(11, "MEC");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(12, "Diffractive");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(13, "EM");
   hTrueIntMode__1->GetXaxis()->SetBinLabel(14, "Weak Mix");
   hTrueIntMode__1->GetXaxis()->SetLabelFont(42);
   hTrueIntMode__1->GetXaxis()->SetTitleOffset(1);
   hTrueIntMode__1->GetXaxis()->SetTitleFont(42);
   hTrueIntMode__1->GetYaxis()->SetLabelFont(42);
   hTrueIntMode__1->GetYaxis()->SetTitleFont(42);
   hTrueIntMode__1->GetZaxis()->SetLabelFont(42);
   hTrueIntMode__1->GetZaxis()->SetTitleOffset(1);
   hTrueIntMode__1->GetZaxis()->SetTitleFont(42);
   hTrueIntMode__1->Draw("HIST");
   c2->Modified();
   c2->SetSelected(c2);
}
