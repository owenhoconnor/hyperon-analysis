#ifdef __CLING__
#pragma cling optimize(0)
#endif
void nTracksShowersSig()
{
//=========Macro generated from canvas: c1/Track-Shower Topology
//=========  (Tue Aug 11 11:58:53 2026) by ROOT version 6.28/12
   TCanvas *c1 = new TCanvas("c1", "Track-Shower Topology",0,0,2000,2000);
   gStyle->SetOptStat(0);
   c1->Range(-2,-1.9,8,7.433333);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.15);
   c1->SetRightMargin(0.15);
   c1->SetBottomMargin(0.15);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH2F *hTracksShowersSig = new TH2F("hTracksShowersSig","",7,-0.5,6.5,7,-0.5,6.5);
   hTracksShowersSig->SetDirectory(nullptr);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTracksShowersSig->SetLineColor(ci);
   hTracksShowersSig->GetXaxis()->SetTitle("Number of tracks");
   hTracksShowersSig->GetXaxis()->CenterTitle(true);
   hTracksShowersSig->GetXaxis()->SetNdivisions(7);
   hTracksShowersSig->GetXaxis()->SetLabelFont(42);
   hTracksShowersSig->GetXaxis()->SetLabelSize(0.04);
   hTracksShowersSig->GetXaxis()->SetTitleSize(0.04);
   hTracksShowersSig->GetXaxis()->SetTitleOffset(1);
   hTracksShowersSig->GetXaxis()->SetTitleFont(42);
   hTracksShowersSig->GetYaxis()->SetTitle("Number of showers");
   hTracksShowersSig->GetYaxis()->CenterTitle(true);
   hTracksShowersSig->GetYaxis()->SetNdivisions(7);
   hTracksShowersSig->GetYaxis()->SetLabelFont(42);
   hTracksShowersSig->GetYaxis()->SetLabelSize(0.04);
   hTracksShowersSig->GetYaxis()->SetTitleSize(0.04);
   hTracksShowersSig->GetYaxis()->SetTitleFont(42);
   hTracksShowersSig->GetZaxis()->SetTitle("Events");
   hTracksShowersSig->GetZaxis()->SetLabelFont(42);
   hTracksShowersSig->GetZaxis()->SetLabelSize(0.04);
   hTracksShowersSig->GetZaxis()->SetTitleSize(0.04);
   hTracksShowersSig->GetZaxis()->SetTitleOffset(1);
   hTracksShowersSig->GetZaxis()->SetTitleFont(42);
   hTracksShowersSig->Draw("COLZ TEXT");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
