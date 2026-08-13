#ifdef __CLING__
#pragma cling optimize(0)
#endif
void nTracksShowersSig()
{
//=========Macro generated from canvas: c1/Track-Shower Topology
//=========  (Thu Aug 13 09:48:51 2026) by ROOT version 6.40.02
   TCanvas *c1 = new TCanvas("c1", "Track-Shower Topology", 0, 0, 2000, 2000);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(63, nullptr);
   c1->Range(-2,-1.9,8,7.433333);
   c1->SetFillColor(0);
   c1->SetFillStyle(1001);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLeftMargin(0.15);
   c1->SetRightMargin(0.15);
   c1->SetBottomMargin(0.15);
   
   TH2F *hTracksShowersSig = new TH2F("hTracksShowersSig", "", 7, -0.5, 6.5, 7, -0.5, 6.5);
   hTracksShowersSig->SetDirectory(nullptr);
   hTracksShowersSig->SetContour(20);
   
   TPaletteAxis *palette = new TPaletteAxis(6.55, -0.5, 7, 6.5, hTracksShowersSig);
   palette->SetNdivisions(510);
   palette->SetAxisColor(1);
   palette->SetLabelColor(1);
   palette->SetLabelFont(42);
   palette->SetLabelOffset(0.004999999888241291);
   palette->SetLabelSize(0.03999999910593033);
   palette->SetMaxDigits(0);
   palette->SetTickLength(0.02999999932944775);
   palette->SetTitleOffset(1);
   palette->SetTitleSize(0.03999999910593033);
   palette->SetTitleColor(1);
   palette->SetTitleFont(42);
   palette->SetTitle("Events");
   palette->SetFillColor(TColor::GetColor("#687557"));
   palette->SetFillStyle(1001);
   hTracksShowersSig->GetListOfFunctions()->Add(palette,"br");
   hTracksShowersSig->SetFillColor(0);
   hTracksShowersSig->SetFillStyle(1001);
   hTracksShowersSig->SetLineColor(TColor::GetColor("#000099"));
   hTracksShowersSig->GetXaxis()->SetTitle("Number of tracks");
   hTracksShowersSig->GetXaxis()->CenterTitle(true);
   hTracksShowersSig->GetXaxis()->SetNdivisions(7);
   hTracksShowersSig->GetXaxis()->SetLabelFont(42);
   hTracksShowersSig->GetXaxis()->SetLabelSize(0.03999999910593033);
   hTracksShowersSig->GetXaxis()->SetTitleSize(0.03999999910593033);
   hTracksShowersSig->GetXaxis()->SetTitleOffset(1);
   hTracksShowersSig->GetXaxis()->SetTitleFont(42);
   hTracksShowersSig->GetYaxis()->SetTitle("Number of showers");
   hTracksShowersSig->GetYaxis()->CenterTitle(true);
   hTracksShowersSig->GetYaxis()->SetNdivisions(7);
   hTracksShowersSig->GetYaxis()->SetLabelFont(42);
   hTracksShowersSig->GetYaxis()->SetLabelSize(0.03999999910593033);
   hTracksShowersSig->GetYaxis()->SetTitleSize(0.03999999910593033);
   hTracksShowersSig->GetYaxis()->SetTitleFont(42);
   hTracksShowersSig->GetZaxis()->SetTitle("Events");
   hTracksShowersSig->GetZaxis()->SetLabelFont(42);
   hTracksShowersSig->GetZaxis()->SetLabelSize(0.03999999910593033);
   hTracksShowersSig->GetZaxis()->SetTitleSize(0.03999999910593033);
   hTracksShowersSig->GetZaxis()->SetTitleOffset(1);
   hTracksShowersSig->GetZaxis()->SetTitleFont(42);
   hTracksShowersSig->Draw("COLZ TEXT");
   c1->Modified();
   c1->SetSelected(c1);
}
