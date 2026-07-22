#ifdef __CLING__
#pragma cling optimize(0)
#endif
void nTracksShowersSig()
{
//=========Macro generated from canvas: c1/Track-Shower Topology
//=========  (Wed Jul 22 15:11:07 2026) by ROOT version 6.40.02
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
   hTracksShowersSig->SetBinContent(11,349);
   hTracksShowersSig->SetBinContent(12,628);
   hTracksShowersSig->SetBinContent(13,589);
   hTracksShowersSig->SetBinContent(14,179);
   hTracksShowersSig->SetBinContent(15,12);
   hTracksShowersSig->SetBinContent(16,1);
   hTracksShowersSig->SetBinContent(20,1352);
   hTracksShowersSig->SetBinContent(21,2728);
   hTracksShowersSig->SetBinContent(22,2189);
   hTracksShowersSig->SetBinContent(23,187);
   hTracksShowersSig->SetBinContent(24,8);
   hTracksShowersSig->SetBinContent(25,1);
   hTracksShowersSig->SetBinContent(29,2437);
   hTracksShowersSig->SetBinContent(30,1946);
   hTracksShowersSig->SetBinContent(31,555);
   hTracksShowersSig->SetBinContent(32,61);
   hTracksShowersSig->SetBinContent(33,1);
   hTracksShowersSig->SetBinContent(38,2165);
   hTracksShowersSig->SetBinContent(39,615);
   hTracksShowersSig->SetBinContent(40,99);
   hTracksShowersSig->SetBinContent(41,10);
   hTracksShowersSig->SetBinContent(42,1);
   hTracksShowersSig->SetBinContent(47,475);
   hTracksShowersSig->SetBinContent(48,123);
   hTracksShowersSig->SetBinContent(49,20);
   hTracksShowersSig->SetBinContent(50,6);
   hTracksShowersSig->SetBinContent(56,73);
   hTracksShowersSig->SetBinContent(57,10);
   hTracksShowersSig->SetBinContent(58,2);
   hTracksShowersSig->SetBinContent(65,5);
   hTracksShowersSig->SetBinContent(66,2);
   hTracksShowersSig->SetBinContent(67,1);
   hTracksShowersSig->SetBinContent(74,1);
   hTracksShowersSig->SetEntries(16831);
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
   palette->SetFillColor(TColor::GetColor("#efefd7"));
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
