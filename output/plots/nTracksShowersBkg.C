#ifdef __CLING__
#pragma cling optimize(0)
#endif
void nTracksShowersBkg()
{
//=========Macro generated from canvas: c2/Background Track-Shower Topology
//=========  (Tue Aug 11 11:58:53 2026) by ROOT version 6.28/12
   TCanvas *c2 = new TCanvas("c2", "Background Track-Shower Topology",0,0,2000,2000);
   gStyle->SetOptStat(0);
   c2->Range(-2,-1.9,8,7.433333);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetLeftMargin(0.15);
   c2->SetRightMargin(0.15);
   c2->SetBottomMargin(0.15);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   
   TH2F *nTracksShowersBkg = new TH2F("nTracksShowersBkg","",7,-0.5,6.5,7,-0.5,6.5);
   nTracksShowersBkg->SetBinContent(11,28);
   nTracksShowersBkg->SetBinContent(12,18);
   nTracksShowersBkg->SetBinContent(13,2);
   nTracksShowersBkg->SetBinContent(14,2);
   nTracksShowersBkg->SetBinContent(20,4);
   nTracksShowersBkg->SetBinContent(21,4);
   nTracksShowersBkg->SetBinContent(22,1);
   nTracksShowersBkg->SetBinContent(25,1);
   nTracksShowersBkg->SetBinContent(29,1);
   nTracksShowersBkg->SetBinContent(30,1);
   nTracksShowersBkg->SetBinContent(38,2);
   nTracksShowersBkg->SetBinContent(39,3);
   nTracksShowersBkg->SetEntries(67);
   nTracksShowersBkg->SetDirectory(nullptr);
   nTracksShowersBkg->SetContour(20);
   nTracksShowersBkg->SetContourLevel(0,0);
   nTracksShowersBkg->SetContourLevel(1,1.4);
   nTracksShowersBkg->SetContourLevel(2,2.8);
   nTracksShowersBkg->SetContourLevel(3,4.2);
   nTracksShowersBkg->SetContourLevel(4,5.6);
   nTracksShowersBkg->SetContourLevel(5,7);
   nTracksShowersBkg->SetContourLevel(6,8.4);
   nTracksShowersBkg->SetContourLevel(7,9.8);
   nTracksShowersBkg->SetContourLevel(8,11.2);
   nTracksShowersBkg->SetContourLevel(9,12.6);
   nTracksShowersBkg->SetContourLevel(10,14);
   nTracksShowersBkg->SetContourLevel(11,15.4);
   nTracksShowersBkg->SetContourLevel(12,16.8);
   nTracksShowersBkg->SetContourLevel(13,18.2);
   nTracksShowersBkg->SetContourLevel(14,19.6);
   nTracksShowersBkg->SetContourLevel(15,21);
   nTracksShowersBkg->SetContourLevel(16,22.4);
   nTracksShowersBkg->SetContourLevel(17,23.8);
   nTracksShowersBkg->SetContourLevel(18,25.2);
   nTracksShowersBkg->SetContourLevel(19,26.6);
   
   TPaletteAxis *palette = new TPaletteAxis(6.55,-0.5,7,6.5,nTracksShowersBkg);
   palette->SetNdivisions(510);
   palette->SetAxisColor(1);
   palette->SetLabelColor(1);
   palette->SetLabelFont(42);
   palette->SetLabelOffset(0.005);
   palette->SetLabelSize(0.04);
   palette->SetMaxDigits(0);
   palette->SetTickLength(0.03);
   palette->SetTitleOffset(1);
   palette->SetTitleSize(0.04);
   palette->SetTitleColor(1);
   palette->SetTitleFont(42);
   palette->SetTitle("Events");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ebef62");
   palette->SetFillColor(ci);
   palette->SetFillStyle(1001);
   nTracksShowersBkg->GetListOfFunctions()->Add(palette,"br");

   ci = TColor::GetColor("#000099");
   nTracksShowersBkg->SetLineColor(ci);
   nTracksShowersBkg->GetXaxis()->SetTitle("Number of tracks");
   nTracksShowersBkg->GetXaxis()->CenterTitle(true);
   nTracksShowersBkg->GetXaxis()->SetNdivisions(7);
   nTracksShowersBkg->GetXaxis()->SetLabelFont(42);
   nTracksShowersBkg->GetXaxis()->SetLabelSize(0.04);
   nTracksShowersBkg->GetXaxis()->SetTitleSize(0.04);
   nTracksShowersBkg->GetXaxis()->SetTitleOffset(1);
   nTracksShowersBkg->GetXaxis()->SetTitleFont(42);
   nTracksShowersBkg->GetYaxis()->SetTitle("Number of showers");
   nTracksShowersBkg->GetYaxis()->CenterTitle(true);
   nTracksShowersBkg->GetYaxis()->SetNdivisions(7);
   nTracksShowersBkg->GetYaxis()->SetLabelFont(42);
   nTracksShowersBkg->GetYaxis()->SetLabelSize(0.04);
   nTracksShowersBkg->GetYaxis()->SetTitleSize(0.04);
   nTracksShowersBkg->GetYaxis()->SetTitleFont(42);
   nTracksShowersBkg->GetZaxis()->SetTitle("Events");
   nTracksShowersBkg->GetZaxis()->SetLabelFont(42);
   nTracksShowersBkg->GetZaxis()->SetLabelSize(0.04);
   nTracksShowersBkg->GetZaxis()->SetTitleSize(0.04);
   nTracksShowersBkg->GetZaxis()->SetTitleOffset(1);
   nTracksShowersBkg->GetZaxis()->SetTitleFont(42);
   nTracksShowersBkg->Draw("COLZ TEXT");
   c2->Modified();
   c2->cd();
   c2->SetSelected(c2);
}
