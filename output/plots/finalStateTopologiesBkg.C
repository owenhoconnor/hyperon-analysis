#ifdef __CLING__
#pragma cling optimize(0)
#endif
void finalStateTopologiesBkg()
{
//=========Macro generated from canvas: cTopology/Final-state topologies
//=========  (Wed Aug 12 17:54:13 2026) by ROOT version 6.40.02
   TCanvas *cTopology = new TCanvas("cTopology", "Final-state topologies", 0, 0, 1800, 1100);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   cTopology->Range(-2.024096,-184.8,14.84337,431.2);
   cTopology->SetFillColor(0);
   cTopology->SetFillStyle(1001);
   cTopology->SetBorderMode(0);
   cTopology->SetBorderSize(2);
   cTopology->SetLeftMargin(0.12);
   cTopology->SetRightMargin(0.05);
   cTopology->SetBottomMargin(0.3);
   
   TH1I *hTopology__4 = new TH1I("hTopology", "Background Final State Topology Distribution After Presel", 14, 0, 14);
   hTopology__4->SetBinContent(1,50);
   hTopology__4->SetBinContent(2,308);
   hTopology__4->SetBinContent(3,110);
   hTopology__4->SetBinContent(4,297);
   hTopology__4->SetBinContent(5,15);
   hTopology__4->SetBinContent(6,84);
   hTopology__4->SetBinContent(7,3);
   hTopology__4->SetBinContent(8,6);
   hTopology__4->SetBinContent(9,10);
   hTopology__4->SetBinContent(10,19);
   hTopology__4->SetBinContent(11,17);
   hTopology__4->SetBinContent(12,3);
   hTopology__4->SetBinContent(13,15);
   hTopology__4->SetBinContent(14,16);
   hTopology__4->SetMinimum(0);
   hTopology__4->SetMaximum(369.6);
   hTopology__4->SetEntries(14);
   hTopology__4->SetStats(0);
   hTopology__4->SetFillColor(TColor::GetColor("#3399ffb2"));
   hTopology__4->SetFillStyle(1001);
   hTopology__4->SetLineWidth(2);
   hTopology__4->GetXaxis()->SetTitle("True final state topology");
   hTopology__4->GetXaxis()->SetBinLabel(1, "CC 0#pi, 1p");
   hTopology__4->GetXaxis()->SetBinLabel(2, "CC 0#pi, #geq2p");
   hTopology__4->GetXaxis()->SetBinLabel(3, "CC 1#pi^{0}");
   hTopology__4->GetXaxis()->SetBinLabel(4, "CC 1#pi^{+}");
   hTopology__4->GetXaxis()->SetBinLabel(5, "CC 1#pi^{-}");
   hTopology__4->GetXaxis()->SetBinLabel(6, "CC multi-#pi");
   hTopology__4->GetXaxis()->SetBinLabel(7, "CC strange");
   hTopology__4->GetXaxis()->SetBinLabel(8, "NC 0#pi, 0p");
   hTopology__4->GetXaxis()->SetBinLabel(9, "NC 0#pi, 1p");
   hTopology__4->GetXaxis()->SetBinLabel(10, "NC 0#pi, #geq2p");
   hTopology__4->GetXaxis()->SetBinLabel(11, "NC 1#pi^{0}");
   hTopology__4->GetXaxis()->SetBinLabel(12, "NC 1#pi^{+}");
   hTopology__4->GetXaxis()->SetBinLabel(13, "NC 1#pi^{-}");
   hTopology__4->GetXaxis()->SetBinLabel(14, "NC multi-#pi");
   hTopology__4->GetXaxis()->SetBit(TAxis::kLabelsVert);
   hTopology__4->GetXaxis()->CenterTitle(true);
   hTopology__4->GetXaxis()->SetLabelFont(42);
   hTopology__4->GetXaxis()->SetTitleOffset(3.5);
   hTopology__4->GetXaxis()->SetTitleFont(42);
   hTopology__4->GetYaxis()->SetTitle("Interactions");
   hTopology__4->GetYaxis()->CenterTitle(true);
   hTopology__4->GetYaxis()->SetLabelFont(42);
   hTopology__4->GetYaxis()->SetTitleOffset(1.299999952316284);
   hTopology__4->GetYaxis()->SetTitleFont(42);
   hTopology__4->GetZaxis()->SetLabelFont(42);
   hTopology__4->GetZaxis()->SetTitleOffset(1);
   hTopology__4->GetZaxis()->SetTitleFont(42);
   hTopology__4->Draw("HIST TEXT0");
   
   TPaveText *pt = new TPaveText(0.15, 0.934552, 0.85, 0.995, "blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_text1 = pt->AddText("Background Final State Topology Distribution After Presel");
   pt->Draw("blNDC");
   cTopology->Modified();
   cTopology->SetSelected(cTopology);
}
