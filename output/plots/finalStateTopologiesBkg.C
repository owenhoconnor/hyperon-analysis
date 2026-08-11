#ifdef __CLING__
#pragma cling optimize(0)
#endif
void finalStateTopologiesBkg()
{
//=========Macro generated from canvas: cTopology/Final-state topologies
//=========  (Tue Aug 11 11:24:54 2026) by ROOT version 6.40.02
   TCanvas *cTopology = new TCanvas("cTopology", "Final-state topologies", 0, 0, 1800, 1100);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(1);
   TColor::SetPalette(57, nullptr);
   cTopology->Range(-0.2891566,-750.6001,2.120482,1751.4);
   cTopology->SetFillColor(0);
   cTopology->SetFillStyle(1001);
   cTopology->SetBorderMode(0);
   cTopology->SetBorderSize(2);
   cTopology->SetLeftMargin(0.12);
   cTopology->SetRightMargin(0.05);
   cTopology->SetBottomMargin(0.3);
   
   TH1I *hTopology__5 = new TH1I("hTopology", "Background Final State Topology Distribution After Presel", 2, 0, 2);
   hTopology__5->SetBinContent(1,1251);
   hTopology__5->SetBinContent(2,257);
   hTopology__5->SetMinimum(0);
   hTopology__5->SetMaximum(1501.2);
   hTopology__5->SetEntries(2);
   hTopology__5->SetStats(0);
   hTopology__5->SetFillColor(TColor::GetColor("#3399ffb2"));
   hTopology__5->SetFillStyle(1001);
   hTopology__5->SetLineWidth(2);
   hTopology__5->GetXaxis()->SetTitle("True final state topology");
   hTopology__5->GetXaxis()->SetBinLabel(1, "CC 0#pi, 0p");
   hTopology__5->GetXaxis()->SetBinLabel(2, "NC 0#pi, 0p");
   hTopology__5->GetXaxis()->SetBit(TAxis::kLabelsVert);
   hTopology__5->GetXaxis()->CenterTitle(true);
   hTopology__5->GetXaxis()->SetLabelFont(42);
   hTopology__5->GetXaxis()->SetTitleOffset(3.5);
   hTopology__5->GetXaxis()->SetTitleFont(42);
   hTopology__5->GetYaxis()->SetTitle("Interactions");
   hTopology__5->GetYaxis()->CenterTitle(true);
   hTopology__5->GetYaxis()->SetLabelFont(42);
   hTopology__5->GetYaxis()->SetTitleOffset(1.299999952316284);
   hTopology__5->GetYaxis()->SetTitleFont(42);
   hTopology__5->GetZaxis()->SetLabelFont(42);
   hTopology__5->GetZaxis()->SetTitleOffset(1);
   hTopology__5->GetZaxis()->SetTitleFont(42);
   hTopology__5->Draw("HIST TEXT0");
   
   TPaveText *pt = new TPaveText(0.15, 0.934552, 0.85, 0.995, "blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_text8 = pt->AddText("Background Final State Topology Distribution After Presel");
   pt->Draw("blNDC");
   cTopology->Modified();
   cTopology->SetSelected(cTopology);
}
