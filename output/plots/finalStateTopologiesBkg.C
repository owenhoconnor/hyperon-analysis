#ifdef __CLING__
#pragma cling optimize(0)
#endif
void finalStateTopologiesBkg()
{
//=========Macro generated from canvas: cTopology/Final-state topologies
//=========  (Mon Aug 10 13:57:24 2026) by ROOT version 6.28/12
   TCanvas *cTopology = new TCanvas("cTopology", "Final-state topologies",0,0,1800,1100);
   gStyle->SetOptStat(0);
   cTopology->Range(-2.168675,-783.0001,15.90361,1827);
   cTopology->SetFillColor(0);
   cTopology->SetBorderMode(0);
   cTopology->SetBorderSize(2);
   cTopology->SetLeftMargin(0.12);
   cTopology->SetRightMargin(0.05);
   cTopology->SetBottomMargin(0.3);
   cTopology->SetFrameBorderMode(0);
   cTopology->SetFrameBorderMode(0);
   
   TH1I *hTopology__5 = new TH1I("hTopology__5","Background Final State Topology Distribution After Presel",15,0,15);
   hTopology__5->SetBinContent(1,35);
   hTopology__5->SetBinContent(2,1305);
   hTopology__5->SetBinContent(3,400);
   hTopology__5->SetBinContent(4,760);
   hTopology__5->SetBinContent(5,40);
   hTopology__5->SetBinContent(6,809);
   hTopology__5->SetBinContent(7,20);
   hTopology__5->SetBinContent(8,6);
   hTopology__5->SetBinContent(9,10);
   hTopology__5->SetBinContent(10,356);
   hTopology__5->SetBinContent(11,116);
   hTopology__5->SetBinContent(12,163);
   hTopology__5->SetBinContent(13,33);
   hTopology__5->SetBinContent(14,242);
   hTopology__5->SetBinContent(15,5);
   hTopology__5->SetMinimum(0);
   hTopology__5->SetMaximum(1566);
   hTopology__5->SetEntries(15);
   hTopology__5->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hTopology__5->SetLineColor(ci);
   hTopology__5->GetXaxis()->SetTitle("True final state topology");
   hTopology__5->GetXaxis()->SetBinLabel(1,"CC 0#pi, 1p");
   hTopology__5->GetXaxis()->SetBinLabel(2,"CC 0#pi, #geq2p");
   hTopology__5->GetXaxis()->SetBinLabel(3,"CC 1#pi^{0}");
   hTopology__5->GetXaxis()->SetBinLabel(4,"CC 1#pi^{+}");
   hTopology__5->GetXaxis()->SetBinLabel(5,"CC 1#pi^{-}");
   hTopology__5->GetXaxis()->SetBinLabel(6,"CC multi-#pi");
   hTopology__5->GetXaxis()->SetBinLabel(7,"CC strange");
   hTopology__5->GetXaxis()->SetBinLabel(8,"NC 0#pi, 0p");
   hTopology__5->GetXaxis()->SetBinLabel(9,"NC 0#pi, 1p");
   hTopology__5->GetXaxis()->SetBinLabel(10,"NC 0#pi, #geq2p");
   hTopology__5->GetXaxis()->SetBinLabel(11,"NC 1#pi^{0}");
   hTopology__5->GetXaxis()->SetBinLabel(12,"NC 1#pi^{+}");
   hTopology__5->GetXaxis()->SetBinLabel(13,"NC 1#pi^{-}");
   hTopology__5->GetXaxis()->SetBinLabel(14,"NC multi-#pi");
   hTopology__5->GetXaxis()->SetBinLabel(15,"NC strange");
   hTopology__5->GetXaxis()->SetBit(TAxis::kLabelsVert);
   hTopology__5->GetXaxis()->CenterTitle(true);
   hTopology__5->GetXaxis()->SetLabelFont(42);
   hTopology__5->GetXaxis()->SetTitleOffset(3.5);
   hTopology__5->GetXaxis()->SetTitleFont(42);
   hTopology__5->GetYaxis()->SetTitle("Interactions");
   hTopology__5->GetYaxis()->CenterTitle(true);
   hTopology__5->GetYaxis()->SetLabelFont(42);
   hTopology__5->GetYaxis()->SetTitleOffset(1.3);
   hTopology__5->GetYaxis()->SetTitleFont(42);
   hTopology__5->GetZaxis()->SetLabelFont(42);
   hTopology__5->GetZaxis()->SetTitleOffset(1);
   hTopology__5->GetZaxis()->SetTitleFont(42);
   hTopology__5->Draw("HIST TEXT0");
   
   TPaveText *pt = new TPaveText(0.15,0.9345522,0.85,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("Background Final State Topology Distribution After Presel");
   pt->Draw();
   cTopology->Modified();
   cTopology->cd();
   cTopology->SetSelected(cTopology);
}
