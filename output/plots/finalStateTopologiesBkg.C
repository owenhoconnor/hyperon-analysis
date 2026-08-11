#ifdef __CLING__
#pragma cling optimize(0)
#endif
void finalStateTopologiesBkg()
{
//=========Macro generated from canvas: cTopology/Final-state topologies
//=========  (Mon Aug 10 18:14:51 2026) by ROOT version 6.28/12
   TCanvas *cTopology = new TCanvas("cTopology", "Final-state topologies",0,0,1800,1100);
   gStyle->SetOptStat(0);
   cTopology->Range(-0.1445783,-0.6000001,1.060241,1.4);
   cTopology->SetFillColor(0);
   cTopology->SetBorderMode(0);
   cTopology->SetBorderSize(2);
   cTopology->SetLeftMargin(0.12);
   cTopology->SetRightMargin(0.05);
   cTopology->SetBottomMargin(0.3);
   cTopology->SetFrameBorderMode(0);
   cTopology->SetFrameBorderMode(0);
   
   TH1I *hTopology__5 = new TH1I("hTopology__5","Background Final State Topology Distribution After Presel",1,0,1);
   hTopology__5->SetBinContent(1,1);
   hTopology__5->SetMinimum(0);
   hTopology__5->SetMaximum(1.2);
   hTopology__5->SetEntries(1);
   hTopology__5->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = 1179;
   color = new TColor(ci, 0.2, 0.6, 1, " ", 0.7);
   hTopology__5->SetFillColor(ci);
   hTopology__5->SetLineWidth(2);
   hTopology__5->GetXaxis()->SetTitle("True final state topology");
   hTopology__5->GetXaxis()->SetBinLabel(1,"CC 0#pi, 0p");
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
