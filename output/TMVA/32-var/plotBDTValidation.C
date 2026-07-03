void plotBDTValidation()
{
    TFile *f = TFile::Open("validTree_BDT.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open validTree_BDT.root" << std::endl;
        return;
    }

    TTree *t = (TTree*)f->Get("validTree_BDT");
    if (!t) {
        std::cerr << "Could not find tree validTree_BDT. File contents:" << std::endl;
        f->ls();
        return;
    }

    gStyle->SetOptStat(0);

    TH1F *hSig = new TH1F("hSig", "BDT response;BDT score;Normalised events", 50, -1.0, 1.0);
    TH1F *hBkg = new TH1F("hBkg", "BDT response;BDT score;Normalised events", 50, -1.0, 1.0);

    t->Draw("BDTScore>>hSig", "IsSignal", "goff");
    t->Draw("BDTScore>>hBkg", "!IsSignal", "goff");

    if (hSig->Integral() > 0) hSig->Scale(1.0 / hSig->Integral());
    if (hBkg->Integral() > 0) hBkg->Scale(1.0 / hBkg->Integral());

    hSig->SetLineColor(kRed+1);
    hSig->SetFillColorAlpha(kRed, 0.35);
    hSig->SetLineWidth(2);

    hBkg->SetLineColor(kBlue+1);
    hBkg->SetFillColorAlpha(kBlue, 0.35);
    hBkg->SetLineWidth(2);

    TCanvas *c = new TCanvas("c", "BDT validation", 900, 700);

    double maxY = std::max(hSig->GetMaximum(), hBkg->GetMaximum());
    hSig->SetMaximum(1.25 * maxY);

    hSig->Draw("HIST");
    hBkg->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hSig, "True signal", "f");
    leg->AddEntry(hBkg, "True background", "f");
    leg->Draw();

    c->SaveAs("BDT_validation_response.png");

    Long64_t nSig      = t->GetEntries("IsSignal");
    Long64_t nSigPass  = t->GetEntries("IsSignal && PassBDT");
    Long64_t nBkg      = t->GetEntries("!IsSignal");
    Long64_t nBkgPass  = t->GetEntries("!IsSignal && PassBDT");
    Long64_t nAll      = t->GetEntries();
    Long64_t nAllPass  = t->GetEntries("PassBDT");

    double sigEff = (nSig > 0) ? double(nSigPass) / nSig : 0.0;
    double bkgEff = (nBkg > 0) ? double(nBkgPass) / nBkg : 0.0;
    double allEff = (nAll > 0) ? double(nAllPass) / nAll : 0.0;

    std::cout << "\n===== BDT validation summary =====" << std::endl;
    std::cout << "Total events:      " << nAll << std::endl;
    std::cout << "Total pass BDT:    " << nAllPass
              << "  efficiency = " << allEff << std::endl;

    std::cout << "\nTrue signal events: " << nSig << std::endl;
    std::cout << "True signal pass:   " << nSigPass
              << "  signal efficiency = " << sigEff << std::endl;

    std::cout << "\nTrue background events: " << nBkg << std::endl;
    std::cout << "True background pass:   " << nBkgPass
              << "  background efficiency = " << bkgEff << std::endl;
    std::cout << "==================================\n" << std::endl;
}
