#define eventClassification_cxx
#include "eventClassification.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

enum eventType {
   Signal,
   Background,
   Dirt,
   Cosmic,
   Hyperon
};

void eventClassification::Loop()
{
//   In a ROOT session, you can do:
//      root> .L eventClassification.C
//      root> eventClassification t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;


   int nEvents[3] = {0};
   nEvents[0] = 155746; // number of events in all hyperon files
   nEvents[1] = nEvents[0]; // any hyperon event
   nEvents[2] = 242298 + nEvents[0]; // num of events in beam files + num of events in hyp files

   // POT Calculation

   TChain *srChainHyperons = new TChain("ana/subRunTree");
   TChain *srChainProd = new TChain("ana/subRunTree");
   srChainHyperons->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_hyperons.root");
   srChainProd->Add("/data/ooconnor/sbnd/hyperons/analyzer_output/merged_anaOut_2026.root");

   double potHyperons = 0;
   double potProd = 0;
   double totPotHyperons = 0;
   double totPotProd = 0;

   srChainHyperons->SetBranchAddress("pot", &potHyperons);
   srChainProd->SetBranchAddress("pot", &potProd);

   const Long64_t nHypSubRuns = srChainHyperons->GetEntries();
   const Long64_t nProdSubRuns = srChainProd->GetEntries();

   for (Long64_t iEntry = 0; iEntry < nHypSubRuns; ++iEntry){
      srChainHyperons->GetEntry(iEntry);
      totPotHyperons += potHyperons;
   }
   for (Long64_t iEntry = 0; iEntry < nProdSubRuns; ++iEntry){
      srChainProd->GetEntry(iEntry);
      totPotProd += potProd;
   }

   double potSBND = 1.0e21;
   double hyperonScaleFactor = potSBND/totPotHyperons;
   double prodScaleFactor = potSBND/totPotProd;


   std::cout<<"Total POT in Hyperon files = "<<totPotHyperons<<std::endl;
   std::cout<<"Total POT in Production files = "<<totPotProd<<std::endl;
   std::cout<<"Scale factor for Hyperon files = "<<hyperonScaleFactor<<std::endl;
   std::cout<<"Scale factor for Production files = "<<prodScaleFactor<<std::endl;

   // Output file creation

   TFile *sigFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output_sig.root", "RECREATE");

   if (!sigFile || sigFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   sigFile->cd();
   TTree *signalTree = fChain->CloneTree(0);
   signalTree->SetName("tree");
   signalTree->SetDirectory(sigFile);

   TFile *bkgFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output_bkg.root", "RECREATE");

   if (!bkgFile || bkgFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   bkgFile->cd();
   TTree *bkgTree = fChain->CloneTree(0);
   bkgTree->SetName("tree");
   bkgTree->SetDirectory(bkgFile);

   TFile *cosmicFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output_cosmic.root", "RECREATE");

   if (!cosmicFile || cosmicFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   cosmicFile->cd();
   TTree *cosmicTree = fChain->CloneTree(0);
   cosmicTree->SetName("tree");
   cosmicTree->SetDirectory(cosmicFile);

   TFile* outFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/eventClassification_output.root", "RECREATE");
   if (!outFile || outFile->IsZombie()){
      std::cerr<<"Could not open file!"<<std::endl;
      return;
   }

   outFile->cd();
   TTree *outTree = fChain->CloneTree(0);
   outTree->SetName("tree");
   outTree->SetDirectory(outFile);

   int sampleType = -1;
   int chosenTruthIdx;
   int nuSliceIdx;
   int nuSliceKey;
   int nTracks = 0;
   int nShowers = 0;
   float eventWeight = 1.0;

   signalTree->Branch("eventWeight", &eventWeight);
   bkgTree->Branch("eventWeight", &eventWeight);
   cosmicTree->Branch("eventWeight", &eventWeight);
   outTree->Branch("eventWeight", &eventWeight);

   signalTree->Branch("sampleType", &sampleType);
   bkgTree->Branch("sampleType", &sampleType);
   cosmicTree->Branch("sampleType", &sampleType);
   outTree->Branch("sampleType", &sampleType);

   signalTree->Branch("chosenTruthIdx", &chosenTruthIdx);
   bkgTree->Branch("chosenTruthIdx", &chosenTruthIdx);
   cosmicTree->Branch("chosenTruthIdx", &chosenTruthIdx);
   outTree->Branch("chosenTruthIdx", &chosenTruthIdx);

   signalTree->Branch("nuSliceIdx", &nuSliceIdx);
   bkgTree->Branch("nuSliceIdx", &nuSliceIdx);
   cosmicTree->Branch("nuSliceIdx", &nuSliceIdx);
   outTree->Branch("nuSliceIdx", &nuSliceIdx);

   signalTree->Branch("nTracks", &nTracks);
   bkgTree->Branch("nTracks", &nTracks);
   cosmicTree->Branch("nTracks", &nTracks);
   outTree->Branch("nTracks", &nTracks);

   signalTree->Branch("nShowers", &nShowers);
   bkgTree->Branch("nShowers", &nShowers);
   cosmicTree->Branch("nShowers", &nShowers);
   outTree->Branch("nShowers", &nShowers);

   float shortestDistTrueToRecoVtx = 0;
   int nSignalHyp = 0;
   int nSignalProd = 0;
   float nSignal = 0.0;
   int nHyperonRaw = 0;
   float nHyperon = 0.0;
   int nBkgRaw = 0;
   float nBkg = 0.0;
   int nDirtHyp = 0; // separate dirt samples as they are from different POT and will need to be scaled differently
   int nDirtProd = 0;
   float nDirt = 0.0;
   float nCosmic = 0.0;
   int nCosmicInTrueFV = 0;
   float nSignalInDentTrueFV = 0;
   int nSignalHypInDentTrueFV = 0;
   int nSignalProdInDentTrueFV = 0;
   float nHyperonInDentTrueFV = 0;
   int nBkgInDentTrueFVRaw = 0;
   float nBkgInDentTrueFV = 0;
   float nCosmicInDentTrueFV = 0;
   float nInRecoFVSig = 0.0;
   float nInRecoFVHyperon = 0.0;
   float nInRecoFVBkg = 0.0;
   float nInRecoFVDirt = 0.0;
   float nInRecoFVCosmic = 0.0;
   float nGoodTopoSig = 0.0;
   float nGoodTopoHyperon = 0.0;
   float nGoodTopoBkg = 0.0;
   float nGoodTopoDirt = 0.0;
   float nGoodTopoCosmic = 0.0;

   int nBeamOrigin = 0;
   int nBeamOriginBkg = 0;
   int nCosmicOrigin = 0;
   int nCosmicOriginBkg = 0;
   int nUnknownOrigin = 0;
   int nUnknownOriginBkg = 0;
   bool isCosmic;

   int nSlicesSkipped = 0;
   int sliceSizeMissmatch = 0;
   int nSingleIntEvents = 0;
   int nMultiIntEvents = 0;
   int nZeroIntEvents = 0;

   Long64_t nentries = fChain->GetEntries();
   std::cout<<"nentries ="<<nentries<<std::endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      const int treeNumber = fChain->GetTreeNumber();
      const bool isHyperonSample = (treeNumber == 0);
      const bool isProdSample = (treeNumber == 1);

      // --------------------------------
      // Define and reset per event variables
      // --------------------------------

      std::cout<<"----------------------------------------------------------------"<<std::endl;
      std::cout<<"Processing event "<<jentry<<" of "<<nentries<<std::endl;

      nuSliceIdx = -1;
      nuSliceKey = -1;
      chosenTruthIdx = -1;
      sampleType = -1;
      eventWeight = 1.0;

      // ---------------------------------------------------
      // Choose slice with highest nuScore (nuSlice) and use to assign event as cosmic or not
      // ---------------------------------------------------

     if (sliceID->size()==0){
      nSlicesSkipped++;
      continue;
     }

     if(sliceID->size() == sliceNuScore->size()){
      float highestNuScore = -1; 
      for (int i = 0; i < sliceID->size(); ++i){
         float nuScore = sliceNuScore->at(i);

         if (nuScore > highestNuScore || i == 0){
            highestNuScore = nuScore;
            nuSliceIdx = i;
            nuSliceKey = sliceKey->at(i);
         }
      }
   }

   if(sliceID->size() != sliceNuScore->size()){
      sliceSizeMissmatch++;
      continue;
   }

     // Assignment of nuSlice PFPs to track or shower based on trackScore -> question: how often does a pfp with <0.5 track score have a shower?

     nTracks = 0;
     nShowers = 0;
     for(int i = 0; i < pfpSelfID->size(); ++i){
      if(pfpSliceKey->at(i) != nuSliceKey){continue;}
       
         if(pfpTrackScore->at(i) > 0.5){
            nTracks++;
         }
         else if (pfpTrackScore->at(i) <= 0.5 && pfpTrackScore->at(i) >= 0.0){
            nShowers++;
         }
      }

   isCosmic = false;
   if (sliceTrueOrigin->at(nuSliceIdx) == 2){
      isCosmic = true;
   }

      // ------------------------------------------------------------
      // Choose the MCTruth index corresponding to the true neutrino interaction vertex 
      // closest to the reconstructed vertex
      // ------------------------------------------------------------

      if(trueNuVtxX->size() == 0){
         nZeroIntEvents++;
         continue;
      }
      if (trueNuVtxX->size() == 1){nSingleIntEvents++;}
      if (trueNuVtxX->size() > 1){nMultiIntEvents++;}

      float recoVtxX = sliceVtxX->at(nuSliceIdx);
      float recoVtxY = sliceVtxY->at(nuSliceIdx);
      float recoVtxZ = sliceVtxZ->at(nuSliceIdx);
      TVector3 recoVtx(recoVtxX, recoVtxY, recoVtxZ);

      for (int i = 0; i < trueNuVtxX->size(); i++){
         TVector3 trueVtx(trueNuVtxX->at(i), trueNuVtxY->at(i), trueNuVtxZ->at(i));
         float distTrueToRecoVtx = (trueVtx - recoVtx).Mag();
         if (i == 0 || distTrueToRecoVtx < shortestDistTrueToRecoVtx){
            shortestDistTrueToRecoVtx = distTrueToRecoVtx;
            chosenTruthIdx = i;
         }

         // Check if origin of MCTruth (cosmic or beam) NOTE: TRUE ORIGIN IS ALWAYS BEAM FOR BOTH HYPERON AND BEAM (use trueSliceOrigin instead), prob should depricate trueOrigin

         //std::cout<<"MCTruth at index "<<i<<" has origin "<<trueOrigin->at(i)<<std::endl;
         //if (trueOrigin->at(i) == 1){nBeamOrigin++;}
         //if (trueOrigin->at(i) == 2){nCosmicOrigin++;}
         //if (trueOrigin->at(i) != 1 && trueOrigin->at(i) != 2){nUnknownOrigin++;}
      }

      //std::cout<<"Chosen MCTruth has index"<<chosenTruthIdx<<" and origin"<<trueOrigin->at(chosenTruthIdx)<<std::endl;

      // ------------------------------------------------------------
      // First loop over primary particles in the event to determine if the event contains a primary Sigma0, primary anti-muon, good Lambda, and good photon
      // ------------------------------------------------------------

      bool isInDentTrueFV = false;
      bool isInTrueFV = false;
      bool isInRecoFV = false;
      bool isInDentRecoFV = false;

      if (std::abs(trueNuVtxX->at(chosenTruthIdx)) < 180 && 
         std::abs(trueNuVtxY->at(chosenTruthIdx)) < 180 && 
         trueNuVtxZ->at(chosenTruthIdx) < 450 &&
         trueNuVtxZ->at(chosenTruthIdx) > 10){
         isInTrueFV = true;
      }

      if (std::abs(trueNuVtxX->at(chosenTruthIdx)) < 180 && 
         trueNuVtxY->at(chosenTruthIdx) < 100 && trueNuVtxY->at(chosenTruthIdx) > -180 && 
         trueNuVtxZ->at(chosenTruthIdx) < 250 &&
         trueNuVtxZ->at(chosenTruthIdx) > 10){
         isInDentTrueFV = true;
      }

      int sigmaTrackID = -1;
      bool hasPrimarySigma0 = false;
      bool hasPrimaryMuPlus = false;
//std::cout<<"got here"<<std::endl;
      for (int i = 0; i < truePDG->size(); i++){

         // skip if the particle is not a primary particle

         if(trueGeneration->at(i) != 0){continue;}

         if (trueMCTruthIndex->at(i) != chosenTruthIdx){continue;}

         if(truePDG->at(i) == 3212){
            hasPrimarySigma0 = true;
            sigmaTrackID = trueTrackID->at(i);
         }

         if(truePDG->at(i) == -13){
            hasPrimaryMuPlus = true;
         }
      }

      // Now loop over second generation particles

      int lambdaTrackID = -1;
      bool hasSigmaLambda = false;  
      bool hasSigmaGamma = false;

      for (int i = 0; i < truePDG->size(); i++){
         if(trueMCTruthIndex->at(i) != chosenTruthIdx){continue;}

         if(trueGeneration->at(i) != 1){continue;}

         if(trueMotherTrackID->at(i) == sigmaTrackID){
            if(truePDG->at(i) == 3122){
               hasSigmaLambda = true;
               lambdaTrackID = trueTrackID->at(i);
            }

            if(truePDG->at(i) == 22){
               hasSigmaGamma = true;
            }
         }
      }

      // Finally follow the Lambda

      bool hasLambdaProton = false;
      bool hasLambdaPionMinus = false;

      for (int i = 0; i < truePDG->size(); i++){
         if(trueMCTruthIndex->at(i) != chosenTruthIdx){continue;}

         if(trueGeneration->at(i) != 2){continue;}

         if(trueMotherTrackID->at(i) == lambdaTrackID){
            if(truePDG->at(i) == 2212){
               hasLambdaProton = true;
            }

            if(truePDG->at(i) == -211){
               hasLambdaPionMinus = true;
            }
         }
      }

      // Signal definition
      bool hasCorrectSigmaDecay = hasSigmaLambda && hasSigmaGamma;
      bool hasCorrectLambdaDecay = hasLambdaProton && hasLambdaPionMinus;
      bool isSignal = isInTrueFV && hasPrimarySigma0 && hasPrimaryMuPlus && hasCorrectSigmaDecay && hasCorrectLambdaDecay && !isCosmic;

      if (isCosmic){ // Cosmics - no cosmics in filtered hyps, so should be fine like this
         sampleType = Cosmic;
         eventWeight = prodScaleFactor;
         nCosmic += eventWeight;
         cosmicTree->Fill();
         /*if(isInTrueFV){
            nCosmicInTrueFV++;
         }
         if(isInDentTrueFV){
            nCosmicInDentTrueFV++;
         }*/
      }
      else if (isHyperonSample && isSignal){ // Sigma0 Signal
         nSignalHyp++;
         sampleType = Signal;
         eventWeight = hyperonScaleFactor; // scale hyperon events to 1e21 POT
         nSignal += eventWeight;
         signalTree->Fill();
         if(isInDentTrueFV){
            nSignalHypInDentTrueFV++;
            nSignalInDentTrueFV += eventWeight;
         }
      }
      else if (isHyperonSample && isInTrueFV){ // Other hyperons
         nHyperonRaw++;
         sampleType = Hyperon;
         eventWeight = hyperonScaleFactor; // scale hyperon events to 1e21 POT
         nHyperon += eventWeight;
         if(isInDentTrueFV){
            nHyperonInDentTrueFV += eventWeight;
         }
      }
      else if (isHyperonSample && !isInTrueFV){ // dirt in hyperon sample
         nDirtHyp++;
         sampleType = Dirt;
         eventWeight = hyperonScaleFactor; // scale hyperon events to 1e21 POT
         nDirt += eventWeight; // add to total dirt count, scaled to 1e21 POT
      }
      else if (isProdSample && isSignal){ // Signal in beam background sample (very rare)
        sampleType = Signal;
        nSignalProd++;
        eventWeight = prodScaleFactor; // scale production events to 1e21 POT
        nSignal += eventWeight; // add to total signal count, scaled to 1e21 POT
        if(isInDentTrueFV){
            nSignalProdInDentTrueFV++;
            nSignalInDentTrueFV += eventWeight;
        }
      }
      else if (isProdSample && isInTrueFV){ // Beam background
         nBkgRaw++;
         sampleType = Background;
         eventWeight = prodScaleFactor; // scale production events to 1e21 POT
         nBkg += eventWeight; // add to total bkg count, scaled to 1e21 POT
         bkgTree->Fill();
         if(isInDentTrueFV){
            nBkgInDentTrueFVRaw++;
            nBkgInDentTrueFV += eventWeight;
         }
      }
      else if (isProdSample && !isInTrueFV){ // Dirt - for now only flag dirt from bnb2026 sample (in future, include from filt hyp too, but needs diff weight)
         sampleType = Dirt; //unsure if this is necessary correct (or best way) to flag dirt, as I believe I may be able to grab dirt property from somewhere, similar to cosmic origin
         nDirtProd++;
         eventWeight = prodScaleFactor; // scale production events to 1e21 POT
         nDirt += eventWeight; // add to total dirt count, scaled to 1e21 POT
      }
      else {continue;} 

      outTree->Fill();

      // Reco FV
      if (std::abs(recoVtxX) < 180 && std::abs(recoVtxY) < 180 && recoVtxZ < 450 && recoVtxZ > 10){
         isInRecoFV = true;
      }
      if (std::abs(recoVtxX) < 180 && recoVtxY < 100 && recoVtxY > -180 && recoVtxZ < 450 && recoVtxZ > 10){
         isInDentRecoFV = true;
      }

      // 3 Track + 1 Shower Topology 
      if(isInRecoFV){
         if(sampleType==Signal){nInRecoFVSig += eventWeight;}
         if(sampleType==Hyperon){nInRecoFVHyperon += eventWeight;}
         if(sampleType==Background){nInRecoFVBkg += eventWeight;}
         if(sampleType==Dirt){nInRecoFVDirt += eventWeight;}
         if(sampleType==Cosmic){nInRecoFVCosmic += eventWeight;}
      
         if (nTracks == 3 && nShowers == 1){
            if(sampleType==Signal){nGoodTopoSig += eventWeight;}
            if(sampleType==Hyperon){nGoodTopoHyperon += eventWeight;}
            if(sampleType==Dirt){nGoodTopoDirt += eventWeight;}
            if(sampleType==Background){nGoodTopoBkg += eventWeight;}
            if(sampleType==Cosmic){nGoodTopoCosmic += eventWeight;}
         }

      }  

   } // End of event loop


   sigFile->cd();
   signalTree->Write("tree"); // Write signal tree and close file
   sigFile->Close();
   delete sigFile;

   bkgFile->cd();
   bkgTree->Write("tree"); // write bkg tree
   bkgFile->Close();
   delete bkgFile;

   cosmicFile->cd();
   cosmicTree->Write("tree");
   cosmicFile->Close();
   delete cosmicFile;

   outFile->cd();
   outTree->Write("tree");
   outFile->Close();
   delete outFile;

   std::cout<<"===================================================================="<<std::endl;
   std::cout<<"number of events where sliceID size != sliceNuScore size = "<<sliceSizeMissmatch<<std::endl;
   std::cout<<"number of events with 1 MCtruth  = "<<nSingleIntEvents<<std::endl;
   std::cout<<"num of events with >1 MCTruth = "<<nMultiIntEvents<<std::endl;
   std::cout<<"num of events wth 0 MCTruth = "<<nZeroIntEvents<<std::endl;
   std::cout<<"------------ Summary of Event Counts ----------------"<<std::endl;
   std::cout<<"Total number of events processed = "<<nentries<<std::endl;
   std::cout<<"Total number of events skipped due to no slices = "<<nSlicesSkipped<<std::endl;
   std::cout<<"# Signal (scaled )=  "<<nSignal<<" from "<<nSignalHyp<<" signal in hyperon files + "<<nSignalProd<<"signal in prod files"<<std::endl;
   std::cout<<"# Hyperons (scaled )= "<<nHyperon<<std::endl;
   std::cout<<"# Background (scaled) = "<<nBkg<<std::endl;
   std::cout<<"# Dirt (scaled) = "<<nDirt<<" from "<<nDirtHyp<<" dirt events in hyperon files + "<<nDirtProd<<" dirt events in prod files"<<std::endl;
   std::cout<<"# Cosmics (scaled) = "<<nCosmic<<std::endl;
   std::cout<<"----------- DENT True FV --------------"<<std::endl;
   std::cout<<"# Signal in DENT true FV = "<<nSignalInDentTrueFV<<std::endl;
   std::cout<<"# Hyperons in DENT true FV = "<<nHyperonInDentTrueFV<<std::endl;
   std::cout<<"# Background in DENT true FV = "<<nBkgInDentTrueFV<<std::endl;
   //std::cout<<"# Cosmic in DENT true FV = "<<nCosmicInTrueFV<<std::endl;
   std::cout<<"----------- Reco FV --------------"<<std::endl;
   std::cout<<"# Signal after reco FV cut = "<<nInRecoFVSig<<std::endl;
   std::cout<<"# Hyperons after reco FV cut = "<<nInRecoFVHyperon<<std::endl;
   std::cout<<"# Background after reco FV cut = "<<nInRecoFVBkg<<std::endl;
   std::cout<<"# Dirt after reco FV cut = "<<nInRecoFVDirt<<std::endl;
   std::cout<<"# Cosmic after reco FV cut = "<<nInRecoFVCosmic<<std::endl;
   std::cout<<"----------- 3+1 Topology --------------"<<std::endl;
   std::cout<<"# Signal after 3+1 topo cut = "<<nGoodTopoSig<<std::endl;
   std::cout<<"# Hyperons after 3+1 topo cut = "<<nGoodTopoHyperon<<std::endl;
   std::cout<<"# Background after 3+1 topo cut = "<<nGoodTopoBkg<<std::endl;
   std::cout<<"# Dirt after 3+1 topo cut = "<<nGoodTopoDirt<<std::endl;
   std::cout<<"# Cosmic after 3+1 topo cut = "<<nGoodTopoCosmic<<std::endl;
   std::cout<<"================================================================"<<std::endl;
   std::cout<<"========================== POT Info ======================================"<<std::endl;
   std::cout<<"Total POT in Hyperon files = "<<totPotHyperons<<std::endl;
   std::cout<<"Total POT in Production files = "<<totPotProd<<std::endl;
   std::cout<<"Scale factor for Hyperon files = "<<hyperonScaleFactor<<std::endl;
   std::cout<<"Scale factor for Production files = "<<prodScaleFactor<<std::endl;


}
