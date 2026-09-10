#define signalDef_cxx
#include "signalDef.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

enum eventType {
   Signal,
   Background,
   Dirt,
   Cosmic
};

void signalDef::Loop()
{
//   In a ROOT session, you can do:
//      root> .L signalDef.C
//      root> signalDef t
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
   nEvents[0] = 146685; // number of events in all hyperon files
   nEvents[1] = nEvents[0]; // any hyperon event
   nEvents[2] = 199499 + nEvents[0]; // num of events in beam files + num of events in hyp files

   TFile *sigFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_sig.root", "RECREATE");

   if (!sigFile || sigFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   sigFile->cd();
   TTree *signalTree = fChain->CloneTree(0);
   signalTree->SetName("tree");
   signalTree->SetDirectory(sigFile);

   TFile *bkgFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_bkg.root", "RECREATE");

   if (!bkgFile || bkgFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   bkgFile->cd();
   TTree *bkgTree = fChain->CloneTree(0);
   bkgTree->SetName("tree");
   bkgTree->SetDirectory(bkgFile);

   TFile *cosmicFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_cosmic.root", "RECREATE");

   if (!cosmicFile || cosmicFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   cosmicFile->cd();
   TTree *cosmicTree = fChain->CloneTree(0);
   cosmicTree->SetName("tree");
   cosmicTree->SetDirectory(cosmicFile);

   TFile* outFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output.root", "RECREATE");
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
   float eventWeight = 1.0;

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

   float shortestDistTrueToRecoVtx = 0;
   int nSignal = 0;
   int nBkg = 0;
   int nDirt = 0;
   int nCosmic = 0;
   int nInRecoFVSig = 0;
   int nInRecoFVBkg = 0;
   int nInRecoFVCosmic = 0;
   int nGoodTopoSig = 0;
   int nGoodTopoBkg = 0;
   int nGoodTopoCosmic = 0;

   int nBeamOrigin = 0;
   int nBeamOriginBkg = 0;
   int nCosmicOrigin = 0;
   int nCosmicOriginBkg = 0;
   int nUnknownOrigin = 0;
   int nUnknownOriginBkg = 0;
   int sliceSizeMissmatch = 0;
   int nSingleIntEvents = 0;
   int nMultiIntEvents = 0;
   int nZeroIntEvents = 0;
   bool isCosmic;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      // --------------------------------
      // Define per event variables
      // --------------------------------

      // ---------------------------------------------------
      // Choose slice with highest nuScore (nuSlice) and use to assign event as cosmic or not
      // ---------------------------------------------------

      if (sliceID->size()==0){continue;}

     std::cout<<"sliceID size = "<<sliceID->size()<<std::endl;
     std::cout<<"sliceNuScore size = "<<sliceNuScore->size()<<std::endl;
     if(sliceID->size() == sliceNuScore->size()){
      float highestNuScore = -1; 
      for (int i = 0; i < sliceID->size(); ++i){
         float nuScore = sliceNuScore->at(i);

         if (nuScore > highestNuScore || i == 0){
            highestNuScore = nuScore;
            nuSliceIdx = i;
         }
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

      if(trueNuVtxX->size() == 0){nZeroIntEvents++;}
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

         // Check if origin of MCTruth (cosmic or beam)

         //std::cout<<"MCTruth at index "<<i<<" has origin "<<trueOrigin->at(i)<<std::endl;
         if (trueOrigin->at(i) == 1){nBeamOrigin++;}
         if (trueOrigin->at(i) == 2){nCosmicOrigin++;}
         if (trueOrigin->at(i) != 1 && trueOrigin->at(i) != 2){nUnknownOrigin++;}
      }

      //std::cout<<"Chosen MCTruth has index"<<chosenTruthIdx<<" and origin"<<trueOrigin->at(chosenTruthIdx)<<std::endl;

      // ------------------------------------------------------------
      // First loop over primary particles in the event to determine if the event contains a primary Sigma0, primary anti-muon, good Lambda, and good photon
      // ------------------------------------------------------------

      bool isInTrueFV = false;
      bool isInRecoFV = false;

      if (std::abs(trueNuVtxX->at(chosenTruthIdx)) < 180 && 
         std::abs(trueNuVtxY->at(chosenTruthIdx)) < 180 && 
         trueNuVtxZ->at(chosenTruthIdx) < 450 &&
         trueNuVtxZ->at(chosenTruthIdx) > 10){
         isInTrueFV = true;
      }

      int sigmaTrackID = -1;
      bool hasPrimarySigma0 = false;
      bool hasPrimaryMuPlus = false;

      for (int i = 0; i < truePDG->size(); i++){

         // skip if the particle is not a primary particle
         if(trueGeneration->at(i) != 0) continue;

         if (trueMCTruthIndex->at(i) != chosenTruthIdx){ continue;}

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

      // Define signal
      bool hasCorrectSigmaDecay = hasSigmaLambda && hasSigmaGamma;
      bool hasCorrectLambdaDecay = hasLambdaProton && hasLambdaPionMinus;
      bool isSignal = isInTrueFV &&hasPrimarySigma0 && hasPrimaryMuPlus && hasCorrectSigmaDecay && hasCorrectLambdaDecay && !isCosmic;

      outTree->Fill();

      if (isCosmic){ // no cosmics in filtered hyps, so should be fine like this
         sampleType = Cosmic;
         cosmicTree->Fill();
         nCosmic++;
      }
      else if (isSignal && jentry < nEvents[0] + 1){
         sampleType = Signal; 
         signalTree->Fill();
         nSignal++;
      }
      else if (!isInTrueFV && jentry > nEvents[0]){ // for now only flag dirt from background sample (in future, include from filt hyp too, but needs diff weight)
         sampleType = Dirt;
         nDirt++;
      }
      else if (!isSignal && isInTrueFV && jentry > nEvents[0]){
         sampleType = Background;
         bkgTree->Fill();
         nBkg++;
         if(trueOrigin->at(chosenTruthIdx) == 0){nUnknownOriginBkg++;}
         if(trueOrigin->at(chosenTruthIdx) == 1){nBeamOriginBkg++;}
         if(trueOrigin->at(chosenTruthIdx) == 2){nCosmicOriginBkg++;}
      }
      else {continue;} // non signal /dirt in filtered hyps, signal in beam background true fv (prob v v rare)


      // Reco FV and 3 Track + 1 Shower Count

      if (std::abs(recoVtxX) < 180 && std::abs(recoVtxY) < 180 && recoVtxZ < 450 && recoVtxZ > 10){
         isInRecoFV = true;
      }


      if(isInRecoFV){
         if(sampleType==Signal){nInRecoFVSig++;}
         if(sampleType==Background){nInRecoFVBkg++;}

        /* if (trackCount == 3 && showerCount == 1){
            if(sampleType==Signal && jentry < nEvents[0] + 1){nGoodTopoSig++;}
            if(sampleType==Background && jentry > nEvents[0]){nGoodTopoBkg++;}
         }*/

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

   std::cout<<"# Signal =  "<<nSignal<<std::endl;
   std::cout<<"# Background = "<<nBkg<<std::endl;
   std::cout<<"# Cosmics = "<<nCosmic<<std::endl;
   std::cout<<"# Signal after reco FV cut = "<<nInRecoFVSig<<std::endl;
   std::cout<<"# Background after reco FV cut = "<<nInRecoFVBkg<<std::endl;
   std::cout<<"# Cosmic after reco FV cut = "<<nInRecoFVCosmic<<std::endl;
   std::cout<<"# Signal after 3+1 topo cut = "<<nGoodTopoSig<<std::endl;
   std::cout<<"# Background after 3+1 topo cut = "<<nGoodTopoBkg<<std::endl;
   std::cout<<"# Cosmic after 3+1 topo cut = "<<nGoodTopoCosmic<<std::endl;
   std::cout<<"================================================================"<<std::endl;
   //std::cout<<"# of MCTruth Objs with Beam Nu Origin = "<<nBeamOrigin<<std::endl;
   //std::cout<<"# of MCTruth Objs with Cosmic Nu Origin = "<<nCosmicOrigin<<std::endl;
   //std::cout<<"# of MCTruths with Unknown/Other Origin = "<<nUnknownOrigin<<std::endl;
   std::cout<<"# of Bkg MCTruths with Beam Neutrino Origin = "<<nBeamOriginBkg<<std::endl;
   std::cout<<"# of Bkg MCTruths with Cosmic Neutrino Origin = "<<nCosmicOriginBkg<<std::endl;
   std::cout<<"# of Bkg MCTruths with Unknown Neutrino Origin = "<<nUnknownOriginBkg<<std::endl;
   std::cout<<"================================================================"<<std::endl;
   std::cout<<"number of events where sliceID size != sliceNuScore size = "<<sliceSizeMissmatch<<std::endl;
   std::cout<<"number of events with 1 MCtruth  = "<<nSingleIntEvents<<std::endl;
   std::cout<<"num of events with >1 MCTruth = "<<nMultiIntEvents<<std::endl;
   std::cout<<"num of events wth 0 MCTruth = "<<nZeroIntEvents<<std::endl;

}
