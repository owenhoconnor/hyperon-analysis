#define newSignalDef_cxx
#include "newSignalDef.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void newSignalDef::Loop()
{
//   In a ROOT session, you can do:
//      root> .L newSignalDef.C
//      root> newSignalDef t
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
   nEvents[0] = 137544; // number of events in all hyperon files
   nEvents[1] = nEvents[0]; // any hyperon event
   nEvents[2] = 209643 + nEvents[0]; // num of events in bkg files + num of events in hyp files

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

   TFile *beamSigFile = TFile::Open("/data/ooconnor/sbnd/hyperons/preselection_output/signalDef_output_beamSig.root", "RECREATE");

   if (!beamSigFile || beamSigFile->IsZombie()){
	   std::cerr<<"Could not open file!"<<std::endl;
	   return;
   }

   beamSigFile->cd();
   TTree *beamSigTree = fChain->CloneTree(0);
   beamSigTree->SetName("tree");
   beamSigTree->SetDirectory(beamSigFile);

   int sampleType;
   int chosenTruthIdx;

   signalTree->Branch("sampleType", &sampleType);
   bkgTree->Branch("sampleType", &sampleType);
   beamSigTree->Branch("sampleType", &sampleType);

   signalTree->Branch("chosenTruthIdx", &chosenTruthIdx);
   bkgTree->Branch("chosenTruthIdx", &chosenTruthIdx);
   beamSigTree->Branch("chosenTruthIdx", &chosenTruthIdx);

   float shortestDistTrueToRecoVtx = 0;
   int nSig = 0;
   int nBkg = 0;
   int nBeamSig = 0;
   int nInRecoFVSig = 0;
   int nInRecoFVBkg = 0;
   int nInRecoFVBeamSig = 0;
   int nGoodTopoSig = 0;
   int nGoodTopoBkg = 0;
   int nGoodTopoBeamSig = 0;

   int nBeamOrigin = 0;
   int nBeamOriginBkg = 0;
   int nCosmicOrigin = 0;
   int nCosmicOriginBkg = 0;
   int nUnknownOrigin = 0;
   int nUnknownOriginBkg = 0;

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


      // ------------------------------------------------------------
      // Choose the MCTruth index corresponding to the true neutrino interaction vertex 
      // closest to the reconstructed vertex
      // ------------------------------------------------------------

      TVector3 recoVtx(RecoVertexX, RecoVertexY, RecoVertexZ);

      for (int i = 0; i < trueNuVtxX->size(); i++){
         TVector3 trueVtx(trueNuVtxX->at(i), trueNuVtxY->at(i), trueNuVtxZ->at(i));
         float distTrueToRecoVtx = (trueVtx - recoVtx).Mag();
         if (i == 0 || distTrueToRecoVtx < shortestDistTrueToRecoVtx){
            shortestDistTrueToRecoVtx = distTrueToRecoVtx;
            chosenTruthIdx = i;
         }

         // Check if origin of MCTruth (cosmic or beam)

         std::cout<<"MCTruth at index "<<i<<" has origin "<<trueOrigin->at(i)<<std::endl;
         if (trueOrigin->at(i) == 1){nBeamOrigin++;}
         if (trueOrigin->at(i) == 2){nCosmicOrigin++;}
         if (trueOrigin->at(i) != 1 && trueOrigin->at(i) != 2){nUnknownOrigin++;}
      }

      std::cout<<"Chosen MCTruth has index"<<chosenTruthIdx<<" and origin"<<trueOrigin->at(chosenTruthIdx)<<std::endl;

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

      bool hasCorrectSigmaDecay = hasSigmaLambda && hasSigmaGamma;
      bool hasCorrectLambdaDecay = hasLambdaProton && hasLambdaPionMinus;
      bool isSignal = isInTrueFV &&hasPrimarySigma0 && hasPrimaryMuPlus && hasCorrectSigmaDecay && hasCorrectLambdaDecay;

      if (isSignal){
         sampleType = 0; 
      }
      else if (!isSignal && isInTrueFV){
         sampleType = 2;
         if(trueOrigin->at(chosenTruthIdx) == 0){nUnknownOriginBkg++;}
         if(trueOrigin->at(chosenTruthIdx) == 1){nBeamOriginBkg++;}
         if(trueOrigin->at(chosenTruthIdx) == 2){nCosmicOriginBkg++;}
      }

      // Reco FV and 3 Track + 1 Shower Cut

      if (std::abs(RecoVertexX) < 180 && std::abs(RecoVertexY) < 180 && RecoVertexZ < 450 && RecoVertexZ > 10){
         isInRecoFV = true;}


      if(isInRecoFV){
         if(sampleType==0){nInRecoFVSig++;}
         if(sampleType==2){nInRecoFVBkg++;}

         if (trackCount == 3 && showerCount == 1){
            if(sampleType==0 && jentry < nEvents[0] + 1){nGoodTopoSig++; signalTree->Fill();}
            if(sampleType==2 && jentry > nEvents[0]){nGoodTopoBkg++; bkgTree->Fill();}
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

   beamSigFile->cd();
   beamSigTree->Write("tree");
   beamSigFile->Close();
   delete beamSigFile;

   std::cout<<"# Signal =  "<<nSig<<std::endl;
   std::cout<<"# Background = "<<nBkg<<std::endl;
   std::cout<<"# Beam Signal = "<<nBeamSig<<std::endl;
   std::cout<<"# Signal after reco FV cut = "<<nInRecoFVSig<<std::endl;
   std::cout<<"# Background after reco FV cut = "<<nInRecoFVBkg<<std::endl;
   std::cout<<"# Beam Signal after reco FV cut = "<<nInRecoFVBeamSig<<std::endl;
   std::cout<<"# Signal after 3+1 topo cut = "<<nGoodTopoSig<<std::endl;
   std::cout<<"# Background after 3+1 topo cut = "<<nGoodTopoBkg<<std::endl;
   std::cout<<"# Beam Signal after 3+1 topo cut = "<<nGoodTopoBeamSig<<std::endl;
   std::cout<<"================================================================"<<std::endl;
   //std::cout<<"# of MCTruth Objs with Beam Nu Origin = "<<nBeamOrigin<<std::endl;
   //std::cout<<"# of MCTruth Objs with Cosmic Nu Origin = "<<nCosmicOrigin<<std::endl;
   //std::cout<<"# of MCTruths with Unknown/Other Origin = "<<nUnknownOrigin<<std::endl;
   std::cout<<"# of Bkg MCTruths with Beam Neutrino Origin = "<<nBeamOriginBkg<<std::endl;
   std::cout<<"# of Bkg MCTruths with Cosmic Neutrino Origin = "<<nCosmicOriginBkg<<std::endl;
   std::cout<<"# of Bkg MCTruths with Unknown Neutrino Origin = "<<nUnknownOriginBkg<<std::endl;
   std::cout<<"================================================================"<<std::endl;

}
