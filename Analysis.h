//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 16 14:37:30 2017 by ROOT version 6.06/01
// from TTree MCTree/MCTree
// found on file: MC_DY.root
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TPaveStats.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           event;
   Float_t         weight;
   vector<float>   *reco_lepton_pt;
   vector<float>   *reco_lepton_energy;
   vector<float>   *reco_lepton_eta;
   vector<float>   *reco_lepton_phi;
   vector<int>     *reco_lepton_charge;
   vector<float>   *reco_lepton_HoverE;
   vector<float>   *reco_lepton_dxy;
   vector<float>   *reco_lepton_dR03TrackSumPt;
   vector<bool>    *reco_lepton_iselectron;
   vector<bool>    *reco_lepton_isGood;
   vector<int>     *gen_lepton_charge;
   vector<float>   *gen_lepton_pt;
   vector<float>   *gen_lepton_eta;
   vector<float>   *gen_lepton_phi;
   vector<float>   *gen_lepton_energy;
   vector<bool>    *gen_lepton_iselectron;
   vector<float>   *gen_quark_pz;
   vector<int>   *gen_quark_pdgId;
   Bool_t          isDYqqOnly;
   Int_t nbOfAdditionalTracksInVertex;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_reco_lepton_pt;   //!
   TBranch        *b_reco_lepton_energy;   //!
   TBranch        *b_reco_lepton_eta;   //!
   TBranch        *b_reco_lepton_phi;   //!
   TBranch        *b_reco_lepton_charge;   //!
   TBranch        *b_reco_lepton_HoverE;   //!
   TBranch        *b_reco_lepton_dxy;   //!
   TBranch        *b_reco_lepton_dR03TrackSumPt;   //!
   TBranch        *b_reco_lepton_iselectron;   //!
   TBranch        *b_reco_lepton_isGood;   //!
   TBranch        *b_gen_lepton_charge;   //!
   TBranch        *b_gen_lepton_pt;   //!
   TBranch        *b_gen_lepton_eta;   //!
   TBranch        *b_gen_lepton_phi;   //!
   TBranch        *b_gen_lepton_energy;   //!
   TBranch        *b_gen_lepton_iselectron;   //!
   TBranch        *b_gen_quark_pz;   //!
   TBranch        *b_gen_quark_pdgId;   //!
   TBranch        *b_isDYqqOnly;   //!
   TBranch        *b_nbOfAdditionalTracksInVertex;
   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     DrawHisto1D(TH1F *Histo, TString name);
   virtual void     Superimpose2Histos(TH1F *Histo1, TH1F *Histo2, string name, string leg1, string leg2);
   virtual Float_t  DeltaR(float eta1, float phi1, float eta2, float phi2);
   virtual TH1F * MakeMyHisto(TString name, TString title, int nbins, double first, double last);
   virtual Int_t QuarkPzSign();
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("tree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("tree","");
      chain->Add("DataEE_FullRun2.root");
      //chain->Add("MC_DY_ee.root");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   reco_lepton_pt = 0;
   reco_lepton_energy = 0;
   reco_lepton_eta = 0;
   reco_lepton_phi = 0;
   reco_lepton_charge = 0;
   reco_lepton_HoverE = 0;
   reco_lepton_dxy = 0;
   reco_lepton_dR03TrackSumPt = 0;
   gen_lepton_charge = 0;
   gen_lepton_pt = 0;
   gen_lepton_eta = 0;
   gen_lepton_phi = 0;
   gen_lepton_energy = 0;
   gen_quark_pz = 0;
   gen_quark_pdgId = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("reco_lepton_pt", &reco_lepton_pt, &b_reco_lepton_pt);
   fChain->SetBranchAddress("reco_lepton_energy", &reco_lepton_energy, &b_reco_lepton_energy);
   fChain->SetBranchAddress("reco_lepton_eta", &reco_lepton_eta, &b_reco_lepton_eta);
   fChain->SetBranchAddress("reco_lepton_phi", &reco_lepton_phi, &b_reco_lepton_phi);
   fChain->SetBranchAddress("reco_lepton_charge", &reco_lepton_charge, &b_reco_lepton_charge);
   fChain->SetBranchAddress("reco_lepton_HoverE", &reco_lepton_HoverE, &b_reco_lepton_HoverE);
   fChain->SetBranchAddress("reco_lepton_dxy", &reco_lepton_dxy, &b_reco_lepton_dxy);
   fChain->SetBranchAddress("reco_lepton_dR03TrackSumPt", &reco_lepton_dR03TrackSumPt, &b_reco_lepton_dR03TrackSumPt);
   fChain->SetBranchAddress("reco_lepton_iselectron", &reco_lepton_iselectron, &b_reco_lepton_iselectron);
   fChain->SetBranchAddress("reco_lepton_isGood", &reco_lepton_isGood, &b_reco_lepton_isGood);
   fChain->SetBranchAddress("gen_lepton_charge", &gen_lepton_charge, &b_gen_lepton_charge);
   fChain->SetBranchAddress("gen_lepton_pt", &gen_lepton_pt, &b_gen_lepton_pt);
   fChain->SetBranchAddress("gen_lepton_eta", &gen_lepton_eta, &b_gen_lepton_eta);
   fChain->SetBranchAddress("gen_lepton_phi", &gen_lepton_phi, &b_gen_lepton_phi);
   fChain->SetBranchAddress("gen_lepton_energy", &gen_lepton_energy, &b_gen_lepton_energy);
   fChain->SetBranchAddress("gen_lepton_iselectron", &gen_lepton_iselectron, &b_gen_lepton_iselectron);
   fChain->SetBranchAddress("gen_quark_pz", &gen_quark_pz, &b_gen_quark_pz);
   fChain->SetBranchAddress("gen_quark_pdgId", &gen_quark_pdgId, &b_gen_quark_pdgId);
   fChain->SetBranchAddress("isDYqqOnly", &isDYqqOnly, &b_isDYqqOnly);
   fChain->SetBranchAddress("nbOfAdditionalTracksInVertex", &nbOfAdditionalTracksInVertex, &b_nbOfAdditionalTracksInVertex);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
