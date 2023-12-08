#define Analysis_cxx
#include "Analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include "AFBhelper.h" 
#include <iostream>

TFile * out;
void Analysis::Loop()
{
   if (fChain == 0) return;
   out = new TFile("out_data.root","recreate");
   Long64_t nentries = fChain->GetEntries();

   bool isMC = false;
   //bool isMC = true;
   
   // Exercice 02 : Etude de la simulation du processus de Drell-Yan
   //// Declare histograms
   // structure de l'argument de la fonction MakeMyHisto: ("nom du fichier", "titre de l'histogramme;titre de l'axe des abcisses; titre de l'axe des ordonnées", nombre de bins, borne min pour l'axe des abcisse, borne max pour l'axe des abcisses)
   
   TH1F *Histo_electron_pt_gen = MakeMyHisto("Histo_electron_pt_gen","GEN-level electron p_{T} ;p_{T}^{e-} [GeV];Events (unitless)",200,0,200);
   TH1F *Histo_electron_phi_gen = MakeMyHisto("Histo_electron_phi_gen","GEN-level electron  #phi ; #phi^{e-} (unitless);Events (unitless)",100,-3.14,3.14);
   TH1F *Histo_electron_eta_gen = MakeMyHisto("Histo_electron_eta_gen","GEN-level electron  #eta; #eta^{e-} (unitless);Events (unitless)",100,-6,6);
   
   // Exercice 02.1 : Cinématique des électrons/positrons au niveau généré 
   
   TH1F *Histo_positron_pt_gen = MakeMyHisto("Histo_positron_pt_gen","GEN-level positron p_{T};p_{T}^{e+} [GeV];Events (unitless)",200,0,200); // Il suffit de rajouter les mêmes types d'histogrammes pour les positrons.
   TH1F *Histo_positron_phi_gen = MakeMyHisto("Histo_positron_phi_gen","GEN-level positron #phi; #phi^{e+} [GeV];Events (unitless) ",100,-3.14,3.14);
   TH1F *Histo_positron_eta_gen = MakeMyHisto("Histo_positron_eta_gen","GEN-level positron  #eta; #eta^{e+} [GeV];Events (unitless)",100,-6,6); // J'ai augmenté le nombre de bins pour que ca soit comparable à l'électron
   
   //Exercice 02.2 : Cinématique des paires electron-positron au niveau généré 
   TH1F *Histo_invMass_gen = MakeMyHisto("Histo_InvMass_gen","GEN-level (e^{-}e^{+}) pairs invariant mass ; M(e^{-}e^{+}) [GeV]; Events (unitless)",200,40,120);
   TH1F *Histo_pz_gen = MakeMyHisto("Histo_pz_gen","GEN-level longitudinal momentum p_{z} of (e^{-}e^{+})  pairs; p_{z}(e^{-}e^{+}) [GeV]; Events (unitless)",1000,-2000,2000);
   TH1F *Histo_y_gen = MakeMyHisto("Histo_y_gen","GEN-level rapidity of (e^{-}e^{+})  pairs; y(e^{-}e^{+}) (unitless); Events (unitless)",100,-6,6);
   
  
   TH1F *Histo_electron_eta_gen_m200 = MakeMyHisto("Histo_electron_eta_gen_m200","GEN-level electron pseudorapidity for M(e^{-}e^{+})>200 [GeV]; #eta^{e} (unitless); Events (unitless)",100,-6,6);
   TH1F *Histo_positron_eta_gen_m200 = MakeMyHisto("Histo_positron_eta_gen_m200","GEN-level positron pseudorapidity for M(e^{-}e^{+})>200 [GeV]; #eta^{e-} (unitless); Events (unitless)",100,-6,6);
   
  //Exercice 02.3 : Distribution angulaire des leptons dans le centre de masse de la paire de leptons 
   TH1F *Histo_electron_costheta_gen_meas= MakeMyHisto("Histo_electron_costheta_gen_meas","GEN-level measured Theta angle in laboratory frame; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_gen_true= MakeMyHisto("Histo_electron_costheta_gen_true","GEN-level true Theta angle in laboratory frame; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_gen_meas_m200= MakeMyHisto("Histo_electron_costheta_gen_meas_m200","GEN-level measured Theta angle in the laboratory frame for M(e^{-}e^{+})>200 [GeV]; M(e^{-}e^{+}) [GeV]; Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_gen_true_m200= MakeMyHisto("Histo_electron_costheta_gen_true"," GEN-level true Theta angle in the laboratory reference frame for M(e^{-}e^{+})>200 [GeV]; M(e^{-}e^{+}) [GeV]; Events (unitless)",100,-1,1);
   
   // Exercice 02.4: Cinématique de la paire électron-positron au niveau généré après coupure d'acceptance
   TH1F *Histo_electron_costheta_gen_meas_acc= MakeMyHisto("Histo_electron_costheta_gen_meas_acc","GEN-level measured Theta angle in laboratory frame for p_{T}> 25 GeV; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_gen_true_acc= MakeMyHisto("Histo_electron_costheta_gen_true_acc","GEN-level measured Theta angle in laboratory frame for p_{T}> 25 GeV; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_gen_meas_m200_acc= MakeMyHisto("Histo_electron_costheta_gen_meas_m200_acc","GEN-level measured Theta angle in laboratory frame for p_{T}> 25 GeV and M(e^{-}e^{+})> 200 GeV; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_gen_true_m200_acc= MakeMyHisto("Histo_electron_costheta_gen_true_m200_acc","GEN-level true Theta angle in laboratory frame for p_{T}> 25 GeV and M(e^{-}e^{+})> 200 GeV; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   
   //Exercice 03.1 : Cinématique de la paire électron-positron au niveau reconstruit
   TH1F *Histo_invMass_reco = MakeMyHisto("Histo_InvMass_reco","RECO-level (e^{-}e^{+}) pairs invariant mass ; M(e^{-}e^{+}) [GeV]; Events (unitless)",200,40,120);
   TH1F *Histo_electron_phi_reco = MakeMyHisto("Histo_electron_phi_reco","RECO-level electron  #phi ; #phi^{e-} (unitless);Events (unitless)",100,-3.14,3.14);
   TH1F *Histo_electron_eta_reco = MakeMyHisto("Histo_electron_eta_reco","RECO-level electron  #eta; #eta^{e-} (unitless);Events (unitless)",100,-6,6);
   TH1F *Histo_electron_costheta_reco_meas_m50= MakeMyHisto("Histo_electron_costheta_reco_meas_m50","RECO-level measured #theta angle in the laboratory frame for M(e^{-}e^{+})>50 [GeV]; cos(#theta) (unitless); Events (unitless)",100,-1,1);
   TH1F *Histo_electron_costheta_reco_meas_m200= MakeMyHisto("Histo_electron_costheta_reco_meas_m200","RECO-level measured #theta angle in the laboratory frame for M(e^{-}e^{+})>200 [GeV]; cos(#theta) (unitless); Events (unitless)",100,-1,1);
      TH1F *Histo_electron_costheta_reco_meas_m500= MakeMyHisto("Histo_electron_costheta_reco_meas_m500","RECO-level measured #theta angle in the laboratory frame for M(e^{-}e^{+})>500 [GeV];cos(#theta) (unitless); Events (unitless)",100,-1,1); // not required but nice to see the evolution
      
    //useful for 03.1 to compare with the GEN-level
    TH1F *Histo_electron_costheta_gen_meas_m50_acc = MakeMyHisto("Histo_electron_costheta_gen_meas_m50_acc","GEN-level measured #theta angle in laboratory frame for p_{T}> 25 GeV and M(e^{-}e^{+})> 50GeV; cos(#theta) (unitless); Events (unitless)",100,-1,1);
    TH1F *Histo_electron_costheta_gen_meas_m500_acc = MakeMyHisto("Histo_electron_costheta_gen_meas_m500_acc","GEN-level measured #theta angle in laboratory frame for p_{T}> 25 GeV and M(e^{-}e^{+})> 500 GeV; cos(#theta) (unitless); Events (unitless)",100,-1,1);
    TH1F *Histo_electron_phi_gen_acc = MakeMyHisto("Histo_electron_phi_gen_acc","GEN-level electron  #phi for p_{T}> 25 GeV; #phi^{e-} (unitless);Events (unitless)",100,-3.14,3.14);
    TH1F *Histo_electron_eta_gen_acc = MakeMyHisto("Histo_positron_eta_gen_acc","GEN-level positron  #eta for p_{T}> 25 GeV ; #eta^{e+} (unitless);Events (unitless)",100,-6,6); 
    
    //Exercice 03.2 : Asymétrie avant-arrière au niveau reconstruit
    TH1F *Histo_Mass_forwardevents = MakeMyHistoVarBins("Histo_Mass_forwardevents","M_{reco}(e^{+}e^{-}) (forward events);M_{reco}(e^{+}e^{-}) (GeV);Events");
    TH1F *Histo_Mass_backwardevents = MakeMyHistoVarBins("Histo_Mass_backwardevents","M_{reco}(e^{+}e^{-}) (backward events);M_{reco}(e^{+}e^{-}) (GeV);Events");
    TH1F *Histo_Mass_forwardevents_largerap = MakeMyHistoVarBins("Histo_Mass_forwardevents_largerap","M_{reco}(e^{+}e^{-}) (forward events) for large rapididty;M_{reco}(e^{+}e^{-}) (GeV);Events");
    TH1F *Histo_Mass_backwardevents_largerap = MakeMyHistoVarBins("Histo_Mass_backwardevents_largerap","M_{reco}(e^{+}e^{-}) (backward events) for large rapidity ;M_{reco}(e^{+}e^{-}) (GeV);Events");


   
   
   // GenGen method
   int nBins_AFB = 20;
   float x_AFB[] = {0,50,60,70,80,90,100,110,120,130,140,150,160,180,200,250,300,350,400,500,800};
   TH1F *Histo_AFB_GenGen_num = new TH1F("Histo_AFB_GenGen_num","Forward-Backward asymetry (Gen-Gen method);M(e^{+}e^{-}) [GeV];A_{FB}",nBins_AFB,x_AFB);
   Histo_AFB_GenGen_num->Sumw2(); // Needed for proper propagation of uncertainties
   TH1F *Histo_AFB_GenGen_den = new TH1F("Histo_AFB_GenGen_den","Forward-Backward asymetry (Gen-Gen method);M(e^{+}e^{-}) [GeV];A_{FB}",nBins_AFB,x_AFB);
   Histo_AFB_GenGen_den->Sumw2();
   TH1F *Histo_forward_GenGen = new TH1F("Histo_forward_GenGen","Invariant e^{+}e^{-} mass distribution (forward events);M(e^{+}e^{-}) [GeV];Events",nBins_AFB,x_AFB);
   Histo_forward_GenGen->Sumw2();
   TH1F *Histo_backward_GenGen = new TH1F("Histo_backward_GenGen","Invariant e^{+}e^{-} mass distribution (backward events);M(e^{+}e^{-}) [GeV];Events",nBins_AFB,x_AFB);
   Histo_backward_GenGen->Sumw2();


   //// Declare 4-vectors
   TLorentzVector v_em_Gen, v_ep_Gen;
   TLorentzVector v_em_Reco, v_ep_Reco;


   //// Declare e+e- variables
   float InvMass_ee_Gen, pZ_ee_Gen;
   float InvMass_ee_Reco, pZ_ee_Reco;

   //// Declare theta angle
   float theta_q;
   float theta_em_Gen;
   float theta_em_Reco;

   Long64_t nbytes = 0, nb = 0;
   //// Loop over events
   cout << "Entries " <<nentries <<endl; 
  

   for (Long64_t jentry=0; jentry<nentries;jentry++){
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%100000==0)cout << "entry number: " <<jentry <<endl;
      ////////// GEN-level analysis //////////
    if(!isMC) weight = 1;
    if(isMC){
        // Loop over GEN electrons
	for(unsigned int i = 0; i< gen_lepton_pt->size(); i++){//First loop to select the electron (not the positron)
	  if((*gen_lepton_charge)[i]>0) continue;
	  if(!(*gen_lepton_iselectron)[i]) continue; //Make sure the lepton is an electron and not a muon
	  for(unsigned int j = 0; j< gen_lepton_pt->size(); j++){//Second loop to select the positron
	    if((*gen_lepton_charge)[j]<0) continue;
	    if(!(*gen_lepton_iselectron)[j]) continue; //Make sure the lepton is an electron and not a muon      
	    Histo_electron_pt_gen->Fill((*gen_lepton_pt)[i],weight);
	    Histo_electron_phi_gen->Fill((*gen_lepton_phi)[i],weight);
	    Histo_electron_eta_gen->Fill((*gen_lepton_eta)[i],weight);
	    
	   ////Exercice 02.2 : Cinématique des paires electron-positron au niveau généré 
	   //i c'est le e- et j c'est e+
	    Histo_positron_eta_gen->Fill((*gen_lepton_eta)[j],weight);
	    Histo_positron_pt_gen->Fill((*gen_lepton_pt)[j],weight);
	    Histo_positron_phi_gen->Fill((*gen_lepton_phi)[j],weight);


	    
	    float px_gen = (*gen_lepton_pt)[i]*cos((*gen_lepton_phi)[i]) + (*gen_lepton_pt)[j]*cos((*gen_lepton_phi)[j]); //impulsion selon l'axe x de la paire electron-positron au niveau généré
	    float py_gen = (*gen_lepton_pt)[i]*sin((*gen_lepton_phi)[i]) + (*gen_lepton_pt)[j]*sin((*gen_lepton_phi)[j]); //impulsion selon l'axe y de la paire electron-positron au niveau généré
	    float pz_gen = (*gen_lepton_pt)[i]*sinh((*gen_lepton_eta)[i]) + (*gen_lepton_pt)[j]*sinh((*gen_lepton_eta)[j]); //impulsion longitudinale de la paire electron-positron au niveau généré
	    float E_gen = (*gen_lepton_energy)[i] +(*gen_lepton_energy)[j]; // Energie de la paire electron-positron
	    float M_gen = sqrt( pow(E_gen,2) - pow(px_gen,2) - pow(py_gen,2) - pow(pz_gen,2)  ); // masse invariante de la paire electron-positron
	    float y_gen = 1.0/2.0*log(( E_gen +pz_gen )/( E_gen -pz_gen));
	    Histo_invMass_gen->Fill(M_gen,weight);
	    Histo_pz_gen ->Fill(pz_gen,weight);
	    Histo_y_gen ->Fill(y_gen,weight);
	    	    
	   if(M_gen>200) Histo_electron_eta_gen_m200->Fill((*gen_lepton_eta)[i],weight);
	   if(M_gen>200) Histo_positron_eta_gen_m200->Fill((*gen_lepton_eta)[j],weight);
	   
	   //Exercice 02.3 : Distribution angulaire des leptons dans le centre de masse de la paire de leptons 
	    
	    TLorentzVector v_ele, v_pos;
	    v_ele.SetPtEtaPhiE((*gen_lepton_pt)[i],(*gen_lepton_eta)[i],(*gen_lepton_phi)[i],(*gen_lepton_energy)[i]); // calcul la vitesse de l'electron 
	    v_pos.SetPtEtaPhiE((*gen_lepton_pt)[j],(*gen_lepton_eta)[j],(*gen_lepton_phi)[j],(*gen_lepton_energy)[j]);  // calcul de la vitesse du positron
	    v_ele.Boost(-(v_ele+v_pos).BoostVector()); //Boost de l'électron dans le referentiel du centre de masse
	    float theta_ele= v_ele.Theta();
	    Histo_electron_costheta_gen_meas->Fill(cos(theta_ele)*abs(pz_gen)/pz_gen,weight);
	    Histo_electron_costheta_gen_true->Fill(cos(theta_ele)*QuarkPzSign(),weight);
	    
	   if(M_gen>200) Histo_electron_costheta_gen_meas_m200->Fill(cos(theta_ele)*abs(pz_gen)/pz_gen,weight);
	   if(M_gen>200) Histo_electron_costheta_gen_true_m200->Fill(cos(theta_ele)*QuarkPzSign(),weight);
	   
	   //Exercice 02.4: Cinématique de la paire électron-positron au niveau généré après coupure d'acceptance
	   if ((*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4 )  Histo_electron_costheta_gen_true_acc->Fill(cos(theta_ele)*QuarkPzSign(),weight);
	   if ((*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4 )  Histo_electron_costheta_gen_meas_acc->Fill(cos(theta_ele)*abs(pz_gen)/pz_gen,weight);
	   
	   if (M_gen>50 and (*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4  )  Histo_electron_costheta_gen_meas_m50_acc->Fill(cos(theta_ele)*abs(pz_gen)/pz_gen,weight);  //Not requested for 02.4 but useful for 03.1 (GEN-level with acceptance conditions)
	   
	   if (M_gen>200 and (*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4 )  Histo_electron_costheta_gen_true_m200_acc->Fill(cos(theta_ele)*QuarkPzSign(),weight);
	   if (M_gen>200 and (*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4 )  Histo_electron_costheta_gen_meas_m200_acc->Fill(cos(theta_ele)*abs(pz_gen)/pz_gen,weight);  //Useful for 03.1
	   

	   if (M_gen>500 and (*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4  )  Histo_electron_costheta_gen_meas_m500_acc->Fill(cos(theta_ele)*abs(pz_gen)/pz_gen,weight); //Not requested for 02.4 but useful for 03.1 (GEN-level with acceptance conditions)
	   if ((*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4 )  Histo_electron_eta_gen_acc->Fill((*gen_lepton_eta)[j],weight); //Not requested for 02.4 but useful for 03.1 (GEN-level with acceptance conditions)
	   if ((*gen_lepton_pt)[j] > 25 and (*gen_lepton_pt)[i] > 25 and abs((*gen_lepton_eta)[i]) <2.4 and abs((*gen_lepton_eta)[j])<2.4  )  Histo_electron_phi_gen_acc->Fill((*gen_lepton_phi)[j],weight); //Not requested for 02.4 but useful for 03.1 (GEN-level with acceptance conditions)




	    
	  }
	}
      }
      // Exercice 3.1: Cinématique de la paire électron-positron au niveau reconstruit
      
      ////////// RECO-level analysis //////////
    
      
      
      for(unsigned int i = 0; i< reco_lepton_pt->size(); i++){//First loop to select the electron (not the positron)
	  if((*reco_lepton_charge)[i]>0) continue;
	  if(!(*reco_lepton_iselectron)[i]) continue; //Make sure the lepton is an electron and not a muon
	  for(unsigned int j = 0; j< reco_lepton_pt->size(); j++){ //Second loop to select the positron
	    if((*reco_lepton_charge)[j]<0) continue;
	    if(!(*reco_lepton_iselectron)[j]) continue; //Make sure the (anti-)lepton is an electron and not a muon 
	    if ((*reco_lepton_pt)[i] < 25 or (*reco_lepton_pt)[j] < 25  or abs((*reco_lepton_eta)[i]) >2.4 or abs((*reco_lepton_eta)[j])> 2.4 ) continue; //Make sure the acceptance conditions are satisfied
	    if (not (*reco_lepton_isGood)[i] or not (*reco_lepton_isGood)[j]) continue; //Make sure these are good leptons     


	    float px_reco = (*reco_lepton_pt)[i]*cos((*reco_lepton_phi)[i]) + (*reco_lepton_pt)[j]*cos((*reco_lepton_phi)[j]); //impulsion selon l'axe x de la paire electron-positron au niveau reconstruit
	    float py_reco = (*reco_lepton_pt)[i]*sin((*reco_lepton_phi)[i]) + (*reco_lepton_pt)[j]*sin((*reco_lepton_phi)[j]); //impulsion selon l'axe y de la paire electron-positron  au niveau reconstruit
	    float pz_reco = (*reco_lepton_pt)[i]*sinh((*reco_lepton_eta)[i]) + (*reco_lepton_pt)[j]*sinh((*reco_lepton_eta)[j]); //impulsion longitudinale de la paire electron-positron  au niveau reconstruit
	    float E_reco = (*reco_lepton_energy)[i] +(*reco_lepton_energy)[j]; // Energie de la paire electron-positron
	    float M_reco = sqrt( pow(E_reco,2) - pow(px_reco,2) - pow(py_reco,2) - pow(pz_reco,2)  ); // masse invariante de la paire electron-positron
	    
	    Histo_invMass_reco->Fill(M_reco,weight);
	    Histo_electron_phi_reco->Fill((*reco_lepton_phi)[i],weight); // distribution eta et phi pour l'électron au niveau reconstruit
	    Histo_electron_eta_reco->Fill((*reco_lepton_eta)[i],weight);
	    
	    //pour les histogrammes en cosinus du Theta mesuré
	    TLorentzVector v_ele_reco, v_pos_reco ; //initialisation des vecteurs 
	    v_ele_reco.SetPtEtaPhiE((*reco_lepton_pt)[i],(*reco_lepton_eta)[i],(*reco_lepton_phi)[i],(*reco_lepton_energy)[i]); //calcul la vitesse de l'electron au niveau reconstruit
	    v_pos_reco.SetPtEtaPhiE((*reco_lepton_pt)[j],(*reco_lepton_eta)[j],(*reco_lepton_phi)[j],(*reco_lepton_energy)[j]);  //calcul de la vitesse du positron au niveau reconstruit
	    v_ele_reco.Boost(-(v_ele_reco+v_pos_reco).BoostVector()); //Boost de l'électron dans le referentiel du centre de masse
	    float theta_ele_reco= v_ele_reco.Theta();
	    if(M_reco>50) Histo_electron_costheta_reco_meas_m50->Fill(cos(theta_ele_reco)*abs(pz_reco)/pz_reco,weight);
	    if(M_reco>200) Histo_electron_costheta_reco_meas_m200->Fill(cos(theta_ele_reco)*abs(pz_reco)/pz_reco,weight);
	    if(M_reco>500) Histo_electron_costheta_reco_meas_m500->Fill(cos(theta_ele_reco)*abs(pz_reco)/pz_reco,weight); //not requested but nice to see the evolution
	    
	    
	    
	    // Exercice 0.3.2: Asymétrie avant-arrière au niveau reconstruit (simulation)
	    float y_reco = 1.0/2.0*log(( E_reco + pz_reco )/( E_reco - pz_reco));
	    
	    if (cos(theta_ele_reco)*abs(pz_reco)/pz_reco > 0) Histo_Mass_forwardevents->Fill(M_reco,weight);
	    if (cos(theta_ele_reco)*abs(pz_reco)/pz_reco < 0) Histo_Mass_backwardevents->Fill(M_reco,weight);
	    
	    if (cos(theta_ele_reco)*abs(pz_reco)/pz_reco > 0 and abs(y_reco)>1) Histo_Mass_forwardevents_largerap->Fill(M_reco,weight);
	    if (cos(theta_ele_reco)*abs(pz_reco)/pz_reco < 0 and abs(y_reco)>1) Histo_Mass_backwardevents_largerap->Fill(M_reco,weight);

      
      }
      }
      }
      
      // Exercice 0.3.2: Asymétrie avant-arrière au niveau reconstruit (simulation)
      ComputeAFB(Histo_Mass_forwardevents,Histo_Mass_backwardevents,"AFB",out);
      ComputeAFB(Histo_Mass_forwardevents_largerap,Histo_Mass_backwardevents_largerap,"AFB_largerap",out);

   // Draw histograms
   // Structure des arguments de DrawHisto1D(ensemble de données de l'histogramme, "titre du fichier de l'histogramme")
   // Structure des arguments de la fonction Superimpose2Histos( premier histogramme, deuxième histogramme, "titre de l'histogramme superposé", "titre du prermier histogramme", "titre du deuxième histogramme")
   
   if(isMC){

     DrawHisto1D(Histo_electron_pt_gen,"Histo_electron_pt_gen");
     DrawHisto1D(Histo_electron_phi_gen,"Histo_electron_phi_gen");
     DrawHisto1D(Histo_electron_eta_gen,"Histo_electron_eta_gen");
     
   // Exercice 02.1 : Cinématique des électrons/positrons au niveau généré 
     DrawHisto1D(Histo_positron_pt_gen,"Histo_positron_pt_gen");
     DrawHisto1D(Histo_positron_phi_gen,"Histo_positron_phi_gen");
     DrawHisto1D(Histo_positron_eta_gen,"Histo_positron_eta_gen");
     
   //Exercice 02.2 : Cinématique des paires électron-positron au niveau généré 
     DrawHisto1D(Histo_invMass_gen,"Histo_invMass_gen");
     DrawHisto1D(Histo_pz_gen,"Histo_pz_gen");
     DrawHisto1D(Histo_y_gen,"Histo_y_gen");
     
     DrawHisto1D(Histo_electron_eta_gen_m200,"Histo_electron_eta_gen_m200");// not requested
     DrawHisto1D(Histo_positron_eta_gen_m200,"Histo_positron_eta_gen_m200");//not requested
     
     
   //Exercice 02.3 : Distribution angulaire des leptons dans le centre de masse de la paire de leptons 
     DrawHisto1D(Histo_electron_costheta_gen_meas,"Histo_electron_costheta_gen_meas");// not requested alone pas sure 
     DrawHisto1D(Histo_electron_costheta_gen_true,"Histo_electron_costheta_gen_true");// not requested alone we just need the superposition
     
   //Exercice 03.1: Cinématique de la paire électron-positron au niveau reconstruit
     DrawHisto1D(Histo_invMass_reco,"Histo_invMass_reco");
     DrawHisto1D(Histo_electron_phi_reco ,"Histo_electron_phi_reco");
     DrawHisto1D(Histo_electron_eta_reco,"Histo_electron_eta_reco");
     DrawHisto1D(Histo_electron_costheta_reco_meas_m50,"Histo_electron_costheta_reco_meas_50");
     DrawHisto1D(Histo_electron_costheta_reco_meas_m200,"Histo_electron_costheta_reco_meas_200");
     DrawHisto1D(Histo_electron_costheta_reco_meas_m500,"Histo_electron_costheta_reco_meas_500");
   
   
     
     
   // Example to superimpose pt(ele) and pt(positron): Create and fill a histo Histo_positron_pt_gen and then uncomment the following line: 
     
   //Exercice 02.1 : Cinématique des électrons/positrons au niveau généré 
     Superimpose2Histos(Histo_electron_pt_gen,Histo_positron_pt_gen,"Histo_positronvselectron_pt_gen", "Electrons", "Positrons");
     Superimpose2Histos(Histo_electron_phi_gen,Histo_positron_phi_gen,"Histo_positronvselectron_phi_gen", "Electrons", "Positrons");
     Superimpose2Histos(Histo_electron_eta_gen,Histo_positron_eta_gen,"Histo_positronvselectron_eta_gen", "Electrons", "Positrons");
     
   //Exercice 02.2 : Cinématique des paires électron-positron au niveau généré  
     Superimpose2Histos(Histo_electron_eta_gen_m200,Histo_positron_eta_gen_m200,"Histo_positronvselectron_eta_gen_m200", "Electrons", "Positrons");
     
   //Exercice 02.3 : Distribution angulaire des leptons dans le centre de masse de la paire de leptons 
     Superimpose2Histos(Histo_electron_costheta_gen_meas,Histo_electron_costheta_gen_true,"Histo_electron_costheta_gen_truevsmeas", "Measured angle", "True angle");// l'angle mesuré est celui du sytème des électrons et le theta true est la valeur dans le système du quark
     Superimpose2Histos(Histo_electron_costheta_gen_meas_m200,Histo_electron_costheta_gen_true_m200,"THisto_electron_costheta_gen_truesvsmeas_m200 ", "Measured angle", "True angle");
     
   // Exercice 02.4: Cinématique de la paire électron-positron au niveau généré après coupure d'acceptance
     Superimpose2Histos(Histo_electron_costheta_gen_meas_acc,Histo_electron_costheta_gen_true_acc,"Histo_electron_costheta_gen_truevsmeas_acc", "Measured angle", "True angle");
     Superimpose2Histos(Histo_electron_costheta_gen_meas_m200_acc,Histo_electron_costheta_gen_true_m200_acc,"Histo_electron_costheta_gen_truevsmeas_m200_acc", "Measured angle", "True angle");
     
   //Exercice 3.1: Cinématique de la paire électron-positron au niveau reconstruit
  
  Superimpose2Histos(Histo_invMass_reco ,Histo_invMass_gen,"Histo_invMass_genvsreco", "Reconstructed level", "Generated level");
  Superimpose2Histos(Histo_electron_phi_reco ,Histo_electron_phi_gen_acc,"Histo_electron_phi_genvsreco_acc", "Reconstructed level", "Generated level");
  Superimpose2Histos(Histo_electron_eta_reco ,Histo_electron_eta_gen_acc,"Histo_electron_eta_genvsreco_acc", "Reconstructed level", "Generated level");
  
  Superimpose2Histos(Histo_electron_costheta_reco_meas_m50 ,Histo_electron_costheta_gen_meas_m50_acc,"Histo_electron_costheta_genvsreco_meas_m50_acc", "Reconstructed level", "Generated level");
  Superimpose2Histos(Histo_electron_costheta_reco_meas_m200,Histo_electron_costheta_gen_meas_m200_acc,"Histo_electron_costheta_genvsreco_meas_m200_acc", "Reconstructed level", "Generated level");
  Superimpose2Histos(Histo_electron_costheta_reco_meas_m500 ,Histo_electron_costheta_gen_meas_m500_acc,"Histo_electron_costheta_genvsreco_meas_m500_acc", "Reconstructed level", "Generated level");
   } 
   }
   
Float_t Analysis::DeltaR(float eta1, float phi1, float eta2, float phi2){

  Float_t DeltaEta = fabs(eta1-eta2);
  Float_t DeltaPhi;
  if(fabs(phi1-phi2)>3.141592) DeltaPhi = (2*3.141592)-fabs(phi1-phi2);
  else DeltaPhi = fabs(phi1-phi2);

  return sqrt(pow(fabs(eta1-eta2),2)+pow(fabs(phi1-phi2),2));
}

void Analysis::DrawHisto1D(TH1F *Histo, TString name){
  TCanvas *c = new TCanvas;
  c->SetName("c_"+name);


   // Choose to show (option==1) or not (option==0) the stat box
   gStyle->SetOptStat(0);
   // Set log scale for the Y axis
   gPad->SetLogy();

   c->cd();

   Histo->GetXaxis()->SetTitleOffset(1.2);
   Histo->GetYaxis()->SetTitleOffset(1.4);

   Histo->Draw("hist E1 X0");


   TString stringname = name+".pdf";

   out->cd();
   Histo->Write();
   c->Write();
   c->SaveAs( stringname);
}

void Analysis::Superimpose2Histos(TH1F *Histo1, TH1F *Histo2, string name, string leg1, string leg2){
   TCanvas *c = new TCanvas ;
   c->SetName(name.c_str());

   
   Histo1->GetXaxis()->SetTitleOffset(1.2);
   Histo1->GetYaxis()->SetTitleOffset(1.4);
   Histo1->GetYaxis()->SetRangeUser(1,2*max(Histo1->GetMaximum(),Histo2->GetMaximum()));
   Histo1->SetTitle("");
   Histo1->Draw();
   Histo1->Draw("hist X0");


   Histo2->SetTitle("");
   Histo2->SetLineColor(kRed);
   Histo2->Draw("sames hist X0");

   const char *charleg1 = leg1.c_str();
   const char *charleg2 = leg2.c_str();

   gStyle->SetLegendBorderSize(0);
   TLegend *leg = new TLegend(0.15,0.7,0.38,0.88);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.038);
   leg->AddEntry(Histo1,charleg1,"lep");
   leg->AddEntry(Histo2,charleg2,"lep");
   leg->Draw();

   std::string stringname = name+".pdf";
   const char *charname = stringname.c_str();
   c->SaveAs(charname);

   out->cd();
   c->Write();
}

TH1F * Analysis::MakeMyHisto(TString name, TString title, int nbins, double first, double last){
  TH1F * h = new TH1F(name, title, nbins, first, last);
  h->Sumw2();
  return h;
}

Int_t Analysis::QuarkPzSign(){

  for(unsigned int i = 0; i< gen_quark_pz->size(); i++){
    if( (*gen_quark_pdgId)[i]  <0) continue;
    if( (*gen_quark_pz)[i] >0) return 1; 
    if( (*gen_quark_pz)[i] <0) return -1; 
  }
  
  for(unsigned int i = 0; i< gen_quark_pz->size(); i++){
    if( (*gen_quark_pdgId)[i]  >0) continue;
    if( (*gen_quark_pz)[i] >0) return -1;
    if( (*gen_quark_pz)[i] <0) return 1;
  }
  return 0;
}

