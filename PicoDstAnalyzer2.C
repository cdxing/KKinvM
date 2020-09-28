/**
 * \brief Example of how to read a file (list of files) using StPicoEvent classes
 *
 * RunPicoDstAnalyzer.C is an example of reading STAR picoDst format.
 * One can use either picoDst file or a list of picoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 *
 * \brief Incoporated Jason's FXT 4.5 Phi analysis code
 *
 * \author Yang Wu
 * \date May 15, 2019
 *
 * \Incortate Yang's  event plane information to reconstruct Phi flow
 *
 * \author Ding Chen
 * \date August 3, 2019
 *
 * \Update to analyze Run18 3p85 GeV 26.5GeV FXT analysis
 *
 * \author Ding Chen
 * \date Feb 22, 2020
 * \Update to analyze Run18 3p85 GeV 26.5GeV FXT analysis
 *
 * \author Ding Chen
 * \date Jul 8, 2020
 * \Update to do Min-bias study
 */

// C++ headers
#include <iostream>
#include <cstdio>
#include <vector>
#include <set>

// ROOT headers
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TSystem.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"

// PicoDst headers
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "StPicoEvent/StPicoDstReader.h"
#include "StPicoEvent/StPicoHelix.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofHit.h"
#include "StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoTrackCovMatrix.h"

// ----------- Class st_k to store useful Kaon track information ---------------
class st_K : public TObject
{
public:
  st_K(){ reset(); }
  virtual ~st_K(){}

  bool b_pos_charge;
  bool b_bad_TOF;

  unsigned int runNumber;
  unsigned int eventNumber;
  unsigned int nGoodTracks;
  unsigned int nFXTMult;

  float px;
  float py;
  float pz;

  float x_vtx;
  float y_vtx;
  float z_vtx;
  float r_vtx;
  float DCA_r;

  float tofBeta;
  float mass2;

  float nHitsDedx;
  float nHitsFit;
  float nHitsRatio;

  float d_TPCnSigmaKaon;
  float d_TOFnSigmaKaon;

  void reset()
  {
    b_pos_charge = 0;
    b_bad_TOF = 0;

    runNumber   = 0;
    eventNumber = 0;
    nGoodTracks = 0;
    nFXTMult  = 0;

    px = 0.0;
    py = 0.0;
    pz = 0.0;

    x_vtx = -9999.0;
    y_vtx = -9999.0;
    z_vtx = -9999.0;
    r_vtx = -9999.0;
    DCA_r = -9999.0;

    tofBeta = -999.0;
    mass2 = -9999.0;

    nHitsDedx = -9999.0;
    nHitsFit = -9999.0;
    nHitsRatio = -9999.0;

    d_TPCnSigmaKaon = -9999.0;
    d_TOFnSigmaKaon = -9999.0;
  }

private:
  ClassDef(st_K,1);
};

const Double_t _massPion     = 0.13957039;
const Double_t _massKaon     = 0.493677;
const Double_t _massProton   = 0.938272081;
// ===================== (1) Main Function to build Kaon Tree ==================
void PicoDstAnalyzer2(
  const Char_t *inFile = "../files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
  TString outFile = "test"
){// -------------------- Read information from PicoDst ------------------------
  std::cout << "Hi! Lets do some physics, my Lord!" << std::endl;
  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);
  std::cout << "Status has been set" << std::endl;
  std::cout << "Now I know what to read, my Lord!" << std::endl;
  if( !picoReader->chain() ){
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Number of events to read: " << events2read << std::endl;
  outFile.Append(".picoDst.result.root");
  TFile *outputFile = new TFile(outFile,"recreate");
  // ------------------- variables and histograms ------------------------------
  Double_t  d_zvtx  = -9999.0;
  Double_t  d_xvtx  = -9999.0;
  Double_t  d_yvtx  = -9999.0;

  unsigned int  i_event     = 0;
  unsigned int  runNumber   = 0;
  unsigned int  eventNumber = 0;
  std::vector <unsigned int> triggerIDs;
  st_K Kaoninfo;
  TTree *t_K = new TTree("t_K","t_K");
  t_K->Branch("Kaoninfo",&Kaoninfo);
  TH2D *hist_pt_y_kaonPlus = new TH2D("hist_pt_y_kaonPlus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_y_kaonMinus = new TH2D("hist_pt_y_kaonMinus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH1D *h_evt       = new TH1D("h_evt","# of events",1,0,1);
  TH2D *h2_dEdxVsPq = new TH2D("h2_dEdxVsPq","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *h2_dEdxVspTq = new TH2D("h2_dEdxVspTq","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH1D *hist_Vz_pri = new TH1D("hist_Vz_pri","V_{Z} [cm]",6000,-300.0,300.0);
  TH2D *hist_VyVx_pri = new TH2D("hist_VyVx_pri","V_{Y} [cm] vs. V_{X} [cm]",500,-5.0,5.0,500,-5.0,5.0);
  TH1D *hist_Vr_pri = new TH1D("hist_Vr_pri","V_{R} [cm]",500,0.0,20.0);
  TH1D *hist_pT       = new TH1D("hist_pT","Transverse momentum pT distribution",64,0.0,32.0);
  TH1D *hist_mom = new TH1D("hist_mom","p_{mom} [GeV/c]",2000,0.0,10.0);
  TH1D *hist_DCA = new TH1D("hist_DCA","hist_DCA",100,0,10.0);
  TH2D *h2_dEdx_All_pq = new TH2D("h2_dEdx_All_pq","dEdx [GeV/cm] vs. q*p [C*GeV/c]",500,-2.0,2.0,2000,0.,10.);
  TH2D *h2_m2_QA_pq = new TH2D("h2_m2_QA_pq","m^{2} vs. q*p [C*GeV/c]",5000,-5.0,5.0,5000,-0.2,1.6);
  TH2D *h2_m2_QA_pT = new TH2D("h2_m2_QA_pT","m^{2} vs. q*pT [C*GeV/c]",5000,-5.0,5.0,5000,-2.0,1.6);
  TH1D *h_mult     = new TH1D("h_mult","# of Good Tracks",1600,0.0,1600.0);
  TH1D *h_K_DCA_r      = new TH1D("h_K_DCA_r","helix DCA",400,0.0,6.0);
  TH1D *h_K_obj_DCA_r  = new TH1D("h_K_obj_DCA_r","gDCA",400,0.0,6.0);
  TH1D *h_K_diff_DCA_r = new TH1D("h_K_diff_DCA_r","Difference between helix DCA and gDCA",800,-6.0,6.0);
  TH1D *hist_SE_mass_Phi     = new TH1D("hist_SE_mass_Phi","Same event invariant mass",200,0.9,1.1);
  // ===================== (2) Analysis on each and every events  ==============
  for(Long64_t iEvent = 0; iEvent < events2read; iEvent++){
    if((iEvent+1)%100 == 0) std::cout << "Working on event #[" << (iEvent+1) << "/" << events2read << "]" << std::endl;
    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ){
      std::cout << "Something went wrong, my Lord! Nothing to analyze..." << std::endl;
      break;
    }
    StPicoDst   *dst   = picoReader->picoDst();
    StPicoEvent *event = dst->event();
    if( !event ){
      std::cout << "Something went wrong, my Lord! Event is hiding from me..." << std::endl;
      break;
    }
    runNumber     = event->runId();
    eventNumber   = event->eventId();
    Int_t nTracks = dst->numberOfTracks();
    const Float_t _f_MagField = event->bField();
    // --------------------------- Trigger Selection ---------------------------
    triggerIDs.clear();
    triggerIDs      = event->triggerIds();
    bool b_bad_trig = true;
    for(unsigned int i=0; i < triggerIDs.size(); i++){
      if(triggerIDs[i] == 630052) b_bad_trig = false; // bbce_tofmult1 production_26p5GeV_fixedTarget_2018
    }
    TVector3 pVtx               = event->primaryVertex();
    Double_t primaryVertex_perp = (Double_t)event->primaryVertex().Perp();
    d_zvtx = pVtx.z();
    d_xvtx = pVtx.x();
    d_yvtx = pVtx.y();
    hist_Vz_pri->Fill(d_zvtx);
    hist_VyVx_pri->Fill(d_xvtx,d_yvtx);
    hist_Vr_pri->Fill(primaryVertex_perp);
    // hist_DCA  ->Fill(picoTrack->gDCA(d_xvtx,d_yvtx,d_zvtx));
    bool b_bad_zvtx   =  ((d_zvtx < 197.6) || (d_zvtx > 202.4)); //FXT_26p5_2018 // loose systematic cut
    bool b_bad_xvtx   =  ((d_xvtx < -1.0) || (d_xvtx > 1.0)); //FXT_26p5_2018
    bool b_bad_yvtx   =  ((d_yvtx < -3.0) || (d_yvtx > -0.5)); //FXT_26p5_2018
    bool b_bad_rvtx   =  sqrt(pow(d_xvtx,2)+pow(d_yvtx+2,2))> 2.4;  // loose systematic cut
    bool b_bad_evt  = b_bad_zvtx || b_bad_trig /*|| b_bad_xvtx || b_bad_yvtx */|| b_bad_rvtx;
    if(b_bad_evt) continue;
    // ====== (3) Track Loop to determine nGoodTracks (centrality)  ============
    int nGoodTracks = 0;
    int nFXTMult = 0;
    std::vector<StPicoTrack *> vGoodTracks; // vector of good tracks for TPC event plane Q-vector loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++){
      StPicoTrack *picoTrack = dst->track(iTrk);
      if(!picoTrack) continue;
      if(!picoTrack->isPrimary()) continue;
      nFXTMult++;
      StPicoBTofPidTraits *trait = NULL;
      double        tofBeta    = -999.;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      if(trait){
        tofBeta   = trait->btofBeta();
        // b_bad_ToF = (trait->btof() <= 0.0);
      }
      double d_px  = picoTrack->pMom().x();
      double d_py  = picoTrack->pMom().y();
      double d_pz  = picoTrack->pMom().z();
      double d_pT  = picoTrack->pPt();
      double d_mom = sqrt(d_pT*d_pT + d_pz*d_pz);
      hist_pT   ->Fill(d_pT);
      hist_mom  ->Fill(d_mom);
      double mass2 = d_mom*d_mom*((1.0/(tofBeta*tofBeta))-1.0);
      bool    b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool    b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.45); // loose cut
      bool b_not_enough_hits = ((double)picoTrack->nHitsFit()) < 10; // loose cut 15;
      bool    b_bad_DCA      = false;// loose cut //(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);
      bool    b_bad_track    = b_bad_dEdx || b_bad_tracking || b_not_enough_hits || b_bad_DCA;
      if(b_bad_track) continue;
      h2_dEdx_All_pq->Fill(d_mom/(picoTrack->charge()),picoTrack->dEdx());
      nGoodTracks++; // nGoodTracks is used to determine centrality later in the event loop
      vGoodTracks.push_back(picoTrack);
      if(tofBeta == -999) continue;
      h2_m2_QA_pq   ->Fill(d_mom/(picoTrack->charge()),mass2);
      h2_m2_QA_pT   ->Fill(d_pT/(picoTrack->charge()),mass2);
    }
    h_mult->Fill(nFXTMult); // use the # of primary tracks to define centrality
    // ====== (4) Track Loop to Get Kaon tracks and build Kaon TTree ===========
    vector<StPicoTrack *> v_pri_tracks;
    vector<StPicoTrack *> v_pri_tracks_pl;
    vector<StPicoTrack *> v_pri_tracks_mi;
    for(unsigned int iTrk=0; iTrk<vGoodTracks.size();iTrk++){
      StPicoTrack* picoTrack = vGoodTracks[iTrk];
      StPicoBTofPidTraits *trait        = NULL;
      Short_t charge;
      Double_t d_pt,d_px,d_py,d_pz,d_pT,ptot;
      Double_t mass2 =-999.0,tofBeta =-999.0;
      Double_t energyProton,energyKaon,energyPion,rapProton,rapKaon,rapPion,mtProton,mtKaon,mtPion;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      bool b_bad_ToF       = false;
      if(trait){
        tofBeta   = trait->btofBeta();
        b_bad_ToF = (trait->btof() <= 0.0);
      }
      charge = picoTrack->charge();
      d_px     = picoTrack->pMom().x();
      d_py     = picoTrack->pMom().y();
      d_pz     = picoTrack->pMom().z();
      d_pT     = picoTrack->pPt();
      ptot   = sqrt(d_pT*d_pT + d_pz*d_pz);

      if(tofBeta != -999.0) mass2 = ptot * ptot *( ( 1.0 / ( tofBeta*tofBeta ) ) - 1.0 );
      Int_t particleType = -999; //default -999. 0,1,2,3,4 indicate p, K+, K-, \Pi+, \Pi-
      energyProton = TMath::Sqrt(ptot*ptot + _massProton*_massProton);
      energyKaon = TMath::Sqrt(ptot*ptot + _massKaon*_massKaon);
      energyPion = TMath::Sqrt(ptot*ptot + _massPion*_massPion);
      rapProton    = 0.5*TMath::Log( (energyProton + d_pz) / (energyProton - d_pz) );
      rapKaon    = 0.5*TMath::Log( (energyKaon + d_pz) / (energyKaon - d_pz) );
      rapPion    = 0.5*TMath::Log( (energyPion + d_pz) / (energyPion - d_pz) );
      mtProton   = TMath::Sqrt(d_pT*d_pT + _massProton*_massProton);
      mtKaon   = TMath::Sqrt(d_pT*d_pT + _massProton*_massProton);
      mtPion   = TMath::Sqrt(d_pT*d_pT + _massProton*_massProton);
      h2_dEdxVsPq->Fill(charge*ptot,picoTrack->dEdx());
      h2_dEdxVspTq->Fill(charge*d_pT,picoTrack->dEdx());
      double d_TPCnSigmaElectron = fabs(picoTrack->nSigmaElectron());
      double d_TPCnSigmaPion   = fabs(picoTrack->nSigmaPion());
      double d_TPCnSigmaProton = fabs(picoTrack->nSigmaProton());
      double d_TPCnSigmaKaon   = fabs(picoTrack->nSigmaKaon());
      double d_dEdx = picoTrack->nHitsDedx();
      Double_t d_nSigmaKaonCut, d_KaonM2low, d_KaonM2high, d_KaonpTlow;
      // default -> loose cuts
      d_nSigmaKaonCut = 4.0; //2.0 // default cuts
      d_KaonM2low     = 0.15; //0.16
      d_KaonM2high    = 0.33; //0.33
      d_KaonpTlow     = 0.0; //0.2
      // ------------------------ Particle identifications ------------------------------
      // Reference: PicoAnalyzer.cxx in Repository EpdAna, brach EpdEpPhiFlow_v2_SysErr
      if( // Proton PID: require both TPC and TOF
        TMath::Abs(picoTrack->nSigmaProton()) < 2.0 &&
        tofBeta != -999.0 && mass2 > 0.7 && mass2 < 1.1 &&
        d_pT > 0.2 &&
        charge > 0
      ){
        particleType=0;// Proton
        // nProtons++;
        // // Fill histograms
        // hist_pt_proton->Fill(d_pT);
        // hist_eta_proton->Fill(eta);
        // hist_y_proton->Fill(rapProton);
        // hist_phi_proton->Fill(phi);
        // hist_rap_eta_proton->Fill(eta,rapProton);
        // hist_pt_y_proton->Fill(rapProton,d_pT,1);
        // hist_pt_eta_proton->Fill(eta,d_pT,1);
        // hist_dEdx_proton->Fill(charge*ptot,picoTrack->dEdx());
        // hist_beta_proton->Fill(charge*ptot,1.0/tofBeta);
        // hist_mass_proton->Fill(charge*ptot,mass2);
      } else if( // Kaons PID: require both TPC and TOF
        TMath::Abs(picoTrack->nSigmaKaon()) < d_nSigmaKaonCut &&
        tofBeta != -999.0 && mass2 > d_KaonM2low && mass2 < d_KaonM2high
        && d_pT > d_KaonpTlow
      ){
        if(charge > 0){
          particleType=1;// K+
          // nKaonPlus++;
          // v_KaonPlus_tracks.push_back(picoTrack); // push back K+ tracks
          // // Fill histograms
          // hist_pt_kaonPlus->Fill(d_pT);
          // hist_eta_kaonPlus->Fill(eta);
          // hist_y_kaonPlus->Fill(rapKaon);
          // hist_phi_kaonPlus->Fill(phi);
          // hist_rap_eta_kaonPlus->Fill(eta,rapKaon);
          // hist_pt_y_kaonPlus->Fill(rapKaon,d_pT,1);
          // hist_pt_eta_kaonPlus->Fill(eta,d_pT,1);
          // hist_dEdx_kaonPlus->Fill(charge*ptot,picoTrack->dEdx());
          // hist_beta_kaonPlus->Fill(charge*ptot,1.0/tofBeta);
          // hist_mass_kaonPlus->Fill(charge*ptot,mass2);
        } else { // charge < 0
          particleType=2;// K-
          // nKaonMinus++;
          // v_KaonMinus_tracks.push_back(picoTrack); // push back K+ tracks
          // // Fill histograms
          // hist_pt_kaonMinus->Fill(d_pT);
          // hist_eta_kaonMinus->Fill(eta);
          // hist_y_kaonMinus->Fill(rapKaon);
          // hist_phi_kaonMinus->Fill(phi);
          // hist_rap_eta_kaonMinus->Fill(eta,rapKaon);
          // hist_pt_y_kaonMinus->Fill(rapKaon,d_pT,1);
          // hist_pt_eta_kaonMinus->Fill(eta,d_pT,1);
          // hist_dEdx_kaonMinus->Fill(charge*ptot,picoTrack->dEdx());
          // hist_beta_kaonMinus->Fill(charge*ptot,1.0/tofBeta);
          // hist_mass_kaonMinus->Fill(charge*ptot,mass2);
        }
      } else if( // Pions PID: require both TPC and TOF
        TMath::Abs(picoTrack->nSigmaPion()) <  2.0 &&
        tofBeta != -999.0 && mass2 > -0.01 && mass2 < 0.05 &&
        d_pT > 0.2  &&
        !(TMath::Abs(mass2)<0.005 && ptot<0.25) // Remove electron influence
      ){
        if(charge > 0){
          particleType=3;// \Pi+
          // nPionPlus++;
          // // Fill histograms
          // hist_pt_pionPlus->Fill(d_pT);
          // hist_eta_pionPlus->Fill(eta);
          // hist_y_pionPlus->Fill(rapPion);
          // hist_phi_pionPlus->Fill(phi);
          // hist_rap_eta_pionPlus->Fill(eta,rapPion);
          // hist_pt_y_pionPlus->Fill(rapPion,d_pT,1);
          // hist_pt_eta_pionPlus->Fill(eta,d_pT,1);
          // hist_dEdx_pionPlus->Fill(charge*ptot,picoTrack->dEdx());
          // hist_beta_pionPlus->Fill(charge*ptot,1.0/tofBeta);
          // hist_mass_pionPlus->Fill(charge*ptot,mass2);
        } else { // charge < 0
          particleType=4;// \Pi-
          // nPionMinus++;
          // // Fill histograms
          // hist_pt_pionMinus->Fill(d_pT);
          // hist_eta_pionMinus->Fill(eta);
          // hist_y_pionMinus->Fill(rapPion);
          // hist_phi_pionMinus->Fill(phi);
          // hist_rap_eta_pionMinus->Fill(eta,rapPion);
          // hist_pt_y_pionMinus->Fill(rapPion,d_pT,1);
          // hist_pt_eta_pionMinus->Fill(eta,d_pT,1);
          // hist_dEdx_pionMinus->Fill(charge*ptot,picoTrack->dEdx());
          // hist_beta_pionMinus->Fill(charge*ptot,1.0/tofBeta);
          // hist_mass_pionMinus->Fill(charge*ptot,mass2);
        }
      }
      // Additional Kaon canditated that there's no TOF -> tofBeta == -999.0
      // # Systematic Analysis
      // sys_cutN == 16; // TPCpid
      if( // Kaons PID: tracks that only have TPC, no TOF
        // sys_cutN == 16 &&
        TMath::Abs(picoTrack->nSigmaKaon()) < d_nSigmaKaonCut &&
        TMath::Abs(picoTrack->nSigmaKaon()) < TMath::Abs(picoTrack->nSigmaElectron()) &&
        TMath::Abs(picoTrack->nSigmaKaon()) < TMath::Abs(picoTrack->nSigmaPion()) &&
        TMath::Abs(picoTrack->nSigmaKaon()) < TMath::Abs(picoTrack->nSigmaProton()) &&
        tofBeta == -999.0
        && d_pT > d_KaonpTlow
      ){
        if(charge > 0){
          particleType=1;// K+
          // nKaonPlus++;
          // v_KaonPlus_tracks.push_back(picoTrack); // push back K+ tracks
          // // Fill histograms
          // hist_pt_kaonPlus->Fill(d_pT);
          // hist_eta_kaonPlus->Fill(eta);
          // hist_y_kaonPlus->Fill(rapKaon);
          // hist_phi_kaonPlus->Fill(phi);
          // hist_rap_eta_kaonPlus->Fill(eta,rapKaon);
          // hist_pt_y_kaonPlus->Fill(rapKaon,d_pT,1);
          // hist_pt_eta_kaonPlus->Fill(eta,d_pT,1);
          // hist_dEdx_kaonPlus->Fill(charge*ptot,picoTrack->dEdx());
          // hist_beta_kaonPlus->Fill(charge*ptot,1.0/tofBeta);
          // hist_mass_kaonPlus->Fill(charge*ptot,mass2);
        } else { // charge < 0
          particleType=2;// K-
          // nKaonMinus++;
          // v_KaonMinus_tracks.push_back(picoTrack); // push back K+ tracks
          // // Fill histograms
          // hist_pt_kaonMinus->Fill(d_pT);
          // hist_eta_kaonMinus->Fill(eta);
          // hist_y_kaonMinus->Fill(rapKaon);
          // hist_phi_kaonMinus->Fill(phi);
          // hist_rap_eta_kaonMinus->Fill(eta,rapKaon);
          // hist_pt_y_kaonMinus->Fill(rapKaon,d_pT,1);
          // hist_pt_eta_kaonMinus->Fill(eta,d_pT,1);
          // hist_dEdx_kaonMinus->Fill(charge*ptot,picoTrack->dEdx());
          // hist_beta_kaonMinus->Fill(charge*ptot,1.0/tofBeta);
          // hist_mass_kaonMinus->Fill(charge*ptot,mass2);
        }
      }

      if(particleType !=1 && particleType !=2) continue; // not a kaon track
      Kaoninfo.reset();
      if(particleType == 1){
        hist_pt_y_kaonPlus->Fill(rapKaon,d_pT);
      } else if(particleType == 2){
        hist_pt_y_kaonMinus->Fill(rapKaon,d_pT);
      }
      if(charge < 0)      Kaoninfo.b_pos_charge = false;
      else                  Kaoninfo.b_pos_charge = true;
      Kaoninfo.b_bad_TOF = b_bad_ToF ;

      Kaoninfo.runNumber   = runNumber;
      Kaoninfo.eventNumber = eventNumber;
      Kaoninfo.nGoodTracks = nGoodTracks;
      Kaoninfo.nFXTMult    = nFXTMult;

      Kaoninfo.px = d_px;
      Kaoninfo.py = d_py;
      Kaoninfo.pz = d_pz;

      Kaoninfo.x_vtx = d_xvtx;
      Kaoninfo.y_vtx = d_yvtx;
      Kaoninfo.z_vtx = d_zvtx;
      Kaoninfo.r_vtx = sqrt(pow(d_xvtx,2)+pow(d_yvtx+2,2));
      Kaoninfo.DCA_r = picoTrack->gDCA(d_xvtx,d_yvtx,d_zvtx);

      Kaoninfo.tofBeta = tofBeta;
      Kaoninfo.mass2 = mass2;

      Kaoninfo.nHitsDedx = d_dEdx;
      Kaoninfo.nHitsFit = (double)picoTrack->nHitsFit();
      Kaoninfo.nHitsRatio = (double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss();

      Kaoninfo.d_TPCnSigmaKaon = d_TPCnSigmaKaon;
      t_K -> Fill();
    }
  }
  hist_pt_y_kaonPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_y_kaonMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_Vz_pri->GetXaxis()->SetTitle("V_{Z} [cm]");
  hist_Vz_pri->GetYaxis()->SetTitle("# of events");
  hist_VyVx_pri->GetXaxis()->SetTitle("V_{X} [cm]");
  hist_VyVx_pri->GetYaxis()->SetTitle("V_{Y} [cm]");
  hist_Vr_pri->GetXaxis()->SetTitle("V_{R} [cm]");
  hist_Vr_pri->GetYaxis()->SetTitle("# of events");
  outputFile->cd();
  outputFile->Write();
}
