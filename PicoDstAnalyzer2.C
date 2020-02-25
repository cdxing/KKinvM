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

/*test
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
    unsigned int nTrkvsCuts;

    float px;
    float py;
    float pz;

    float x_vtx;
    float y_vtx;
    float z_vtx;

    float tofBeta;
    float dEdxFit;
    float d_dEdxTru;
    float DCA_r;

    float d_TPCnSigmaKaon;
    float d_TOFnSigmaKaon;

    void reset()
    {
      b_pos_charge = 0;
      b_bad_TOF = 0;

      runNumber   = 0;
      eventNumber = 0;
      nGoodTracks = 0;
      nTrkvsCuts  = 0;

      px = 0.0;
      py = 0.0;
      pz = 0.0;

      x_vtx = 0.0;
      y_vtx = 0.0;
      z_vtx = 0.0;

      d_TPCnSigmaKaon = 0.0;
      d_TOFnSigmaKaon = 0.0;

      DCA_r = -9999.0;
    }

private:
    ClassDef(st_K,1);
};
*/
//=============================== Main Function ==============================================
void PicoDstAnalyzer2(const Char_t *inFile = "../files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test")
{
  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();

  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);

  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !picoReader->chain() ) {
      std::cout << "No chain has been found." << std::endl;
  }

  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;

  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Number of events to read: " << events2read << std::endl;

  outFile.Append(".picoDst.result.root");

  TFile *outputFile = new TFile(outFile,"recreate");

  // Histogramning
  double  d_zvtx  = -9999.0;
  double  d_xvtx  = -9999.0;
  double  d_yvtx  = -9999.0;

  unsigned int  i_event = 0;

  unsigned int  runNumber        = 0;
  unsigned int  eventNumber      = 0;

  std::vector <unsigned int> triggerIDs;
  /* test
  st_K Kaoninfo;
  // st_track trackinfo;

  TTree *  t_K        = new TTree("t_K","t_K");
  t_K         -> Branch("Kaoninfo",&Kaoninfo);
  */
  double d_PI_m2   = 0.019479955;
  double d_K_m2    = 0.24371698;
  double d_PRO_m2  = 0.8803545;

  TH1D *  h_evt       = new TH1D("h_evt","h_evt",1,0,1);
  TH1D *  h_zvtx      = new TH1D("h_zvtx","h_zvtx",100,200,220);
  TH1D *  h_pT       = new TH1D("h_pT","h_pT",64,0.0,32.0);
  TH2D *  h2_dEdx_PI_pq = new TH2D("h2_dEdx_PI_pq","h2_dEdx_PI_pq",500,-2.0,2.0,500,0.6,3);
  TH2D *  h2_dEdx_PRO_pq = new TH2D("h2_dEdx_PRO_pq","h2_dEdx_PRO_pq",500,-2.0,2.0,500,0.6,3);
  TH2D *  h2_dEdx_K_pq = new TH2D("h2_dEdx_K_pq","h2_dEdx_K_pq",500,-2.0,2.0,500,0.6,3);

  TH2D *hist_dEdx = new TH2D("hist_dEdx","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  hist_dEdx->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx->GetYaxis()->SetTitle("dE/dx (keV/cm)");

  TH2D *  h2_dEdx_All_pq = new TH2D("h2_dEdx_All_pq","h2_dEdx_All_pq",500,-2.0,2.0,2000,0.,10.);
  TH2D *  h2_dEdx_All_pq_1 = new TH2D("h2_dEdx_All_pq_1","h2_dEdx_All_pq_1",500,-2.0,2.0,2000,0.,10.);
  TH2D *  h2_m2_QA_pq = new TH2D("h2_m2_QA_pq","h2_m2_QA_pq",5000,-5.0,5.0,5000,-0.2,1.6);
  TH2D *  h2_m2_QA_pq_1 = new TH2D("h2_m2_QA_pq_1","h2_m2_QA_pq_1",5000,-5.0,5.0,5000,-0.2,1.6);
  TH2D *  h2_m2_QA_pT = new TH2D("h2_m2_QA_pT","h2_m2_QA_pT",5000,-5.0,5.0,5000,-2.0,1.6);
  TH2D *  h2_m2_QA_pT_1 = new TH2D("h2_m2_QA_pT_1","h2_m2_QA_pT_1",5000,-5.0,5.0,5000,-2.0,1.6);

  TH1D *  h_mult     = new TH1D("h_mult","h_mult",1600,0.0,1600.0);
  // TH2D *  h2_mult     = new TH2D("h2_mult","h2_mult",400,0.0,400.0,100000,0.,100000.);//test

  TH1D *  h_DCA_r    = new TH1D("h_DCA_r","h_DCA_r",100,0.0,2.0);
  TH1D *  h_DCA_PRO    = new TH1D("h_DCA_PRO","h_DCA_PRO",100,0.0,2.0);
  TH1D *  h_DCA_PIM    = new TH1D("h_DCA_PIM","h_DCA_PIM",100,0.0,2.0);
  TH1D *  h_DCA_mother_r    = new TH1D("h_DCA_mother_r","h_DCA_mother_r",100,0.0,2.0);
  TH1D *  h_decay_length    = new TH1D("h_decay_length","h_decay_length",1000,0.0,20.0);

  TH1D *  h_K_DCA_r      = new TH1D("h_K_DCA_r","h_K_DCA_r",400,0.0,6.0);
  TH1D *  h_K_obj_DCA_r  = new TH1D("h_K_obj_DCA_r","h_K_obj_DCA_r",400,0.0,6.0);
  TH1D *  h_K_diff_DCA_r = new TH1D("h_K_diff_DCA_r","h_K_diff_DCA_r",800,-6.0,6.0);

  TH1D *  h_rapidity_pm = new TH1D("h_rapidity_pm","h_rapidity_pm",100,-2.0,2.0);
  TH2D *  h2_mT_rapidity_pm = new TH2D("h2_mT_rapidity_pm","h2_mT_rapidity_pm",100,0.0,1.0,20,-2.05,-0.05);

  TH1D *  h_nodecaylength_cut_inv_m_PHI    = new TH1D("h_nodecaylength_cut_inv_m_PHI","h_nodecaylength_cut_inv_m_PHI",1000,0.9,1.1);//.,2.);
  TH1D *  h_nodipangle_cut_inv_m_PHI    = new TH1D("h_nodipangle_cut_inv_m_PHI","h_nodipangle_cut_inv_m_PHI",1000,0.9,1.1);//.,2.);

  TH1D *  h_inv_m_PHI         = new TH1D("h_inv_m_PHI","h_inv_m_PHI",1000,0.9,1.1);
  TH1D *  h_prim_inv_m_PHI    = new TH1D("h_prim_inv_m_PHI","h_prim_inv_m_PHI",1000,0.9,1.1);

  TH1D *  h_inv_m_K0S         = new TH1D("h_inv_m_K0S","h_inv_m_K0S",1000,0.0,1.0);
  TH1D *  h_inv_m_RHO         = new TH1D("h_inv_m_RHO","h_inv_m_RHO",2000,0.45,1.2);
  TH1D *  h_inv_m_LAMBDA = new TH1D("h_inv_m_LAMBDA","h_inv_m_LAMBDA",1000,1.0,1.18);

  TH2D *  h2_pidprob_vs_nsigma_PI  = new TH2D("h2_pidprob_vs_nsigma_PI","h2_pidprob_vs_nsigma_PI",100,0.0,1.0,100,0.0,10.0);
  TH2D *  h2_pidprob_vs_nsigma_PRO = new TH2D("h2_pidprob_vs_nsigma_PRO","h2_pidprob_vs_nsigma_PRO",100,0.0,1.0,100,0.0,10.0);
  TH2D *  h2_pidprob_vs_nsigma_K   = new TH2D("h2_pidprob_vs_nsigma_K","h2_pidprob_vs_nsigma_K",100,0.0,1.0,100,0.0,10.0);

  TH2D *  h2_TOF_nsigma_vs_TPC_nsigma_PI  = new TH2D("h2_TOF_nsigma_vs_TPC_nsigma_PI", "h2_TOF_nsigma_vs_TPC_nsigma_PI",100,0.0,10.0,100,0.0,10.0);
  TH2D *  h2_TOF_nsigma_vs_TPC_nsigma_PRO = new TH2D("h2_TOF_nsigma_vs_TPC_nsigma_PRO","h2_TOF_nsigma_vs_TPC_nsigma_PRO",100,0.0,10.0,100,0.0,10.0);
  TH2D *  h2_TOF_nsigma_vs_TPC_nsigma_K   = new TH2D("h2_TOF_nsigma_vs_TPC_nsigma_K",  "h2_TOF_nsigma_vs_TPC_nsigma_K",100,0.0,10.0,100,0.0,10.0);

  TH2D *  h2_TOF_nsigma_K_vs_PI  = new TH2D("h2_TOF_nsigma_K_vs_PI",   "h2_TOF_nsigma_K_vs_PI",100,0.0,10.0,100,0.0,10.0);
  TH2D *  h2_TPC_nsigma_K_vs_PI  = new TH2D("h2_TPC_nsigma_K_vs_PI",   "h2_TPC_nsigma_K_vs_PI",100,0.0,10.0,100,0.0,10.0);
  TH2D *  h2_TOF_nsigma_K_vs_PRO = new TH2D("h2_TOF_nsigma_K_vs_PRO",  "h2_TOF_nsigma_K_vs_PRO",100,0.0,10.0,100,0.0,10.0);
  TH2D *  h2_TPC_nsigma_K_vs_PRO = new TH2D("h2_TPC_nsigma_K_vs_PRO",  "h2_TPC_nsigma_K_vs_PRO",100,0.0,10.0,100,0.0,10.0);

  // TH1D *  h_PHI_decay_length     = new TH1D("h_PHI_decay_length","h_PHI_decay_length",100,0,10.0);
  // TH1D *  h_dip_angle            = new TH1D("h_dip_angle","h_dip_angle",1000,-4,4);

  // Histogramming End

  //=============================== Event Loop ===================================================
  for(Long64_t iEvent = 0; iEvent < events2read; iEvent++)//events2read
  {
    if((iEvent+1)%100 == 0) std::cout << "Working on event #[" << (iEvent+1) << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);

    if( !readEvent )
    {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    StPicoDst *dst = picoReader->picoDst();
    // Retrieve picoDst

    StPicoEvent *event = dst->event();
    // Retrieve event information

    Int_t runId_Dst   = event->runId();
    Int_t eventId_Dst = event->eventId();
    // std::cout << " runId from picoDst = "   <<  runId_Dst << std::endl;
    // std::cout << " eventId from picoDst = " <<  eventId_Dst << std::endl;

    if( !event )
    {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    const Float_t B = event->bField();

    int i_muEvent_num = event -> eventId();

    runNumber        = event->runId();
    eventNumber      = event->eventId();

    double d_MagField = event->bField();

    //============================ Trigger Selection ==============================================
    triggerIDs.clear();
    triggerIDs       = event->triggerIds();
    bool b_bad_trig = true;

    // loop for the trigger ids and see if any == 1
    for(unsigned int i=0; i < triggerIDs.size(); i++)
      {
        if(triggerIDs[i] == 630802) b_bad_trig = false; // hlt_fixedTargetGood 7.2GeV
      }

    //=========================== End Trigger Slection ===========================================

    TVector3 pVtx     = event->primaryVertex();
    Double_t primaryVertex_perp = (Double_t)event->primaryVertex().Perp();

    // Primary Vertex

    //=========================== Z-VTX Selection =================================================
    TVector3 v3D_vtx  = event->primaryVertex();
    d_zvtx = pVtx.z();
    d_xvtx = pVtx.x();
    d_yvtx = pVtx.y();

    h_zvtx -> Fill(d_zvtx);
    bool b_bad_zvtx   = ((d_zvtx < 199.0) || (d_zvtx > 202.0)); //FXT_26p5_2018
    bool b_bad_xvtx   =  ((d_xvtx < -1.0) || (d_xvtx > 1.0)); //FXT_26p5_2018
    bool b_bad_yvtx   =  ((d_yvtx < -3.0) || (d_yvtx > -0.5)); //FXT_26p5_2018
    bool b_bad_rvtx   =  primaryVertex_perp > 3.0;

    //======================== END Z-VTX Selection =================================================

    bool b_bad_evt  = b_bad_zvtx || b_bad_trig || b_bad_xvtx || b_bad_yvtx || b_bad_rvtx;
    if(b_bad_evt) continue;
    // Bad Event Cut

    //======================== Track Analysis =====================================================

    Int_t nTracks = dst->numberOfTracks();
    //Number of tracks in this event

    bool b_cent_01  = false; // 0  < centrality <= 5%
    bool b_cent_02  = false; // 5  < centrality <= 10%
    bool b_cent_03  = false; // 10 < centrality <= 15%
    bool b_cent_04  = false; // 15 < centrality <= 20%
    bool b_cent_05  = false; // 20 < centrality <= 25%
    bool b_cent_06  = false; // 25 < centrality <= 30%
    bool b_cent_07  = false; // centrality > 30%

    //========================= Track loop to define good tracks ======================================================
    int nGoodTracks = 0;
    int nTrkvsCuts  = 0;

    /*test
    h2_dEdx_All_pq_1->Reset();
    h2_m2_QA_pq_1   ->Reset();
    h2_m2_QA_pT_1   ->Reset();
    */
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    {
      StPicoTrack *picoTrack = dst->track(iTrk);
      // Retrieve i-th Pico Track

      if(!picoTrack) continue;
      // PicoTrack Cut

      bool b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      bool b_bad_track    = b_bad_dEdx || b_bad_tracking;
      if(b_bad_track) continue;
      // Track-Level Cut

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut

      // To get beta from Btof to calculate mass2
      StPicoBTofPidTraits *trait        = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      double d_tofBeta0                 = -999;
      if(trait) d_tofBeta0              = trait->btofBeta();

      double d_px0      = picoTrack->gMom().x();
      double d_py0      = picoTrack->gMom().y();
      double d_pz0      = picoTrack->gMom().z();
      double d_pT0      = picoTrack->gPt();
      double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
      double mass2      = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);

      /*test
      h2_dEdx_All_pq_1->Fill(d_mom0/(picoTrack->charge()),picoTrack->dEdx());
`     */
      h2_dEdx_All_pq->Fill(d_mom0/(picoTrack->charge()),picoTrack->dEdx());

      nGoodTracks++;

      if(d_tofBeta0 == -999) continue;
      // TOF Beta Cut

      /*test
      h2_m2_QA_pq_1   ->Fill(d_mom0/(picoTrack->charge()),mass2);
      h2_m2_QA_pT_1   ->Fill(d_pT0/(picoTrack->charge()),mass2);
      */
      h2_m2_QA_pq   ->Fill(d_mom0/(picoTrack->charge()),mass2);
      h2_m2_QA_pT   ->Fill(d_pT0/(picoTrack->charge()),mass2);

      h_pT            ->Fill(d_pT0);

      nTrkvsCuts++;

    }

    h_mult->Fill(nGoodTracks);
    /*test
    h2_dEdx_All_pq->Add(h2_dEdx_All_pq_1);
    h2_m2_QA_pq   ->Add(h2_m2_QA_pq_1);
    h2_m2_QA_pT   ->Add(h2_m2_QA_pT_1);
    */
    //======================== END Track loop to define good tracks ====================================================

    //========================= Track Cut Settings ===============================================
    int i_Nhits_min = 9;  //minimum number of hits

    //Mon Jul  3 09:17:54 EDT 2017 from phi paper
    double d_pT_min = 0.1;//0.05;
    double d_pT_max = 10.0;//0.05;
    double d_mom_min = 0.1;//0.05;
    double d_mom_max = 10.0;// 1.0;//0.05;
    double d_SigmaCutLevel = 4.0;//3.0;//4.0; // PID Sigma Cut // Mon Jul  3 09:16:46 EDT 2017

    //--- Not being used Fri Jun 16 10:48:30 EDT 2017
    // double d_PRO_daughter_DCA = 0.3; // trk must be greater than this dca
    // double d_PI_daughter_DCA  = 1.0; // trk must be greater than this dca

    // Lambda: 1.0, K0s 1.2
    //  double d_cut_dca_daughters    = 2.0; //distance between daughters, must be less than this
    // double d_cut_dca_daughters_lam    = 1.0;//2.0; //distance between daughters, must be less than this
    // double d_cut_dca_daughters_k0s    = 1.2; //distance between daughters, must be less than this


    // Lambda: 1.5 Anti Lambda: 2.0 K0s: 1.0
    //double d_cut_dca_mother       = 5.0; // must be less than this
    // double d_cut_dca_mother_lam       = 1.5;//5.0; // must be less than this
    // double d_cut_dca_mother_k0s       = 1.0; // must be less than this

    // Lambda: 3.0, K0s: 2.0
    //  double d_cut_mother_decay_length = 3.0; // must be greater than this
    // double d_cut_mother_decay_length_lam = 3.0; // must be greater than this
    // double d_cut_mother_decay_length_k0s = 2.0; // must be greater than this
    // double d_cut_mother_decay_length_RHO = 0.5; // must be LESS    than this
    // double d_cut_mother_decay_length_PHI = 0.5; // must be LESS    than this
    //======================== END Track Settings ================================================

    //======================= Primary Track Loop  to build kaon TTree =================================================
    vector<StPicoTrack *> v_pri_tracks;
    vector<StPicoTrack *> v_pri_tracks_pl;
    vector<StPicoTrack *> v_pri_tracks_mi;

    int index = 0;

    double d_PI_m    = 0.13957018;
    double d_PRO_m   = 0.9382720813;
    double d_K_m     = 0.493677;

    for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    {
      StPicoTrack *picoTrack = dst->track(iTrk);
      // Retrieve i-th pico track

      if(picoTrack == NULL)       continue;

      if(!picoTrack->isPrimary()) continue;
      // Primary Track Cut

      // To get beta from Btof to calculate mass2
      StPicoBTofPidTraits *trait = NULL;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );

      //nHits minimum cut
      unsigned short nHits = picoTrack->nHits();
      if(nHits < i_Nhits_min)     continue;

      bool b_bad_dEdx      = false;
      bool b_bad_tracking  = false;
      bool b_bad_ToF       = false;

      b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);

      if(trait) b_bad_ToF       = (trait->btof() <= 0.0);
      bool b_bad_track          = b_bad_dEdx || b_bad_tracking;
      if(b_bad_track)              continue;
      // Bad Track Cut

      double d_TPCnSigmaElectron = fabs(picoTrack->nSigmaElectron());
      double d_TPCnSigmaPion   = fabs(picoTrack->nSigmaPion());
      double d_TPCnSigmaProton = fabs(picoTrack->nSigmaProton());
      double d_TPCnSigmaKaon   = fabs(picoTrack->nSigmaKaon());

      double tofBeta       = -999;
      if(trait) tofBeta    = trait->btofBeta();

      double d_tofBeta0    = -999;
      if(trait) d_tofBeta0 = trait->btofBeta();

      // test new sets of PID cuts
      bool b_PI  = false; // fabs(d_TPCnSigmaPion)   < d_SigmaCutLevel;
      bool b_PRO = false; // fabs(d_TPCnSigmaProton) < d_SigmaCutLevel;
      bool b_K   = false; //fabs(d_TPCnSigmaKaon)   < d_SigmaCutLevel;
      bool b_E   = false; //fabs(d_TPCnSigmaElectron)< d_SigmaCutLevel;

      /********* Remove this part of PID cut
      if( b_E
         && (d_TPCnSigmaElectron < d_TPCnSigmaPion)
         && (d_TPCnSigmaElectron < d_TPCnSigmaProton)
         && (d_TPCnSigmaElectron < d_TPCnSigmaKaon) )
        { b_E = true; b_PI = false; b_PRO = false; b_K = false;}

      if( b_PI
         && (d_TPCnSigmaPion < d_TPCnSigmaElectron)
         && (d_TPCnSigmaPion < d_TPCnSigmaKaon)
         && (d_TPCnSigmaPion < d_TPCnSigmaProton) )
         { b_E = false; b_PI = true; b_PRO = false; b_K = false;}

      if( b_PRO
         && (d_TPCnSigmaProton < d_TPCnSigmaElectron)
         && (d_TPCnSigmaProton < d_TPCnSigmaPion)
         && (d_TPCnSigmaProton < d_TPCnSigmaKaon) )
         { b_E = false; b_PI = false; b_PRO = true; b_K = false;}

      if( b_K
         && (d_TPCnSigmaKaon < d_TPCnSigmaElectron)
         && (d_TPCnSigmaKaon < d_TPCnSigmaProton)
         && (d_TPCnSigmaKaon < d_TPCnSigmaPion) )
         { b_E = false; b_PI = false; b_PRO = false; b_K = true;}
      */

      double d_charge     = picoTrack->charge();
      double d_px0        = picoTrack->pMom().x();
      double d_py0        = picoTrack->pMom().y();
      double d_pz0        = picoTrack->pMom().z();
      double d_pT0        = picoTrack->pPt();
      double d_mom0       = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);

      double mass2        = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);
      double d_E_PRO      = sqrt(d_PRO_m*d_PRO_m + d_mom0*d_mom0);
      double d_E_PI       = sqrt(d_PI_m*d_PI_m + d_mom0*d_mom0);
      double d_E_K        = sqrt(d_K_m*d_K_m + d_mom0*d_mom0);

      double d_y_PRO      = 0.5*TMath::Log((d_E_PRO + d_pz0)/(d_E_PRO - d_pz0));
      double d_y_PI       = 0.5*TMath::Log((d_E_PI + d_pz0)/(d_E_PI - d_pz0));
      double d_y_K        = 0.5*TMath::Log((d_E_K + d_pz0)/(d_E_K - d_pz0));

      if(d_charge > 0.0)
        {
          if( b_E
            && (fabs(d_TPCnSigmaElectron) < 3.0)
            && (mass2 < 0.005)
            && (mass2 > -0.008)
            && (d_pT0 < 0.3)){
              b_E = true; b_PI = false; b_PRO = false; b_K = false;
            // h_E_plus_pT -> Fill(d_pT0);
            // h_E_plus_y  -> Fill(d_y_K);
            // h2_E_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
            // h_E_plus_mT_Diff -> Fill((d_mT_E - d_E_m));
          }
          // if( b_K
          //   && (fabs(d_TPCnSigmaKaon) < 3.0)
          //   && (mass2 > 0.15)
          //   && (mass2 < 0.35))
          if( (fabs(d_TPCnSigmaKaon) < 2.0)
             && ( d_tofBeta0 != -999.0
                 && mass2 > 0.19
                 && mass2 < 0.3
                )
             && d_pT0 > 0.2
             && d_pT0 < 1.6
          )
          {
               b_E = false; b_PI = false; b_PRO = false; b_K = true;
            // h_K_plus_pT -> Fill(d_pT0);
            // h_K_plus_y  -> Fill(d_y_K);
            // h2_K_plus_pT_vs_y -> Fill(d_y_K,d_pT0);
            // h_K_plus_mT_Diff -> Fill((d_mT_K - d_K_m));
          }
          if( b_PI
            && (fabs(d_TPCnSigmaPion) < 3.0)
            && (mass2 > -0.03)
            && (mass2 < 0.06)
            && ((d_pT0 > 0.4) || (mass2 > 0.008))
          ){
            b_E = false; b_PI = true; b_PRO = false; b_K = false;
            // h_PI_plus_pT -> Fill(d_pT0);
            // h_PI_plus_y  -> Fill(d_y_PI);
            // h2_PI_plus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            // h_PI_plus_mT_Diff -> Fill((d_mT_PI - d_PI_m));
          }
          if( b_PRO
            && (fabs(d_TPCnSigmaProton) < 3.0)
            && (mass2 > 0.7)
            && (mass2 < 1.1))
            {
              b_E = false; b_PI = false; b_PRO = true; b_K = false;
            // h_PRO_plus_pT -> Fill(d_pT0);
            // h_PRO_plus_y  -> Fill(d_y_PRO);
            // h2_PRO_plus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
            // h_PRO_plus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));
          }
        }
      else
        {
          if( b_E
            && (fabs(d_TPCnSigmaElectron) < 3.0)
            && (mass2 < 0.005)
            && (mass2 > -0.008)
            && (d_pT0 < 0.3)){
              b_E = true; b_PI = false; b_PRO = false; b_K = false;
            // h_E_minus_pT -> Fill(d_pT0);
            // h_E_minus_y  -> Fill(d_y_K);
            // h2_E_minus_pT_vs_y -> Fill(d_y_K,d_pT0);
            // h_E_minus_mT_Diff -> Fill((d_mT_E - d_E_m));
          }
          // if( b_K
          //   && (fabs(d_TPCnSigmaKaon) < 3.0)
          //   && (mass2 > 0.15)
          //   && (mass2 < 0.35))
          if( (fabs(d_TPCnSigmaKaon) < 2.0)
             && ( d_tofBeta0 != -999.0
                 && mass2 > 0.19
                 && mass2 < 0.3
                )
             && d_pT0 > 0.2
             && d_pT0 < 1.6
          )
          {
               b_E = false; b_PI = false; b_PRO = false; b_K = true;
            // h_K_minus_pT -> Fill(d_pT0);
            // h_K_minus_y  -> Fill(d_y_K);
            // h2_K_minus_pT_vs_y -> Fill(d_y_K,d_pT0);
            // h_K_minus_mT_Diff -> Fill((d_mT_K - d_K_m));
          }
          if( b_PI
            && (fabs(d_TPCnSigmaPion) < 3.0)
            && (mass2 > -0.03)
            && (mass2 < 0.06)
            && ((d_pT0 > 0.4) || (mass2 > 0.008))
          ){
            b_E = false; b_PI = true; b_PRO = false; b_K = false;
            // h_PI_minus_pT -> Fill(d_pT0);
            // h_PI_minus_y  -> Fill(d_y_PI);
            // h2_PI_minus_pT_vs_y -> Fill(d_y_PI,d_pT0);
            // h_PI_minus_mT_Diff -> Fill((d_mT_PI - d_PI_m));
          }
          if( b_PRO
            && (fabs(d_TPCnSigmaProton) < 3.0)
            && (mass2 > 0.7)
            && (mass2 < 1.1))
            {
              b_E = false; b_PI = false; b_PRO = true; b_K = false;
            // h_PRO_minus_pT -> Fill(d_pT0);
            // h_PRO_minus_y  -> Fill(d_y_PRO);
            // h2_PRO_minus_pT_vs_y -> Fill(d_y_PRO,d_pT0);
            // h_PRO_minus_mT_Diff -> Fill((d_mT_PRO - d_PRO_m));
          }
        }

        if(d_charge > 0.0)
          {
            if(b_PI||b_PRO||b_E)    continue;
            if(!b_K)                continue;
          }
        else
          {
            if(b_E||b_PI||b_PRO)    continue;
            if(!b_K)                continue;
          }
        // Kaon Cut

        if( (d_pT0<d_pT_min)   || (d_pT0 > d_pT_max))   continue;
        if( (d_mom0<d_mom_min) || (d_mom0 > d_mom_max)) continue;
        //pT Min, Mom Min Cut

        bool b_bad_DCA       = false;
        StPicoPhysicalHelix trackhelix = picoTrack->helix(B);
        double helixpathl              = trackhelix.pathLength(v3D_vtx, false);
        TVector3 v3D_dca               = trackhelix.at(helixpathl)-v3D_vtx;
        double d_helix_DCA_r           = v3D_dca.Mag();
        double d_DCA_r_cut             = 3.0;
        TVector3 v3D_obj_DCA           = picoTrack->gDCA(pVtx);
        double d_obj_DCA               = v3D_obj_DCA.Mag();

        h_K_DCA_r      -> Fill(d_helix_DCA_r);
        h_K_obj_DCA_r  -> Fill(d_obj_DCA);
        h_K_diff_DCA_r -> Fill(d_helix_DCA_r-d_obj_DCA);

        if(d_helix_DCA_r > d_DCA_r_cut) b_bad_DCA == true;
        if(b_bad_DCA) continue;
        //Kaon DCA Cut
        /*test
        Kaoninfo.reset();
        */
        if(d_charge > 0.0)
        {
          hist_pt_y_kaonPlus->Fill(d_y_K,d_pT0);
          v_pri_tracks_pl.push_back(picoTrack);
        }
        else if(d_charge < 0.0)
        {
          hist_pt_y_kaonMinus->Fill(d_y_K,d_pT0);
          v_pri_tracks_mi.push_back(picoTrack);
        }
        v_pri_tracks.push_back(picoTrack);
        /* test
        if(d_charge < 0)      Kaoninfo.b_pos_charge = false;
        else                  Kaoninfo.b_pos_charge = true;
        Kaoninfo.runNumber   = runNumber;
        Kaoninfo.eventNumber = eventNumber;
        Kaoninfo.nGoodTracks = nGoodTracks;

        Kaoninfo.px = d_px0;
        Kaoninfo.py = d_py0;
        Kaoninfo.pz = d_pz0;

        Kaoninfo.tofBeta = tofBeta;

        Kaoninfo.x_vtx = v3D_vtx.x();
        Kaoninfo.y_vtx = v3D_vtx.y();
        Kaoninfo.z_vtx = v3D_vtx.z();

        Kaoninfo.d_TPCnSigmaKaon = d_TPCnSigmaKaon;

        Kaoninfo.DCA_r = d_obj_DCA;

        Kaoninfo.b_bad_TOF = b_bad_ToF ;//|| b_bad_TOF_match;
        t_K -> Fill();
        */
        // Fill Tree
        index++;

    }
    //=================== END Primary Track Loop =================================================

    i_event++;
    h_evt -> Fill(0.5);
  }
  //=============================== End Event Loop ===============================================

  outputFile->cd();

  // t_K ->Write();
  hist_pt_y_kaonPlus->Write();
  hist_pt_y_kaonMinus->Write();

  h_evt      -> Write();
  h_zvtx      -> Write();
  h_pT       -> Write();
  // h2_dEdx_PI_pq -> Write();
  // h2_dEdx_PRO_pq -> Write();
  // h2_dEdx_K_pq -> Write();
  h2_dEdx_All_pq->Write();
  h2_m2_QA_pq->Write();
  h2_m2_QA_pT->Write();
  h_mult     -> Write();
  h_DCA_r    -> Write();
  h_DCA_PRO    -> Write();
  h_DCA_PIM    -> Write();
  h_DCA_mother_r    -> Write();
  // h_rapidity_pm -> Write();
  // h2_mT_rapidity_pm -> Write();

  h_nodecaylength_cut_inv_m_PHI -> Write();
  h_nodipangle_cut_inv_m_PHI -> Write();
  h_inv_m_PHI    -> Write();
  h_prim_inv_m_PHI -> Write();
  h_inv_m_RHO    -> Write();
  h_inv_m_K0S    -> Write();
  h_inv_m_LAMBDA -> Write();

  h_decay_length -> Write();

  h_K_DCA_r       -> Write();
  h_K_obj_DCA_r   -> Write();
  h_K_diff_DCA_r  -> Write();

  // h2_TOF_nsigma_K_vs_PI  -> Write();
  h2_TPC_nsigma_K_vs_PI  -> Write();
  // h2_TOF_nsigma_K_vs_PRO -> Write();
  h2_TPC_nsigma_K_vs_PRO -> Write();

  // h_PHI_decay_length     -> Write();
  // h_dip_angle   -> Write();
}
//============================= End Main Function  ===============================================
