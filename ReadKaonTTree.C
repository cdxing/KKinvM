/**
  * \Update to analyze 3p85_FXT_2018
  * \author Ding Chen
  * \date Feb 24, 2020
  */
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStreamerElement.h>
#include <TStyle.h>

#include "Riostream.h"
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TLeaf.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGaxis.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TChain.h"
#include "TLegend.h"
#include "TFractionFitter.h"
#include "TVirtualFitter.h"
#include "TCut.h"
#include "TObject.h"
#include "TLine.h"
#include "TSpline.h"
#include "TPaveText.h"
#include "Fit/FitResult.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"
#include "TDatime.h"
#include "Math/MinimizerOptions.h"

using namespace std;

bool b_K_min_check = false;
bool b_K_pl_check  = false;
// Kaon + - check?
struct st_track
{
  bool  b_pos_charge;
  bool  b_bad_TOF;

  int   runNumber;
  int   eventNumber;
  int   nGoodTracks;

  double  px;
  double  py;
  double  pz;

  double  x_vtx;
  double  y_vtx;
  double  z_vtx;

  double  tofBeta;
  double  DCA_r;
  double  d_TPCnSigmaKaon;

};
//st_track structure

struct st_event
{
  vector<st_track> v_trk_pl;
  vector<st_track> v_trk_mi;
};

TFile * tf_in = 0;
TChain * t_K   = new TChain("t_K");
// TTree * t_K   = 0;

/////////////////////////////// Main Function //////////////////////////////////
void ReadKaonTTree( string FileName,
TString outFile = "test"//,
// double inputParameter1 = 0
)
{
  // Systematic analysis parameters
  // Double_t d_entryRange = inputParameter1;
  const double d_K_m        = 0.493696;

  outFile.Append(".readKTree.result.root");
  TFile * tf_KPlusEffTable_in   = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/KPlusEffTable.root","READ");
  TFile * tf_KMinusEffTable_in  = new TFile("/star/data01/pwg/dchen/offline/paper/psn0716/p_pi_v1/fixed_target_efficiency_acceptance_factors_from_embedding/KMinusEffTable.root","READ");
  TFile * tf_out = new TFile(outFile,"RECREATE");

  TH3D * h3_KPlusEffTable_input  = (TH3D*) tf_KPlusEffTable_in->Get("KPlusEffTable");
  TH3D * h3_KMinusEffTable_input = (TH3D*) tf_KMinusEffTable_in->Get("KMinusEffTable");

  TH1D * h_prim_inv_m_PHI     = new TH1D("h_prim_inv_m_PHI","h_prim_inv_m_PHI",100,0.9,1.1);
  TH1D * h_mx_prim_inv_m_PHI = new TH1D("h_mx_prim_inv_m_PHI","h_mx_prim_inv_m_PHI",100,0.9,1.1);
  TH1D * h_dip_angle = new TH1D("h_dip_angle","h_dip_angle",1000,-1,1.0);
  TH1D * h_Kmin_TPC2_TOF3_pT = new TH1D("h_Kmin_TPC2_TOF3_pT","h_Kmin_TPC2_TOF3_pT",200,0.0,10);
  TH1D * h_Kpl_TPC2_TOF3_pT  = new TH1D("h_Kpl_TPC2_TOF3_pT","h_Kpl_TPC2_TOF3_pT",200,0.0,10);
  TH1D * h_phi_TPC2_TOF3_pT  = new TH1D("h_phi_TPC2_TOF3_pT","h_phi_TPC2_TOF3_pT",200,0.0,10);
  TH1D * h_phi_TPC2_TOF3_mT  = new TH1D("h_phi_TPC2_TOF3_mT","h_phi_TPC2_TOF3_mT",200,0.0,10);
  TH1D * h_phi_TPC2_TOF3_y  = new TH1D("h_phi_TPC2_TOF3_y","h_phi_TPC2_TOF3_y",200,-10.,10);
  TH1D * h_phi_TPC2_TOF3_tight_y  = new TH1D("h_phi_TPC2_TOF3_tight_y","h_phi_TPC2_TOF3_tight_y",200,-10.0,10.0);

  // QA plots
  TH2D *hist_pt_y_kaonPlus = new TH2D("hist_pt_y_kaonPlus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_kaonPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH2D *hist_pt_y_kaonMinus = new TH2D("hist_pt_y_kaonMinus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  hist_pt_y_kaonMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");

  TH1D * h_eff_corr0                      = new TH1D("h_eff_corr0","h_eff_corr0",200,-10,10);
  TH1D * h_eff_corr1                      = new TH1D("h_eff_corr1","h_eff_corr1",200,-10,10);

  TH2D * h2_TOF_beta_pq                  = new TH2D("h2_TOF_beta_pq","h2_TOF_beta_pq",500,-2,2,500,0,3);
  TH2D * h2_dEdx_pq                      = new TH2D("h2_dEdx_pq","h2_dEdx_pq",500,-2,2,500,2e-6,10e-6);

  TH3D * h3_TOF_beta_pq_vs_nsig          = new TH3D("h3_TOF_beta_pq_vs_nsig","h3_TOF_beta_pq_vs_nsig",500,-2,2,500,0,3,100,0.0,10.0);
  TH3D * h3_dEdx_pq_vs_nsig              = new TH3D("h3_dEdx_pq_vs_nsig","h3_dEdx_pq_vs_nsig",500,-2,2,500,2e-6,10e-6,100,0.0,10.0);

  TH2D * h2_TOF_beta_TPC2sig_TOF_3sig_pq = new TH2D("h2_TOF_beta_TPC2sig_TOF_3sig_pq","h2_TOF_beta_TPC2sig_TOF_3sig_pq",500,-2,2,500,0,3);
  TH2D * h2_dEdx_TPC2sig_TOF_3sig_pq     = new TH2D("h2_dEdx_TPC2sig_TOF_3sig_pq","h2_dEdx_TPC2sig_TOF_3sig_pq",500,-2,2,500,2e-6,10e-6);

  TH1D * h_K_pl_pT  = new TH1D("h_K_pl_pT","h_K_pl_pT",200,0.0,10);
  TH1D * h_K_min_pT = new TH1D("h_K_min_pT","h_K_min_pT",200,0.0,10);

  TH1D * h_K_pl_eta  = new TH1D("h_K_pl_eta","h_K_pl_eta",200,-10,10);
  TH1D * h_K_min_eta = new TH1D("h_K_min_eta","h_K_min_eta",200,-10,10);
  // END QA plots
  TH1D * h_TPC2_TOF3_invM = new TH1D("h_TPC2_TOF3_invM","M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}",100,0.9,1.1);
  TH1D * h_TPC2_TOF3_invM_y_bin1 = new TH1D("h_TPC2_TOF3_invM_y_bin1","M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -2.0 <= y < -1.5",100,0.9,1.1);
  TH1D * h_TPC2_TOF3_invM_y_bin2 = new TH1D("h_TPC2_TOF3_invM_y_bin2","M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -1.5 <= y < -1.0",100,0.9,1.1);
  TH1D * h_TPC2_TOF3_invM_y_bin3 = new TH1D("h_TPC2_TOF3_invM_y_bin3","M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -1.0 <= y < -0.5",100,0.9,1.1);
  TH1D * h_TPC2_TOF3_invM_y_bin4 = new TH1D("h_TPC2_TOF3_invM_y_bin4","M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -0.5 <= y <= 0.0",100,0.9,1.1);

  TH1D * h_mx_TPC2_TOF3_invM = new TH1D("h_mx_TPC2_TOF3_invM","Mixed Events M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}",100,0.9,1.1);
  TH1D * h_mx_TPC2_TOF3_invM_y_bin1 = new TH1D("h_mx_TPC2_TOF3_invM_y_bin1","Mixed Events M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -2.0 <= y < -1.5",100,0.9,1.1);
  TH1D * h_mx_TPC2_TOF3_invM_y_bin2 = new TH1D("h_mx_TPC2_TOF3_invM_y_bin2","Mixed Events M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -1.5 <= y < -1.0",100,0.9,1.1);
  TH1D * h_mx_TPC2_TOF3_invM_y_bin3 = new TH1D("h_mx_TPC2_TOF3_invM_y_bin3","Mixed Events M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -1.0 <= y < -0.5",100,0.9,1.1);
  TH1D * h_mx_TPC2_TOF3_invM_y_bin4 = new TH1D("h_mx_TPC2_TOF3_invM_y_bin4","Mixed Events M_{K+K-} TPC < 2#sigma, 0.15 < M_{K}^{2} < 0.35(Gev/c^{2})^{2}, -0.5 <= y <= 0.0",100,0.9,1.1);

  TH1D *mHistKKInvMpT[12];
  for(int pt=0; pt<12; pt++)
  {
    mHistKKInvMpT[pt] = new TH1D(Form("histKKInvMpT%d",pt),Form("histKKInvMpT%d",pt),100,0.9,1.1);
  }

  TH1D *mHistKKInvMpTMixed[12];
  for(int pt=0; pt<12; pt++)
  {
    mHistKKInvMpTMixed[pt] = new TH1D(Form("histKKInvMpT%d_Mixed",pt),Form("histKKInvMpT%d_Mixed",pt),100,0.9,1.1);
  }
  // tf_in = TFile::Open( FileName );
  // if( !tf_in || tf_in->IsZombie() ){cout << "Error:Could not open, "<<FileName << endl;return;}
  // else{ cout << "Opened " << FileName << endl; }
  // // t_K = (TTree*) tf_in->Get("t_K");
  // string const dirFile = FileName.data();
  if( FileName.find(".list") != string::npos ||
      FileName.find(".lis") != string::npos ) {

    std::ifstream inputStream( FileName.c_str() );

    if(!inputStream) {
      cout << "ERROR: Cannot open list file " << FileName << endl;
    }

    Int_t nFile = 0;
    string file;
    while(getline(inputStream, file)) {
      if(file.find(".picoDst.result.root") != string::npos) {
        TFile* ftmp = TFile::Open(file.c_str());
        if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) {
          cout << " Read in picoDst file " << file << endl;
          t_K->Add(file.c_str());
          ++nFile;
        } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

        if (ftmp) {
          ftmp->Close();
        } //if (ftmp)
      } //if(file.find(".picoDst.root") != string::npos)
    } //while (getline(inputStream, file))

    cout << " Total " << nFile << " files have been read in. " << endl;
  } //if(FileName.find(".list") != string::npos || FileName.find(".lis" != string::npos))
  else if(FileName.find(".picoDst.result.root") != string::npos) {
    t_K->Add(FileName.c_str());
  }
  else {
    cout << " No good input file to read ... " << endl;
  }
  // t_K->Add(FileName.Data());

  if( t_K->GetEntries() == 0 ){cout << "Error:Could not find 't_K' in file: " << FileName << endl; return;}

  // if( t_K==0 ){cout << "Error:Could not find 't_K' in file: " << FileName << endl; return;}

  // int N_entries = t_K -> GetEntries();
  // Get the number of entries in the TTree

  int N_entries = /*(((d_entryRange+1)*5550) < (t_K -> GetEntries())) ? ((d_entryRange+1)*5550) : */ (t_K -> GetEntries());
  // multi jobs

  int pre_runNumber   = -9999;
  int pre_eventNumber = -9999;
  // To make sure run the next event
  vector<st_track> v_trk_pl;
  vector<st_track> v_trk_mi;

  vector<st_event> v_evt;
  ////////////////////////////// Read Kaon Info ////////////////////////////////
  int N_events = 0; //count # of events
  // unsigned long int N_max_events = 10;
  for( int i_entries = 0 /*(d_entryRange*5550)*/; i_entries< N_entries; i_entries++)
  {
    // if (i_entries > N_max_events) break;

    t_K -> GetEntry(i_entries);
    bool    b_pos_charge     = t_K-> GetLeaf("b_pos_charge")->GetValue(0);
    bool    b_bad_TOF        = t_K-> GetLeaf("b_bad_TOF")->GetValue(0);
    int     runNumber        = t_K-> GetLeaf("runNumber")->GetValue(0);
    int     eventNumber      = t_K-> GetLeaf("eventNumber")->GetValue(0);
    int     nGoodTracks      = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
    double  px               = t_K-> GetLeaf("px")->GetValue(0);
    double  py               = t_K-> GetLeaf("py")->GetValue(0);
    double  pz               = t_K-> GetLeaf("pz")->GetValue(0);
    double  x_vtx            = t_K-> GetLeaf("x_vtx")->GetValue(0);
    double  y_vtx            = t_K-> GetLeaf("y_vtx")->GetValue(0);
    double  z_vtx            = t_K-> GetLeaf("z_vtx")->GetValue(0);
    double  DCA_r            = t_K-> GetLeaf("DCA_r")->GetValue(0);
    double  d_TPCnSigmaKaon  = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);

    double d_tofBeta         = t_K -> GetLeaf("tofBeta")->GetValue(0);
    double dEdxFit           = t_K -> GetLeaf("dEdxFit")->GetValue(0);
    double d_dEdxTru         = t_K -> GetLeaf("d_dEdxTru")->GetValue(0);

    double d_inv_tofBeta     = -9999.0;

    double pT  = sqrt(px*px+py*py);
    double mom = sqrt(px*px+py*py+pz*pz);
    double E   = sqrt((px*px+py*py+pz*pz)+d_K_m*d_K_m);
    double y   = ((E-pz) != 0.0) ? 0.5*TMath::Log( (E + pz) / (E - pz) ) : -9999;
    double eta = ((mom - pz) != 0.0) ? 0.5*TMath::Log( (mom + pz) / (mom - pz) ) : -9999.0;
    double mass2      = mom*mom*((1.0/(d_tofBeta*d_tofBeta))-1.0);
    // Get useful values

    if(d_tofBeta != 0.0) d_inv_tofBeta = 1.0 / d_tofBeta;

    bool b_new_event = false;
    if((pre_runNumber != runNumber)||(pre_eventNumber != eventNumber))
    {
      N_events++;
      b_new_event     = true;
      pre_runNumber   = runNumber;
      pre_eventNumber = eventNumber;

      st_event evt;
      evt.v_trk_pl = v_trk_pl;
      evt.v_trk_mi = v_trk_mi;

      v_evt.push_back(evt);
      v_trk_pl.clear();
      v_trk_mi.clear();
    }
    // push back k+ k- tracks info into v_evt

    // PID QA histograms
    double d_q            = 1.;
    if(!b_pos_charge) d_q = -1.;
    double d_p            = sqrt(px*px+py*py+pz*pz);
    double d_pq           = fabs(d_p) * d_q;

    h2_TOF_beta_pq  -> Fill(d_pq,d_inv_tofBeta);
    h2_dEdx_pq      -> Fill(d_pq,d_dEdxTru);

    h3_dEdx_pq_vs_nsig      -> Fill(d_pq,d_dEdxTru,d_TPCnSigmaKaon);

    if((fabs(d_TPCnSigmaKaon) < 3.0)
    && (mass2 > 0.15)
    && (mass2 < 0.35))
    {
      h2_TOF_beta_TPC2sig_TOF_3sig_pq -> Fill(d_pq,d_inv_tofBeta);
      h2_dEdx_TPC2sig_TOF_3sig_pq     -> Fill(d_pq,d_dEdxTru);
    }

    st_track trk;
    trk.b_pos_charge    = b_pos_charge;
    trk.b_bad_TOF       = b_bad_TOF;
    trk.runNumber       = runNumber;
    trk.eventNumber     = eventNumber;
    trk.nGoodTracks     = nGoodTracks;
    trk.px              = px;
    trk.py              = py;
    trk.pz              = pz;
    trk.x_vtx           = x_vtx;
    trk.y_vtx           = y_vtx;
    trk.z_vtx           = z_vtx;
    trk.DCA_r           = DCA_r;
    trk.d_TPCnSigmaKaon = d_TPCnSigmaKaon;
    trk.tofBeta         = d_tofBeta;

    if(b_pos_charge) v_trk_pl.push_back(trk);
    else             v_trk_mi.push_back(trk);
    // push back track vector by trk

    bool b_K_eta   = (eta >= -1.47) && (eta <= 0.0);
    // if(!b_K_eta)   continue;
    // We want the psedorapidity region covered by TOF

    // if((y<-2.0)||(y>0.0)) continue;
    // rapidity cut

    // if(pT > 3.0) continue;

    if((fabs(d_TPCnSigmaKaon) < 3.0)
    && (mass2 > 0.15)
    && (mass2 < 0.35))
    {
      if(b_pos_charge) h_K_pl_pT -> Fill(pT);
      else  h_K_min_pT -> Fill(pT);

      if(b_pos_charge) h_K_pl_eta -> Fill(eta);
      else  h_K_min_eta -> Fill(eta);
    }
  }
  // fill track vectors and event vectors, Number of events draw some QA plots
  ////////////////////////// END Read Kaon Info ////////////////////////////////

    cout<<endl<<"  LOOPING OVER EVENTS "<< v_evt.size()<<endl<<endl;

    //////////////////////// Normal Invariant Mass Event Loop /////////////////////////
    for(unsigned int i = 0; i < v_evt.size(); i++)
    {

      st_event         evt      = v_evt[i];
      vector<st_track> v_trk_pl = v_evt[i].v_trk_pl;
      vector<st_track> v_trk_mi = v_evt[i].v_trk_mi;
      ////////////////////// Positive Track Loop ///////////////////////////////
      for(unsigned int j = 0; j < v_trk_pl.size(); j++)
      {
        st_track trk0 = v_trk_pl[j];

        bool    b_pos_charge0    = trk0.b_pos_charge;
        bool    b_bad_TOF0       = trk0.b_bad_TOF;
        int     runNumber0       = trk0.runNumber;
        int     eventNumber0     = trk0.eventNumber;
        int     nGoodTracks0     = trk0.nGoodTracks;
        double  px0              = trk0.px;
        double  py0              = trk0.py;
        double  pz0              = trk0.pz;
        double  x_vtx0           = trk0.x_vtx;
        double  y_vtx0           = trk0.y_vtx;
        double  z_vtx0           = trk0.z_vtx;
        double  DCA_r0           = trk0.DCA_r;
        double  d_TPCnSigmaKaon0 = trk0.d_TPCnSigmaKaon;
        double  d_tofBeta0       = trk0.tofBeta;

        double E_0   = sqrt((px0*px0+py0*py0+pz0*pz0)+d_K_m*d_K_m);
        double y_0   = ((E_0-pz0) != 0.0) ? 0.5*TMath::Log( (E_0 + pz0) / (E_0 - pz0) ) : -9999;

        double d_pT0 = sqrt(px0*px0+py0*py0);
        double d_mom0 = sqrt(px0*px0+py0*py0+pz0*pz0);
        double eta0 = ((d_mom0 - pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + pz0) / (d_mom0 - pz0) ) : 1.0;

        double mass2_0      = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);

        hist_pt_y_kaonPlus->Fill(y_0,d_pT0);
        ////////////////////// Negative Track Loop ///////////////////////////////
        for(unsigned int k = 0; k < v_trk_mi.size(); k++)
        {
          st_track trk1 = v_trk_mi[k];

          bool    b_pos_charge1    = trk1.b_pos_charge;
  	      bool    b_bad_TOF1       = trk1.b_bad_TOF;
  	      int     runNumber1       = trk1.runNumber;
  	      int     eventNumber1     = trk1.eventNumber;
  	      int     nGoodTracks1     = trk1.nGoodTracks;
  	      double  px1              = trk1.px;
  	      double  py1              = trk1.py;
  	      double  pz1              = trk1.pz;
  	      double  x_vtx1           = trk1.x_vtx;
  	      double  y_vtx1           = trk1.y_vtx;
  	      double  z_vtx1           = trk1.z_vtx;
  	      double  DCA_r1           = trk1.DCA_r;
  	      double  d_TPCnSigmaKaon1 = trk1.d_TPCnSigmaKaon;
          double  d_tofBeta1       = trk1.tofBeta;

          if(b_K_min_check) d_TPCnSigmaKaon1 = 0.0;
          //turn off pid for k minus to get more statistics

          double d_pT1 = sqrt(px1*px1+py1*py1);
          double d_mom1 = sqrt(px1*px1+py1*py1+pz1*pz1);
          double eta1 = ((d_mom1 - pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + pz1) / (d_mom1 - pz1) ) : 1.0;

          if(runNumber0 != runNumber1) continue;
  	      if(eventNumber0 != eventNumber1) continue;
          //double check events
          bool b_PHI = true;

          double d_dip_angle = TMath::ACos((d_pT0*d_pT1+pz0*pz1) / (d_mom0*d_mom1) );
  	      h_dip_angle        -> Fill(d_dip_angle);

          Int_t centrality0 = 0;
          Int_t centrality1 = 0;

          double d_M1 = d_K_m;
          double d_M0 = d_K_m;

          double d_E0 = sqrt((px0*px0+py0*py0+pz0*pz0)+d_M0*d_M0);
          double d_E1 = sqrt((px1*px1+py1*py1+pz1*pz1)+d_M1*d_M1);

          double d_y0 = 0.5*TMath::Log((d_E0 + pz0)/(d_E0 - pz0));
          double d_y1 = 0.5*TMath::Log((d_E1 + pz1)/(d_E1 - pz1));

          double d_mT0        = sqrt(d_pT0*d_pT0 + d_M0*d_M0);
          double d_mT1        = sqrt(d_pT1*d_pT1 + d_M1*d_M1);

          hist_pt_y_kaonMinus->Fill(d_y1,d_pT1);

          Int_t i_ybin0 = h3_KPlusEffTable_input->GetYaxis()->FindBin(d_y0);
          Int_t i_zbin0 = h3_KPlusEffTable_input->GetZaxis()->FindBin(d_mT0-d_M0);

          Int_t i_ybin1 = h3_KMinusEffTable_input->GetYaxis()->FindBin(d_y1);
          Int_t i_zbin1 = h3_KMinusEffTable_input->GetZaxis()->FindBin(d_mT1-d_M1);

          // double d_eff_corr0 = h3_KPlusEffTable_input ->GetBinContent(centrality0,i_ybin0,i_zbin0);
          // double d_eff_corr1 = h3_KMinusEffTable_input->GetBinContent(centrality1,i_ybin1,i_zbin1);

          // efficiency corrections
          // d_eff_corr0 = (d_eff_corr0 <= 0.01 || d_eff_corr0 >= 1) ? 1 : d_eff_corr0;
          // d_eff_corr1 = (d_eff_corr1 <= 0.01 || d_eff_corr1 >= 1) ? 1 : d_eff_corr1;

          // Primary normal events
          double d_eff_corr0 = 1.0;
          double d_eff_corr1 = 1.0;


          h_eff_corr0->Fill(d_eff_corr0);
          h_eff_corr1->Fill(d_eff_corr1);

          double d_inv_m = sqrt(  d_M0*d_M0
                                + d_M1*d_M1
                                + 2.0 *d_E0*d_E1
                                - 2.0 *(px0*px1+py0*py1+pz0*pz1) );

          if(b_PHI) h_prim_inv_m_PHI    -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
          // Raw distribution of Phi invariant mass
          double d_min_DCA     = (DCA_r0 < DCA_r1) ? DCA_r0 : DCA_r1;
          double d_max_TPCnsig = (fabs(d_TPCnSigmaKaon0) > fabs(d_TPCnSigmaKaon1)) ? d_TPCnSigmaKaon0 : d_TPCnSigmaKaon1;

          int i_max_multi      = (nGoodTracks0 > nGoodTracks1) ? nGoodTracks0 : nGoodTracks1;
          double d_max_multi   = (double) i_max_multi;

          bool b_pileup0   = (nGoodTracks0 > 270);
          bool b_pileup1   = (nGoodTracks1 > 270);
          bool b_low_mult0 = (nGoodTracks0 <= 10);
          bool b_low_mult1 = (nGoodTracks1 <= 10);

          if(b_pileup0||b_pileup1) continue;

          bool a_b_cent0[10];

  	      a_b_cent0[0] = ( nGoodTracks0 >= 200)   && ( nGoodTracks0 <= 270);
  	      a_b_cent0[1] = ( nGoodTracks0 >= 150) && ( nGoodTracks0 < 200);
  	      a_b_cent0[2] = ( nGoodTracks0 >= 124)  && ( nGoodTracks0 < 150);
  	      a_b_cent0[3] = ( nGoodTracks0 >= 100)   && ( nGoodTracks0 < 124);
          a_b_cent0[4] = ( nGoodTracks0 >= 72)   && ( nGoodTracks0 < 100);
          a_b_cent0[5] = ( nGoodTracks0 >= 50) && ( nGoodTracks0 < 72);
          a_b_cent0[6] = ( nGoodTracks0 >= 40)  && ( nGoodTracks0 < 50);
          a_b_cent0[7] = ( nGoodTracks0 >= 30)   && ( nGoodTracks0 < 40);
          a_b_cent0[8] = ( nGoodTracks0 >= 20)  && ( nGoodTracks0 < 30);
          a_b_cent0[9] = ( nGoodTracks0 >= 10)   && ( nGoodTracks0 < 20);

          bool a_b_cent1[10];
          a_b_cent1[0] = ( nGoodTracks1 >= 200)   && ( nGoodTracks1 <= 270);
  	      a_b_cent1[1] = ( nGoodTracks1 >= 150) && ( nGoodTracks1 < 200);
  	      a_b_cent1[2] = ( nGoodTracks1 >= 124)  && ( nGoodTracks1 < 150);
  	      a_b_cent1[3] = ( nGoodTracks1 >= 100)   && ( nGoodTracks1 < 124);
          a_b_cent1[4] = ( nGoodTracks1 >= 72)   && ( nGoodTracks1 < 100);
          a_b_cent1[5] = ( nGoodTracks1 >= 50) && ( nGoodTracks1 < 72);
          a_b_cent1[6] = ( nGoodTracks1 >= 40)  && ( nGoodTracks1 < 50);
          a_b_cent1[7] = ( nGoodTracks1 >= 30)   && ( nGoodTracks1 < 40);
          a_b_cent1[8] = ( nGoodTracks1 >= 20)  && ( nGoodTracks1 < 30);
          a_b_cent1[9] = ( nGoodTracks1 >= 10)   && ( nGoodTracks1 < 20);

          // 10-70% centrality
          if( nGoodTracks0>=200 || nGoodTracks1>=200||nGoodTracks0<40||nGoodTracks1<40)

          bool a_b_pT0[12];
          a_b_pT0[0]    = (d_pT0 >= 0.0) && (d_pT0 < 0.3);
          a_b_pT0[1]    = (d_pT0 >= 0.3) && (d_pT0 < 0.6);
          a_b_pT0[2]    = (d_pT0 >= 0.6) && (d_pT0 < 0.9);
          a_b_pT0[3]    = (d_pT0 >= 0.9) && (d_pT0 < 1.2);
          a_b_pT0[4]    = (d_pT0 >= 1.2) && (d_pT0 < 1.5);
          a_b_pT0[5]    = (d_pT0 >= 1.5) && (d_pT0 < 1.8);
          a_b_pT0[6]    = (d_pT0 >= 1.8) && (d_pT0 < 2.1);
          a_b_pT0[7]    = (d_pT0 >= 2.1) && (d_pT0 < 2.4);
          a_b_pT0[8]    = (d_pT0 >= 2.4) && (d_pT0 < 2.7);
          a_b_pT0[9]    = (d_pT0 >= 2.7) && (d_pT0 < 3.0);
          a_b_pT0[10]    = (d_pT0 >= 3.0) && (d_pT0 < 3.3);
          a_b_pT0[11]    = (d_pT0 >= 3.3) && (d_pT0 < 3.6);

          bool a_b_pT1[12];
          a_b_pT1[0]    = (d_pT1 >= 0.0) && (d_pT1 < 0.3);
          a_b_pT1[1]    = (d_pT1 >= 0.3) && (d_pT1 < 0.6);
          a_b_pT1[2]    = (d_pT1 >= 0.6) && (d_pT1 < 0.9);
          a_b_pT1[3]    = (d_pT1 >= 0.9) && (d_pT1 < 1.2);
          a_b_pT1[4]    = (d_pT1 >= 1.2) && (d_pT1 < 1.5);
          a_b_pT1[5]    = (d_pT1 >= 1.5) && (d_pT1 < 1.8);
          a_b_pT1[6]    = (d_pT1 >= 1.8) && (d_pT1 < 2.1);
          a_b_pT1[7]    = (d_pT1 >= 2.1) && (d_pT1 < 2.4);
          a_b_pT1[8]    = (d_pT1 >= 2.4) && (d_pT1 < 2.7);
          a_b_pT1[9]    = (d_pT1 >= 2.7) && (d_pT1 < 3.0);
          a_b_pT1[10]    = (d_pT1 >= 3.0) && (d_pT1 < 3.3);
          a_b_pT1[11]    = (d_pT1 >= 3.3) && (d_pT1 < 3.6);

          bool b_K0_eta   = (eta0 >= -1.47) && (eta0 <= 0.0);
  	      bool b_K1_eta   = (eta1 >= -1.47) && (eta1 <= 0.0);
          // eta fiducial cut

          double d_Phi_pT = sqrt(px0*px0 + py0*py0 +px1*px1 +py1+py1 + 2.*px0*px1 + 2.*py0*py1);
          bool a_b_Phi_pT[12]={};
          a_b_Phi_pT[0]    = (d_Phi_pT >= 0.0) && (d_Phi_pT < 0.3);
          a_b_Phi_pT[1]    = (d_Phi_pT >= 0.3) && (d_Phi_pT < 0.6);
          a_b_Phi_pT[2]    = (d_Phi_pT >= 0.6) && (d_Phi_pT < 0.9);
          a_b_Phi_pT[3]    = (d_Phi_pT >= 0.9) && (d_Phi_pT < 1.2);
          a_b_Phi_pT[4]    = (d_Phi_pT >= 1.2) && (d_Phi_pT < 1.5);
          a_b_Phi_pT[5]    = (d_Phi_pT >= 1.5) && (d_Phi_pT < 1.8);
          a_b_Phi_pT[6]    = (d_Phi_pT >= 1.8) && (d_Phi_pT < 2.1);
          a_b_Phi_pT[7]    = (d_Phi_pT >= 2.1) && (d_Phi_pT < 2.4);
          a_b_Phi_pT[8]    = (d_Phi_pT >= 2.4) && (d_Phi_pT < 2.7);
          a_b_Phi_pT[9]    = (d_Phi_pT >= 2.7) && (d_Phi_pT < 3.0);
          a_b_Phi_pT[10]    = (d_Phi_pT >= 3.0) && (d_Phi_pT < 3.3);
          a_b_Phi_pT[11]    = (d_Phi_pT >= 3.3) && (d_Phi_pT < 3.6);

          double m_phi = 1.019455;
          double d_mT_phi = sqrt(d_Phi_pT*d_Phi_pT + m_phi*m_phi );
          double d_phi_pz = pz0+pz1;
  	      double d_phi_E  = d_E0+d_E1;
  	      double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;

          h_phi_TPC2_TOF3_y -> Fill(d_phi_y);

          // if(d_mT_phi > 2.0) continue;

          if((d_inv_m <= 0.9) || (d_inv_m >= 1.1)) continue;

          if(nGoodTracks0 != nGoodTracks1) continue;

          if(b_K_pl_check) d_max_TPCnsig = d_TPCnSigmaKaon1;
          // Turn off the pid for K minus to get more statistics

          /* Delete PID in read kaon tree test
          if( (fabs(d_max_TPCnsig) <= 2.0) && (mass2_1 > 0.15) && (mass2_1 < 0.35) && (mass2_0 > 0.15) && (mass2_0 < 0.35) )
             {*/
          if(j==0) h_Kmin_TPC2_TOF3_pT -> Fill( d_pT1);
          if(k==0) h_Kpl_TPC2_TOF3_pT  -> Fill( d_pT0);
          if(b_PHI)
          {
            for(int pt=0; pt<12; pt++)
            {
              if(a_b_Phi_pT[pt]) mHistKKInvMpT[pt]->Fill(d_inv_m);
              h_phi_TPC2_TOF3_pT -> Fill(d_Phi_pT);
              h_phi_TPC2_TOF3_mT -> Fill(d_mT_phi-m_phi);
            }
          }

            if(b_PHI&&(a_b_cent0[2]||a_b_cent0[1]||a_b_cent0[0])&&(a_b_cent1[2]||a_b_cent1[1]||a_b_cent1[0]))
            {
              if((d_inv_m > 1.00138) && (d_inv_m < 1.03956)) h_phi_TPC2_TOF3_tight_y -> Fill(d_phi_y);

              h_TPC2_TOF3_invM -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));

              if(d_phi_y>=-2.0 && d_phi_y<-1.5) h_TPC2_TOF3_invM_y_bin1 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
              if(d_phi_y>=-1.5 && d_phi_y<-1.0) h_TPC2_TOF3_invM_y_bin2 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
              if(d_phi_y>=-1.0 && d_phi_y<-0.5) h_TPC2_TOF3_invM_y_bin3 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
              if(d_phi_y>=-0.5 && d_phi_y<=0.0) h_TPC2_TOF3_invM_y_bin4 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));

            }
          // }
        }
        ////////////////// End Negative Track Loop ///////////////////////////////
      }
      ////////////////// End Positive Track Loop ///////////////////////////////
    }
    //////////////////// END Normal Invariant mass Event Loop /////////////////////////

    int j_start = 0;
    unsigned long long pre_event_pl_runnumber = 9999;
    unsigned long long this_mixed_event_pl_runnumber = 9999;

    // make a list of events used in mixed events
    bool b_found_mixed_evt = false;
    unsigned int i_found_mixed_evt = 0;

    bool b_half  = false;
    bool b_quart = false;
    set<unsigned long long> s_used_mx_events;
    // ////////////////// Mixed Invariant mass Event Loop ///////////////////////
    // for(int i = (d_entryRange*5550); i < N_entries; i++)
    // {
    //   // if(i > N_max_events) break;
    //   if(( i > ((double) (N_entries-d_entryRange*5550)/4.0))&&(!b_quart)&&(!b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
    //   if(( i > ((double) (3.0*(N_entries-d_entryRange*5550))/4.0))&&(!b_quart)&&(b_half)) {cout<<" 0.25 done"<<endl; b_quart = true;}
    //   if(( i > ((double) (N_entries-d_entryRange*5550)/2.0))&&(!b_half)) {cout<<" halfway"<<endl; b_half = true; b_quart = false;}
    //   t_K -> GetEntry(i);
    //   bool b_pos_charge0 = t_K -> GetLeaf("b_pos_charge")->GetValue(0);
    //   unsigned int i_runNumber0   = t_K -> GetLeaf("runNumber")->GetValue(0);
    //   unsigned int i_eventNumber0 = t_K -> GetLeaf("eventNumber")->GetValue(0);
    //   unsigned long long event_pl_runnumber0 = (i_runNumber0-19150000)*10000000+i_eventNumber0;
    //   if(event_pl_runnumber0 != pre_event_pl_runnumber) {j_start = i; pre_event_pl_runnumber = event_pl_runnumber0; b_found_mixed_evt = false; }
    //
    //   int    nGoodTracks0  = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
    //
    //   double d_xvtx0 = t_K -> GetLeaf("x_vtx")->GetValue(0);
    //   double d_yvtx0 = t_K -> GetLeaf("y_vtx")->GetValue(0);
    //   double d_zvtx0 = t_K -> GetLeaf("z_vtx")->GetValue(0);
    //
    //   double px0 = t_K -> GetLeaf("px")->GetValue(0);
    //   double py0 = t_K -> GetLeaf("py")->GetValue(0);
    //   double pz0 = t_K -> GetLeaf("pz")->GetValue(0);
    //
    //   double  d_tofBeta0       = t_K -> GetLeaf("tofBeta")->GetValue(0);
    //   double  DCA_r0           = t_K-> GetLeaf("DCA_r")->GetValue(0);
    //   double  d_TPCnSigmaKaon0 = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
    //
    //   double d_pT0    = sqrt(px0*px0+py0*py0);
    //   double d_mom0   = sqrt(px0*px0+py0*py0+pz0*pz0);
    //   double eta0     = ((d_mom0 - pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + pz0) / (d_mom0 - pz0) ) : 1.0;
    //   double mass2_0  = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);
    //
    //   ///////////////////////// Mixed Track Loop ///////////////////////////////
    //   for(int j = j_start; j < N_entries; j++)
    //   {
    //     t_K -> GetEntry(j);
    //     unsigned int i_runNumber1   = t_K -> GetLeaf("runNumber")->GetValue(0);
    //     unsigned int i_eventNumber1 = t_K -> GetLeaf("eventNumber")->GetValue(0);
    //     unsigned long long event_pl_runnumber1 = (i_runNumber1-19150000)*10000000+i_eventNumber1;
    //
    //     int     nGoodTracks1     = t_K-> GetLeaf("nGoodTracks")->GetValue(0);
    //     double d_xvtx1 = t_K -> GetLeaf("x_vtx")->GetValue(0);
    //     double d_yvtx1 = t_K -> GetLeaf("y_vtx")->GetValue(0);
    //     double d_zvtx1 = t_K -> GetLeaf("z_vtx")->GetValue(0);
    //
    //
    //
    //     if(event_pl_runnumber0 == event_pl_runnumber1)
    // 	    {
    // 	      //cout<<"                            same event continueing"<<endl;
    // 	      continue;
    // 	    }
    // 	  if(b_found_mixed_evt && (event_pl_runnumber1 != this_mixed_event_pl_runnumber) && i_found_mixed_evt == 5 )
    // 	    {
    // 	      //cout<<" hit next mixed event breaking "<<event_pl_runnumber1 <<endl;
    // 	      break;
    // 	    }
    // 	  if(fabs(d_zvtx1 - d_zvtx0)>0.5)
    // 	    {
    // 	      // cout<<"                            bad zvtx continuing "<<event_pl_runnumber1<<endl;
    // 	      continue;
    // 	    }
    //
    //     if(!b_found_mixed_evt) s_used_mx_events.insert(event_pl_runnumber1);
    //     b_found_mixed_evt = true;
    //     i_found_mixed_evt++;
    //     this_mixed_event_pl_runnumber = event_pl_runnumber1;
    //
    //     bool b_pos_charge1 = t_K -> GetLeaf("b_pos_charge")->GetValue(0);
    //
    //     if(b_pos_charge0 == b_pos_charge1) continue;
    //
    //
    //     double px1 = t_K -> GetLeaf("px")->GetValue(0);
    // 	  double py1 = t_K -> GetLeaf("py")->GetValue(0);
    // 	  double pz1 = t_K -> GetLeaf("pz")->GetValue(0);
    //
    //     double  DCA_r1           = t_K-> GetLeaf("DCA_r")->GetValue(0);
    //     double  d_TPCnSigmaKaon1 = t_K-> GetLeaf("d_TPCnSigmaKaon")->GetValue(0);
    //     double  d_tofBeta1       = t_K -> GetLeaf("tofBeta")->GetValue(0);
    //
    //     double d_M1 = d_K_m;
    //     double d_M0 = d_K_m;
    //     double d_E0 = sqrt((px0*px0+py0*py0+pz0*pz0)+d_M0*d_M0);
    //     double d_E1 = sqrt((px1*px1+py1*py1+pz1*pz1)+d_M1*d_M1);
    //
    //     double d_pT1   = sqrt(px1*px1+py1*py1);
    //     double d_mom1  = sqrt(px1*px1+py1*py1+pz1*pz1);
    //     double eta1    = ((d_mom1 - pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + pz1) / (d_mom1 - pz1) ) : 1.0;
    //     double mass2_1 = d_mom1*d_mom1*((1.0/(d_tofBeta1*d_tofBeta1))-1.0);
    //
    //     double d_y0 = 0.5*TMath::Log((d_E0 + pz0)/(d_E0 - pz0));
    //     double d_y1 = 0.5*TMath::Log((d_E1 + pz1)/(d_E1 - pz1));
    //
    //     double d_mT0        = sqrt(d_pT0*d_pT0 + d_M0*d_M0);
    //     double d_mT1        = sqrt(d_pT1*d_pT1 + d_M1*d_M1);
    //
    //     double d_inv_m = sqrt( d_M0*d_M0
    //                           +d_M1*d_M1
    //                           +2.0*d_E0*d_E1
    //                           -2.0*(px0*px1+py0*py1+pz0*pz1) );
    //
    //
    //     bool   b_PHI       = true;
    //     double d_dip_angle = TMath::ACos((d_pT0*d_pT1+pz0*pz1) / (d_mom0*d_mom1) );
    //
    //     Int_t centrality0 = 0;
    //     Int_t centrality1 = 0;
    //
    //     bool b_pileup0   = (nGoodTracks0 > 270);
    //     bool b_pileup1   = (nGoodTracks1 > 270);
    //     bool b_low_mult0 = (nGoodTracks0 <= 10);
    //     bool b_low_mult1 = (nGoodTracks1 <= 10);
    //
    //     if(b_pileup0||b_pileup1) continue;
    //
    //     bool a_b_cent0[10];
    //
    //     a_b_cent0[0] = ( nGoodTracks0 >= 200)   && ( nGoodTracks0 <= 270);
    //     a_b_cent0[1] = ( nGoodTracks0 >= 150) && ( nGoodTracks0 < 200);
    //     a_b_cent0[2] = ( nGoodTracks0 >= 124)  && ( nGoodTracks0 < 150);
    //     a_b_cent0[3] = ( nGoodTracks0 >= 100)   && ( nGoodTracks0 < 124);
    //     a_b_cent0[4] = ( nGoodTracks0 >= 72)   && ( nGoodTracks0 < 100);
    //     a_b_cent0[5] = ( nGoodTracks0 >= 50) && ( nGoodTracks0 < 72);
    //     a_b_cent0[6] = ( nGoodTracks0 >= 40)  && ( nGoodTracks0 < 50);
    //     a_b_cent0[7] = ( nGoodTracks0 >= 30)   && ( nGoodTracks0 < 40);
    //     a_b_cent0[8] = ( nGoodTracks0 >= 20)  && ( nGoodTracks0 < 30);
    //     a_b_cent0[9] = ( nGoodTracks0 >= 10)   && ( nGoodTracks0 < 20);
    //
    //     bool a_b_cent1[10];
    //     a_b_cent1[0] = ( nGoodTracks1 >= 200)   && ( nGoodTracks1 <= 270);
    //     a_b_cent1[1] = ( nGoodTracks1 >= 150) && ( nGoodTracks1 < 200);
    //     a_b_cent1[2] = ( nGoodTracks1 >= 124)  && ( nGoodTracks1 < 150);
    //     a_b_cent1[3] = ( nGoodTracks1 >= 100)   && ( nGoodTracks1 < 124);
    //     a_b_cent1[4] = ( nGoodTracks1 >= 72)   && ( nGoodTracks1 < 100);
    //     a_b_cent1[5] = ( nGoodTracks1 >= 50) && ( nGoodTracks1 < 72);
    //     a_b_cent1[6] = ( nGoodTracks1 >= 40)  && ( nGoodTracks1 < 50);
    //     a_b_cent1[7] = ( nGoodTracks1 >= 30)   && ( nGoodTracks1 < 40);
    //     a_b_cent1[8] = ( nGoodTracks1 >= 20)  && ( nGoodTracks1 < 30);
    //     a_b_cent1[9] = ( nGoodTracks1 >= 10)   && ( nGoodTracks1 < 20);
    //
    //     bool a_b_pT0[12];
    //     a_b_pT0[0]    = (d_pT0 >= 0.0) && (d_pT0 < 0.3);
    //     a_b_pT0[1]    = (d_pT0 >= 0.3) && (d_pT0 < 0.6);
    //     a_b_pT0[2]    = (d_pT0 >= 0.6) && (d_pT0 < 0.9);
    //     a_b_pT0[3]    = (d_pT0 >= 0.9) && (d_pT0 < 1.2);
    //     a_b_pT0[4]    = (d_pT0 >= 1.2) && (d_pT0 < 1.5);
    //     a_b_pT0[5]    = (d_pT0 >= 1.5) && (d_pT0 < 1.8);
    //     a_b_pT0[6]    = (d_pT0 >= 1.8) && (d_pT0 < 2.1);
    //     a_b_pT0[7]    = (d_pT0 >= 2.1) && (d_pT0 < 2.4);
    //     a_b_pT0[8]    = (d_pT0 >= 2.4) && (d_pT0 < 2.7);
    //     a_b_pT0[9]    = (d_pT0 >= 2.7) && (d_pT0 < 3.0);
    //     a_b_pT0[10]    = (d_pT0 >= 3.0) && (d_pT0 < 3.3);
    //     a_b_pT0[11]    = (d_pT0 >= 3.3) && (d_pT0 < 3.6);
    //
    //     bool a_b_pT1[12];
    //     a_b_pT1[0]    = (d_pT1 >= 0.0) && (d_pT1 < 0.3);
    //     a_b_pT1[1]    = (d_pT1 >= 0.3) && (d_pT1 < 0.6);
    //     a_b_pT1[2]    = (d_pT1 >= 0.6) && (d_pT1 < 0.9);
    //     a_b_pT1[3]    = (d_pT1 >= 0.9) && (d_pT1 < 1.2);
    //     a_b_pT1[4]    = (d_pT1 >= 1.2) && (d_pT1 < 1.5);
    //     a_b_pT1[5]    = (d_pT1 >= 1.5) && (d_pT1 < 1.8);
    //     a_b_pT1[6]    = (d_pT1 >= 1.8) && (d_pT1 < 2.1);
    //     a_b_pT1[7]    = (d_pT1 >= 2.1) && (d_pT1 < 2.4);
    //     a_b_pT1[8]    = (d_pT1 >= 2.4) && (d_pT1 < 2.7);
    //     a_b_pT1[9]    = (d_pT1 >= 2.7) && (d_pT1 < 3.0);
    //     a_b_pT1[10]    = (d_pT1 >= 3.0) && (d_pT1 < 3.3);
    //     a_b_pT1[11]    = (d_pT1 >= 3.3) && (d_pT1 < 3.6);
    //
    //
    //     double d_Phi_pT = sqrt(px0*px0 + py0*py0 +px1*px1 +py1+py1 + 2.*px0*px1 + 2.*py0*py1);
    //     bool a_b_Phi_pT[12]={};
    //     a_b_Phi_pT[0]    = (d_Phi_pT >= 0.0) && (d_Phi_pT < 0.3);
    //     a_b_Phi_pT[1]    = (d_Phi_pT >= 0.3) && (d_Phi_pT < 0.6);
    //     a_b_Phi_pT[2]    = (d_Phi_pT >= 0.6) && (d_Phi_pT < 0.9);
    //     a_b_Phi_pT[3]    = (d_Phi_pT >= 0.9) && (d_Phi_pT < 1.2);
    //     a_b_Phi_pT[4]    = (d_Phi_pT >= 1.2) && (d_Phi_pT < 1.5);
    //     a_b_Phi_pT[5]    = (d_Phi_pT >= 1.5) && (d_Phi_pT < 1.8);
    //     a_b_Phi_pT[6]    = (d_Phi_pT >= 1.8) && (d_Phi_pT < 2.1);
    //     a_b_Phi_pT[7]    = (d_Phi_pT >= 2.1) && (d_Phi_pT < 2.4);
    //     a_b_Phi_pT[8]    = (d_Phi_pT >= 2.4) && (d_Phi_pT < 2.7);
    //     a_b_Phi_pT[9]    = (d_Phi_pT >= 2.7) && (d_Phi_pT < 3.0);
    //     a_b_Phi_pT[10]    = (d_Phi_pT >= 3.0) && (d_Phi_pT < 3.3);
    //     a_b_Phi_pT[11]    = (d_Phi_pT >= 3.3) && (d_Phi_pT < 3.6);
    //
    //     Int_t i_ybin0 = h3_KPlusEffTable_input->GetYaxis()->FindBin(d_y0);
    //     Int_t i_zbin0 = h3_KPlusEffTable_input->GetZaxis()->FindBin(d_mT0-d_M0);
    //
    //     Int_t i_ybin1 = h3_KMinusEffTable_input->GetYaxis()->FindBin(d_y1);
    //     Int_t i_zbin1 = h3_KMinusEffTable_input->GetZaxis()->FindBin(d_mT1-d_M1);
    //
    //     double d_eff_corr0 = h3_KPlusEffTable_input ->GetBinContent(centrality0,i_ybin0,i_zbin0);
    //     double d_eff_corr1 = h3_KMinusEffTable_input->GetBinContent(centrality1,i_ybin1,i_zbin1);
    //
    //     // efficiency corrections
    //     d_eff_corr0 = (d_eff_corr0 <= 0.01 || d_eff_corr0 >= 1) ? 1 : d_eff_corr0;
    //     d_eff_corr1 = (d_eff_corr1 <= 0.01 || d_eff_corr1 >= 1) ? 1 : d_eff_corr1;
    //
    //     h_eff_corr0->Fill(d_eff_corr0);
    //     h_eff_corr1->Fill(d_eff_corr1);
    //
    //     double d_min_DCA     = (DCA_r0 < DCA_r1) ? DCA_r0 : DCA_r1;
    //     double d_max_TPCnsig = (fabs(d_TPCnSigmaKaon0) > fabs(d_TPCnSigmaKaon1)) ? d_TPCnSigmaKaon0 : d_TPCnSigmaKaon1;
    //
    //     int i_max_multi      = (nGoodTracks0 > nGoodTracks1) ? nGoodTracks0 : nGoodTracks1;
    // 	  double d_max_multi   = (double) i_max_multi;
    //
    //     if(b_PHI) h_mx_prim_inv_m_PHI    -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
    //
    //     // centrality cut
    //
    //     //eta fiducial cut
    // 	  bool b_K0_eta   = (eta0 >= -1.47) && (eta0 <= 0.0);
    // 	  bool b_K1_eta   = (eta1 >= -1.47) && (eta1 <= 0.0);
    //
    // 	  double m_phi = 1.019455;
    // 	  double d_mT_phi = sqrt(d_Phi_pT*d_Phi_pT + m_phi*m_phi );
    //
    // 	  double d_phi_pz = pz0+pz1;
    // 	  double d_phi_E  = d_E0+d_E1;
    // 	  double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;
    //
    //     // if(d_Phi_pT > 3.0) continue;
    //
    //     // if((d_phi_y<-2.0)||(d_phi_y>0.0)) continue;
    //
    //     if((d_inv_m <= 0.9) || (d_inv_m >= 1.1)) continue;
    //
    //     if(b_K_pl_check)  d_max_TPCnsig = d_TPCnSigmaKaon1;
    //
    //     if(b_PHI)
    //     {
    //       for(int pt=0; pt<12; pt++)
    //       {
    //         if(a_b_Phi_pT[pt]) mHistKKInvMpTMixed[pt]->Fill(d_inv_m);
    //       }
    //
    //       h_mx_TPC2_TOF3_invM -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
    //       if(d_phi_y>=-2.0 && d_phi_y<-1.5) h_mx_TPC2_TOF3_invM_y_bin1 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
    //       if(d_phi_y>=-1.5 && d_phi_y<-1.0) h_mx_TPC2_TOF3_invM_y_bin2 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
    //       if(d_phi_y>=-1.0 && d_phi_y<-0.5) h_mx_TPC2_TOF3_invM_y_bin3 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
    //       if(d_phi_y>=-0.5 && d_phi_y<=0.0) h_mx_TPC2_TOF3_invM_y_bin4 -> Fill(d_inv_m,1/(d_eff_corr0*d_eff_corr1));
    //
    //     } //test
    //   }
    //   ///////////////////// End Mixed Track Loop ///////////////////////////////
    // }
    // //////////// END Mixed Invariant mass Event Loop ///////////////////////
    // h_mx_prim_inv_m_PHI    -> SaveAs("./result/h_mx_prim_inv_m_PHI_5mx_tot.root");

    // for( int k = 0; k < 8; k++)
    //   {
    //     if(k >= 8) break;
    //     cout<<"saving: "<<k<<endl;
    //     a_h3_mx_prim_inv_nsig[k] -> SaveAs(Form("h3_mx_prim_inv_nsig_%s.root",TS_ctr[k].Data()));
    //   }

    h_phi_TPC2_TOF3_tight_y -> SetTitle("#phi Uncorrected dN/dy");
    h_phi_TPC2_TOF3_tight_y -> SetXTitle("y");
    h_phi_TPC2_TOF3_tight_y -> SetYTitle("Counts");
    // h_mx_TPC2_TOF3_invM -> SaveAs("./result/h_mx_TPC2_TOF3_invM_out_5mx_tot.root");//tfile
    // h_mx_TPC2_TOF3_invM_y_bin1 -> SaveAs("./result/h_mx_TPC2_TOF3_invM_y_bin1_out_5mx_tot.root");
    // h_mx_TPC2_TOF3_invM_y_bin2 -> SaveAs("./result/h_mx_TPC2_TOF3_invM_y_bin2_out_5mx_tot.root");
    // h_mx_TPC2_TOF3_invM_y_bin3 -> SaveAs("./result/h_mx_TPC2_TOF3_invM_y_bin3_out_5mx_tot.root");
    // h_mx_TPC2_TOF3_invM_y_bin4 -> SaveAs("./result/h_mx_TPC2_TOF3_invM_y_bin4_out_5mx_tot.root");

    // new TCanvas();
    // h_TPC2_TOF3_invM    -> Draw();
    // h_TPC2_TOF3_invM_y_bin1    -> Draw();
    // h_TPC2_TOF3_invM_y_bin2    -> Draw();
    // h_TPC2_TOF3_invM_y_bin3    -> Draw();
    // h_TPC2_TOF3_invM_y_bin4    -> Draw();
    //
    // h_mx_TPC2_TOF3_invM -> SetFillColor(kYellow);
    // h_mx_TPC2_TOF3_invM_y_bin1 -> SetFillColor(kYellow);
    // h_mx_TPC2_TOF3_invM_y_bin2 -> SetFillColor(kYellow);
    // h_mx_TPC2_TOF3_invM_y_bin3 -> SetFillColor(kYellow);
    // h_mx_TPC2_TOF3_invM_y_bin4 -> SetFillColor(kYellow);

    tf_out->Write();
    return;
}
/////////////////////////// END Main Function //////////////////////////////////
