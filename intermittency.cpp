#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TPad.h"
#include "math.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TLegend.h"
#include "TText.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TAttAxis.h"
#include "TString.h"

using namespace ROOT::Math;

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace std;

#ifdef __MAKECINT__
    #pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >++;
#endif

vector<TString>* getListOfFiles(TString strfiles){

  vector<TString>* vfiles = new vector<TString>;

  if(strfiles.Contains(".root")){
    TChain chain("tree/tree","");
    chain.Add(strfiles);
    TObjArray* fileElements=chain.GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl=0;
    while (( chEl=(TChainElement*)next() )) {
      vfiles->push_back(TString(chEl->GetTitle()));
    }
  }
  else if(strfiles.Contains(".txt")){
    ifstream txtfile;
    txtfile.open(strfiles);
    if(!txtfile) {
      cout<<"Unable to read the txt file where the rootfiles are." << endl ;
      cout << strfiles << " doesn't exist." << endl << "Aborting ...";
      exit(0);
    }
    string filename;
    while(txtfile>>filename && filename!="EOF")
      vfiles->push_back(TString(filename));
    txtfile.close();
  }
  else {
    cout << "Unknown type of input to get files. Must contain either .root or .txt extension." << endl << "Aborting ..." << endl;
    exit(0);
  }
  cout << "[getListOfFiles] Will run on " << vfiles->size() << " files" << endl;
  return vfiles;
}

void intermittency()
{
	TH1::SetDefaultSumw2(1);
	TH2::SetDefaultSumw2(1);
	
    //===========================implement cuts================================
    const double eta_cut = 2.4;
    const double vtx_number_cut = 1.0;
    const float vtxxysize = 0.2;
    const float vtxzsize = 10;
    const double pt_cut = 0.5;
    const int lumi_cut = 90;
    const float ptErr_pt_cut = 0.05;
    const float dz_dzErr_cut = 3;
    const float d0_d0Err_cut = 3;
    const float dof_cut = 4;
	const float fb_eta_upper = 0.8;
	const float fb_eta_lower = 0.2;

    //==============================================Variables==============================================================

    double fdata_evt = 0;
    double fdata_trkphi = 0;
    double fdata_trketa = 0;
    double fdata_trkd0 = 0;
    double fdata_trkpt = 0;
    double fdata_trkdz = 0;
    double fdata_trkdpt = 0;
    double fdata_dz_sigmadz, fdata_d0_sigmad0, fdata_sigmapt_pt, fdata_dz_sigmadzcalc, fdata_d0_sigmad0calc, fdata_sigmad0, fdata_sigmad0run1, fdata_d0, fdata_d0_sigmad0run1, fdata_trkdx, fdata_trkdy;
    double fdata_wx, fdata_wy,fdata_wz, fdata_vtxxysize, fdata_vtxzsize, fdata_sigmad0calc, fdata_dz, fdata_sigmadz, fdata_dz_sigmadzrun1, fdata_d0calc, fdata_dzcalc, fdata_d0calc_sigmad0, fdata_dzcalc_sigmadz;
    double fdata_dzBS_dzErr, fdata_dzvtxBS_dzErr, fdata_dxyBS_d0Err, fdata_dxyvtxBS_d0Err;
	double fdata_numberofvtxxBS, fdata_vtxxBSvalue, fdata_vtxxBSlower, fdata_vtxxBSupper;
    double fdata_numberofvtxyBS, fdata_vtxyBSvalue, fdata_vtxyBSlower, fdata_vtxyBSupper;
    double fdata_numberofvtxzBS, fdata_vtxzBSvalue, fdata_vtxzBSlower, fdata_vtxzBSupper;
    double fdata_vtxzminusvtxz, fdata_multiplicity, fdata_forward_multiplicity, fdata_backward_multiplicity;
    double fdata_multiplicity_norm = 0;
    double fdata_sqvtxx = 0;
    double fdata_sqvtxxnumber = 0;
    double fdata_sqvtxy = 0;
    double fdata_sqvtxynumber = 0;
    double fdata_sqvtxz = 0;
    double fdata_sqvtxznumber = 0;
    double fdata_numberoftrkdx = 0;
    double fdata_numberoftrkdy = 0;
    double fdata_numselectedvtxz = 0;
    double fdata_numvtxzminusvtxz = 0;
    double fdata_trkvalidhits = 0;
    double fdata_trkchi2n = 0;
	float fsum_forward_multiplicity = 0;
	float fsum_backward_multiplicity = 0;
	float fsum_fb_evts = 0;
    int ndata_numberofvtxx, ndata_numberofvtxy, ndata_numberofvtxz, ndata_numberofvtxxBS, ndata_numberofvtxyBS, ndata_numberofvtxzBS;
    int ndata_totaltrk;
	float eta;


//============================================== Histos for pT, eta, phi ==============================================================

    TH1D *data_pt_histo = new TH1D ("data pT", "Data p_{T}", 200, 0, 10);
    const int eta_bins = 480;
	TH1D *data_eta_histo = new TH1D("data eta", "Data #eta", 500, -2.5, 2.5);
    TH1D *data_fulleta_histo = new TH1D("data_fulleta", "Full #eta", 500, -2.5, 2.5);
	TH1D *data_phi_histo = new TH1D ("data phi", "Data #phi", 160, -4, 4);

//================================================== Histos for Multiplicity =========================================================================

	TH1D *data_multiplicity = new TH1D("Multiplicity", "Multiplicity", 200, 0, 200);
	TH2D *multiplicity_eta = new TH2D("multiplicity_eta", "multiplicity_eta", 200, 0, 200, 480, -2.4, 2.4);
	TH2D *eta_multiplicity = new TH2D("eta_multiplicity", "eta_multiplicity", 480, -2.4, 2.4, 200, 0, 200);
	
//===========================retrieve ROOT file for loop over multiple ROOT files that requires CMSSW============================

	vector<TString> *datafiles = new vector<TString>();
    cout << "Getting list of files..." << endl;

	//FBlist_of_files.txt contains only 10 files to loop over to test code.
    datafiles = getListOfFiles("/afs/cern.ch/user/h/hang/Multiplicity/Forward_Backward/FBlist_of_files.txt");
    cout << "File list stored" << endl;

        TFile *datafile;
        TTree *datatree;

	for(vector<TString>::iterator itlistdatafiles = datafiles->begin() ; itlistdatafiles != datafiles->end(); ++itlistdatafiles)
	{
        cout << "\nOpening new file " << *itlistdatafiles << endl;

        datafile = new TFile(*itlistdatafiles, "READ");

        cout << "Opened " << *itlistdatafiles << endl;
        datatree = (TTree*)datafile->Get("tree/tree");
	
		cout << "Reading ROOT File..." << endl;
		
		//===========================define variables to read TTree================

		vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *data_tracks = 0;
		datatree->SetBranchAddress("trP4", &data_tracks);

		int ndata_run = 0;
		datatree->SetBranchAddress("Run", &ndata_run);

		int ndata_lumi;
		datatree->SetBranchAddress("Lumi", &ndata_lumi);

		// int ndata_zerobias;
		// datatree->SetBranchAddress("trgZeroBias", &ndata_zerobias);

		vector<double> *fvecdata_vtxx = 0;
		datatree->SetBranchAddress("vtxx", &fvecdata_vtxx);

		vector<double> *fvecdata_vtxy = 0;
		datatree->SetBranchAddress("vtxy", &fvecdata_vtxy);

		vector<double> *fvecdata_vtxz = 0;
		datatree->SetBranchAddress("vtxz", &fvecdata_vtxz);

		Double_t fvecdata_BSx = 0;
		datatree->SetBranchAddress("BSx", &fvecdata_BSx);

		Double_t fvecdata_BSy = 0;
		datatree->SetBranchAddress("BSy", &fvecdata_BSy);

		Double_t fvecdata_BSz = 0;
		datatree->SetBranchAddress("BSz", &fvecdata_BSz);

		vector<double> *fvecdata_vtxxBS = 0;
		datatree->SetBranchAddress("vtxxBS", &fvecdata_vtxxBS);

		vector<double> *fvecdata_vtxyBS = 0;
		datatree->SetBranchAddress("vtxyBS", &fvecdata_vtxyBS);

		vector<double> *fvecdata_vtxzBS = 0;
		datatree->SetBranchAddress("vtxzBS", &fvecdata_vtxzBS);

		vector<int> *nvecdata_highpurity = 0;
		datatree->SetBranchAddress("highPurity", &nvecdata_highpurity);

		vector<int> *nvecdata_vtxndof = 0;
		datatree->SetBranchAddress("vtxndof", &nvecdata_vtxndof);

		vector<double> *fvecdata_dzvtxBS = 0;
		datatree->SetBranchAddress("dzvtxBS", &fvecdata_dzvtxBS);

		vector<double> *fvecdata_dxyvtxBS = 0;
		datatree->SetBranchAddress("dxyvtxBS", &fvecdata_dxyvtxBS);

		vector<double> *fvecdata_dzerr = 0;
		datatree->SetBranchAddress("dzerr", &fvecdata_dzerr);

		vector<double> *fvecdata_d0err = 0;
		datatree->SetBranchAddress("d0err", &fvecdata_d0err);

		vector<double> *fvecdata_pterr = 0;
		datatree->SetBranchAddress("ptErr", &fvecdata_pterr);

		vector<double> *fvecdata_vtxxerr = 0;
		datatree->SetBranchAddress("vtxxErr", &fvecdata_vtxxerr);

		vector<double> *fvecdata_vtxyerr = 0;
		datatree->SetBranchAddress("vtxyErr", &fvecdata_vtxyerr);

		vector<double> *fvecdata_vtxzerr = 0;
		datatree->SetBranchAddress("vtxzErr", &fvecdata_vtxzerr);

		vector<double> *fvecdata_vtxxErrBS = 0;
		datatree->SetBranchAddress("vtxxErrBS", &fvecdata_vtxxErrBS);

		vector<double> *fvecdata_vtxyErrBS = 0;
		datatree->SetBranchAddress("vtxyErrBS", &fvecdata_vtxyErrBS);

		vector<double> *fvecdata_vtxzErrBS = 0;
		datatree->SetBranchAddress("vtxzErrBS", &fvecdata_vtxzErrBS);

		vector<int> *nvecdata_validhits = 0;
		datatree->SetBranchAddress("nValidHits", &nvecdata_validhits);

		vector<double> *fvecdata_trackschi2n = 0;
		datatree->SetBranchAddress("chi2n", &fvecdata_trackschi2n);

		//triggers not needed for FB correlations
		vector<int> *triggerZB = 0;
		datatree->SetBranchAddress("trigger", &triggerZB);

		vector<int> *triggerMB = 0;
		datatree->SetBranchAddress("triggerMB", &triggerMB);

		Int_t ndata_totalEvt = (Int_t)datatree->GetEntries();
		cout << "There is a total of " << ndata_totalEvt << " events." << endl;

//========================================================= Start of Evt Loop ================================================================
		float eta_bincontent = 0;
		for (Int_t i = 0; i < ndata_totalEvt; ++i)
		{
			if (i%10000 == 0)
			{	
				cout << "Looping through event..." << endl;
			}
			datatree->GetEntry(i);

			if
			(
				//good lumisections for runs
				((ndata_lumi >= 1 && ndata_lumi <=103) && ndata_run == 254987) ||
				(((ndata_lumi >= 1 && ndata_lumi <=120) || (ndata_lumi >= 153 && ndata_lumi <=214) ) && ndata_run == 254989) ||
				((ndata_lumi >= 1 && ndata_lumi <=28) && ndata_run == 254993) ||
				((ndata_lumi >= 4 && ndata_lumi <=414) && ndata_run == 255019) ||
				(((ndata_lumi >= 1 && ndata_lumi <=34) || (ndata_lumi >= 36 && ndata_lumi <=209) || (ndata_lumi >= 306 && ndata_lumi <=325) || (ndata_lumi >= 327 && ndata_lumi <=343) ) && ndata_run == 255029) ||
				(((ndata_lumi >= 31 && ndata_lumi <=101) || (ndata_lumi >= 125 && ndata_lumi <=231) || (ndata_lumi >= 233 && ndata_lumi <=493) || (ndata_lumi >= 586 && ndata_lumi <=1054) || (ndata_lumi >= 1096 && ndata_lumi <=1199)) && ndata_run == 255031)  //255031 lumicuts
			)
			{
				//if (triggerMB[i] == 1)
				{
					ndata_totaltrk = data_tracks->size();
					++fdata_evt;

					int vtxdof = 0;

					ndata_numberofvtxx = fvecdata_vtxx->size();
					ndata_numberofvtxy = fvecdata_vtxy->size();
					ndata_numberofvtxz = fvecdata_vtxz->size();
					fdata_backward_multiplicity = 0;
					fdata_forward_multiplicity = 0;
					
	//========================================================= Start of Vertex Loop ================================================================

					for (int vtxnumber = 0; vtxnumber != ndata_numberofvtxx; ++vtxnumber)
					{
						//if (fdata_vtxxysize <= vtxxysize && fdata_vtxzsize <= vtxzsize)				
						if ((ndata_numberofvtxz == vtx_number_cut) && ((*nvecdata_vtxndof)[vtxdof] > dof_cut))
						{

	//========================================================= Start of Trk Loop ================================================================

							for (int t = 0; t != ndata_totaltrk; ++t)
							{
								XYZTVector data_vec = (*data_tracks)[t];
								
								if ((*nvecdata_highpurity)[t] == 1) //only high purity tracks are selected
								{								
									data_eta_histo->Fill(data_vec.Eta());
									data_fulleta_histo->Fill(data_vec.Eta());
								}
									
							}
							
	//========================================================= End of Trk Loop ================================================================
						}					
						
					}
					
	//========================================================= End of Vertex Loop ================================================================
					
				}
				//write the loop to extract data_eta_histo bin content here.
				
				eta = 0;
				for(int netabin = 0; netabin <= eta_bins; ++netabin)
				{
					eta = netabin/100.00 - 2.40;
					multiplicity_eta->Fill((data_eta_histo->GetBinContent(netabin+1)), eta);
					eta_multiplicity->Fill(eta,(data_eta_histo->GetBinContent(netabin+1)));
				}
				
				data_eta_histo->Reset();
			}
	//========================================================= End of Lumisection Loop ================================================================
		}
//========================================================= End of Evt Loop ================================================================
	}
//========================================================= End of File Loop ===========================================================================
	
	cout << "Creating root file to write to..." << endl;
    TFile data_plot("/afs/cern.ch/user/h/hang/Multiplicity/Intermittency/Plots/multiplicity_eta.root", "RECREATE");
	cout << "Created root file to write to..." << endl;
	//gPad->SetLogy();
	
	/*fb_profile->GetXaxis()->SetTitle("N_{f}");
	fb_profile->GetYaxis()->SetTitle("<N_{b}>_{N_{f}}");
	fb_profile->GetYaxis()->SetTitleOffset(1.3);
	fb_profile->SetErrorOption("i");
	fb_profile->Fit("pol1", "", "SAME", 0, 50); //y-intercept = 0.0455878, b_{corr} = 0.807957
	fb_profile->Write();*/
	
	data_fulleta_histo->GetXaxis()->SetTitle("#eta");
	data_fulleta_histo->GetYaxis()->SetTitleOffset(1.3);
	data_fulleta_histo->GetYaxis()->SetTitle("Frequency");
	data_fulleta_histo->Write();

    multiplicity_eta->GetXaxis()->SetTitle("N_{ch}");
    multiplicity_eta->GetYaxis()->SetTitleOffset(1.3);
    multiplicity_eta->GetYaxis()->SetTitle("#eta");
    multiplicity_eta->Write();
	
	eta_multiplicity->GetXaxis()->SetTitle("#eta");
	eta_multiplicity->GetYaxis()->SetTitleOffset(1.3);
	eta_multiplicity->GetYaxis()->SetTitle("N_{ch}");
	eta_multiplicity->Write();

    data_plot.Write();
	data_plot.Close();
}
