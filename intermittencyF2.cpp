#ifndef HISTO_PLOTTING_CPP_INCLUDED
#define HISTO_PLOTTING_CPP_INCLUDED



#endif // HISTO_PLOTTING_CPP_INCLUDED

#include "TROOT.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <algorithm>
#include "math.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TLine.h"
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
#include "TChain.h"
#include "TChainElement.h"
#include "TPaveText.h"


using namespace std;
using namespace ROOT::Math;

void intermittencyF2()
{
    TFile *data_plot = TFile::Open("Plots/multiplicity_eta.root");
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //gStyle->SetTitleAlign(33);

    //============================================= Eta ===============================================================

    TH2D *eta_multiplicity = (TH2D*)data_plot->Get("eta_multiplicity");
	TH1D *F_2 = new TH1D("F_2_vs_etabinsize", "F_{2} vs #delta#eta", 130, 0, 1.3);
	
    TCanvas *eta_canvas = new TCanvas ("eta_canvas", "multiplicity_eta_working_plots", 2);
    TLegend *legend_eta = new TLegend (0.6, 0.6, 0.8, 0.8);
	const float Neta_bin = 480;
	float f_NSum, Nf_NSum, frequency_N, N_ave, N_ave_sum, N_ave_sumsq, N_Nminus1_Sum, N_Nminus1_ave_sum, F_culmulant;
	int Nbins = 0;
		
	eta_canvas->cd();
	
	
	
		for(int data_point = 0; data_point <= 7; ++data_point){ // maximum value of 8 for data_point as I want a total of 8 points
																// per cumulant plot, with each data point having half the number 
																// of bins as the previous point, starting from 480 bins 
			//eta_multiplicity->RebinX(2^data_point);
			N_ave_sum = 0;
			N_ave_sumsq = 0;
			N_Nminus1_ave_sum = 0;
			F_culmulant = 0;
			
			for (int eta_bins = 1; eta_bins <= (480/(2^data_point)); ++eta_bins){
				
				f_NSum = 0;
				Nf_NSum = 0;
				N_Nminus1_Sum = 0;
				
				for (int multiplicity = 0; multiplicity <= 200; ++multiplicity){ //multiplicity denote the multiplicity axis of eta_multiplicity
					
					frequency_N = 0;
					frequency_N = eta_multiplicity->GetBinContent(eta_multiplicity->GetBin(eta_bins, (multiplicity+1)));
					
					Nf_NSum += multiplicity*frequency_N;								// correspond to summ_N*frequency_N
					N_Nminus1_Sum += multiplicity*(multiplicity - 1)*frequency_N;		// correspond to summ_N*(N-1)*frequency_N
					
					/*for (int i = 0; i <= data_point; ++i){					//automate the calculation of <n(n-1)...(n-p+1)>

						variable = (multiplicity - i)

					}*/
					
					f_NSum += frequency_N;												// correspond to summ_f_N
				}
				
				N_ave = 0;						
				N_ave = Nf_NSum/f_NSum;							// N_ave correspond to a particular eta bin
				
				N_Nminus1_Sum = 0;
				N_Nminus1_ave = N_Nminus1_Sum/f_NSum;			// correspond to <N(N-1)>
				
				N_ave_sum += N_ave;								// correspond to summation_<N_m>
				N_ave_sumsq += (N_ave)^2						// correspond to summation_<N_m>^2
				N_Nminus1_ave_sum += N_Nminus1_ave;
			}
			
			F_culmulant = (0.01*(2^data_point)/4.8)*(N_Nminus1_ave_sum) / (N_ave_sumsq);

			F_2->SetBinContent((0.01*(2^data_point)),F_culmulant);
		}
	gPad->SetLogX();
	
	F_2->GetXaxis()->SetTitle("#delta#eta");
	F_2->GetYaxis()->SetTitleOffset(1.3);
	F_2->GetYaxis()->SetTitle("F_{2}");
	F_2->Write();	
	
    eta_canvas->Update();
    eta_canvas->Draw();
	
	data_plot.Write();
	data_plot.Close();
}

/*data_eta_histo->GetXaxis()->SetRangeUser(-2.4, 2.4);
    data_eta_histo->SetMinimum(0);
    data_eta_histo->SetMaximum(0.04);
    //data_eta_histo->SetTitle("Ongoing Analysis #sqrt{s} = 13TeV");
    data_eta_histo->GetXaxis()->SetLabelSize(0.02);
    data_eta_histo->GetYaxis()->SetTitleOffset(2.0);
    data_eta_histo->GetYaxis()->SetTitleSize(0.02);
    data_eta_histo->GetYaxis()->SetLabelOffset(0.001);
    data_eta_histo->GetYaxis()->SetLabelSize(0.02);
    data_eta_histo->SetLineStyle(0);
    data_eta_histo->SetLineColorAlpha(kBlack,1);
    data_eta_histo->SetMarkerStyle(20);
    data_eta_histo->Draw("E1");

    legend_eta->SetFillColor(0);
    legend_eta->SetFillStyle(0);
    legend_eta->SetBorderSize(0);
    legend_eta->SetTextSize(0.02);
    legend_eta->AddEntry(data_eta_histo, "Uncorrected Data", "pe");
    legend_eta->AddEntry(herwig_eta_histo, "Herwig", "l");
    legend_eta->Draw("SAME");

    text_eta_title->SetFillColor(0);
    text_eta_title->SetFillStyle(0);
    text_eta_title->AddText("Ongoing Analysis #eqrt{s} = 13TeV");
    text_eta_title->Draw("SAME");*/
