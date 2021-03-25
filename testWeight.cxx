#include <numeric>      // std::iota
#include <vector>
#include <algorithm>
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TH1.h"
#include "TCut.h"
#include "TChain.h"
#include "TApplication.h"
#include "math.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TLatex.h"

#include "TH2.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TVirtualFitter.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>

#include <stdio.h>
#include <stdlib.h>
using namespace std;

//------PARAMETERS-----//
//const int N_Nu = 5;  // number of Nu bins
const double E_max = 2.0;
const double E_min = 0.5;
const double limit_xf = 0.1;
const int nshift_E = 99;
const double step_E = 1.0/1000.0;
const int nbins = 100;

int Nu_bin, N_Nu, N_Q2;      //number of Nu bin
TString Nuclei_Type;
Double_t Q2_min;
Double_t Q2_max;
Double_t Q2_bin;



//-------MODIFIED FEYNMAN X--------//
double Calculate_Modified_Xf(Float_t Shift, Float_t Nu, Float_t P, Float_t Pt, Float_t Q2, Float_t W, Float_t Zh)
{
	//cout << "Shift " << Shift << endl;
  	//double xf = ((Nu + 0.9385)*(TMath::Sqrt(P*P-Pt*Pt) - TMath::Sqrt(Q2+Nu*Nu)*Zh*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0  .1395*0.1395*W*W)/2./W));
  	double xf = ((Nu + 0.9385)*(TMath::Sqrt((P+Shift)*(P+Shift)-(P+Shift)*(Pt/P)*(P+Shift)*(Pt/P))-TMath::Sqrt(Q2+Nu*Nu)*(Zh+Shift/Nu)*Nu/(Nu+0.9385))/W)/((TMath::Sqrt(TMath::Power(W*W-0.9392*0.9392+0.1395*0.1395,2)-4.*0.1395*0.1395*W*W)/2./W));
  	//cout<< "Xf " << xf << endl;
  	return xf;
}


int main(int argc, char *argv[]){

	Nuclei_Type = (TString) argv[1];
  	N_Nu = atoi(argv[2]);
  	int n = atoi(argv[3]);  // ORDER OF THE CHEBYSHEV FUNC

	const Double_t Q2_down = 1.;
	const Double_t Q2_up = 4.;
	const Int_t N_Q2= 6.;//6

	//double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
	//double Nu_max = Nu_min + (4.2-3.2)/N_Nu;


	cout << "The Nuclei type studied is " << Nuclei_Type << endl;
	//cout << "Nu interval studied : " << Nu_min << " - " << Nu_max << endl;
	ostringstream nubin_label;
	//nubin_label << Nu_min  <<" < #nu < " << Nu_max << " GeV ";
	cout<< "The cut on Xf is " << limit_xf << endl;
	//------Opening data files-----//
	//TFile *file = new TFile(Form("/Users/gbibim/Documents/data/" + Nuclei_Type + "_data.root"));
	TFile *file = new TFile(Form("/user/h/hamilton/ThesisProj/data/" + Nuclei_Type + "_data.root"));

	//TFile* fout = new TFile(Form("OUTPUT/Xfmod_Ehist_"+Nuclei_Type+"_nubin%d_.root",Nu_bin),"RECREATE");

	//-----Opening TTree----//
	TTree* tree = (TTree*)file->Get("ntuple_data");
	//Reading Branches with appropiate variables. 
	Float_t TargType;
	Float_t Q2;
	Float_t Nu;
	Float_t Xb;
	Float_t W;
	Float_t SectorEl;
	Float_t ThetaPQ;
	Float_t PhiPQ;
	Float_t Zh;
	Float_t Pt;
	Float_t W2p;
	Float_t Xf;
	Float_t T;
	Float_t P;
	Float_t T4;
	Float_t deltaZ;
	Float_t NmbPion;
	tree->SetBranchAddress("TargType",&TargType);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("Nu",&Nu);
	tree->SetBranchAddress("Xb",&Xb);
	tree->SetBranchAddress("W",&W);
	tree->SetBranchAddress("SectorEl",&SectorEl);
	tree->SetBranchAddress("ThetaPQ",&ThetaPQ);
	tree->SetBranchAddress("PhiPQ",&PhiPQ);
	tree->SetBranchAddress("Zh",&Zh);
	tree->SetBranchAddress("Pt",&Pt);
	tree->SetBranchAddress("W2p",&W2p);
	tree->SetBranchAddress("Xf",&Xf);
	tree->SetBranchAddress("T",&T);
	tree->SetBranchAddress("P",&P);
	tree->SetBranchAddress("T4",&T4);
	tree->SetBranchAddress("deltaZ",&deltaZ);
	tree->SetBranchAddress("NmbPion",&NmbPion);

	Int_t nentries = tree->GetEntries();
	//Int_t nentries = 100000;

	//-----Creating output file-----//	
	TFile *fout = new TFile(Form("ROOT_"+Nuclei_Type+"_xFmod_Ecut%d_%dbins_cheb%d_Ebins%d.root", int(E_max), N_Nu, n, nbins), "RECREATE");

	//-----Creating the Graphs for the Eloss Shift Values------//
	TGraph *gElossKS = new TGraph();  //  Graph for Eloss values for the KS test
	TGraph *gElossWKS = new TGraph(); //  Graph for Eloss values fot the Weighted KS test
	TGraph *gElossKSb = new TGraph(); //  Graph for Eloss values for the Binned KS test
	TGraph *gElossWKSb = new TGraph(); // Graph for Eloss values for the Binned Weighted KS test


	for(Nu_bin = 0; Nu_bin < N_Nu; Nu_bin++){

		double Nu_min = 3.2 + Nu_bin*((4.2-3.2)/N_Nu); 
		double Nu_max = Nu_min + (4.2-3.2)/N_Nu;

		cout << "Nu interval studied : " << Nu_min << " - " << Nu_max << endl;

		//--------Starting Histos and Graphs-------//

		TGraph *gpKS = new TGraph();      //Graph for the Unbinned KS Test
		TGraph *gpWKS = new TGraph();     //Graph for the Unbinned Weighted KS Test
		TGraph *gpKSb = new TGraph();     // Graph for the Binned KS Test
		TGraph *gpWKSb = new TGraph();    // Graph for the Binned Weighted KS Test

		//TGraph *gpCSb = new TGraph();     // Graph for the ChiSquared Test

		Double_t parD[6][n+1];
		Double_t parS[6][n+1];


		TF1 * funcD = (TF1*) gROOT->GetFunction(Form("chebyshev%d", n));
		funcD->SetRange(E_min,E_max);
		TF1 * funcS = (TF1*) gROOT->GetFunction(Form("chebyshev%d", n));
		funcS->SetRange(E_min,E_max);

		TFile *acc = new TFile(Form("fout_"+Nuclei_Type+"_nubin%d.root", Nu_bin));

		//-----Histograms with Energy distribution-----//
		
		TH1F *D = new TH1F("D","D",nbins,E_min,E_max);
		D->Sumw2();

		std::map<int,TH1F*> histograms;
		for(int i = 0; i<= nshift_E; i++){
			histograms[i] = new TH1F(Form("Nuclei_bin%d",i),Form("Nuclei_bin%d",i),nbins,E_min,E_max);
			histograms[i]->Sumw2();
		}

		TH1F *DW = new TH1F("D","D",nbins,E_min,E_max);
		DW->Sumw2();

		std::map<int,TH1F*> histogramsW;
		for(int i = 0; i<= nshift_E; i++){
			histogramsW[i] = new TH1F(Form("Nuclei_bin%d",i),Form("Nuclei_bin%d",i),nbins,E_min,E_max);
			histogramsW[i]->Sumw2();
		}

		for (int Q2_bin = 0; Q2_bin < N_Q2; ++Q2_bin){

			//-------Fitting Deuterium:--------//
			cout << "Fitting Deuterium" << endl;

			//TVirtualFitter::SetMaxIterations(100000);

			TH1F *DAcc = (TH1F*)acc->Get(Form("acceptance_histoD%d%d", Nu_bin, Q2_bin));
			for (int i = 0; i <=n; ++i) funcD->SetParameter(i,1); 

			gStyle->SetOptFit(1);
			DAcc->Fit(funcD, "R");
			DAcc->Draw();
		
			fout->cd();
			DAcc->SetName(Form("fitD_"+Nuclei_Type+"_nubin%d_Q2bin%d", Nu_bin, Q2_bin));
			DAcc->Write();

			funcD->GetParameters(&parD[Q2_bin][0]);
			//cout << parD[Q2_bin][0] << "  ,  " << parD[Q2_bin][1] << endl;

		}

		//--------LOOP OVER THE SHIFTS---------//

		cout << "Starting loop over Shifts" << endl;

		for (int i = 0; i <= nshift_E; ++i){ 

			for (int Q2_bin = 0; Q2_bin < N_Q2; ++Q2_bin){

	   			cout << "Fitting" << endl;
	   			//Fitting the Solid Target
    
    			TH1F *h = (TH1F*)acc->Get(Form("acceptance_histo"+Nuclei_Type+"_shift%d_%d%d", i, Nu_bin, Q2_bin));
				for (int i = 0; i <=n; ++i) funcS->SetParameter(i,1); 

	    		h->Fit(funcS, "R");

				gStyle->SetOptFit(1);

				h->Draw();
			
				fout->cd();
				h->SetName(Form("fit_"+Nuclei_Type+"_nubin%d_Q2bin%d_shift%d", Nu_bin, Q2_bin, i));
				h->Write();

	    		funcS->GetParameters(&parS[Q2_bin][0]);	    		
	    	}  //END LOOP OVER Q2

		   	std::cout<<"Entering the loop over entries "<< nentries <<std::endl;
		  	vector<double> dataS;
		   	vector<double> dataD;
		   	vector<double> weightD;
		   	vector<double> weightS;

     		double energy_shift = step_E*i;      // is doing the right shift

		    //-------Loop over the ENTRIES: Filling the histos and vectors--------//
		    for(Int_t j=1; j<= nentries; j++){
		     	if(j%5000000==0) std::cout<< "Processing event number " << j/1000000.0 << "M "<< std::endl;
		      	tree->GetEntry(j);
		      	//Apply Cuts bin in Nu
		      	if(Nu > Nu_max || Nu < Nu_min) continue; 
		      	//Apply cuts for Deuterium, fill with Energy = Zh*Nu
		    	int K;
		    	if ( Q2>1. && Q2<1.5 ) K = 0;
		    	if ( Q2>1.5 && Q2<2. ) K = 1;
		    	if ( Q2>2. && Q2<2.5 ) K = 2;
		    	if ( Q2>2.5 && Q2<3. ) K = 3;
		    	if ( Q2>3. && Q2<3.5 ) K = 4;
		    	if ( Q2>3.5 && Q2<4. ) K = 5;

		      	if(TargType==1 && Xf> limit_xf && Zh*Nu < E_max && Zh*Nu > E_min){
		      		dataD.push_back(Zh*Nu);    //for the unbinned KS
		        	funcD->SetParameters(parD[K]);
		        	double w = funcD->Eval(Zh*Nu, 0, 0);
		      		weightD.push_back(1./w);
		      		DW->Fill(Zh*Nu, 1./w);
		        	D->Fill(Zh*Nu);

		      	}

		      	if(TargType==2 && Zh*Nu + energy_shift < E_max && Zh*Nu + energy_shift > E_min){ 
		        	double Xf_Nuclei = Calculate_Modified_Xf(energy_shift, Nu, P, Pt,  Q2,  W, Zh);
		        	if(Xf_Nuclei>limit_xf){
		          		dataS.push_back((Zh*Nu)+energy_shift);    // Storing Energy value in vector
		          		histograms[i]->Fill(Zh*Nu+energy_shift);  // Filling Energy Histogram
		          		funcS->SetParameters(parS[K]);            //
		        		double w = funcS->Eval(Zh*Nu+energy_shift, 0, 0);  //weight for E value
		      			weightS.push_back(1./w);
		      			histogramsW[i]->Fill(Zh*Nu+energy_shift, 1./w);
		        	}                                                
		      	}
	    	}

	    	std::cout<< "END LOOP OVER ENTRIES" << std::endl;
	    	//--------Sorting vectors and filling energy Graphs-------//
	    	int nD = dataD.size();
	    	int nS = dataS.size();
 			vector<int> idxD(nD);
 			vector<int> idxS(nS);

 			int x=0, y=0;
 			std::iota(idxD.begin(), idxD.end(),x++); //Initializing
 			std::iota(idxS.begin(), idxS.end(),x++); //Initializing

 			sort(idxD.begin(), idxD.end(), [&](int i,int j){return dataD[i]<dataD[j];} );
 			sort(idxS.begin(), idxS.end(), [&](int i,int j){return dataS[i]<dataS[j];} );

	    	sort(dataD.begin(), dataD.end());
	    	sort(dataS.begin(), dataS.end());

	    	//double* dataD_array = new double[nD];
	    	//double* dataS_array = new double[nS];

	    	//double* weightD_array = new double[nD];
	    	//double* weightS_array = new double[nS];

	   		//copy(dataD.begin(), dataD.end(), dataD_array);
		    //copy(dataS.begin(), dataS.end(), dataS_array);

		    //copy(weightD.begin(), weightD.end(), weightD_array);
		    //copy(weightS.begin(), weightS.end(), weightS_array);

		    //-----Unbinned Kolmogorov Test-----//
		    double pKS = TMath::KolmogorovTest(nD, &dataD[0], nS, &dataS[0], "D");
		    //gpKS->SetPoint(i, i, -1*TMath::Log10(pKS));
		    gpKS->SetPoint(i, i, pKS);


		    //WEIGHTED KOLMOGOROV TEST

		    int j1=0, j2=0;
		    Double_t d=0;
		    Double_t w1sum = std::accumulate(weightD.begin(), weightD.end(), 0);
		    Double_t w2sum = std::accumulate(weightS.begin(), weightS.end(), 0);
			Double_t j1w=0., j2w=0., fn1=0., fn2=0.;
			/*
		    while (j1<nD && j2<nS)
		    {
		    	Double_t d1 = dataD_array[j1];
		    	Double_t d2 = dataS_array[j2];
		    	int a = idxD[j1];
		    	int b = idxS[j2];
		    	Double_t w1 = weightD_array[a]/w1sum;
		    	//cout << d1 << "   " << idxD[j1] << "   " << weightD[j1] << "   " << weightD[idxD[j1]] << endl;
		    	Double_t w2 = weightS_array[b]/w2sum;
		    	
		    	cout << nD << "   " << w1sum << endl;

		    	if (d1<=d2){
            		j1+=1;
            		j1w+=w1;
            		fn1=(j1w);
            	}
		        if (d2<=d1){
		            j2+=1;
		            j2w+=w2;
		            fn2=(j2w);
		        }
		        if (TMath::Abs(fn2-fn1)>d){
		            d=TMath::Abs(fn2-fn1);
		        }
		    }
*/

		    //ROOT TEST
	        Double_t rdiff = 0;
		    Bool_t ok = kFALSE;
		    //Double_t sa  = 1./nD;
    		//Double_t sb  = 1./nS;
		    for (int l = 0; l < nD+nS; ++l){
		    	int a = idxD[j1];
		    	int b = idxS[j2];
		    	Double_t w1 = (weightD[a])/w1sum;
		    	//cout << d1 << "   " << idxD[j1] << "   " << weightD[j1] << "   " << weightD[idxD[j1]] << endl;
		    	Double_t w2 = (weightS[b])/w2sum;
		    	//cout << w1 << "  " << sa << endl;
				if (dataD[j1] < dataS[j2]) {
				    rdiff -= w1;
				    j1++;
				    if (j1 >= nD) {ok = kTRUE; break;}
				} else if (dataD[j1] > dataS[j2]) {
				    rdiff += w2;
				    j2++;
				    if (j2 >= nS) {ok = kTRUE; break;}
				} else {
				    // special cases for the ties
				    double x = dataD[j1];
				    while(j1 < nD && dataD[j1] == x) {
				        rdiff -= w1;
				        j1++;
				    }
				    while(j2 < nS && dataS[j2] == x) {
				        rdiff += w2;
				        j2++;
				    }
				    if (j1 >= nD) {ok = kTRUE; break;}
				    if (j2 >= nS) {ok = kTRUE; break;}
				}
				d = TMath::Max(d,TMath::Abs(rdiff));
			}
		    
		    //R__ASSERT(ok);
		    Double_t z = d * TMath::Sqrt(nD*nS/(nD+nS));
       		double pWKS = TMath::KolmogorovProb(z);
       		//gpKS->SetPoint(i, i, -1*TMath::Log10(pKS));
		    gpWKS->SetPoint(i, i, pWKS);


		    //-----Binned Chi2 and Kolmogorov-Smirnov Test-----//
		    D->Scale(1.0/D->Integral());
		    histograms[i]->Scale(1.0/histograms[i]->Integral());
		    //double pCSbinned = D->Chi2Test(histograms[i], "NORM");
		    double pKSbinned = D->KolmogorovTest(histograms[i], "D");
		    //gpCSb->SetPoint(i, i, -1*TMath::Log10(pCSbinned));
		    //gpKSb->SetPoint(i, i, -1*TMath::Log10(pKSbinned));
		    //gpCSb->SetPoint(i, i, pCSbinned);
		    gpKSb->SetPoint(i, i, pKSbinned);

		    //-----Binned Weighted KS Test-----//
		    DW->Scale(1.0/DW->Integral());
		    histogramsW[i]->Scale(1.0/histogramsW[i]->Integral());
		    double pWKSbinned = DW->KolmogorovTest(histogramsW[i], "D");
		    //gWpKSb->SetPoint(i, i, -1*TMath::Log10(WpKSbinned));
		    gpWKSb->SetPoint(i, i, pWKSbinned);	
	  	}//-------End of the Loop over shifts------//

	  	std::cout<< "END LOOP OVER SHIFTS" << std::endl;

		//--------Statistical Tests Graphs--------//
		gpKS->SetName(Form("pKS_"+Nuclei_Type+"_nubin%d", Nu_bin));
		gpWKS->SetName(Form("pWKS_"+Nuclei_Type+"_nubin%d", Nu_bin));
		gpKSb->SetName(Form("pKSb_"+Nuclei_Type+"_nubin%d", Nu_bin));
		gpWKSb->SetName(Form("pWKSb_"+Nuclei_Type+"_nubin%d", Nu_bin));


		gpKS->SetLineColor(4);   
		gpWKS->SetLineColor(1);
		gpKSb->SetLineColor(6); 
		gpWKSb->SetLineColor(2);

		gpKS->SetLineWidth(3);
		gpWKS->SetLineWidth(3);
		gpKSb->SetLineWidth(3);
		gpWKSb->SetLineWidth(3);

		gpKS->SetTitle("Unbinned KS");
		gpWKS->SetTitle("Weighted Unbinned KS");
		gpKSb->SetTitle("Binned KS");
		gpWKSb->SetTitle("Weighted Binned KS"); //Unbinned

		TCanvas *canvas = new TCanvas();
		TMultiGraph *multi = new TMultiGraph();
		multi->Add(gpKS, "L");
		multi->Add(gpWKS,"L");
		multi->Add(gpKSb,"L");
		multi->Add(gpWKSb,"L");

		multi->Draw("AL");
		multi->SetMaximum(1.2);
		multi->GetYaxis()->SetNdivisions(2);
		multi->GetXaxis()->SetTitle("dE [MeV]");
		multi->GetYaxis()->SetTitle("p_{0}"); //"-Log(p_{0})"
		
		canvas->BuildLegend();
		canvas->SaveAs(Form("ROOTW_Xfmod_"+Nuclei_Type+"_nubin%d_%d_%dEcut%d_cheb%d_Ebins%d.pdf", Nu_bin, nentries, int(E_min), int(E_max*100), n, nbins));

	
		//-----ELOSS HISTOGRAMS PLOTS-----//
		
		
		Double_t elossKS=0;
		Double_t elossWKS=0;
		Double_t elossKSb=0;
		Double_t elossWKSb=0;
		Int_t i_KS=0;
		Int_t i_WKS=0;
		Int_t i_KSb=0;
		Int_t i_WKSb=0;

		for (int i = 0; i < nshift_E; ++i){
			Double_t x, y;

			gpKS->GetPoint(i, x, y);
			if(y>elossKS){
				elossKS=y;
				i_KS = i;
			}

			gpWKS->GetPoint(i, x, y);
			if(y>elossWKS){
				elossWKS=y;
				i_WKS = i;
			}

			gpKSb->GetPoint(i, x, y);
			if(y>elossKSb){
				elossKSb=y;
				i_KSb = i;
			}

			gpWKSb->GetPoint(i, x, y); 
			if(y>elossWKSb){
				elossWKSb=y;
				i_WKSb = i;
			}
		}
		cout << "ELOSS VALUE FOR KS: " << i_KS << "   PROB: " << elossKS << endl;
		cout << "ELOSS VALUE FOR WKS: " << i_WKS << "   PROB: " << elossWKS << endl;
		cout << "ELOSS VALUE FOR KSb: " << i_KSb << "   PROB: " << elossKSb << endl;
		cout << "ELOSS VALUE FOR WKS: " << i_WKSb << "   PROB: " << elossWKSb << endl;


		//-----Energy spectra distributions-----//
		/*
		TCanvas *c1 = new TCanvas();
		c1->SetTitle("Energy Distribution");
	    histograms[i_KS]->Scale(1.0/histograms[i_KS]->Integral());
	    histograms[0]->Scale(1.0/histograms[0]->Integral());
	    D->Scale(1.0/D->Integral());
	       
	    histograms[i_KS]->SetLineColor(6);
	    histograms[i_KS]->SetMarkerColor(6);
	    histograms[i_KS]->SetMarkerSize(1);
	    histograms[i_KS]->SetMarkerStyle(2);
	    histograms[i_KS]->SetLineStyle(1);
	    histograms[i_KS]->SetStats(0);

	    histograms[0]->SetLineColor(2);
	    histograms[0]->SetMarkerColor(2);
	    histograms[0]->SetMarkerSize(1);
	    histograms[0]->SetMarkerStyle(1);
	    histograms[0]->SetStats(0);

	    D->SetStats(0);
	    D->SetName("Energy Distribution");
	    D->Draw("Ehist");
	    histograms[i_KS]->Draw("Esame");
	    histograms[0]->Draw("Ehistsame");

	    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
	  	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	  	legend->AddEntry(D,"Deuterium","l");
	   	legend->AddEntry(histograms[0],"Carbon","l");
	   	legend->AddEntry(histograms[i_KS],Form("Carbon shift %d", i_KS),"l");
	   	legend->Draw();
	    c1->BuildLegend();
		c1->SaveAs(Form("ROOTW_PART%d_Spectrum_"+Nuclei_Type+"_nubin%d_%d_Ecut_%d.pdf", F, Nu_bin, nentries, int(E_max)));
*/
		fout->cd();

		histograms[i_KS]->SetName(Form("KS_%d", Nu_bin));
		histograms[i_KS]->Write();

		histogramsW[i_WKS]->SetName(Form("WKS_%d", Nu_bin));
		histogramsW[i_WKS]->Write();

		histograms[i_KSb]->SetName(Form("KSb_%d", Nu_bin));
		histograms[i_KSb]->Write();

		histogramsW[i_WKSb]->SetName(Form("WKSb_%d", Nu_bin));
		histogramsW[i_WKSb]->Write();

		D->SetName(Form("D_%d", Nu_bin));
		D->Write();

		DW->SetName(Form("Weighted_D_%d", Nu_bin));
		DW->Write();

	    //------- ELOSS GRAPHS -------//
	    gElossKS->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_KS);
	    gElossWKS->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_WKS);
	    gElossKSb->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_KSb);
	    gElossWKSb->SetPoint(Nu_bin, (Nu_min+Nu_max)/2, i_WKSb);
	}//  END OF LOOP OVER NU BINS


	fout->cd();
	
	//  ELOSS PLOTS  //
	gElossKS->SetName(Form("ElossKS_"+Nuclei_Type));
	gElossWKS->SetName(Form("ElossWKS_"+Nuclei_Type));
	gElossKSb->SetName(Form("ElossKSb_"+Nuclei_Type));
	gElossWKSb->SetName(Form("ElossWKS_"+Nuclei_Type));

	gElossKS->Write();
	gElossWKS->Write();
	gElossKSb->Write();
	gElossWKSb->Write();

	gElossKS->SetMarkerColor(4);    
	gElossWKS->SetMarkerColor(1);
	gElossKSb->SetMarkerColor(6);   
	gElossWKSb->SetMarkerColor(2);  

	gElossKS->SetMarkerStyle(20);
	gElossWKS->SetMarkerStyle(21);
	gElossKSb->SetMarkerStyle(22);
	gElossWKSb->SetMarkerStyle(23);

	gElossKS->SetTitle("Unbinned KS");
	gElossWKS->SetTitle("Weighted Unbinned KS");
	gElossKSb->SetTitle("Binned KS");
	gElossWKSb->SetTitle("Weighted Binned KS");

	TCanvas *canvas = new TCanvas();
	TMultiGraph *multi = new TMultiGraph();
	multi->Add(gElossKS, "");
	multi->Add(gElossWKS,"");
	multi->Add(gElossKSb,"");
	multi->Add(gElossWKSb, "");

	multi->Draw("AP");
	multi->Write();
	multi->SetMaximum(100);
	multi->GetYaxis()->SetNdivisions(2);
	multi->GetXaxis()->SetTitle("Nu [GeV]");
	multi->GetYaxis()->SetTitle("dE [MeV]"); 

	canvas->BuildLegend();
	canvas->SaveAs(Form("ROOTW_Eloss_Xfmod_"+Nuclei_Type+"_nubin%d_%d_Ecut%d_cheb%d_Ebins%d.pdf", Nu_bin, nentries, int(E_max), n, nbins));

	std::cout<<" ABOUT TO CLOSE " << std::endl;
	fout->Close();
	std::cout<< " BYE BYE " << std::endl;
	return 1;
}
