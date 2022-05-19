#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <random>
#include <boost/math/special_functions/bessel.hpp>
#include <TStyle.h>
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TRandom.h"
#endif

using namespace std;


double F(double y) {

	double rho=0; 
	double k=0;
	double I_y=0;
	rho=y;
	//funzione di BesselJ1 per l'Airy Pattern
	k=TMath::BesselJ1(rho);
	I_y=pow((2*k/rho),2);

	return I_y ; //restituisco valore della funzione

}

//macro per incompatibilità con il pacchetto RooUnfold
void test6(){
	
	gSystem->Load("libMathMore.so");

	gStyle -> SetOptFit(0111);
	
	TApplication * grafica = new TApplication("grafica", 0, NULL);
	TCanvas * c1 = new TCanvas("histo1", "histo1");
	
	vector <double> real ;
	vector <double> smear ;
	
	int nEventi = 50000;
	int nBin = 50;
	TH1F *h2 = new TH1F("histo1","real", nBin, -10, 10);
	
	int iEv = 0;
	random_device rd;  

    mt19937 thetaGen(rd());
	uniform_real_distribution<> dis(0, 6.28);

	mt19937 gen3(rd()+1);
	//Random numbers between 0 and 0.18 (the maximum of F)
	uniform_real_distribution<> disMax(0., 1.0);

	//Genero il primo istogramma utilizzando un meccanismo HIT MISS
	while(iEv < nEventi) {

		//Random theta 
		double y = 10*sin(dis(thetaGen));
		double function = F(y);
		double check1 = disMax(gen3);
		
		if(check1 < function) {

			h2 -> Fill(y);
			real.push_back (y);
			iEv += 1;
			}
		}
		
	//SMEARING GAUSSIANO	
	double c=0.6;
	double dx=c*((double)20/50);
	
	double k =real.size ();
	TH1F *h3 = new TH1F("histo2","smear", nBin, -10, 10);
	for(int i=0;i<k;i++){
	double z= real.at(i) + gRandom->Gaus(0,dx); //sommo il contributo gaussiano per lo smearing (uso gaus della Trandom)
	smear.push_back (z);
	h3 -> Fill(z);
	}

	//applico i metodi del pacchetto RooUnfold 
	//response mi restituisce la matrice di risposta
	RooUnfoldResponse response1 (nBin ,-10, 10);
	RooUnfoldResponse response2 (nBin , -10, 10);
	for (int i = 0; i < real.size (); i++) {
	response1 .Fill( real [i], real [i]); //matrice di risposta tra la real e la real
	response2 .Fill( smear [i], real [i]);//matrice di risposta tra la distribuzione con lo smearing e la real
	}

	RooUnfoldBayes unfold1 (& response1 , h2 ,4); //con metodo Bayesiano inverto matrice di risposta
	RooUnfoldBayes unfold2 (& response2 , h3 ,4);
	TH1D* hReco1 = (TH1D *) unfold1.Hreco (); //concludo con la ricostruzione della funzione originale
	TH1D* hReco2 = (TH1D *) unfold2.Hreco ();
	
	h2 -> Draw();
	TCanvas * c2 = new TCanvas("histo2", "histo2");
	h3 -> Draw();
	TCanvas * c3 = new TCanvas("histo3", "histo3");
	hReco1->Draw();
	TCanvas * c4 = new TCanvas("histo4", "histo4");
	hReco2->Draw();

	grafica -> Run();
}

#ifndef __CINT__
int main () { test6(); return 0; }  // Main program when run stand-alone
#endif


		
