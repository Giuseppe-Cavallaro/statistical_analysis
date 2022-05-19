#include <iostream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <TStyle.h>
#include "TH1F.h"
#include "TH2F.h"
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
#include <cstdlib>


using namespace std;

int main() {

	gStyle->SetOptStat(0000); 
	TApplication * grafica = new TApplication("grafica", 0, NULL);
	srand(time(NULL));

	int nBin = 100;
		
	//modulo momento tra 0 e 10 Gev/c e angolo tra 0 e 2pi
	double min_p = 0.;  
	double max_p = 10.; 
	double min_teta = 0; 
	double max_teta = 6.28; 

	//vettore per contenere pt e pl, somma servirà per calcolare la media di pt
	double pt[5000]; 
	double pl[5000]; 
	double sum[nBin]; 
	int counter[nBin] = {0}; //contatore per la media

	for(int i = 0; i < 5000; i++){
		
		//generazione momento e angolo
		double p =  10*(double)rand()/(double) RAND_MAX; 
		double teta = 6.28*(double)rand()/(double) RAND_MAX; 
		pt[i] = abs(p * sin(teta)); 
		pl[i] = abs(p * cos(teta)); 

		//media del momento trasverso
		for(int k = 0; k < nBin; k++){
			if(pl[i] > 10/(double)nBin*k && pl[i] <= 10/(double)nBin*(k+1)){
				sum[k] += pt[i];
				counter[k]++;
				k == nBin; 
			}
		}
	}
	
	TH1F * histo = new TH1F("histo", "", nBin, 0, 10);
	//Riempimento Istogramma con valori pesati
	for (int j = 0; j < nBin; j++){
		histo -> Fill(10/(double)nBin * j, sum[j] / counter[j]); 
	}

	TCanvas * c = necounter TCanvas("histo", "histo");
	histo -> GetXaxis() -> SetTitle("p_{L}");
	histo -> GetYaxis() -> SetTitle("<p_{T}>");
	histo -> Draw();
	grafica -> Run();
	return 0;
}
