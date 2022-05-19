#include <iostream>
#include <ctime>
#include <sstream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "TStyle.h"
#include "TH1F.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TApplication.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TNtuple.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Types.h"

//definisco parametri
#define PI 3.14159
#define sigma_s 0.3  
#define rho_s 0.5 
#define mu_s 0. 
#define sigma_f 1. 
#define rho_f 0.4 
#define mu_f 4. 

using namespace std;


//Segnale
double segnale (double x, double y){ 
	return (1/(2 * PI * sigma_s * sigma_s * sqrt(1-pow(rho_s, 2)))) * exp((-1)/(2 * (1-pow(rho_s, 2))) * (pow((x-mu_s)/sigma_s, 2) 
		+ pow((y-mu_s)/sigma_s, 2) - 2 * rho_s * (x-mu_s)/sigma_s * (y-mu_s)/sigma_s));
}

//Fondo
double fondo (double x, double y){ 
	return (1/(2 * PI * sigma_f * sigma_f * sqrt(1-pow(rho_f, 2)))) * exp((-1)/(2 * (1-pow(rho_f, 2))) * (pow((x-mu_f)/sigma_f, 2) 
		+ pow((y-mu_f)/sigma_f, 2) - 2 * rho_f * (x-mu_f)/sigma_f * (y-mu_f)/sigma_f));
}

//Genera numeri casuali distributi su un intervallo
double rand_range(double min, double max){
	return min + (max - min) * rand() / RAND_MAX;
}

//Genera numeri casuali distributi secondo una certa funzione func
tuple<double, double> rand_TAC(double xMin, double xMax, double yMin, double yMax, double fMin, double fMax, bool t){
	double x = 0., y = 0., f = 0.;
	x = rand_range(xMin, xMax);
	y = rand_range(yMin, yMax);
	f = rand_range(fMin, fMax);
	if (t == 1){ //segnale
		if (f < segnale(x,y)) return make_tuple (x,y);
		if (f > segnale(x,y)) return make_tuple (0,0);
	}
	if (t == 0){ //fondo
		if (f < fondo(x,y)) return make_tuple (x,y);
		if (f > fondo(x,y)) return make_tuple (0,0);
	}
	return make_tuple (0,0);
}

void test7(){

	srand(520); //Inizializzazione seme per la generazione MC

	int N = 500;

	double iniz = -1;
	double fin = 7; 


	TF1 * funz = new TF1("funz", "1.98 - 0.99*x", -3, 7);
	funz -> SetLineColor(kBlack);


	TCanvas * c1 = new TCanvas ("c1", "funzioni"); //Canvas per stampare le funzioni
	c1 -> Divide(1,2);
	TF2 *f_s = new TF2("f2","(1/(2 * 3.14 * 0.3  * 0.3 * sqrt(1-0.5^2))) * exp((-1)/(2 * (1-0.5^2)) * ((x/0.3)^2 + (y/0.3)^2 - 2 * 0.5 * x/0.3 * y/0.3))",-1,1,-1,1);
	TF2 *f_b = new TF2("f2","(1/(2 * 3.14 * sqrt(1-0.4^2))) * exp((-1)/(2 * (1-0.4^2)) * ((x-4)^2 + (y-4)^2 - 2 * 0.4 * (x-4) * (y-4)))",1,7,1,7);
	c1 -> cd(1);
	f_s -> Draw("surf1");
	c1 -> cd(2);
	f_b -> Draw("surf1");
	c1 -> Print("Es_7_funzioni.eps", "eps");

	TNtuple * sgl = new TNtuple ("sgl", "sgl", "x:y"); //TNtuple per il segnale che verrà usato da TMVA
	TNtuple * bkg = new TNtuple ("bkg", "bkg", "x:y"); //TNtuple per il fodno che verrà usato da TMVA

	int i = 0;
	double x_s[N], y_s[N] ; //Variabili per i grafici
	while(i < N){ //Segnale
		auto signal = rand_TAC(iniz, fin, iniz, fin, 0., 2.2, 1); //TRUE per segnale
		if(get<0>(signal) == 0 && get<1>(signal) == 0) continue;
		x_s[i] = get<0>(signal);
		y_s[i] = get<1>(signal);
		sgl -> Fill(get<0>(signal), get<1>(signal));
		i++;
	}

	i = 0;
	double x_f[N], y_f[N] ; //Variabili per i grafici
	while(i < N){ //Fondo
		auto background = rand_TAC(iniz, fin, iniz, fin, 0., 0.18, 0); //FALSE per fondo
		if(get<0>(background) == 0 && get<1>(background) == 0) continue;
		x_f[i] = get<0>(background);
		y_f[i] = get<1>(background);
		bkg -> Fill(get<0>(background), get<1>(background));
		i++;
	}

	TCanvas * c2 = new TCanvas ("c2", "grafici"); //Canvas per i grafici

	TGraph * g_s = new TGraph (N, x_s, y_s);
	g_s -> SetMarkerColor(kBlue);
	TGraph * g_f = new TGraph (N, x_f, y_f);
	g_f -> SetMarkerColor(kRed);

	TMultiGraph * mg = new TMultiGraph();

	mg -> Add(g_s, "AP");
	mg -> Add(g_f, "AP");
	mg -> Draw("A");
	funz -> Draw("same");

	c2 -> Print("discrimine.jpg", "jpg");

	TFile f("fisher.root", "RECREATE"); // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
	TMVA::Factory * factory = new TMVA::Factory("TMVAanalysis", &f, "");
	TMVA::DataLoader * dataloader = new TMVA::DataLoader ("data");

	dataloader -> AddSignalTree(sgl);
	dataloader -> AddBackgroundTree(bkg);

	dataloader -> AddVariable("x", 'F');
	dataloader -> AddVariable("y", 'F');

	factory -> BookMethod(dataloader, TMVA::Types::kLD, "LD", "");
	factory -> TrainAllMethods();
	factory -> TestAllMethods();
	factory -> EvaluateAllMethods();	
}

#ifndef __CINT__
int main () { test7(); return 0; }  // Main program when run stand-alone
#endif

