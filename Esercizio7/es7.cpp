#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
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
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Types.h"
#include "TMVA/MethodCategory.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"

using namespace std;

//segnale
double segnale (double x, double y){
	
	double sigmaS=0.3;
	double rhoS=0.5;
	double muS=0;
	
	return (1/(2*3.14*sigmaS*sigmaS*sqrt(1-pow(rhoS,2))))*exp((-1)/(2*(1-pow(rhoS,2)))*(pow((x-muS)/sigmaS,2)+pow((y-muS)/sigmaS,2)-2*rhoS*(x-muS)/sigmaS*(y-muS)/sigmaS));
}

//fondo
double fondo (double x, double y){ 

	double sigmaF=1.0;
	double rhoF=0.4;
	double muF=4.0;

	return (1/(2*3.14*sigmaF*sigmaF*sqrt(1-pow(rhoF,2))))*exp((-1)/(2*(1-pow(rhoF,2)))*(pow((x-muF)/sigmaF,2)+pow((y-muF)/sigmaF,2)-2*rhoF*(x-muF)/sigmaF*(y-muF)/sigmaF));
}

//generatore numeri casuali
double rand_range(double min, double max){
	return min+(max-min)*rand()/RAND_MAX;
}

//generatore di coppie di punti
tuple<double, double> rand_TAC(double xMin, double xMax, double yMin, double yMax, double fMin, double fMax, bool t){
	
	double x = 0.0;
	double y = 0.0;
	double	f = 0.0;
	x = rand_range(xMin, xMax);
	y = rand_range(yMin, yMax);
	f = rand_range(fMin, fMax);
	
	if (t == 1){ //segnale
		if (f < segnale(x,y)){
		return make_tuple (x,y);}
		if (f > segnale(x,y)){
		return make_tuple (0,0);}
	}
	if (t == 0){ //fondo
		if (f < fondo(x,y)){ 
		return make_tuple (x,y);}
		if (f > fondo(x,y)){
		return make_tuple (0,0);}
	}
	return make_tuple (0,0);
}

void test7(){

	srand(NULL); 

	int N = 50000;

	double iniz = -10;
	double fin = 10; 

	TF2 *f_s = new TF2("f1","(1/(2 * 3.14 * 0.3  * 0.3 * sqrt(1-0.5^2))) * exp((-1)/(2 * (1-0.5^2)) * ((x/0.3)^2 + (y/0.3)^2 - 2 * 0.5 * x/0.3 * y/0.3))",-5,5,-5,5);
	TF2 *f_b = new TF2("f2","(1/(2 * 3.14 * sqrt(1-0.4^2))) * exp((-1)/(2 * (1-0.4^2)) * ((x-4)^2 + (y-4)^2 - 2 * 0.4 * (x-4) * (y-4)))",0,10,0,10);
	
	TCanvas *c1 = new TCanvas ("c1", "funzione1");
	f_s -> Draw("surf");
	TCanvas *c2 = new TCanvas ("c2", "funzione2");
	f_b -> Draw("surf");
	c1 -> Print("Es_7_funzioni1.jpg", "jpg");
	c2 -> Print("Es_7_funzioni2.jpg", "jpg");

	//tnuple per segnale e fondo 
	TNtuple *s = new TNtuple ("segnale", "segnale", "x:y");
	TNtuple *b = new TNtuple ("fondo", "fondo", "x:y"); 

	int i = 0;
	//vettori per le coppie del segnale
	double x_s[N]; 
	double y_s[N]; 
	
	//genero le coppie per il segnale
	while(i < N){ 
	
		auto signal = rand_TAC(iniz, fin, iniz, fin, 0., 2.2, 1); 
		if(get<0>(signal) == 0 && get<1>(signal) == 0) continue;
		x_s[i] = get<0>(signal);
		y_s[i] = get<1>(signal);
		s -> Fill(get<0>(signal), get<1>(signal));
		i++;
	}

	i = 0;
	//vettori per le coppie del fondo
	double x_f[N]; 
	double y_f[N]; 
	
	//genero le coppie per il fondo
	while(i < N){ 
	
		auto background = rand_TAC(iniz, fin, iniz, fin, 0., 0.18, 0); //FALSE per fondo
		if(get<0>(background) == 0 && get<1>(background) == 0) continue;
		x_f[i] = get<0>(background);
		y_f[i] = get<1>(background);
		b -> Fill(get<0>(background), get<1>(background));
		i++;
	}
	
	//rappresento i punti del segnale e del fondo, separati dalla retta
	TCanvas * c3 = new TCanvas ("c3", "grafici"); 

	TGraph * g_s = new TGraph (N, x_s, y_s);
	g_s -> SetMarkerColor(kRed);
	TGraph * g_f = new TGraph (N, x_f, y_f);
	
	TMultiGraph * mg = new TMultiGraph();
	TF1 *retta = new TF1("retta", "2 - x", -5, 10); //funzione scritta dopo aver fatto il test con TMVA

	mg -> Add(g_s, "AP");
	mg -> Add(g_f, "AP");
	mg -> Draw("A");
	retta -> Draw("same");
	c3 -> Print("Es_7.jpg", "jpg");

	//TMVA USO IL METODO LD
	TFile f("tmva.root", "RECREATE"); // file .root con tutti i dati che posso aprire con la GUI
	TMVA::Factory *factory = new TMVA::Factory("TMVAanalysis", &f, "");
	TMVA::DataLoader *dataloader = new TMVA::DataLoader ("data"); //dati che verranno poi valutati col metodo LD

	dataloader -> AddSignalTree(s);
	dataloader -> AddBackgroundTree(b);
	dataloader -> AddVariable("x", 'F');
	dataloader -> AddVariable("y", 'F');

	factory -> BookMethod(dataloader, TMVA::Types::kLD, "LD", ""); //analisi dei dati col metodo LD
	
	//altri metodi che ho testato
	//factory -> BookMethod(dataloader, TMVA::Types::kFisher, "Fisher", "");
	//factory -> BookMethod(dataloader, TMVA::Types::kCFMlpANN, "CF_ANN", "");
	//factory -> BookMethod(dataloader, TMVA::Types::kLikelihood, "Likelihood", "");
	factory -> TrainAllMethods(); 
	factory -> TestAllMethods();
	factory -> EvaluateAllMethods();
	
	//serve per aprire la gui di TMVA
	/*if (not gROOT->IsBatch()) {
	TMVA::TMVAGui("/home/heisenberg/tmva.root");
	}
	
	//serve per estrarre la funzione ROC 
	auto roc = factory->GetROCCurve(dataloader);
	TCanvas * c4 = new TCanvas ("c4", "roc"); //Canvas per stampare le funzioni
	roc->Draw();*/
}

#ifndef __CINT__
int main () { test7(); return 0; }  // Main program when run stand-alone
#endif

