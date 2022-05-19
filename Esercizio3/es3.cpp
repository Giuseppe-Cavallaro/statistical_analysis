#include<iostream>
#include<fstream>
#include<cstring>
#include<math.h>
#include<cmath>
#include<TFitter.h>
#include<TLatex.h>
#include<TStyle.h>
#include<random>
#include<algorithm>
#include<TRandom3.h>
#include<TGraph2D.h>
#include<TMinuit.h>
#include<TFumili.h>
#include<TMath.h>
#include<TCanvas.h>
#include<TGraph.h>
#include<TH2.h>
#include<TF1.h>
#include<TF2.h>
#include<TApplication.h>
#include<TGraphErrors.h>
#include<TMultiGraph.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDatah3t.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChi2Var.h"
#include "RooMinimizer.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAbsReal.h"
using namespace std;

//non ho potuto utilizzare la funzione per problemi di incompatibilità della macro col compilatore (la allego comunque)
/*double F(double j, double k) {

	double alpha = 0.65;
	double beta = 0.06;
	double gamma = -0.18;

	double firstParam = 0.5 * (1-alpha);
	double secondParam = 0.5 * (3*alpha - 1) * pow(cos(j), 2);
	double thirdParam = beta * pow(sin(j), 2) * cos(2*k);
	double fourthParam = sqrt(2) * gamma * sin(2*j) * cos(k);

	return (3/(4*3.14)) * (firstParam + secondParam - thirdParam - fourthParam);
}*/
//per il pacchetto Roofit sonos stato costretto a convertire il programma in una macro di ROOT
void bis(){

	int eventi = 2000; //basso numero di eventi per evitare un p-value di zero
	TH1F * p1 = new TH1F("p-value", "", 100, 0., 1.);

	int nBin = 50;
	double min = 0;
	double max = 6.28;
	double alpha = 0.65;
	double beta = 0.06;
	double gamma = -0.18;

	//per il fit utilizzo Roofit, andando a scrivere le variabili apposite
	RooRealVar x("x","x",min,max) ;
	RooRealVar y("y","y",min,max) ;
	RooRealVar a("a","a",0.65,-1,1) ;
	RooRealVar b("b","b",0.06,-1,1) ;
	RooRealVar c("c","c",-0.18,-1,1) ;
	
	RooGenericPdf gp("gp","Generic PDF","(3/(4*3.14))*(0.5*(1-a)+0.5*(3*a-1)*pow(cos(x),2)-b*pow(sin(x),2)*cos(2*y)-sqrt(2)*c*sin(2*x)*cos(y))",RooArgSet(x,y,a,b,c));
	RooDataSet* data = gp.generate(RooArgSet(x,y),50000); //50000 eventi

	gp.fitTo(*data); //fit col metodo della maximum likelihood, i dati sono unbinned
	
	RooDatah3t* binnedData = data->binnedClone() ;
	gp.chi2FitTo(*binnedData); //fit col metodo del Chi Quadro (i dati sono binnati)
	
	//generazione dei numeri casuali seguendo la funzione indicata
	random_device rd;  

    mt19937 thetaGen(rd());
	mt19937 phiGen(rd()+1);

	//Numeri casuali tra 0 e 2pi
    uniform_real_distribution<> dis(0., 6.28);

	mt19937 yGen(rd()+2);
	//Numeri casuali tra 0 e il max della funzione
	uniform_real_distribution<> disMax(0., 0.18);
	
	TH2F *h2 = new TH2F("histo2","Distribuzione Angolare F", nBin, min, max, nBin, min, max);

	//Genero i dati della distribuzione con decadimento isotropo
	mt19937 thetaIs_gen(rd()+3);
	mt19937 phiIs_gen(rd()+4);

	TH2F *h3 = new TH2F("histo3","Distribuzione Angolare Istropa", nBin, min, max, nBin, min, max);
	
	//ciclo per il test di ipotesi	
	for (int i = 0; i < 10; ++i) {

		h2 -> Reset();
		h3 -> Reset(); 

		int j = 0;

		//primo istogramma (funzione F)
		while(j < eventi) {

			//genero theta and phi
			double theta = dis(thetaGen);
			double phi = dis(phiGen);

			
			double function = (3/(4*3.14))*(0.5*(1-alpha)+0.5*(3*alpha-1)*pow(cos(theta),2)-beta*pow(sin(theta),2)*cos(2*phi)-sqrt(2)*gamma*sin(2*theta)*cos(phi));
			double check = disMax(yGen);
			
			//metodo hit-miss
			if(check < function) {

				h2 -> Fill(theta, phi);
				j += 1;
			}
		}	

 		//Generate the second h3togram (isotropic decay)
		for(int k = 0; k < eventi; ++k) {

			double thetaIs = dis(thetaIs_gen);
			double phiIs = dis(phiIs_gen);

			h3->Fill(thetaIs, phiIs);
		}
		
		//istogramma dei p-value
		double chiResult = h2->Chi2Test(h3, "p");
		//cout << "Test number: " << i << " -> p-value: " << chiResult << endl;
		p1 -> Fill(chiResult);

	}

		TCanvas *c1 = new TCanvas("histo1", "histo1");
		h2 -> Draw("");
		TCanvas *c2 = new TCanvas("histo2", "histo2");
		h3 -> Draw("");
		TCanvas *c3 = new TCanvas("histo3", "histo3");
		p1 -> Draw("");

	
}

#ifndef __CINT__
int main () { bis(); return 0; }  //programma eseguito dalla macro
#endif

