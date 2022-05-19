#include<iostream>
#include<fstream>
#include<cstring>
#include<math.h>
#include<cmath>
#include<TMath.h>
#include<TCanvas.h>
#include<TGraph.h>
#include<TH2.h>
#include<TF1.h>
#include<TApplication.h>
#include<TGraphErrors.h>
#include<TMultiGraph.h>
#include<TLegend.h>
#include<TStyle.h>
#include<random>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(){

	gStyle -> SetOptFit(1111);

	TApplication *gr = new TApplication("grafica", 0, NULL);
 	//TCanvas *c = new TCanvas("histo", "histo");

	//Range della Distribuzione di Landau
	double min = -5;
	double max = 20;
	
	//istogrammi che raccoglieranno i dati delle p-value
	TH1F * h1 = new TH1F("histo1", "", 100, 0., 1.);
	TH1F * h2 = new TH1F("histo2", "", 100, 0., 1.);
	/*TH1F * Histo1 = new TH1F("first sample", "", 1000, min, max);
	TH1F * Histo2 = new TH1F("second sample", "", 1000, min, max);*/

	//uso random_device e mt19337 per il numero casuale
	random_device rd;  

    mt19937 gen(rd()); 
	

	//estraggo numero casuale da una distribuzione uniforme tra 0-1  
    uniform_real_distribution<> dis(0., 1.);

	//Massimo della distribuzione di Landau con i parametri di Default
	double func_max = TMath::Landau(-0.22278, 0, 1);

	//vettori per il test K-S
	vector<double> vettore1;
	vector<double> vettore2;

	//Numeri casuali 
	int N = 1000000;
	srand(time(0)); 

	//ciclo for per calcolare 1000 p-value
	for (int index = 0; index < 1000; ++index) { 

		//Generatore di numeri casuali che voglio testare con Rand()/(Rand_MAX)
		for (int i=0; i < N; ++i) {
			
			//Numeri casuali tra 0 e 1
			double num = ((double) rand() / (RAND_MAX));
			
			//Tra il minimo e il massimo
			double var = min + num * (max - min);

			//Distribuzione di Landau usando la funzione contenuta in Tmath
			double dist = TMath::Landau(var, 0, 1);
			

			//Per l'Hit-miss di Von Neumann estraggo un numero casuale compreso tra 0 e il massimo della funzione
			double num2 = ((double) rand() / (RAND_MAX))*func_max;

			
			if (num2 < dist) {

				vettore1.push_back(var);
				Histo1 -> Fill(var);
			} 
		}

		//Generatore di numeri casuali di riferimento estratti da una distr.uniforme tra 0-1
		for (int j=0; j < N; ++j) {

			double num = dis(gen);
			double var = min + num * (max - min);

			double dist = TMath::Landau(var, 0, 1);
			double num2 = dis(gen)*func_max;

			if (num2 < dist) {
				vettore2.push_back(var);
				Histo2 -> Fill(var);
			} 
		}

        //Test di Kolmogorov-Smirnov (dati non binnati)
		sort(vettore1.begin(), vettore1.end());
		sort(vettore2.begin(), vettore2.end());

		int m = vettore1.size();
		int n = vettore2.size();

		double elenco1[m];
		double elenco2[n];

		for (int index1=0; index1<m; index1++) {
		    elenco1[index1] = vettore1[index1];}

		for (int index2=0; index2<n; index2++) {
			elenco2[index2] = vettore2[index2];}

		
		double ks = TMath::KolmogorovTest(m, elenco1, n, elenco2, "p");
		//cout << "KS " << test << endl;
		h1 -> Fill(ks); 

		vettore1.clear();
		vettore2.clear();

		//Test del Chi quadrato (Dati Binnati)
		double chi = Histo1 -> Chi2Test(Histo2, "p");
		//cout << "CHI" << chi << endl;
		h2 -> Fill(chi);

		Histo1 -> Reset();
		Histo2 -> Reset(); }
	

	/*TCanvas *c1 = new TCanvas("histo1", "histo1");
	Histo1 -> Draw();
	TCanvas *c2 = new TCanvas("histo2", "histo2");
	Histo2 -> Draw();*/
	TCanvas *c1 = new TCanvas("Distribuzione P-Value (KS)", "Distribuzione P-Value (KS)");
	h1->Draw();
	TCanvas *c2 = new TCanvas("Distribuzione P-Value (CHI)", "Distribuzione P-Value (CHI)");
	h2->Draw();
	
	gr -> Run();

	return 0;
}

