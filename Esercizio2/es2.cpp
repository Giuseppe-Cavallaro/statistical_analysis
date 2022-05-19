#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include "TGraph.h"
#include "TStyle.h"
#include "TLorentzVector.h"
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


using namespace std;

//queste funzioni servono per il calcolo di Gamma col metodo Monte Carlo (useremo questo metodo per calcolare l'Integrale)
double f (double x){
  return -log(log(pow(x,-1)));
}

double integral (double (*f)(double)){
	int count;
	double val=0;
  for (count=0; count < 100; count++){
    double u1 = (double)rand()/(double)RAND_MAX;
    val += f(u1);
    }
  
  double ris=val;
  val=0;
  return ris/count;
}

int main(){
	
	double g_t=0.577215664901532; //valore della Gamma di Eulero Mascheroni da utilizzare come riferimento
	
	//calcolo diretto della Gamma con la definizione
   	double tot=0;
	double gamma1=0;
	double error[50];
	double iter[50];
	int n=0;
	int z=0;
	
	for (int i=2; i<50; i++){ //questo ciclo serve per aumentare via via il numero di termini della serie
		for (double k=1; k<i; k++){ //applico definizione
			tot += (1/k);
			z=k;}
		gamma1=tot-log(z);
		error[i-2]=fabs(gamma1 - g_t);
		iter[i-2]=i;
		n=i;
		//cout <<"gamma con la definizione"<<gamma1<<" con errore "<<fabs(gamma1 - g_t)<<endl;
		tot=0;}
	
	//rappresento l'errore assoluto in un grafico con il numero di iterazioni
	TCanvas *c1 = new TCanvas("c1","Gamma con la definizione",200,10,500,300);
	TGraph* gr = new TGraph(n, iter, error);
   	gr->Draw("AC*");
	c1 -> Print("2_definizione.jpg", "jpg");	
	
	//Metodo di Eulero 
	double error_2[50];
	double iter_2[50];
	for (int i=2; i<50; i++){
        	for (double k=1; k<i; k++){ 
			tot += (1/k) - log(1 + (1/k)) ;
			z=k;}
		error_2[i-2]=fabs(tot - g_t);
		iter_2[i-2]=i;
		//cout <<"gamma con metodo di Eulero "<<tot<<" con errore "<<fabs(tot - g_t)<<endl;
		n=i;
		tot=0;}
		
	TCanvas *c2 = new TCanvas("c2","Gamma con metodo di Eulero",200,10,500,300);
	TGraph* gr2 = new TGraph(n, iter_2, error_2);
   	gr2->Draw("AC*");
	c2 -> Print("2_eulero.jpg", "jpg");
	
	
	//Metodo con le Z di Riemann 
	tot=0;
	double p=0;
	double error_3[50];
	double iter_3[50];
	for (int j=3; j<50; j++){
		for (double i=2; i<j; i++){ 
			for (double k=1; k<j;k++){
				p+=pow((1/k),i);}
			tot += pow(-1,i)*(1/i)*p ;
			p=0;}
		error_3[j-3]=fabs(tot - g_t);
		iter_3[j-3]=j;
		//cout <<"gamma con metodo delle Z "<<tot<<" con errore "<<fabs(tot - g_t)<<endl;
		tot=0;
		n=j;}
	
	TCanvas *c3 = new TCanvas("c3","Gamme con metodo delle Z",200,10,500,300);
	TGraph* gr3 = new TGraph(n, iter_3, error_3);
   	gr3->Draw("AC*");
	c3 -> Print("2_Riemann.jpg", "jpg");
      
	//Metodo di Gerst 
	tot=0;
	int m=0; //int per averne la parte intera di log2(k)
	double error_4[50];
	double iter_4[50];
	for (int i=2; i<50; i++){
        	for (double k=1; k<i; k++){ 
			m=log2(k);
			tot += pow(-1,k)*(1/k)*m ;}
		error_4[i-2]=fabs(tot - g_t);
		iter_4[i-2]=i;
		//cout <<"gamma con metodo di Gerst "<<tot<<" con errore "<<fabs(tot - g_t)<<endl;
		tot=0;
		n=i;}
		
	TCanvas *c4 = new TCanvas("c4","Gamma con metodo di Gerst",200,10,500,300);
	TGraph* gr4 = new TGraph(n, iter_4, error_4);
   	gr4->Draw("AC*");
	c4 -> Print("2_Gerst.jpg", "jpg");
	
	//uso un metodo Monte Carlo per calcolare l'integrale 
	int size=1000;
    	double values[size];
	double error_5[size];
	double iter_5[size];
  	for (int j=1; j<size; j++){
        	for (int i=0; i<j; i++){
			values[i]=integral(f);}
			double sum = 0;
			for (int i = 0; i < j; i++){
            			sum += values[i];}
        		double avg=sum / size;
			//cout << "Media gamma con mc:\t" << avg <<" con errore "<<fabs(avg - g_t)<<endl;
			error_5[j-1]=fabs(avg - g_t);
			iter_5[j-1]=j;
		sum=0;
		n=j;}
	
	TCanvas *c5 = new TCanvas("c5","Gamma con Monte Carlo",200,10,500,300);
	TGraph* gr5 = new TGraph(n, iter_5, error_5);
   	gr5->Draw("AC*");
	c5 -> Print("2_MC.jpg", "jpg");

	//Metodo di Havil 
	tot=0;
	p=0;
	double g=0;
	double error_6[50];
	double iter_6[50];
	double b=0;
    	for (int i=1; i<11; i++){ //per questo metodo c'è un limite, provo 10 casi con multipli di 50 iterazioni
		for (double k=1; k<50*i;k++){
			b=i*50;
			g=ceil(b/k); //estraggo la mantissa
			p+=g - (b/k);}
		tot= (1/b)*p ;
        	p=0;
		error_6[i-1]=fabs(tot - g_t);
		iter_6[i-1]=i;
		//cout <<"Metodo di Havil"<<tot<<" con errore "<<fabs(tot - g_t)<<endl;
		tot=0;}
		
	TCanvas *c6 = new TCanvas("c6","Gamma con Metodo di Havil",200,10,500,300);
	TGraph* gr6 = new TGraph(10, iter_6, error_6);
   	gr6->Draw("AC*");
	c6 -> Print("2_Havil.jpg", "jpg");
      
	return 0;}
        
