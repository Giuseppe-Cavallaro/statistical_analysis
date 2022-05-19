#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
using namespace std;


// funzione da integrare
float f (float x){
  return exp(x);
}


//MC CRUDE
int count;
float val=0;
float integral (float (*f)(float)){
  for (count=0; count < 100; count++){
    float u1 = (float)rand()/(float)RAND_MAX;
    val += f(u1);
    }
  
  float ris=val;
  val=0;
  return ris/count;
}

//Sampling stratificato
int count1;
float val1=0;
float val2=0;
float integral1 (float (*f)(float)){
  for (count1=0; count1 < 50; count1++){
    float u1 = (0.5)*(float)rand()/(float)RAND_MAX;
    val1 += f(u1);
	}
  float ris1=val1;
  val1=0;
 
  for (count1=0; count1 < 50; count1++){
    float u2 =0.5 + (0.5)*(float)rand()/(float)RAND_MAX;
    val2 += f(u2);
    }  
  
  float ris2=val2;
  val2=0;
  
  
  return (ris1+ris2)/(2*count1);
}

// F normalizzata per importance sampling
float f1 (float x){
  return exp(x)*(pow(exp(1),2)-exp(1))/exp(x+1); //normalizzazione
}


//Importance Sampling
float integral2 (float (*f)(float)){
int count2=0;
float val3=0;
  for (count2=0; count2 < 100; count2++){
    float u1 = (float)rand()/(float)RAND_MAX;
    float u2 = log(exp(1)+pow(exp(1),2)*u1) - 1;
    
    val3 += f1(u2);
    }
  
  float ris=val3;
  val=0;
  return ris/count2;
}

//Variabili antitetiche
float integral3 (float (*f)(float)){
int count3=0;
float val4=0;
float v1[50];
  for (count3=0; count3 < 50; count3++){
    float u1 = (float)rand()/(float)RAND_MAX;
    v1[count3]=u1;
    val4 += f(u1);
    }

  for (count3=0; count3 < 50; count3++){
    float u2=1-v1[count3];
    val4 += f(u2);
    }
  
  float ris=val4;
  val=0;
  return ris/(2*count3);
}


int main(){

	int size=1000;
	float values[size];
  
  //metodo MC crude, no riduzione della varianza e no hit miss
	
	for (int i=0; i<1000; i++){
	values[i]=integral(f);}

	float sum = 0;
	for (int i = 0; i < size; i++){
            sum += values[i];
        }
	float avg=sum / size;
	float var = 0;

        for (int i = 0; i < size; i++) {
            var +=pow((values[i] - avg),2);
        }
	float sigma=0;
	sigma=var/size;
	cout << "Media del Caso MC CRUDE\t" << avg << " con varianza:\t" << sigma << endl;

  //Sampling stratificato, divido l'intervallo in due sottointervalli

	float values1[size];

	for (int i=0; i<1000; i++){
	values1[i]=integral1(f);}

	 sum = 0;
	for (int i = 0; i < size; i++) {
            sum += values1[i];
        }
  	avg=sum / size;
  	var = 0;

    for (int i = 0; i < size; i++) {
			var +=pow((values1[i] - avg),2);
        }
	sigma=0;
	sigma=var/size;
	cout << "Media con Sampling Stratificato\t" << avg << " con varianza:\t" << sigma << endl;	

  //Importance Sampling effettuo un cambio di variabile con exp(x+1)

	float values2[size];
  
	for (int i=0; i<1000; i++){
	values2[i]=integral2(f);}
  
	sum = 0;

	for (int i = 0; i < size; i++) {
		sum += values2[i];
		}
	avg=sum / size;
	var = 0;

    for (int i = 0; i < size; i++) {
        var +=pow((values2[i] - avg),2);
		}
	sigma=0;
	sigma=var/size;
	cout << "Media con Importance Sampling\t" << avg << "con varianza:\t" << sigma << endl;

 //Metodo delle variabili antitetiche uso la covarianza negativa di f(x) e f(1-x) per funzioni monotone crescenti

	float values3[size];
	for (int i=0; i<1000; i++){
	values3[i]=integral3(f);}
	
	sum = 0;

    for (int i = 0; i < size; i++) {
            sum += values3[i];
        }
	avg=sum / size;
	var = 0;

    for (int i = 0; i < size; i++) {
            var +=pow((values3[i] - avg),2);
        }
	sigma=0;
	sigma=var/size;
	cout << "Media con Variabili antitetiche\t" << avg << " con varianza:\t" << sigma << endl;

return 0;
}




