#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sys/time.h>
#include "eliminacion_gaussiana.h"
#include <unistd.h>
using namespace std;

#define temperatura_isoterma 500
#define pi 3.14159265358979323846
timeval start, end;

void init_time()
{
	gettimeofday(&start,NULL);
}

double get_time()
{
    gettimeofday(&end,NULL);
    return (1000000*(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec))/1000000.0;
}

char * appendCharToCharArray(char * array, char a)
{
	char * ret = (char*)malloc(sizeof(array) + 1 + 1);
	strcpy(ret,array);
	ret[strlen(ret)] = a;
	ret[sizeof(ret)] = '\0';
	return ret;
}


void guardar_resultado(vector< vector<double> >result,char *file){
  ofstream myfile(file);
  if (myfile.good()){
  	for(int j=0;j<result.size();j++){
	    for (int i=0;i<result[j].size();i++){
			myfile<< std::fixed<<std::setprecision(6)<<result[j][i];
			myfile<<endl;
	    }
	}
    myfile.close();
  }
  else{
  	cout<<"error al escribir el archivo de salida"<<endl;
  	exit(-1);
  }
}

void mostrar_resultado(vector<double> v){
	for (int i=0;i<v.size();i++){
		cout<<v[i];
		if(i%10!=9){
			cout<<"\t";
		}
		else{
			if((i+1<v.size())&&i>0)
				cout<<endl;
		}
	}
	cout<<endl;
}

int main(int argc, char *argv[]) {

	// argc es la cantidad de argumentos
	// argv son los argumentos. argv[0] es el nombre del programa compilado, a partir de 1 es el primer parametro y asi.

	// Abro el archivo y parseo los datos...
	string line;
	ifstream myfile(argv[1]);

	if(argc<=3){
		cout<<"Error en el ingreso de parametros"<<endl;
		cout<<"Modo de ejecucion: "<<endl;
		cout<<"./nombre <archivo de entrada> <archivo de salida> <modo> <guardar isoterma>"<<endl;
		exit(-1);
	}
	int guardar_isoterma=0;
	if(argc>4){
		//definimos un parametro extra para no hacer operaciones innecesarias
		guardar_isoterma=atoi(argv[4]);
	}
	int modo=atoi(argv[3]);
	double datos[8];
	//cout<<"Abriendo el archivo: "<<argv[1]<<endl;
	if (myfile.is_open()) {
		getline (myfile,line);
		stringstream ss(line);
		int ind = 1;
		while( ss.good() ) {
		   string substr;
		   getline( ss, substr, ' ' );
		   datos[ind] = atof(substr.c_str());
		   ind++;
		}			
		myfile.close();
	}
	else{
		cout<<"fallo al abrir"<<endl;
		exit(-1);
	}
	
	// Datos Iniciales
	double rinterno = datos[1];
	double rexterno = datos[2];
	int m = datos[3];
	int n = datos[4];
	double isoterma = datos[5];
	int ninstancias = datos[6];
//	cout<<"Parametros de entrada"<<endl;
//	cout<<"Ri: "<<rinterno<<endl;
//	cout<<"Re: "<<rexterno<<endl;
//	cout<<"M +1: "<<m<<endl;	//particion radial
//	cout<<"N: "<<n<<endl;		//particion angular
//	cout<<"Temperatura de la isoterma: "<<isoterma<<endl;
//	cout<<"# Instancias: "<<ninstancias<<endl;

	vector< vector <double> > constante;
	constante.resize(ninstancias);
	for(int i=0;i< ninstancias;i++){
		constante[i].resize(m*n,0);
	}

	ifstream myfile_b(argv[1]);
	if (myfile_b.is_open()) {
		getline (myfile_b,line);
		int inst=0;
		int numerito = 0;
		while ( getline (myfile_b,line) )	{
			stringstream ss(line);
			int j=0;
			int i=0;
			bool primeros_n=true;
			bool salir=false;
			while( ss.good() && !salir) {
				string substr;
				getline( ss, substr, ' ' );
				while(strcmp(substr.c_str(),"")==0){
					getline( ss, substr, ' ' );
				}	
				constante[inst][i]=(atof(substr.c_str()));
				i++;
				if(primeros_n&&i>=n){
					i=n*(m-1);
					primeros_n=false;
				}	
				if(i>=n*m){
					salir=true;
				}	
			}	
			inst++;
		}
		myfile_b.close();
	}
	else{
		cout<<"fallo al leer el archivo"<<endl;
		exit(-1);
	}
	/**
	Indices
	T (r,o)
	T (j-1,k)  	= i-n
	T (j+1,k) 	= i+n
	T (j,k)		= i
	T (j,k-1)  	= i+1
	T (j,k+1)	= i-1
	*/
	matriz a(n*m,n*m);

	double rActual=rinterno;
	double deltaR=(double )(rexterno-rinterno)/(m-1);
	double deltar2=pow(deltaR,2);
	double angulo=0;
	double deltaO2=(double)(pow(2*(pi/n),2));

	for(int radio=1;radio<m-1;radio++){
		rActual=rActual+deltaR;
		for(int angulo=0;angulo<n;angulo++){
		 //coeficientes originales
			int diagonal= radio*n + angulo; //conversion(radio,angulo,n);	
			a.cargar_posicion(diagonal,diagonal-n,((1/(deltar2))-(1/(deltaR*rActual))));		//posicion j-1,k
			a.cargar_posicion(diagonal,diagonal+n,1/(deltar2));									//posicion j+1,k
			a.cargar_posicion(diagonal,diagonal,((1/(deltaR*rActual)) - (2/(deltar2)) -(2/(deltaO2*pow(rActual,2)))));//valor del elemento central
			if(angulo +1 == n)
			{
				a.cargar_posicion(diagonal,diagonal-n+1,(1/(deltaO2*pow(rActual,2)))); 	//posicion j,k+1
				a.cargar_posicion(diagonal,diagonal-1,(1/(deltaO2*pow(rActual,2))));	//posicion j,k-1
			}else
				if(angulo==0){
				a.cargar_posicion(diagonal,diagonal+n-1,(1/(deltaO2*pow(rActual,2))));
				a.cargar_posicion(diagonal,diagonal+1,(1/(deltaO2*pow(rActual,2))));
			}
			else{
				a.cargar_posicion(diagonal,diagonal-1,(1/(deltaO2*pow(rActual,2))));
				a.cargar_posicion(diagonal,diagonal+1,(1/(deltaO2*pow(rActual,2))));
			}
		}
	}
	for(int i=0;i<n;i++){
		a.cargar_posicion(i,i,1);
		a.cargar_posicion(m*n-i-1,m*n-i-1,1);
	}
	//a.guardar_matlab("matlab.txt");
	vector< vector <double> > resultados;
	resultados.resize(ninstancias);
	double tiempo_total =0;
	double tiempo_i=0;
	clock_t ti;
	clock_t tf;
	ti=clock();
	time_t timer;
	time_t timer2;
	double seconds;
	time(&timer);
   	if( modo == 0 ) {
    	// Resuelvo por Eliminacion Gaussiana
    	cout<<"Resolviendo por EG"<<endl;
		for(int i=0;i<ninstancias;i++){
			//cout<<"Instancia: "<<i+1<<endl;
			//init_time();
			resultados[i]=a.eliminacion_gausiana(constante[i]);
			//cout<<"Resultado:"<<endl;
			tiempo_i =get_time();
			tiempo_total +=tiempo_i;
			//mostrar_resultado(constante[i]);
		}	
	}
	else if ( modo == 1 ) {
		cout<<"Resolviendo por LU"<<endl;
		double tiempo_total =0;
		double tiempo_i=0;
		a.calcular_LU();
		for(int i=0;i<ninstancias;i++){
			//init_time();
			//cout<<"Istancia: "<<i+1<<endl;	
			resultados[i]=a.resolver_LU(constante[i]);
			//cout<<"Resultado:"<<endl;
			tiempo_i =get_time();
			tiempo_total +=tiempo_i;
		}
	}
	else{
		cout<<"Modo de resolucion Invalido"<<endl;
		cout<<"Elegir 0 para resolver por EG"<<endl;
		cout<<"Elegir 1 para resolver por LU"<<endl;
	}
	time(&timer2);
	seconds = difftime(timer2,timer);
	tf=clock();
	double milisec=(double)((clock()-ti)%CLOCKS_PER_SEC)/CLOCKS_PER_SEC;
	//esto hay que revisar por el tema de la precision, con uno no se puede tener en ms si el tiempo es mayor a 40 min
	if(guardar_isoterma==0){
		//si no esta deshabilitado el calculo de isotermas
		string nombre(argv[1]);
 		size_t found= nombre.find_last_of("/\\");
 		size_t fin 	= nombre.find_last_of(".");
  		string directorio=nombre.substr(0,found+1);
  		string arch=nombre.substr(found+1,fin-found);
  		cout<<arch<<endl;
  		directorio+=string("Tiempo ")+arch+string("txt");
		//ofstream myfileTiempos(directorio.c_str(),std::ofstream::app);
		ofstream myfileTiempos(directorio.c_str());
		if(!myfileTiempos.is_open()){
			cout<<"Error al abrir el archivo para guardar tiempos"<<endl;
			exit(-1);
		}
		//cout<<"Planteando sistema"<<endl;
		myfileTiempos<<"Archivo: "<<argv[1]<<" Modo: "<<modo<<endl;
		myfileTiempos<<"Tiempo total: "<<seconds<<" Tiempo Resol(ms): "<<(double)(tf-ti)/CLOCKS_PER_SEC<<endl;
		myfileTiempos<<"Radio interno: "<<rinterno<<" Radio Externo: "<<rexterno<<" Temp Isoterma: "<<isoterma<<endl;
		myfileTiempos<<"M: "<<m<<" N: "<<n<<" Instancias: "<<ninstancias<<endl;	
		//Calculo de isotermas para cada instancia
		for(int i=0;i<ninstancias;i++){
			myfileTiempos<<"Posicion de la isoterma en funcion de tita en la instancia "<<i+1<<endl;
			double radio_iso=0;
			double dist_max_iso=0;
			for(int angulo=0;angulo<n;angulo++){
				bool encontrado=false;
				int radio=0;
				while((radio<m)&&(!encontrado)){
					if(resultados[i][radio*n+angulo]<isoterma){
						encontrado=true;
					}
					else{
						radio++;
					}
				}
				if(radio==m){//si no estoy en los extremos
					radio_iso=rexterno;
				}
				else if(radio==0){
					radio_iso=rinterno;
				}
				else{
					double ta=resultados[i][(radio-1)*n+angulo];
					double tb =resultados[i][radio*n+angulo];
					if(son_iguales(ta-tb,0)){//y esta es una isoterma, en este intervalo la temperatura es la misma
						radio_iso=rinterno+deltaR*radio;
					}
					else{//calculo la interpolacion lineal entre esos dos puntos
						radio_iso=rinterno+deltaR*(radio-((double) isoterma-ta)/(tb-ta));
					}
					
				}
				dist_max_iso=max(dist_max_iso,radio_iso);
				myfileTiempos<<"Tita "<<angulo<<": "<<std::fixed<<std::setprecision(3)<<radio_iso<<endl;
			}
			myfileTiempos<<"Distancia maxima de la Isoterma "<<i+1<<": "<<std::fixed<<std::setprecision(3)<<dist_max_iso<<endl<<endl;

		}
		myfileTiempos.close();
	}
	if(argc>2){
		//cout<<"Guardando"<<endl;
		guardar_resultado(resultados,argv[2]);
	}	
	return 0;
}
