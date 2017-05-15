#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <vector>    
#include <iostream>
#include <fstream>


using namespace std;

bool son_iguales(double a, double b, double tolerancia = 1.0e-10){
  return abs(a - b) < tolerancia;
}

class matriz{
    public:
    
    matriz();
    matriz(int a,int b);
    //funciones normales
    void cargar_a_mano();
    void const mostrar_matriz();
    void cargar_posicion(int fila,int columna,double valor); 
    void guardar_matriz(char *file);
    //para resolver matrices completas
    vector <double> eliminacion_gausiana(vector <double> independiente) ; 
    vector <double> resolver_triangular_superior(vector <double> independiente);
    void calcular_LU() ; 
    vector <double> resolver_LU(vector <double> independiente);
    void guardar_matlab(char *file);

    //variables
    int filas;
    int columnas;
    
    vector< vector < double > > posicion;
  private:
    bool transformadaLU;
    vector< vector < double > > LU;
};

matriz::matriz(){
  filas=0;
  columnas=0;
  transformadaLU=false;
}

matriz::matriz(int fil, int col){
  filas=fil;
  columnas=col;
  posicion.resize(filas); 
  transformadaLU=false;
  //redimensiono para que cada fila tenga m columnas
  for (int i = 0; i < filas; i++){
    posicion[i].resize(columnas);
  }
}

void matriz::cargar_a_mano(){
  cout<<"Especificar tamanio Fila: "<<endl;
  cin>>filas;	//filas
  cout<<"Especificar tamanio Columna: "<<endl;
  transformadaLU=false;
  cin>>columnas;	//columnas
  //redimensiono para que tenga n filas
  posicion.resize(filas);  //puede llegar a ser hasta O(cant_filas)
  //redimensiono para que cada fila tenga m columnas
  for (int i = 0; i< filas; i++){
    posicion[i].resize(columnas);
  }
  for(int i = 0; i < filas; i++)
    for(int j = 0; j< columnas; j++){
        cout<<"cargar elemento ["<<i<<"] ["<<j<<"]"<<endl;
        cin>>posicion[i][j];
    }
  
}
void matriz::cargar_posicion(int fila,int columna,double valor){
  posicion[fila][columna]=valor;
  transformadaLU=false;
}

void const matriz::mostrar_matriz(){
    for(int i = 0; i < filas; i++){
        for(int j = 0; j < columnas; j++){
            cout<<posicion[i][j];
            if(j+1<columnas)
                cout<<"\t";
        }
        cout<<endl;
    }
}

void matriz::guardar_matriz(char *file){
  ofstream myfile(file);
  if (myfile.good()){
    for (int i=0;i < filas;i++){
      for(int j=0;j < columnas;j++){
      	double valor=posicion[i][j];
      	myfile<<valor;
      	if(j<columnas-1)
      	  myfile<<"\t";
      }
      myfile<<endl;
    }
    myfile.close();
  }
}

void matriz::guardar_matlab(char *file){
  ofstream myfile(file);
  if (myfile.good()){
    for (int i=0;i < filas;i++){
      for(int j=0;j < columnas;j++){
        double valor=posicion[i][j];
        myfile<<valor;
        if(j<columnas-1)
          myfile<<" ";
      }
      myfile<<endl;
    }
    myfile.close();
  }
}

/************************************************************************************************
 *************************************************************************************************
 Seccion de codigo de eliminacion gausiana
 ************************************************************************************************
/************************************************************************************************/
vector <double> matriz::resolver_triangular_superior(vector <double> independiente){
  vector <double> res;
  res.resize(independiente.size());
  //voy restando la columna a cada fila
  for(int i=columnas-1;0<=i;i--){
      //recorro la fila y le resto el valor conocido
    	for(int fila=0;fila<i;fila++){
        double coeficiente=LU[fila][i]/LU[i][i];
        independiente[fila]+=-coeficiente*independiente[i];
      }
      res[i]=independiente[i]/LU[i][i];
  }
  //devuelvo el sistema triangular resuelto
  return res;
}

vector <double> matriz::eliminacion_gausiana(vector <double> independiente){
    //voy diagonalizando la matriz
  LU = posicion;
  transformadaLU=false;//Porque el espacio reservado es utilizado como auxiliar para la EG
  for(int diagonal=0;diagonal<columnas;diagonal++){
    //si el valor de la diagonal es cero paso a la siguiente columna
    if(! son_iguales(LU[diagonal][diagonal], 0) ){
      //si no calculo el coeficiente aij/aii para todo i<j<n y resto a la fila j - coef * fila i
      for(int fila=diagonal+1;fila<filas;fila++){
        double coeficiente=LU[fila][diagonal]/LU[diagonal][diagonal];
        if(! son_iguales(LU[fila][diagonal], 0) ){
          
          for(int col=diagonal+1;col<columnas;col++){
            LU[fila][col]+= -coeficiente*LU[diagonal][col];
          }
          independiente[fila]+= -coeficiente*independiente[diagonal];
        }
        LU[fila][diagonal]=0;
      }    
    }
  }
  return resolver_triangular_superior(independiente);
}

void matriz::calcular_LU(){
//voy diagonalizando la matriz desde la columna 0 a n
  if(!transformadaLU){
    LU=posicion;
    for(int diagonal=0;diagonal<columnas;diagonal++){
        //si el valor de la diagonal es cero paso a la siguiente columna
        if(! son_iguales(LU[diagonal][diagonal], 0) ){
            //si guardo los coeficientes en la columna debajo de la diagonal y voy guardando la U arriba
            for(int fila=diagonal+1;fila<filas;fila++){
              LU[fila][diagonal]=LU[fila][diagonal]/LU[diagonal][diagonal];
              if(! son_iguales(LU[fila][diagonal], 0) ){
                for(int col=diagonal+1;col<columnas;col++){
                  LU[fila][col]+= -LU[fila][diagonal]*LU[diagonal][col];
                }
              }
              else{
                LU[fila][diagonal]=0;
              }  
            }
        }
    }
    transformadaLU=true;
  }
}

vector <double> matriz::resolver_LU(vector <double> independiente){
  if(!transformadaLU){
    cout<<"Te falto calcular la LU antes"<<endl;
    calcular_LU();
    //Ax=LUx=b
  }
  //calculo la inversa de Ly=banda
  //multiplico por L-1 a indep para tener  el vector y
  for(int fila=0;fila<filas;fila++){
    for(int columna=0;columna<fila;columna++){
      independiente[fila]+= -LU[fila][columna]*independiente[columna];
    }
  }
  //resuelvo Ux=y
  return resolver_triangular_superior(independiente);

}
