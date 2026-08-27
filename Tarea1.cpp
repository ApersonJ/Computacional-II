#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

const double a=0, b=2*(3.1415926), tol=0.000001;

void SolicitarDatos(int &N, int &n, double &R, double &rMin, double &rMax, double &zMin, double &zMax);

void GenerarMalla(double r[], double z[], double rMin, double rMax, double zMin, double zMax, int N);

double Funcion(double r[], double z[], double R, double phi, int k);

double Simpson(double r[], double z[], double R, int n, int k);

void Comprobar(double R, int n);

void generarArchivoSalida(const string &archivoSalida,
                          int n, double R, double N,
                          double r[], double z[], double Resultados[]);
                          
void graficarDatos(const string &archivoSalida);

int main(){
	
	double R, rMin, rMax, zMin, zMax;
	int n, N;
	SolicitarDatos(N, n, R, rMin, rMax, zMin, zMax);
	
	double r[N*N], z[N*N];
	GenerarMalla(r, z, rMin, rMax, zMin, zMax, N);
	
	double Resultados[N*N];
	for (int k=0; k<=N*N-1; k++){
    if (fabs(r[k] - R) < tol && fabs(z[k]) < tol){
      Resultados[k] = 20;
    }
    else {
		  Resultados[k]=Simpson(r, z, R, n, k);
    }
	}
	
	string archivoSalida = "../resultados/resultados_integracion.dat";
	
	generarArchivoSalida(archivoSalida, n, R, N, r, z, Resultados);
	graficarDatos(archivoSalida);
    Comprobar(R, n);
	
return 0;
}

void SolicitarDatos(int &N, int &n, double &R, double &rMin, double &rMax, double &zMin, double &zMax){
	cout<<"Introduzca la resolucion de la malla (# de puntos por eje): "; cin>>N;
	
	cout<<"Ingrese el numero de intervalos: "; cin>>n;
	
	cout<<"Ingrese el radio de circunferencia: "; cin>>R;
	
	// La regla compuesta de Simpson requiere un numero par de intervalos.
	if (n % 2 != 0)
    {
        cout << "Para Simpson, n debe ser par." << endl;
        cout << "Se usara n = " << n + 1 << "." << endl;
        n++;
    }
  
  cout<<"Ingrese un radio minimo: "; cin>>rMin;
  cout<<"Ingrese un radio Maximo: "; cin>>rMax;
  while(rMin < 0){
    cout<<"El radio minimo no puede ser negativo: "; cin>>rMin;
  }
  
  while(rMax < rMin){
    cout<<"Error, el radio maximo no puede ser menor "; cin>>rMax;
  }
  
  cout<<"Ingrese una altura minima: "; cin>>zMin;
  cout<<"Ingrese una altura Maxima: "; cin>>zMax;
  while(zMax < zMin){
    cout<<"Error, la altura maxima no puede ser menor: "; cin >> zMax;
  }
}

void GenerarMalla(double r[], double z[], double rMin, double rMax, double zMin, double zMax, int N)
{
  
  double dr = (rMax - rMin)/(N- 1);
  double dz = (zMax - zMin)/(N - 1);
  
    int k = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            r[k] = rMin + i*dr;
            z[k] = zMin + j*dz;

            k++;
        }
    }
}

double Funcion(double r[], double z[], double R, double phi, int k){
	double D2 = r[k]*r[k] + z[k]*z[k] + R*R - 2*R*r[k]*cos(phi);
		return R/sqrt(D2);
	}


double Simpson(double r[], double z[], double R, int n, int k){
	
    double dphi = (b - a) / n;

    // Extremos: f(x_0) + f(x_n)
    double suma = Funcion(r, z, R, a, k) + Funcion(r, z, R, b, k);

    // T�rminos pares interiores: 2 f(x_{2j})
    for (int j = 1; j <= (n / 2) - 1; j++)
    {
        int i = 2 * j;
        double phi_2j = a + i * dphi;

        suma += 2.0 * Funcion(r, z, R, phi_2j, k);
    }

    // T�rminos impares interiores: 4 f(x_{2j-1})
    for (int j = 1; j <= n / 2; j++)
    {
        int i = 2 * j - 1;
        double phi_2j_1 = a + i * dphi;

        suma += 4.0 * Funcion(r, z, R, phi_2j_1, k);
    }

    return (dphi / 3.0) * suma;
}


void generarArchivoSalida(const string &archivoSalida,
                          int n, double R, double N,
                          double r[], double z[], double Resultados[])
{
    ofstream archivo(archivoSalida.c_str());

    archivo << setprecision(8);

    archivo << "# Potencial electrico de una espira circular de radio "<<R<<endl;
    archivo << "# Numero de intervalos: " << n << endl;
    
    archivo <<"#los resultados son:"<<endl;
    archivo << setw(10) << "#r"
         << setw(10) << "z"
         << setw(10) << "V(r,z)" << endl;
    
    for (int i = 0; i < N*N; i++){
        archivo << r[i] << "\t"
                << z[i] << "\t"
                << Resultados[i] << endl;
    }
    archivo.close();

    cout << "Archivo generado: " << archivoSalida << endl;
}

void graficarDatos(const string &archivoSalida)
{
    // Para obtener una vista 3D de la superficie
    string script =
        "set terminal pngcairo size 800,600\n"
        "set output '../grafica/potencial_grafica.png'\n"
        "set xlabel 'r'\n"
        "set ylabel 'z'\n"
        "set zlabel 'V(r,z)'\n"
        "set title 'Potencial electrico'\n"
        "set dgrid3d 50,50\n"
        "set pm3d\n"
        "set hidden3d\n"
        "set view 60, 30, 1, 1\n"
        "splot '" + archivoSalida + "' using 1:2:3 with pm3d\n";

    ofstream scriptFile("../scripts/potencial_grafica.gnuplot");

    if (scriptFile.is_open())
    {
        scriptFile << script;
        scriptFile.close();

        system("gnuplot ../scripts/potencial_grafica.gnuplot");

        cout << "Grafico generado: potencial_grafica.png" << endl;
    }
    else
    {
        cerr << "Error: No se pudo crear el archivo de script gnuplot" << endl;
    }
}

void Comprobar(double R, int n){

    cout << fixed << setprecision(8);
    cout << "\n Validacion de Simpson\n";
    cout << setw(12) << "z"
         << setw(16) << "Numerico"
         << setw(16) << "Analitico"
         << setw(16) << "ErrorRel" << endl;

    double reje[1]={0.0};
    double zeje[]= {0.0, 0.5*R, R, 2*R, 5*R, 10*R, 100*R };

    for(int i = 0; i < 7; i++){
        double zval = zeje[i];
        
        double z[1] = {zval};
        double Vnum = Simpson(reje, z, R, n, 0);
        double Vana = b*R / sqrt(R*R + zval*zval);
        double errorRel = fabs(Vnum - Vana) / fabs(Vana);

        cout << setw(12) << zval
             << setw(16) << Vnum
             << setw(16) << Vana
             << setw(16) << errorRel << fixed << endl;
        }
}
