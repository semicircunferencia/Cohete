/*SIMULACIÓN DE UN COHETE CON RUNGE-KUTTA 4

MÉTODO DE RUNGE-KUTTA
Tenemos un conjunto de ecuaciones diferenciales ypunto=f(y(t);t), donde y es un vector P-dimensional, y
f es una función de R_(P+1) a R^P. El método de Runge-Kutta consiste de los siguientes pasos:

1) Dar condición inicial y(0)

2) Evaluar: 
k_1j=h*fj(y;t) para j=1,...,P (la y de dentro es un vector)
k_2j=h*fj(y+k_1/2;t) para j=1,...,P (y+k_1/2 es un vector)
k_3j=h*fj(y+k_2/2;t) para j=1,...,P (igual)
k_4j=h*fj(y+k_3;t) para j=1,...,P (igual)

3) y(t+h)=y+(k_1+2*k_2+2*k_3+k_4)/6 (suma de vectores de R^P)

4) t=t+h. Volver a 2)

VARIABLES REESCALADAS:
rgorro=r/dTL
(p_r)gorro=p_r/(m*dTL)
(p_phi)gorro=p_phi/(m*dTL^2)
t en meses (tgorro=t/(60^2*24*30))

La f en función de las variables, tanto reescalada como sin reescalar, puede encontrarse en el
guion del ejercicio.

OBJETIVO: COMPROBAR QUE H'=H-omega*p_phi es una constante del movimiento. Usar variables reescaladas


*/

#include <iostream>
#include <fstream> // Ficheros
#include <cmath>
#include <iomanip> // Formato en los ficheros

#define P 4
#define h 2e-7 // Paso temporal (en meses)
#define itertemp 5e6

// Constantes
#define w 6.8991264 // En rad/mes
#define delta 47.12830601 // En mes^-2
#define mu 0.01230246418 // Adimensional

using namespace std;

void leercondiniciales(string nombre, double yCohete[]);
void RungeKutta(double y[], double t);
double funcionypunto(int j, double y[], double t);
double HamiltonianoModificado(double y[], double t);

/**************************************************** FUNCIÓN PRINCIPAL ****************************************************/
int main(void) {

    // Declaro variables
    double rLuna[2]; // Posición de la luna. Es (cos(wt),sin(wt)) en variables reescaladas
    double yCohete[4]; // Son (r,phi,pr,pphi), reescaladas, claro

    // Condiciones iniciales
    rLuna[0]=1; rLuna[1]=0;
    leercondiniciales("Condiniciales.txt", yCohete);

    // Abro el fichero donde escribiré las posiciones y el del H modificado
    ofstream fich_posiciones, fich_H;

    fich_posiciones.open("posiciones.dat");
    fich_H.open("H_modificado.dat");

    // Itero en el tiempo
    for(int k=0; k<itertemp; k++) {

        // Imprimo los datos
        if(k%1000==0) {
            fich_posiciones << 0 << "," << 0 << "\n"; // Tierra
            fich_posiciones << rLuna[0] << "," << rLuna[1] << "\n"; // Luna
            fich_posiciones << yCohete[0]*cos(yCohete[1]) << "," <<  yCohete[0]*sin(yCohete[1]) << "\n\n"; // Cohete

            fich_H << HamiltonianoModificado(yCohete, k*h) << "\n";
        }
        

        // Llamo a la función RungeKutta para calcular la siguiente yCohete
        RungeKutta(yCohete, k*h);

        // Actualizo la Luna
        rLuna[0]=cos(w*h*(k+1)); rLuna[1]=sin(w*h*(k+1));

    }

    // Imprimo una última vez
    fich_posiciones << 0 << "," << 0 << "\n"; // Tierra
    fich_posiciones << rLuna[0] << "," << rLuna[1] << "\n"; // Luna
    fich_posiciones << yCohete[0]*cos(yCohete[1]) << "," <<  yCohete[0]*sin(yCohete[1]); // Cohete

    fich_H << HamiltonianoModificado(yCohete, h*itertemp);

    fich_posiciones.close();
    fich_H.close();

    return 0;
}
/***************************************************************************************************************************/

/*Función leercondiniciales. Lee las condiciones iniciales desde un fichero*/
void leercondiniciales(string nombre, double yCohete[]) {

    // Condiciones modificables
    ifstream condiniciales;
    condiniciales.open(nombre);

    if(condiniciales.is_open()) {
        for(int j=0; j<4; j++) condiniciales >> yCohete[j];
    }

    return;
}

/*Función RungeKutta. Aplica el método de Runge-Kutta 4 para determinar el vector y(t+h) a partir
del vector y(t)*/
void RungeKutta(double y[], double t) {

    // Contador
    int j;

    // Vectores
    double k[4][P];
    double ymask[3][P];

    // Calculo k1
    for(j=0; j<P; j++) k[0][j]=h*funcionypunto(j,y,t);

    // Calculo y+k_1/2 y k2
    for(j=0; j<P; j++) ymask[0][j]=y[j]+k[0][j]/2.0;
    for(j=0; j<P; j++) k[1][j]=h*funcionypunto(j,ymask[0],t+h/2.0);

    // Calculo y+k_2/2 y k3
    for(j=0; j<P; j++) ymask[1][j]=y[j]+k[1][j]/2.0;
    for(j=0; j<P; j++) k[2][j]=h*funcionypunto(j,ymask[1],t+h/2.0);

    // Calculo y+k_3 y k4
    for(j=0; j<P; j++) ymask[2][j]=y[j]+k[2][j];
    for(j=0; j<P; j++) k[3][j]=h*funcionypunto(j,ymask[2],t+h);

    // El nuevo y
    for(j=0; j<P; j++) y[j]+=(k[0][j]+2*k[1][j]+2*k[2][j]+k[3][j])/6.0;

    return;
}

/*Función funcionypunto. Devuelve la coordenada j-ésima de la función f de la ecuación diferencial
ypunto=f(y(t),t). Modificar según el modelo, claro*/
double funcionypunto(int j, double y[], double t) {

    // Componente 0
    if(j==0) return y[2];

    // Componente 1
    else if(j==1) return y[3]/(y[0]*y[0]);

    // Componente 2
    else if(j==2) {

        // Variables auxiliares
        double r2=y[0]*y[0], rprimacubo=pow(1+r2-2*y[0]*cos(y[1]-w*t),1.5);

        return y[3]*y[3]/(r2*y[0])-delta*(1.0/r2+mu*(y[0]-cos(y[1]-w*t))/rprimacubo);
    }

    // Componente 3
    else return -delta*mu*y[0]*sin(y[1]-w*t)/pow(1+y[0]*y[0]-2*y[0]*cos(y[1]-w*t),1.5);
}

/*Función HamiltonianoModificado. Evalúa el hamiltoniano modificado en un determinado instante de tiempo.
Lo hago todo en variables reescaladas. Por tanto, el resultado está en unidades de m*d_TL^2/meses^2 (unidad de
energía, eso está bien, por lo menos)*/
double HamiltonianoModificado(double y[], double t) {

    // Calculo el rprima
    double rprima=sqrt(1+y[0]*y[0]-2*y[0]*cos(y[1]-w*t));

    return y[2]*y[2]/2.0+y[3]*y[3]/(2*y[0]*y[0])-delta/y[0]-delta*mu/rprima-w*y[3];
}