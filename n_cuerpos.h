#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>

const double Mt = 5.972*10e24;
const double Rt = 149597870700.0;
const double Ms_Mt = 332946;
const double Ms = Ms_Mt*Mt;
const double Gt = 6.67*10e-11;
const double Tt = 2*M_PI*Rt*sqrt(Rt/(Gt*Ms));



//const double G = Gt*Mt*Tt*Tt/(Rt*Rt*Rt);
const double G = 1;
const std::vector<std::vector<double>> parameters = {
  {0.2},
  {0.075, 0.225},
  {44.0/45.0, -1.0*56.0/15.0, 32.0/9.0},
  {19372.0/6561.0, -1*25360.0/2187.0,64448.0/6561.0, -1.0*212.0/729.0},
  {9017.0/3168.0, -1.0*355.0/33.0, 46732.0/5247.0, 49.0/176.0, -1*5103.0/18656.0},
  {35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -1*2187.0/6784.0, 11.0/84.0},
  {5179.0/57600.0, 0.0, 7571.0/16695.0, 393.0/640.0,	-92097.0/339200.0,	187.0/2100.0,	1.0/40.0}
};

struct Punto
{
  Punto(double x1, double y1): x(x1), y(y1) {}
  Punto(): x(0), y(0) {}  
  double x, y;
};

class n_body_2d
{
  public: 
    n_body_2d(int N, double h, double T, std::vector<double> masa, std::vector<Punto> i, std::vector<Punto> v); 
    void RK_4() { __RK_4(); } 
    void DP(double er_rel, double er_abs) { __DP(er_rel,er_abs); }
    
  private: 
    int m_N; 
    double m_h;  
    int M;
    double m_T; 
    std::vector<double> m_masas; 
    std::vector<Punto> m_i;
    std::vector<Punto> m_v; 
    Punto acceleration(Punto p, std::vector<Punto> v, int j);
    Punto RK_aux(Punto v, double n, Punto p);
    Punto DP_aux(std::vector<Punto> p, int n, double h, int index, std::vector<std::vector<Punto>> v);
    Punto escalar(std::vector<Punto> v, double h, int index, int orden, std::vector<std::vector<Punto>> K);
    void __RK_4();
    void __DP(double er_rel, double er_abs); 
};


