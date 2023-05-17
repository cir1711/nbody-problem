#include "n_cuerpos.h"
n_body_2d::n_body_2d(int N, double h, double T, std::vector<double> masa, std::vector<Punto> i, std::vector<Punto> v)
{
  m_T = T; 
  m_N = N; 
  m_h = h; 
  M = int(T/h+1);
  m_masas = masa; 
  m_i = i;
  m_v = v;
}

Punto n_body_2d::acceleration(Punto p, std::vector<Punto> v, int j)
{
  Punto res = Punto(); 
  for(int i = 0; i < m_N; i++)
  {
    if(i != j)
    {
      double modx = std::sqrt((v[i].x - p.x)*(v[i].x - p.x)+(v[i].y - v[j].y)*(v[i].y - v[j].y));
		  res.x = res.x + G*m_masas[i]*(v[i].x-p.x)/(modx*modx*modx);
      double mody = std::sqrt((v[i].y - p.y)*(v[i].y - p.y)+(v[i].x - v[j].x)*(v[i].x - v[j].x));
		  res.y = res.y + G*m_masas[i]*(v[i].y-p.y)/(mody*mody*mody);
    }
  }
  return res;
}

Punto n_body_2d::RK_aux(Punto v, double n, Punto p)
{
  Punto res = Punto();
  res.x = v.x + n * m_h * p.x; 
  res.y = v.y + n * m_h * p.y;
  return res;
}

void n_body_2d::__RK_4()
{
  std::vector<std::vector<double>> X;
	std::vector<std::vector<double>> Y;
	std::vector<std::vector<double>> VX;
	std::vector<std::vector<double>> VY;
	
	for(int i = 0; i < m_N; i++){
		X.push_back(std::vector<double>(M,0));
		Y.push_back(std::vector<double>(M,0));
		VX.push_back(std::vector<double>(M,0));
		VY.push_back(std::vector<double>(M,0));
    X[i][0] = m_i[i].x;
		Y[i][0] = m_i[i].y;
		VX[i][0] = m_v[i].x;
		VY[i][0] = m_v[i].y;
	}

  std::vector<std::vector<Punto>> K = std::vector<std::vector<Punto>>(4, std::vector<Punto>(m_N, Punto()));
  std::vector<std::vector<Punto>> L = std::vector<std::vector<Punto>>(4, std::vector<Punto>(m_N, Punto()));
  
  std::vector<Punto> old = std::vector<Punto>(m_N, Punto());
  std::vector<Punto> v_old = std::vector<Punto>(m_N, Punto()); 
  std::vector<Punto> nw = std::vector<Punto>(m_N, Punto());
  std::vector<Punto> v_new = std::vector<Punto>(m_N, Punto()); 
  
  old = m_i;
  v_old = m_v; 

  FILE* data1;
  FILE* data2;
  FILE* data3;
  FILE* data4;
  data1 = fopen("posiciones_RK.dat","w");
  data2 = fopen("velocidades_RK.dat","w");
  data3 = fopen("tiempo_RK.dat","w");
  data4 = fopen("energias_RK.dat","w");

  fprintf(data3,"%f \n", 0.0);
  
  double tiempo = 0.0;

  for(int i = 0; i < M; i++)
  {
    std::vector<Punto> aux(m_N, Punto());
    for(int j = 0; j < m_N; j++)
    {
      L[0][j] = acceleration(old[j], old, j);
      K[0][j] = v_old[j];
    }
    for(int j = 0; j < m_N; j++)
     	aux[j] = RK_aux(old[j], 0.5, K[0][j]);
    
    for(int j = 0; j < m_N; j++)
    {
      L[1][j] = acceleration(aux[j], aux, j);
      K[1][j] = RK_aux(v_old[j], 0.5, L[0][j]);
    }
    
    for(int j = 0; j < m_N; j++)
     	aux[j] = RK_aux(old[j], 0.5, K[1][j]); 
    
    for(int j = 0; j < m_N; j++)
    {
      L[2][j] = acceleration(aux[j], aux, j);
      K[2][j] = RK_aux(v_old[j], 0.5, L[1][j]);
    }

     for(int j = 0; j < m_N; j++)
     	aux[j] = RK_aux(old[j], 1.0, K[2][j]);
     
     for(int j = 0; j < m_N; j++)
    {      
      L[3][j] = acceleration(aux[j], aux, j);
      K[3][j] = RK_aux(v_old[j], 1.0, L[2][j]);

      nw[j].x = old[j].x + ((K[0][j].x + K[1][j].x*2 + K[2][j].x*2 + K[3][j].x)*m_h)/6;
      nw[j].y = old[j].y + ((K[0][j].y + K[1][j].y*2 + K[2][j].y*2 + K[3][j].y)*m_h)/6;
      v_new[j].x = v_old[j].x +((L[0][j].x + L[1][j].x*2 + L[2][j].x*2 + L[3][j].x)*m_h)/6;
      v_new[j].y = v_old[j].y +((L[0][j].y + L[1][j].y*2 + L[2][j].y*2 + L[3][j].y)*m_h)/6;

      X[j][i] = nw[j].x;
      Y[j][i] = nw[j].y;
      VX[j][i] = v_new[j].x;
      VY[j][i] = v_new[j].y;
    }
    old = nw; 
    v_old = v_new;
    tiempo += m_h;
    fprintf(data3,"%f \n", tiempo);
    
  }

  std::vector<double> E_c(M,0);
  std::vector<double> E_p(M,0);

  for(int i=0; i<M; i++){
    for(int j=0; j<m_N; j++){
      E_c[i] += 0.5*m_masas[j]*(pow(VX[j][i],2)+pow(VY[j][i],2));
      for(int k=0; k<m_N; k++){
        if(k!=j){
          E_p[i] += -G*m_masas[j]*m_masas[k]/sqrt(
            pow(X[k][i]-X[j][i],2)+pow(Y[k][i]-Y[j][i],2)
          );
        }
      } 
    }
  }	
  
  
	
	for(int t=0; t<M; t++){
		for(int l=0; l<m_N; l++){
			fprintf(data1,"%f \t %f \t", X[l][t], Y[l][t] );
			fprintf(data2,"%f \t %f \t", VX[l][t], VY[l][t] );
		}
    fprintf(data4, "%f \n", E_c[t]+E_p[t]);
		fprintf(data1,"\n");
		fprintf(data2,"\n");
	}
}
