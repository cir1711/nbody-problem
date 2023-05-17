#include "n_cuerpos.h"

Punto n_body_2d::DP_aux(std::vector<Punto> p, int n, double h, int index, std::vector<std::vector<Punto>> v)
{
  Punto res = Punto(p[index].x, p[index].y);
  for(int i = 0; i < n; i++)
  {
    res.x += parameters[n-1][i] * h * v[i][index].x;
    res.y += parameters[n-1][i] * h * v[i][index].y;
  }
  return res;
}

Punto n_body_2d::escalar(std::vector<Punto> v, double h, int index, int orden, std::vector<std::vector<Punto>> K)
{
  Punto res(v[index].x, v[index].y);
  for(int i = 0; i < 7; i++)
  {
    if(i != 1){
      res.x += h*parameters[orden][i]*K[i][index].x;
      res.y += h*parameters[orden][i]*K[i][index].y;
    }
  }
  return res;
}

void n_body_2d::__DP(double er_rel, double er_abs)
{
  double h = m_h, i = 0.0, anterior = 0; 
  int index = 1, counter = 0; 
  std::vector<std::vector<double>> X;
  std::vector<std::vector<double>> Y;
  std::vector<std::vector<double>> VX;
  std::vector<std::vector<double>> VY;

  for(int i = 0; i < m_N; i++){
	X.push_back(std::vector<double>(M*20,0));
	Y.push_back(std::vector<double>(M*20,0));
	VX.push_back(std::vector<double>(M*20,0));
	VY.push_back(std::vector<double>(M*20,0));
	X[i][0] = m_i[i].x;
	Y[i][0] = m_i[i].y;
	VX[i][0] = m_v[i].x;
	VY[i][0] = m_v[i].y;
  }
  std::vector<std::vector<Punto>> K = std::vector<std::vector<Punto>>(7, std::vector<Punto>(m_N, Punto()));
  std::vector<std::vector<Punto>> L = std::vector<std::vector<Punto>>(7, std::vector<Punto>(m_N, Punto()));

  std::vector<Punto> old = std::vector<Punto>(m_N, Punto());
  std::vector<Punto> v_old = std::vector<Punto>(m_N, Punto()); 
  std::vector<Punto> nw = std::vector<Punto>(m_N, Punto());
  std::vector<Punto> v_new = std::vector<Punto>(m_N, Punto()); 
  std::vector<Punto> nw4 = std::vector<Punto>(m_N, Punto());
  std::vector<Punto> v_new4 = std::vector<Punto>(m_N, Punto());
  
  old = m_i;
  v_old = m_v; 
  K[0] = v_old; 

  FILE* data1;
  FILE* data2;
  FILE* data3;
  FILE* data4;
  data1 = fopen("posiciones_DP.dat","w");
  data2 = fopen("velocidades_DP.dat","w");
  data3 = fopen("tiempo_DP.dat","w");
  data4 = fopen("energias_DP.dat","w");

  fprintf(data3,"%f \n", 0.0);

  while(i <= m_T)
  {
    double error_aux = 0, error = 0; 
    std::vector<Punto> aux(m_N); 
    for(int j = 0; j < m_N; j++)
    {
      K[0][j] = v_old[j];
      L[0][j] = acceleration(old[j], old, j);
    }
    
    for(int j = 0; j < m_N; j++)
    {
      aux[j] = DP_aux(old, 1, h, j, K);
    }
    for(int j = 0; j < m_N; j++)
    {
        K[1][j] = DP_aux(v_old, 1, h, j, L); 
        L[1][j] = acceleration(aux[j] , aux, j);
    }
    
    for(int j = 0; j < m_N; j++)
    {
      aux[j] = DP_aux(old, 2, h, j, K);
    }
    for(int j = 0; j < m_N; j++)
    {
        K[2][j] = DP_aux(v_old, 2, h, j, L); 
        L[2][j] = acceleration(aux[j], aux, j);
    }
    
    for(int j = 0; j < m_N; j++)
    {
      aux[j] = DP_aux(old, 3, h, j, K);
    }
    for(int j = 0; j < m_N; j++)
    {
        K[3][j] = DP_aux(v_old, 3, h, j, L); 
        L[3][j] = acceleration(aux[j], aux, j);
    }
    
    for(int j = 0; j < m_N; j++)
    {
      aux[j] = DP_aux(old, 4, h, j, K);
    }
    for(int j = 0; j < m_N; j++)
    {
        K[4][j] = DP_aux(v_old, 4, h, j, L); 
        L[4][j] = acceleration(aux[j] , aux, j);
    }
    
    for(int j = 0; j < m_N; j++)
    {
      aux[j] = DP_aux(old, 5, h, j, K);
    }
    for(int j = 0; j < m_N; j++)
    {
        K[5][j] = DP_aux(v_old, 5, h, j, L); 
        L[5][j] = acceleration(aux[j] , aux, j);
    }

    for(int j = 0; j < m_N; j++)
    {
      aux[j] = DP_aux(old, 5, h, j, K);
    }
    for(int j = 0; j < m_N; j++)
    {
      K[6][j] = DP_aux(v_old, 5, h, j, L); 
      L[6][j] = acceleration(aux[j], aux, j);
      
      nw[j] = escalar(old, h, j, 6, K);
      v_new[j] = escalar(v_old, h, j, 6, L); 

      nw4[j] = escalar(old, h, j, 5, K);
      v_new4[j] = escalar(v_old, h, j, 5, L); 


      error_aux += std::pow((nw4[j].x-nw[j].x)/(er_rel*std::max(nw4[j].x,nw[j].x)+er_abs),2)/(4*m_N);
      error_aux += std::pow((nw4[j].y-nw[j].y)/(er_rel*std::max(nw4[j].y,nw[j].y)+er_abs),2)/(4*m_N);
      error_aux += std::pow((v_new4[j].x-v_new[j].x)/(er_rel*std::max(v_new4[j].x,v_new[j].x)+er_abs),2)/(4*m_N);
      error_aux += std::pow((v_new4[j].y-v_new[j].y)/(er_rel*std::max(v_new4[j].y,v_new[j].y)+er_abs),2)/(4*m_N);
      
      X[j][index] = nw[j].x;
      Y[j][index] = nw[j].y;
      VX[j][index] = v_new[j].x;
      VY[j][index] = v_new[j].y;
    } 
    error = std::sqrt(error_aux);

    if(error != 0)
      h = 0.5*std::pow(1/error,0.2) * h;
    else
      h = er_rel;
      	
    if(h==anterior)
    {
	counter++;
    }
    else
    {
	counter = 0; 
	anterior = h; 
    }
    if(counter > 10000)
    {
	std::cout<<"\n\n--------------------------------------------------\nPossible error con esos er_abs y er_r\n--------------------------------------------------\n";
	break;
    } 
    //std::cout<<h<<"\t";
    i += h;
    old = nw;
    v_old = v_new; 
    index++;
    fprintf(data3,"%f\n",i);
  }

  std::vector<double> E_c(index,0);
  std::vector<double> E_p(index,0);

  for(int i=0; i<index; i++){
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
  
  for(int t=0; t<index; t++){
      for(int l=0; l<m_N; l++){
          fprintf(data1,"%f \t %f \t", X[l][t], Y[l][t]);
          fprintf(data2,"%f \t %f \t", VX[l][t], VY[l][t]);

      }
      fprintf(data4,"%f \n", E_c[t]+E_p[t]);
      fprintf(data1,"\n");
      fprintf(data2,"\n");
  }
}

