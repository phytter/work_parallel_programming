#include <stdio.h>
#include <omp.h>

double start, end_point;

int main(){
  int dim=1000000;
  int vec[dim];
  int soma=0;
  start = omp_get_wtime();
  #pragma omp parallel
  {
    #pragma omp for
      for (int i=0; i<dim; i++){
        vec[i] = 1;
      }
    
    #pragma omp for reduction(+: soma)
      for (int i=0; i<dim;i++)
        soma += vec[i];
  }
  end_point = omp_get_wtime();
  printf("Levou %lf segundos\n", end_point-start);
  printf("Soma = %d\n", soma);
}