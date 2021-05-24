#include <stdio.h>
#include <math.h>
#include <omp.h>

omp_lock_t lock_area;
double start, end_point;

int main(){
  omp_init_lock(&lock_area);
  int dim=1000000;
  int vec[dim];
  int vec2[dim];
  int soma=0;
  int soma_total=0;

  start = omp_get_wtime();
  #pragma omp parallel private(soma)
  {
    #pragma omp for
      for (int i=0; i<dim; i++){
        vec[i] = 1;
        vec2[i] = 1;
      }
    
    #pragma omp for
      for (int i=0; i<dim;i++)
        soma += vec[i] * vec2[i];
    omp_set_lock(&lock_area);
      soma_total += soma;
    omp_unset_lock(&lock_area);
  }
  end_point = omp_get_wtime();
  printf("Soma = %d\n", soma_total);
	printf("Levou %lf segundos\n", end_point-start);
  omp_destroy_lock(&lock_area);

}