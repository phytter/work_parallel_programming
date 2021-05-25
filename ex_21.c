#include <stdio.h>
#include <math.h>
#include <omp.h>

omp_lock_t lock_area;
double start, eend;

int main()
{
  int TAM = 800, soma = 0, total = 0;
  int matrizA[TAM][TAM], matrizB[TAM][TAM], matrizC[TAM][TAM];
  int i, j, k, aux = 0;

  for (i = 0; i < TAM; i++)
  {
    for (j = 0; j < TAM; j++)
    {
      matrizA[i][i] = 1;
      matrizB[i][i] = 2;
    }
  }

  start = omp_get_wtime();
#pragma omp parallel num_threads(16)
{

  #pragma omp for schedule(static) collapse(2)
    for(i = 0; i < TAM; i++) {
      for(j = 0; j < TAM; j++) {
        matrizC[i][j] = 0;
        for(int x = 0; x < TAM; x++) {
          aux +=  matrizA[i][x] * matrizB[x][j];
        }
        matrizC[i][j] = aux;
        aux = 0;
      }
    }
}

  eend = omp_get_wtime();

  printf("Levou %lf segundos\n", eend - start);

  return 0;
}
