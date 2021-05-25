#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <map>
#include <limits>
#include <cmath>
#include <cstdio>
#include <string.h>

using namespace std;
// Warn about use of deprecated functions.
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"
#define dim 51
#define NUM_THREADS	4

double start, end_point;

void mostrar_todos_historico(char dir_base[], int n);

void mostrar_historico(char dir[]);

void salvar_matriz(char dir[], float U[dim][dim]);

void copiar(float copiado[dim][dim], float original[dim][dim]);

void show_matriz(float original[dim][dim]);

void pause_if_needed();

void demo_image(float original[dim][dim]);

void operacao_maxvalue_matriz(void);

void* insert_matriz(void *intervalo);

typedef struct 
{
    long long inicio;
    long long fim;
} Intervalo;

typedef struct 
{
  long long value;
} MaxValue;

void operacao();

void *calc_mat(void *intervalo);

void operacao_insert_matriz(void);

void* maxvalue_matriz(void *intervalo);

float U[dim][dim];
float U_old[dim][dim];
float ERR = 0.0;
float DX = 0.1;
float DY = 0.1;
int alpha = 5;
float DT = pow(DX, 2) / float(2.0 * alpha);
pthread_mutex_t mutexerror;
MaxValue maxvalue[NUM_THREADS];

int main()
{
  int fram = 0;
  int Ncount = 0;
  int loop = 1;
  bool save_breakpoint = false;
  char base_dir_save[] = "./save/heatmap_save_";
  int M = 2000;

  // for (int i = 0; i < dim; i++)
  // {
  //   for (int j = 0; j < dim; j++)
  //   {
  //     U[i][j] = 0.0;
  //   }
  // }

start = omp_get_wtime();

operacao_insert_matriz();

for (int j = 0; j < dim; j++)
{
  U[0][j] = 100.0;
}

for (int i = 23; i < 29; i++)
{
  for (int j = 23; j < 29; j++)
  {
    U[i][j] = 1000.0;
  }
}

float maxValue = 0.0;
for (int i = 0; i < dim; i++)
{
  for (int j = 0; j < dim; j++)
  {
    if (U[i][j] > maxValue)
    {
      maxValue = U[i][j];
    }
  }
}

while (loop)
{
  copiar(U_old, U);
  ERR = 0.0;
  operacao();

    // salvar visualizacao
    if (ERR >= 0.01 * maxValue) // allowed error limit is 1% of maximum temperature
    {
      Ncount = Ncount + 1;
      if (Ncount % 50 == 0 && save_breakpoint)
      { // displays movie frame every 50 time steps
        fram = fram + 1;
        char dir[50];
        sprintf(dir, "%s%d.txt", base_dir_save, fram);
        salvar_matriz(dir, U);
      }

      if (Ncount > M)
      {
        loop = 0;
        printf("solution do not reach steady state in %d time steps\n", M);
        printf(" %.2f error\n", ERR);
      }
    }
    else
    {
      loop = 0;
      printf("solution reach steady state in %d time steps\n", Ncount);
    }
  }
  end_point = omp_get_wtime();
  printf("Levou %lf segundos\n", end_point-start);
  // mostrar_todos_historico(base_dir_save, fram);
  // demo_image(U);
}

void* insert_matriz(void *intervalo){
  Intervalo *args = ( Intervalo *)intervalo;
  for (long long x = args->inicio; x < args->fim; x++){
    for (long long y = 0; y < dim; y++){
      U[x][y] = 0.0;
    }
  }
  pthread_exit(NULL);
}

void* maxvalue_matriz(void *intervalo){
  Intervalo *args = ( Intervalo *)intervalo;
  for (long long x = args->inicio; x < args->fim; x++){
    for (long long y = 0; y < dim; y++){
      U[x][y] = 0.0;
    }
  }
  pthread_exit(NULL);
}

void operacao_maxvalue_matriz(void){
  pthread_t thread[NUM_THREADS];
  Intervalo intervalo[NUM_THREADS];
  Intervalo intervalo[NUM_THREADS];
  long long peace_col = dim / NUM_THREADS;
  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    long long pos_vector = n_thread * peace_col;
    long long pos_right = pos_vector + peace_col;
    if (n_thread + 1 == NUM_THREADS) {
      pos_right = dim;
    }
    intervalo[n_thread].inicio = pos_vector;
    intervalo[n_thread].fim = pos_right;

    pthread_create(&thread[n_thread], NULL, insert_matriz, &intervalo[n_thread]);
  }
  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }
}

void operacao_insert_matriz(void){
  pthread_t thread[NUM_THREADS];
  Intervalo intervalo[NUM_THREADS];
  
  long long peace_col = dim / NUM_THREADS;
  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    long long pos_vector = n_thread * peace_col;
    long long pos_right = pos_vector + peace_col;
    if (n_thread + 1 == NUM_THREADS) {
      pos_right = dim;
    }
    intervalo[n_thread].inicio = pos_vector;
    intervalo[n_thread].fim = pos_right;

    pthread_create(&thread[n_thread], NULL, insert_matriz, &intervalo[n_thread]);
  }
  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }
}

void *calc_mat(void *intervalo){
  Intervalo *args = ( Intervalo *)intervalo;
  float residuo = 0.0;
  for (long long i = args->inicio; i < args->fim; i++){
    for (long long j = 1; j < dim - 1; j++){
      residuo = (DT * ((U_old[i + 1][j] - 2.0 * U_old[i][j] + U_old[i - 1][j]) / pow(DX, 2) + (U_old[i][j + 1] - 2.0 * U_old[i][j] + U_old[i][j - 1]) / pow(DY, 2)) + U_old[i][j]) - U[i][j];
        U[i][j] = U[i][j] + residuo;
        pthread_mutex_lock (&mutexerror);
        ERR = ERR + fabs(residuo);
        pthread_mutex_unlock (&mutexerror);
    }
  }
  pthread_exit(NULL);
}

void operacao(){
  pthread_t thread[NUM_THREADS];
  Intervalo intervalo[NUM_THREADS];
  long long peace_col = (dim - 1) / NUM_THREADS;
  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    long long pos_vector = n_thread * peace_col + 1;
    long long pos_right = pos_vector + peace_col;
    if (n_thread + 1 == NUM_THREADS) {
      pos_right = dim;
    }
    intervalo[n_thread].inicio = pos_vector;
    intervalo[n_thread].fim = pos_right;

    pthread_create(&thread[n_thread], NULL, calc_mat, &intervalo[n_thread]);
  }
  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }
}

void mostrar_todos_historico(char dir_base[], int n)
{
  for (int i; i < n; i++)
  {
    int size_str = strlen(dir_base) + 4 + 1;
    if (i > 9)
      size_str += 1;
    char n_dir[size_str];
    sprintf(n_dir, "%s%d.txt", dir_base, i);
    mostrar_historico(n_dir);
  }
}

void mostrar_historico(char dir[])
{
  FILE *fp;
  float U[dim][dim];
  int am, an;
  fp = fopen(dir, "r");
  if (fp == NULL)
    return;
  fscanf(fp, "%d %d", &am, &an);
  for (int i = 0; i < am; i++)
    for (int j = 0; j < an; j++)
      fscanf(fp, "%f", &U[i][j]);
  fclose(fp);
  demo_image(U);
}

void salvar_matriz(char dir[], float U[dim][dim])
{
  FILE *fp;
  fp = fopen(dir, "w");
  if (fp == NULL)
    return;
  fprintf(fp, "%d %d\n", dim, dim);
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
      fprintf(fp, " %.3f", U[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

void copiar(float copiado[dim][dim], float original[dim][dim])
{
  int count;

  for (count = 0; count < dim; count++)
    for (int y = 0; y < dim; y++)
      copiado[count][y] = original[count][y];
}

void show_matriz(float original[dim][dim])
{

  for (int count = 0; count < dim; count++)
  {
    for (int y = 0; y < dim; y++)
      printf("%.2f ", original[count][y]);
    printf("\n");
  }
}

void pause_if_needed()
{
#ifdef _WIN32
  // For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
  // the gnuplot window doesn't get closed.
  std::cout << "Press enter to exit." << std::endl;
  std::cin.get();
#endif
}

void demo_image(float original[dim][dim])
{
  // Example of plotting an image.  Of course you are free (and encouraged) to
  // use Blitz or Armadillo rather than std::vector in these situations.

  Gnuplot gp;

  std::vector<std::vector<double>> image;
  for (int j = 0; j < dim; j++)
  {
    std::vector<double> row;
    for (int i = 0; i < dim; i++)
    {
      // double x = (i - 50.0) / 5.0;
      // double y = (j - 50.0) / 5.0;
      // double z = std::cos(sqrt(x * x + y * y));
      row.push_back(original[j][i]);
    }
    image.push_back(row);
  }

  // It may seem counterintuitive that send1d should be used rather than
  // send2d.  The explanation is as follows.  The "send2d" method puts each
  // value on its own line, with blank lines between rows.  This is what is
  // expected by the splot command.  The two "dimensions" here are the lines
  // and the blank-line-delimited blocks.  The "send1d" method doesn't group
  // things into blocks.  So the elements of each row are printed as columns,
  // as expected by Gnuplot's "matrix with image" command.  But images
  // typically have lots of pixels, so sending as text is not the most
  // efficient (although, it's not really that bad in the case of this
  // example).  See the binary version below.
  //
  //gp << "plot '-' matrix with image\n";
  //gp.send1d(image);

  // To be honest, Gnuplot's documentation for "binary" and for "image" are
  // both unclear to me.  The following example comes by trial-and-error.
  gp << "plot '-' binary" << gp.binFmt2d(image, "array") << "with image\n";
  gp.sendBinary2d(image);

  pause_if_needed();
}
