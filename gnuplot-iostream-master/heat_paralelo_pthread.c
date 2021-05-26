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
#include "gnuplot-iostream.h"

using namespace std;
// Warn about use of deprecated functions.
#define GNUPLOT_DEPRECATE_WARN
#define dim 51
#define NUM_THREADS	4

double start, end_point;

void show_all_history(char dir_base[], int n);

void show_history(char dir[]);

void save_matriz(char dir[], float U[dim][dim]);

void copy_matriz(float copied[dim][dim], float original[dim][dim]);

void show_matriz(float original[dim][dim]);

void pause_if_needed();

void demo_image(float original[dim][dim]);

float operation_maxvalue_matriz(void);

void* insert_matriz(void *interval);

typedef struct 
{
    long long inicio;
    long long fim;
} Interval;

struct ThreadParameters
{
    int start;
    int end;

    float largest;
};

void operation();

void *calc_mat(void *interval);

void operation_insert_matriz(void);

void* find_max(void* args);

float U[dim][dim];
float U_old[dim][dim];
float ERR = 0.0;
float DX = 0.1;
float DY = 0.1;
int alpha = 5;
float DT = pow(DX, 2) / float(2.0 * alpha);
pthread_mutex_t mutexerror;

int main()
{
  int fram = 0;
  int Ncount = 0;
  int loop = 1;
  bool save_breakpoint = false;
  char base_dir_save[] = "./save/heatmap_save_";
  int M = 5000;

start = omp_get_wtime();

operation_insert_matriz();

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

float maxValue = operation_maxvalue_matriz();

while (loop)
{
  copy_matriz(U_old, U);
  ERR = 0.0;
  operation();

  // salvar visualizacao
  if (ERR >= 0.01 * maxValue) // allowed error limit is 1% of maximum temperature
  {
    Ncount = Ncount + 1;
    if (Ncount % 50 == 0 && save_breakpoint)
    { // displays movie frame every 50 time steps
      fram = fram + 1;
      char dir[50];
      sprintf(dir, "%s%d.txt", base_dir_save, fram);
      save_matriz(dir, U);
    }

    if (Ncount > M)
    {
      loop = 0;
      printf("Solução não alcançou o estado em %d passos\n", M);
      printf(" %.2f error\n", ERR);
    }
  }
  else
  {
    loop = 0;
    printf("Solução alcançou o estado em %d passos\n", Ncount);
  }
}
  end_point = omp_get_wtime();
  printf("Levou %lf segundos\n", end_point-start);
  // show_all_history(base_dir_save, fram);
  demo_image(U);
}

void* insert_matriz(void *interval){
  Interval *args = ( Interval *)interval;
  for (long long x = args->inicio; x < args->fim; x++){
    for (long long y = 0; y < dim; y++){
      U[x][y] = 0.0;
    }
  }
  pthread_exit(NULL);
}

void* find_max(void* args)
{
    struct ThreadParameters* params = (struct ThreadParameters*)args;
    int start = params->start;
    int end = params->end;
    int largest = 0;

    for (int i = start; i < end; i++)
    {
      for (int j = 0; j < dim; j++)
        if (U[i][j] > largest)
        {
            largest = U[i][j];
        }
    }

    // write the result back to the parameter structure
    params->largest = largest;

    return NULL;
}

float operation_maxvalue_matriz(void){
    float largest;

    pthread_t threads[NUM_THREADS] = {0};
    struct ThreadParameters thread_parameters[NUM_THREADS]  = {0};

    largest = 0;

    for (int i = 0; i < NUM_THREADS; i++)
    {
        thread_parameters[i].start = i * (dim / NUM_THREADS);
        thread_parameters[i].end = (i+1) * (dim / NUM_THREADS);
        thread_parameters[i].largest = 0;
        if (i == NUM_THREADS-1)
          thread_parameters[i].end = dim;
        pthread_create(&threads[i], NULL, find_max, &thread_parameters[i]);
    }

    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], NULL);
    }
 
    for (int i = 0; i < NUM_THREADS; i++)
    {
      if (thread_parameters[i].largest > largest)
      {
          largest = thread_parameters[i].largest;
      }
    }

  return largest;
}

void operation_insert_matriz(void){
  pthread_t thread[NUM_THREADS];
  Interval interval[NUM_THREADS];
  
  long long peace_col = dim / NUM_THREADS;
  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    interval[n_thread].inicio = n_thread * peace_col;
    interval[n_thread].fim = (n_thread + 1) * peace_col;

    if (n_thread == NUM_THREADS-1)
      interval[n_thread].fim = dim;

    pthread_create(&thread[n_thread], NULL, insert_matriz, &interval[n_thread]);
  }
  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }
}

void *calc_mat(void *interval){
  Interval *args = ( Interval *)interval;
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

void operation(){
  pthread_t thread[NUM_THREADS];
  Interval interval[NUM_THREADS];
  long long peace_col = (dim - 1) / NUM_THREADS;
  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    interval[n_thread].inicio = (n_thread * peace_col) + 1;
    interval[n_thread].fim = (n_thread + 1) * peace_col + 1;

    if (n_thread == NUM_THREADS-1)
      interval[n_thread].fim = dim - 1;

    pthread_create(&thread[n_thread], NULL, calc_mat, &interval[n_thread]);
  }
  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }
}

void show_all_history(char dir_base[], int n)
{
  for (int i; i < n; i++)
  {
    int size_str = strlen(dir_base) + 4 + 1;
    if (i > 9)
      size_str += 1;
    char n_dir[size_str];
    sprintf(n_dir, "%s%d.txt", dir_base, i);
    show_history(n_dir);
  }
}

void show_history(char dir[])
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

void save_matriz(char dir[], float U[dim][dim])
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

void copy_matriz(float copied[dim][dim], float original[dim][dim])
{
  int count;

  for (count = 0; count < dim; count++)
    for (int y = 0; y < dim; y++)
      copied[count][y] = original[count][y];
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
  Gnuplot gp;

  std::vector<std::vector<double>> image;
  for (int j = 0; j < dim; j++)
  {
    std::vector<double> row;
    for (int i = 0; i < dim; i++)
    {
      row.push_back(original[j][i]);
    }
    image.push_back(row);
  }

  gp << "plot '-' binary" << gp.binFmt2d(image, "array") << "with image\n";
  gp.sendBinary2d(image);

  pause_if_needed();
}
