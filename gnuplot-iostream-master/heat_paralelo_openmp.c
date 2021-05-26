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

// Warn about use of deprecated functions.
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"
#define dim 100

double start, end_point;

void show_all_history(char dir_base[], int n);

void show_history(char dir[]);

void save_matriz(char dir[], float U[dim][dim]);

void copy_matriz(float copiado[dim][dim], float original[dim][dim]);

void show_matriz(float original[dim][dim]);

void pause_if_needed();

void demo_image(float original[dim][dim]);

int main()
{
  float U[dim][dim];
  float U_old[dim][dim];
  float DX = 0.1;
  float DY = 0.1;

  int alpha = 5;
  float DT = pow(DX, 2) / float(2.0 * alpha);
  int M = 5000;
  int fram = 0;
  int Ncount = 0;
  int loop = 1;
  float ERR = 0.0;
  bool save_breakpoint = false;
  char base_dir_save[] = "./save/heatmap_save_";
  float residuo = 0.0;
  struct Compare { float val; };    
  #pragma omp declare reduction(maximum : struct Compare : omp_out = omp_in.val > omp_out.val ? omp_in : omp_out)

  struct Compare max; 
  max.val = 0; 
  float maxValue = 0.0;

  int num_t = 4;

  start = omp_get_wtime();

  #pragma omp parallel num_threads(num_t) shared(U)
  {
    #pragma omp for collapse(2)
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          U[i][j] = 0.0;

    #pragma omp for
    for (int j = 0; j < dim; j++)
    {
      U[0][j] = 100.0;
    }

    #pragma omp for collapse(2)
    for (int i = 23; i < 29; i++)
    {
      for (int j = 23; j < 29; j++)
      {
        U[i][j] = 1000.0;
      }
    }

    #pragma omp for reduction(maximum:max) collapse(2)
      for (int i = 0; i < dim; i++)
      {
        for (int j = 0; j < dim; j++)
        {
          if (U[i][j] > max.val)
          {
            max.val = U[i][j];
          }
        }
      }
  }
  maxValue = max.val;

  while (loop)
  {
    copy_matriz(U_old, U);
    ERR = 0.0;
    #pragma omp parallel num_threads(num_t) shared(U, U_old, residuo, ERR)
    {
      #pragma omp for reduction(+: ERR) collapse(2)
      for (int i = 1; i < dim - 1; i++)
      {
        for (int j = 1; j < dim - 1; j++)
        {
          residuo = (DT * ((U_old[i + 1][j] - 2.0 * U_old[i][j] + U_old[i - 1][j]) / pow(DX, 2) + (U_old[i][j + 1] - 2.0 * U_old[i][j] + U_old[i][j - 1]) / pow(DY, 2)) + U_old[i][j]) - U[i][j];
          U[i][j] = U[i][j] + residuo;
          ERR = ERR + fabs(residuo);
        }
      }
    }

    if (ERR >= 0.01 * maxValue) // Permite o limite do erro ate 1% do valor maximo
    {
      Ncount = Ncount + 1;
      if (Ncount % 50 == 0 && save_breakpoint)
      {
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
  // show_history("./save/heatmap_save_1.txt");
  demo_image(U);
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

void copy_matriz(float copiado[dim][dim], float original[dim][dim])
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
