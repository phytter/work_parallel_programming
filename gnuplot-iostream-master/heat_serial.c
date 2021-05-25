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

void copy_matriz(float copied[dim][dim], float original[dim][dim]);

void show_matriz(float original[dim][dim]);

void pause_if_needed();

void demo_image(float original[dim][dim]);

int main()
{
  int vec[dim];
  float DX = 0.1;
  float DY = 0.1;

  float U[dim][dim];

  start = omp_get_wtime();

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      U[i][j] = 0.0;
    }
  }

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

  int alpha = 5;
  float DT = pow(DX, 2) / float(2.0 * alpha);
  int M = 5000;

  // finite difference scheme
  int fram = 0;
  int Ncount = 0;
  int loop = 1;
  float ERR = 0.0;
  char base_dir_save[] = "./save/heatmap_save_";

  float U_old[dim][dim];

  float residuo = 0.0;

  while (loop)
  {
    copy_matriz(U_old, U);
    ERR = 0.0;
    for (int i = 1; i < dim - 1; i++)
    {
      for (int j = 1; j < dim - 1; j++)
      {
        residuo = (DT * ((U_old[i + 1][j] - 2.0 * U_old[i][j] + U_old[i - 1][j]) / pow(DX, 2) + (U_old[i][j + 1] - 2.0 * U_old[i][j] + U_old[i][j - 1]) / pow(DY, 2)) + U_old[i][j]) - U[i][j];
        ERR = ERR + fabs(residuo);
        U[i][j] = U[i][j] + residuo;
      }
    }

    // salvar visualizacao
    if (ERR >= 0.01 * maxValue) // allowed error limit is 1% of maximum temperature
    {
      Ncount = Ncount + 1;
      if (Ncount % 50 == 0)
      { // displays movie frame every 50 time steps
        fram = fram + 1;
        char dir[50];
        // sprintf(dir, "%s%d.txt", base_dir_save, fram);
        // save_matriz(dir, U);
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

  // show_matriz(U_old);
  // show_all_history(base_dir_save, fram);

  // show_history("./save/heatmap_save_5.txt");
  end_point = omp_get_wtime();
  printf("Levou %lf segundos\n", end_point-start);
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
    // printf("%s", n_dir);
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
