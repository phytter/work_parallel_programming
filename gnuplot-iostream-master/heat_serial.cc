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
#define dim 51

double start, end_point;

void mostrar_todos_historico(char dir_base[], int n);

void mostrar_historico(char dir[]);

void salvar_matriz(char dir[], float U[dim][dim]);

void copiar(float copiado[dim][dim], float original[dim][dim]);

void show_matriz(float original[dim][dim]);

void pause_if_needed();

void demo_image(float original[dim][dim]);

int main()
{
  int vec[dim];
  float DX = 0.1;
  float DY = 0.1;
  int Nx = 5;
  int Ny = 5;

  float U[dim][dim];

  int tamX = int(Nx / DX) + 1;
  int tamY = int(Ny / DY) + 1;
  float X[tamX];
  float Y[tamY];
  float sum = 0.0;
  X[0] = 0.0;
  Y[0] = 0.0;

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

  for (int i = 1; i < tamX; i++)
  {
    sum += DX;
    X[i] = sum;
  }

  sum = 0.0;
  for (int i = 1; i < tamY; i++)
  {
    sum = sum + DY;
    Y[i] = sum;
  }

  int alpha = 5;
  float DT = pow(DX, 2) / float(2.0 * alpha);
  int M = 2000;

  // finite difference scheme
  int fram = 0;
  int Ncount = 0;
  int loop = 1;
  float ERR = 0.0;
  char base_dir_save[] = "./save/heatmap_save_";

  float U_old[dim][dim];

  float residuo = 0.0;
  start = omp_get_wtime();

  while (loop)
  {
    copiar(U_old, U);
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
        // salvar_matriz(dir, U);
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
  // mostrar_todos_historico(base_dir_save, fram);

  // mostrar_historico("./save/heatmap_save_5.txt");
  end_point = omp_get_wtime();
  printf("Levou %lf segundos\n", end_point-start);
  demo_image(U);

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
    // printf("%s", n_dir);
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
