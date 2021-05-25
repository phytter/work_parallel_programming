#include <pthread.h>
#include <stdio.h>
#include <unistd.h> // para o sleep
#include <fstream>
#include <omp.h>
#include <math.h>

using namespace std;

#define NUM_THREADS	4
#define N_VETOR 50000

double start, end_point;

typedef struct 
{
    long long inicio;
    long long fim;
} Intervalo;

long long vetor[N_VETOR];
float soma_l2 = 0;
pthread_mutex_t mutexsum;

void* insert_vetor(void *intervalo){
    Intervalo *args = ( Intervalo *)intervalo;
    for (long long i = args->inicio; i < args->fim; i++){
      vetor[i] = 1;
    }
    pthread_exit(NULL);
}

void* calc_l2(void *intervalo){
    Intervalo *args = ( Intervalo *)intervalo;
    int soma = 0;
    for (long long i = args->inicio; i < args->fim; i++){
      soma += vetor[i] * vetor[i];
    }
    pthread_mutex_lock (&mutexsum);
    soma_l2 += soma;
    pthread_mutex_unlock (&mutexsum);
    pthread_exit(NULL);
}

void operacao(void){
  pthread_t thread[NUM_THREADS];
  Intervalo intervalo[NUM_THREADS];
  long long peace_col = N_VETOR / NUM_THREADS;
  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    long long pos_vector = n_thread * peace_col;
    long long pos_right = pos_vector + peace_col;
    if (n_thread + 1 == NUM_THREADS) {
      pos_right = N_VETOR;
    }
    intervalo[n_thread].inicio = pos_vector;
    intervalo[n_thread].fim = pos_right;

    pthread_create(&thread[n_thread], NULL, insert_vetor, &intervalo[n_thread]);
  }
  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }

  Intervalo intervalo_l2[NUM_THREADS];

  for (int n_thread=0; n_thread < NUM_THREADS; n_thread++){
    long long pos_vector = n_thread * peace_col;
    long long pos_right = pos_vector + peace_col;
    if (n_thread + 1 == NUM_THREADS) {
      pos_right = N_VETOR;
    }
    intervalo_l2[n_thread].inicio = pos_vector;
    intervalo_l2[n_thread].fim = pos_right;
    pthread_create(&thread[n_thread], NULL, calc_l2, &intervalo_l2);
  }

  for(int i=0; i<NUM_THREADS; i++) {
    pthread_join(thread[i], NULL);
  }
  printf("L2 == %lf\n", sqrt(soma_l2));
  
}

void salva_arquivo(long long * vetor, char* nome)
{
	ofstream file;
	file.open(nome);

	for (int i = 0; i < N_VETOR; i++)
		file << vetor[i] << " ";
	
	file.close();		
}

int main()
{
  pthread_mutex_init(&mutexsum, NULL);
    start = omp_get_wtime();
	operacao();
    end_point = omp_get_wtime();
	printf("Levou %lf segundos\n", end_point-start);

	//salva_arquivo(vetor, "resultado_paralelo");

	return 0;
}