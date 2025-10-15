#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define PI_O 3.141592653589

double epoch(int i){
	
	printf("====EXECUTING EPOCH %d CURRENTLY====\n", i);
    long long int n = 1e8;
    
    double x, y;
    long long int count = 0;
    #pragma omp parallel
    {
      //unsigned int seed = time(NULL) ^ omp_get_thread_num();
      //unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)(12345 + 10007 * omp_get_thread_num());
      unsigned int seed = (unsigned int)(time(NULL) + i * 1009 + omp_get_thread_num() * 7919);
      long long int c = 0;
      struct drand48_data buffer;
	srand48_r(seed, &buffer);
      #pragma omp for schedule(static)
      for (long long int j = 0;j<n;j++){
        //printf("LOOP %lld : ", i);
        //x = (double)rand_r(&seed) / RAND_MAX;
        //y = (double)rand_r(&seed) / RAND_MAX;
        drand48_r(&buffer, &x);
    	drand48_r(&buffer, &y);
        //printf("X: %.2f, Y: %.2f\n", x, y);
        if(x*x + y*y <= 1) c++;
      }
      
      #pragma omp atomic
      count += c;
     }
    
   
    printf("C: %lld, N: %lld\n", count, n);
    double pi = 4.0 * count / n;
    printf("Intermediate Estimated Value of PI: %.12f\n", pi);
    
    return pi;
    
}
int main (int argc, char *argv[]) {
	printf("Estimating Value of PI\n");
	int n = 10;
    printf("----------EXECUTING FOR %d EPOCHS----------\n", n);
    unsigned int start = time(NULL);
    printf("-----START: TIME = %u------\n", start);
    double pi = 0;
    for(int i = 0;i<n;i++){
        pi += epoch(i+1);
    }
    pi /= n;
    unsigned int stop = time(NULL);
    printf("-----STOP: TIME = %u------\n", stop);
    printf("Time Taken: %u SECONDS\n", stop-start);
    double accuracy = (pi-PI_O < 0)?((1
    +((pi - PI_O) / PI_O)) * 100):((1
    -((pi - PI_O) / PI_O)) * 100);
    printf("\nEstimated PI Value: %.12f\n", pi);
    printf("\nACCURACY: %.4f%\n", accuracy);
    printf("Total NUMBER of THREADS: %d", omp_get_num_threads());
    return 0;
}

