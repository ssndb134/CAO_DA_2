#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define PI_O 3.141592653589

double epoch(int i){
    printf("====EXECUTING %d EPOCH CURRENTLY====\n", i);
    long long int n = 1e8;
    srand(time(NULL));
    double x, y;
    long long int c = 0;
   
    for (long long int i = 0;i<n;i++){
        //printf("LOOP %lld : ", i);
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;
        //printf("X: %.2f, Y: %.2f\n", x, y);
        if(x*x + y*y <= 1) c++;
    }
    printf("C: %lld, N: %lld\n", c, n);
    double pi = 4.0 * c / n;
    printf("Intermediate Estimated Value of PI: %.12f\n", pi);
   
    return pi;
}
int main() {
    // Write C code here
    printf("Estimating Value of PI\n");
    int n = 10;
    printf("----------EXECUTING FOR %d EPOCHS----------\n", n);
    unsigned int start = time(NULL);
    double start_time = omp_get_wtime();
    printf("-----START: TIME = %u------\n", start);
    double pi = 0;
    for(int i = 0;i<n;i++){
        pi += epoch(i+1);
    }
    pi /= n;
    double stop_time = omp_get_wtime();
    unsigned int stop = time(NULL);
    printf("-----STOP: TIME = %u------\n", stop);
    printf("Time Taken: %u SECONDS", stop-start);
    printf("Time Taken (OMP): %.2f SECONDS", stop_time-start_time);
    double accuracy = (pi-PI_O < 0)?((1
    +((pi - PI_O) / PI_O)) * 100):((1
    -((pi - PI_O) / PI_O)) * 100);
    printf("\nEstimated PI Value: %.12f\n", pi);
    printf("\nACCURACY: %.4f%\n", accuracy);
    return 0;
}
