#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#define PI_O 3.141592653589
/*
double epoch_parallel(int i, int t){
	
	printf("====EXECUTING EPOCH %d CURRENTLY====\n", i);
    long long int n = 1e8;
    
    double x, y;
    long long int count = 0;
    //omp_set_num_threads(t);
    #pragma omp parallel
    {
      //unsigned int seed = time(NULL) ^ omp_get_thread_num();
      //unsigned int seed = (unsigned int)time(NULL) ^ (unsigned int)(12345 + 10007 * omp_get_thread_num());
      unsigned int seed = (unsigned int)(time(NULL) + i * 1009 + omp_get_thread_num() * 7919);
      long long int c = 0;
      struct drand48_data buffer;
	srand48_r(seed, &buffer);
      #pragma omp for reduction(+:count)
      for (long long int j = 0;j<n;j++){
        //printf("LOOP %lld : ", i);
        //x = (double)rand_r(&seed) / RAND_MAX;
        //y = (double)rand_r(&seed) / RAND_MAX;
        drand48_r(&buffer, &x);
    	drand48_r(&buffer, &y);
        //printf("X: %.2f, Y: %.2f\n", x, y);
        if(x*x + y*y <= 1) count++;
      }
      /*
      #pragma omp atomic
      count += c;
      
     }
    
   
    printf("C: %lld, N: %lld\n", count, n);
    double pi = 4.0 * count / n;
    printf("Intermediate Estimated Value of PI: %.12f\n", pi);
    
    return pi;
    
}

*/

double epoch_parallel(int i, int t){
    printf("====EXECUTING EPOCH %d CURRENTLY====\n", i);
    long long int n = 1e8;
    
    double x, y;
    long long int count = 0;

    // ---- NEW: Pre-generate random numbers ----
    double *x_arr = malloc(n * sizeof(double));
    double *y_arr = malloc(n * sizeof(double));
    for (long long j = 0; j < n; j++) {
        x_arr[j] = (double)rand() / RAND_MAX;
        y_arr[j] = (double)rand() / RAND_MAX;
    }

    // ---- Parallel counting loop ----
    #pragma omp parallel for reduction(+:count)
    for (long long j = 0; j < n; j++) {
        if(x_arr[j]*x_arr[j] + y_arr[j]*y_arr[j] <= 1) count++;
    }

    // ---- Free memory ----
    free(x_arr);
    free(y_arr);

    printf("C: %lld, N: %lld\n", count, n);
    double pi = 4.0 * count / n;
    printf("Intermediate Estimated Value of PI: %.12f\n", pi);
    
    return pi;
}

double epoch_seq(int i){
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
double pi_estimate_par(int t){
	printf("Estimating Value of PI\n");
	int n = 10;
    printf("----------EXECUTING FOR %d EPOCHS----------\n", n);
    unsigned int start = time(NULL);
    printf("-----START: TIME = %u------\n", start);
    //omp_set_num_threads(t);
    double pi = 0;
    for(int i = 0;i<n;i++){
    
        pi += epoch_parallel(i+1, t);
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
    return pi;
}

double pi_estimate_seq(){
	printf("Estimating Value of PI\n");
    int n = 10;
    printf("----------EXECUTING FOR %d EPOCHS----------\n", n);
    unsigned int start = time(NULL);
    double start_time = omp_get_wtime();
    printf("-----START: TIME = %u------\n", start);
    double pi = 0;
    for(int i = 0;i<n;i++){
        pi += epoch_seq(i+1);
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
    return pi;
}

double f(double x) {
    return x * x;  // Example: f(x) = x^2
}

// Riemann sum using left endpoint method
double riemann_sum_left_par(double a, double b, long n) {
    double delta = (b - a) / n;  // Width of each rectangle
    double sum = 0.0;
    
    #pragma omp parallel for reduction(+:sum)
    for (long i = 0; i < n; i++) {
        double x = a + i * delta;
        sum += f(x);
    }
    
    return sum * delta;
}

// Riemann sum using right endpoint method
double riemann_sum_right_par(double a, double b, long n) {
    double delta = (b - a) / n;
    double sum = 0.0;
    
    #pragma omp parallel for reduction(+:sum)
    for (long i = 1; i <= n; i++) {
        double x = a + i * delta;
        sum += f(x);
    }
    
    return sum * delta;
}

// Riemann sum using midpoint method
double riemann_sum_midpoint_par(double a, double b, long n) {
    double delta = (b - a) / n;
    double sum = 0.0;
    
    #pragma omp parallel for reduction(+:sum)
    for (long i = 0; i < n; i++) {
        double x = a + (i + 0.5) * delta;
        sum += f(x);
    }
    
    return sum * delta;
}

double riemann_sums_par(int t){
	double a = 0.0;      // Lower bound
    double b = 1.0;      // Upper bound
    long n = 10000000;   // Number of rectangles
    
    // Set number of threads (optional)
    //omp_set_num_threads(8);
    //omp_set_num_threads(t);
    printf("Computing Riemann Sum for f(x) = x^2 from %.2f to %.2f\n", a, b);
    printf("Number of rectangles: %ld\n", n);
    printf("Number of threads: %d\n\n", omp_get_max_threads());
    
    // Left endpoint method
    double start = omp_get_wtime();
    double result_left = riemann_sum_left_par(a, b, n);
    double end = omp_get_wtime();
    printf("Left Riemann Sum: %.10f (Time: %.6f seconds)\n", 
           result_left, end - start);
    
    // Right endpoint method
    start = omp_get_wtime();
    double result_right = riemann_sum_right_par(a, b, n);
    end = omp_get_wtime();
    printf("Right Riemann Sum: %.10f (Time: %.6f seconds)\n", 
           result_right, end - start);
    
    // Midpoint method
    start = omp_get_wtime();
    double result_mid = riemann_sum_midpoint_par(a, b, n);
    end = omp_get_wtime();
    printf("Midpoint Riemann Sum: %.10f (Time: %.6f seconds)\n", 
           result_mid, end - start);
    
    // Exact integral for f(x) = x^2 from 0 to 1 is 1/3
    printf("\nExact value (for x^2 from 0 to 1): %.10f\n", 1.0/3.0);
    return ((result_left + result_right + result_mid )/ 3.0);
}


double riemann_sum_left(double a, double b, long n) {
    double delta = (b - a) / n;  // Width of each rectangle
    double sum = 0.0;
    
    
    for (long i = 0; i < n; i++) {
        double x = a + i * delta;
        sum += f(x);
    }
    
    return sum * delta;
}

// Riemann sum using right endpoint method
double riemann_sum_right(double a, double b, long long n) {
    double delta = (b - a) / n;
    double sum = 0.0;
    
    //#pragma omp parallel for reduction(+:sum)
    for (long i = 1; i <= n; i++) {
        double x = a + i * delta;
        sum += f(x);
    }
    
    return sum * delta;
}

// Riemann sum using midpoint method
double riemann_sum_midpoint(double a, double b, long n) {
    double delta = (b - a) / n;
    double sum = 0.0;
    
    for (long i = 0; i < n; i++) {
        double x = a + (i + 0.5) * delta;
        sum += f(x);
    }
    
    return sum * delta;
}

double riemann_sums_seq(){
	double a = 0.0;      // Lower bound
    double b = 1.0;      // Upper bound
    long long n = 1e10;   // Number of rectangles
    
    // Set number of threads (optional)
    omp_set_num_threads(1);
    
    printf("Computing Riemann Sum for f(x) = x^2 from %.2f to %.2f\n", a, b);
    printf("Number of rectangles: %lld\n", n);
    printf("Number of threads: %d\n\n", omp_get_max_threads());
    
    // Left endpoint method
    double start = omp_get_wtime();
    double result_left = riemann_sum_left(a, b, n);
    double end = omp_get_wtime();
    printf("Left Riemann Sum: %.10f (Time: %.6f seconds)\n", 
           result_left, end - start);
    
    // Right endpoint method
    start = omp_get_wtime();
    double result_right = riemann_sum_right(a, b, n);
    end = omp_get_wtime();
    printf("Right Riemann Sum: %.10f (Time: %.6f seconds)\n", 
           result_right, end - start);
    
    // Midpoint method
    start = omp_get_wtime();
    double result_mid = riemann_sum_midpoint(a, b, n);
    end = omp_get_wtime();
    printf("Midpoint Riemann Sum: %.10f (Time: %.6f seconds)\n", 
           result_mid, end - start);
    
    // Exact integral for f(x) = x^2 from 0 to 1 is 1/3
    printf("\nExact value (for x^2 from 0 to 1): %.10f\n", 1.0/3.0);
    return ((result_left + result_right + result_mid )/ 3.0); 
}

int main (int argc, char *argv[]) {
	double time_arr[5];
	int c = 0;
	printf("Total Threads Available: %d\n", omp_get_num_procs());
	printf("------Estimating value of PI and using it to find Integration of function (PI * x ^2)--------\n");
	printf("======SEQUENTIAL Execution=====\n");
	double start = omp_get_wtime();
	double piv = pi_estimate_seq();
	double intg = riemann_sums_seq();
	double stop = omp_get_wtime();
	printf("-------RESULTS for SEQUENTIAL EXECUTION--------\n");
	printf("Estimated PI Value: %.12f\n", piv);
	printf("Estimated Integral Value: %.12f\n", intg);
	printf("Total TIME taken: %.4f\n", (stop - start));
	time_arr[c++] = (stop-start);
	
	
	int threads[4] = {4, 8, 12, 16};
	
	for(int i = 0;i<4;i++){
		printf("======PARALLEL Execution with %d THREADS=====\n", threads[i]);
		omp_set_num_threads(threads[i]);
  omp_set_proc_bind(omp_proc_bind_true);
		double start = omp_get_wtime();
		double piv = pi_estimate_par(threads[i]);
		double intg = riemann_sums_par(threads[i]);
		double stop = omp_get_wtime();
		printf("-------RESULTS for PARALLEL EXECUTION with %d THREADS--------\n", threads[i]);
		printf("Estimated PI Value: %.12f\n", piv);
		printf("Estimated Integral Value: %.12f\n", intg);
		printf("Total TIME taken: %.4f\n", (stop - start));
		time_arr[c++] = (stop-start);
	}
	
	printf("\n============FINAL RESULTS===============\n");
	printf("MODE\tTHREADS\tTIME TAKEN\n");
	for(int i = 0;i<5;i++){
		if(i == 0){
		printf("SEQUENTIAL\t1\t%.4f\n", time_arr[0]);
		}
		else{
		printf("PARALLEL\t%d\t%.4f\n", threads[i-1], time_arr[i]);
		}
	}
    return 0;
}

