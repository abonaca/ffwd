#include "ffwd.h" 

double compare(double **x_obs, double **x_mod, double **err_obs, double ***sigma_inv, double *logdet, int N, int k, int K)
{//Compare two pointlike datasets
	
	double deltap, log_pdf, temp, temp2, limit, std, mean, med;
	int i, l, m1, m2;
	
	// 1d array
	double *log_pdfs = (double*) malloc(N*sizeof(double));
	double *delta = (double*) malloc(k*sizeof(double));

	// calculate model likelihood for each observed star
	for(i=0; i<N; i++){
		
		// log likelihood when we only care about diagonal elements
		deltap = 1e-100;
		for(l=0; l<K; l++){
			temp2 = 0.;
			for(m1=0; m1<k; m1++){
				if(err_obs[m1][i]>0){
					temp = 0.;
					for(m2=0; m2<k; m2++){
						delta[m2] = x_obs[m2][i] - x_mod[m2][l];
						temp += sigma_inv[i][m1][m2]*delta[m2];
					}
					temp2 += delta[m1]*temp;
					//printf("%lf %lf %lf \t", delta[m1], x_obs[m1][i], x_mod[m1][l]);
				}
			}
			deltap += exp(-0.5*temp2);
		}
		//printf("%e\t%e\n", deltap, logdet[i]);
		
		// logdet = logdet - 2*log(1./K * pow(2*pi, -k/2.)), where k number of data dimensions for a given star
		log_pdfs[i] = log(deltap) - 0.5*logdet[i];
	}
	
	// mean
	mean = 0.;
	for(i=0; i<N; i++)
		if(isfinite(log_pdfs[i]))
			mean += log_pdfs[i];
	mean /= (double)N;

	// std
	std = 0.;
	for(i=0; i<N; i++){
 		//printf("%lf %lf %lf %lf\n", x_obs[0][i], log_pdfs[i], mean, pow((log_pdfs[i] - mean),2) );
		if(isfinite(log_pdfs[i]))
			std += pow((log_pdfs[i] - mean),2);
// 		printf("%lf %lf\n", std, pow((log_pdfs[i] - mean),2));
	}
	std = sqrt(std/(double)N);
	
	// median
	med = median(N, log_pdfs);
	
	// remove outliers
	limit = med - 3.*std;
	
	log_pdf = 0.;
	for(i=0; i<N; i++){
		if((log_pdfs[i] >= limit) && (isfinite(log_pdfs[i])))
			log_pdf += log_pdfs[i];
	}
// 	printf("%lf %lf %lf %lf %lf\n", log_pdfs[0], log_pdf, med, std, limit);
	
	// Free memory
	free(log_pdfs);
	free(delta);
	
	return log_pdf;
}

double median(int n, double x[]) {
    double temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

// void sort(unsigned long n, double arr[])
// {
// 	unsigned long i,ir=n,j,k,l=0,*istack;
// 	int jstack=0;
// 	double a,temp;
// 
// 	istack=lvector(0,NSTACK-1);
// 	for (;;) {
// 		if (ir-l < M) {
// 			for (j=l;j<ir;j++) {
// 				a=arr[j];
// 				for (i=j-1;i>=l;i--) {
// 					if (arr[i] <= a) break;
// 					arr[i+1]=arr[i];
// 				}
// 				arr[i+1]=a;
// 			}
// 			if (jstack == 0) break;
// 			ir=istack[jstack--];
// 			l=istack[jstack--];
// 		} else {
// 			k=(l+ir) >> 1;
// 			SWAP(arr[k],arr[l+1])
// 			if (arr[l] > arr[ir]) {
// 				SWAP(arr[l],arr[ir])
// 			}
// 			if (arr[l+1] > arr[ir]) {
// 				SWAP(arr[l+1],arr[ir])
// 			}
// 			if (arr[l] > arr[l+1]) {
// 				SWAP(arr[l],arr[l+1])
// 			}
// 			i=l+1;
// 			j=ir;
// 			a=arr[l+1];
// 			for (;;) {
// 				do i++; while (arr[i] < a);
// 				do j--; while (arr[j] > a);
// 				if (j < i) break;
// 				SWAP(arr[i],arr[j]);
// 			}
// 			arr[l+1]=arr[j];
// 			arr[j]=a;
// 			jstack += 2;
// 			if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
// 			if (ir-i+1 >= j-l) {
// 				istack[jstack]=ir;
// 				istack[jstack-1]=i;
// 				ir=j-1;
// 			} else {
// 				istack[jstack]=j-1;
// 				istack[jstack-1]=l;
// 				l=i;
// 			}
// 		}
// 	}
// 	free_lvector(istack,0,NSTACK-1);
// }
