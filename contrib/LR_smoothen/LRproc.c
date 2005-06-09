#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
				
void cubgcv_(double *x, double *f, double *df, int *N, 
	double *y, double *c, int *ic, double *var, int *job, 
	double *se, double *wk, int *ier);

int getnocol(int argc, char *argv[]){
	int nocol;
	if (argc!=2) {
		printf("Usage: %s <number_of_columns>\n", argv[0]);
		exit(EXIT_FAILURE);
	} else {
		nocol = strtol(argv[1], (char **)NULL, 10);
	}

	if (nocol<=0){
		printf("illegal value for number of columns: %s\n", 
				argv[1]);
		exit(EXIT_FAILURE);
	} else {
		return nocol;
	}
}	

void spsmoothen(double **cols, int nocol, int nlinecount){
	int i, j;
	int N, ic, job, ier;
	double *y, *df, *c, *wk, var, *se;

	N = nlinecount;
	ic = nlinecount-1;
	
	y = malloc(N*sizeof(double)); df = malloc(N*sizeof(double));
	se = malloc(N*sizeof(double));
	wk = malloc((7*(N+2))*sizeof(double));
	c = malloc((ic*3)*sizeof(double));
	
	for(j=1; j<nocol; j++){
		job = 1; var = -1;
		for (i=0; i<N; i++) {
			df[i] = 1.0;
		}
		cubgcv_(cols[0], cols[j], df, &N, y, c, 
			&ic, &var, &job, se, wk, &ier);
		for (i=0; i<N; i++) {
			cols[j][i] = y[i];
		}
	}

	free(y); free(df); free(se); free(wk); free(c);
}	

int main(int argc, char *argv[]){
	int nocol;
	char line[65536];
	int nline, nlinecount, nlineinc;
	double **cols;
	int i, j;

	nocol = getnocol(argc, argv);
	cols = malloc(nocol*sizeof(double *));
	nline = 10; nlineinc = 10;
	for (i=0; i<nocol; i++){
		cols[i] = malloc(nline*sizeof(double));
	}

	
	nlinecount = 0;
	while(  fgets(line, 65536, stdin) != NULL ){
		if (line[0] == '#') continue;
		
		nlinecount++;
		if (nlinecount>nline) {
			nline += nlineinc;
			for (i=0; i<nocol; i++){
				cols[i] = realloc(cols[i],
						  nline*sizeof(double));
			}
		}
		

		cols[0][nlinecount-1] = strtod(strtok(line, " "), NULL);
		for (i=1; i<nocol; i++){
			cols[i][nlinecount-1] = strtod(strtok(NULL, " "), NULL);
		}
	}
	
	spsmoothen(cols, nocol, nlinecount); 
	
	for (i=0; i<nlinecount; i++){
		for (j=0; j<nocol; j++){
			printf("%e ", cols[j][i]);
		}
		printf("\n");
	}

	return 0;
}

