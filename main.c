#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void sobel(unsigned char **m, unsigned char **matrizFiltrada, int i, int col) {
    // Declaramos y llenamos matrices Sobel
    int c[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int f[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    int sumC = 0;
    int sumF = 0;
	
	for(int j=1; j<col-1; j++){
		sumC = ((c[0][0]*m[i-1][j-1]) + (c[0][1]*m[i-1][j]) + (c[0][2]*m[i-1][j+1]) + (c[1][0]*m[i][j-1]) + (c[1][1]*m[i][j]) + (c[1][2]*m[i][j+1]) + (c[2][0]*m[i+1][j-1]) + (c[2][1]*m[i+1][j]) + (c[2][2]*m[i+1][j+1]));
		sumF = ((f[0][0]*m[i-1][j-1]) + (f[0][1]*m[i-1][j]) + (f[0][2]*m[i-1][j+1]) + (f[1][0]*m[i][j-1]) + (f[1][1]*m[i][j]) + (f[1][2]*m[i][j+1]) + (f[2][0]*m[i+1][j-1]) + (f[2][1]*m[i+1][j]) + (f[2][2]*m[i+1][j+1]));
		matrizFiltrada[i][j] = (unsigned char) sqrt(pow(sumC,2)+pow(sumF,2));
	}
    
    
}

void exSimetricaSobel(unsigned char **m, int altura, int anchura){
    for(int i=1;i<=anchura;i++){
        m[0][i] = m[2][i];
    }
    for(int i=1;i<=anchura;i++){
        m[altura+1][i] = m[altura-1][i];
    }
    for(int i=1;i<=altura;i++){
        m[i][0] = m[i][2];
    }
    for(int i=1;i<=altura;i++){
        m[i][anchura+1] = m[i][anchura-1];
    }
    
    m[0][0] = m[2][2];
    m[0][anchura+1] = m[2][anchura-1];
    m[altura+1][0] = m[altura-1][2];
    m[altura+1][anchura+1] = m[altura-1][anchura-1];
}

void extensionSimetrica(unsigned char **m, int filas, int columnas){
    memcpy(m[1],m[3],columnas+4);
    memcpy(m[0],m[4],columnas+4);
    memcpy(m[filas+2],m[filas],columnas+4);
    memcpy(m[filas+3],m[filas-1],columnas+4);
    for(int i=0;i<=filas+3;i++){
        memcpy(&m[i][1],&m[i][3],1);
    }
    for(int i=0;i<=filas+3;i++){
        memcpy(&m[i][0],&m[i][4],1);
    }
    for(int i=0;i<=filas+3;i++){
        memcpy(&m[i][columnas+2],&m[i][columnas],1);
    }
    for(int i=0;i<=filas+3;i++){
        memcpy(&m[i][columnas+3],&m[i][columnas-1],1);
    }
}

unsigned char calculo(unsigned char *v){
    unsigned char min = 255;
    unsigned char max = 0;
    int minpos=0;
    int maxpos=0;

    for(int i=0; i<25 ; i++){
        if(v[i]<min){
            min = v[i];
            minpos = i;
        }
    }
    for(int i=0; i<25 ; i++){
        if(v[i]>max){
            max = v[i];
            maxpos = i;
        }
    }

    int sum = 0;
    // Sumar todos los valores excepto el mínimo y el máximo
    for (int i = 0; i < 25; i++) {
        if (i != minpos && i != maxpos)
        {
            sum += v[i];
        }
    }
    return sum/23;
}

void media(unsigned char **sim, unsigned char **resultado,int fil, int col, int inicio){
    unsigned char *vector = (unsigned char*) malloc (sizeof(unsigned char)*25);
	
    for(int i=2+inicio; i<fil+2+inicio ; i++){
        for(int j=2; j<col+2 ; j++){
            vector[0] = sim[i - 2][j - 2];
            vector[1] = sim[i - 2][j - 1];
            vector[2] = sim[i - 2][j];
            vector[3] = sim[i - 2][j + 1];
            vector[4] = sim[i - 2][j + 2];
            vector[5] = sim[i - 1][j - 2];
            vector[6] = sim[i - 1][j - 1];
            vector[7] = sim[i - 1][j];
            vector[8] = sim[i - 1][j + 1];
            vector[9] = sim[i - 1][j + 2];
            vector[10] = sim[i][j - 2];
            vector[11] = sim[i][j - 1];
            vector[12] = sim[i][j];
            vector[13] = sim[i][j + 1];
            vector[14] = sim[i][j + 2];
            vector[15] = sim[i + 1][j - 2];
            vector[16] = sim[i + 1][j - 1];
            vector[17] = sim[i + 1][j];
            vector[18] = sim[i + 1][j + 1];
            vector[19] = sim[i + 1][j + 2];
            vector[20] = sim[i + 2][j - 2];
            vector[21] = sim[i + 2][j - 1];
            vector[22] = sim[i + 2][j];
            vector[23] = sim[i + 2][j + 1];
            vector[24] = sim[i + 2][j + 2];
            resultado[i-2][j-2] = calculo(vector);
        }
    }
}


// main.exe filas columnas fichero opcion numero_hilos
int main(int argc, char *argv[]){
	int np, iam, filas, columnas, opcion;
	unsigned char **matriz;
	unsigned char **matrizFiltrada;
	FILE *f, *fout;
	clock_t inicio, fin;
	
	filas = atoi(argv[1]);
	columnas = atoi(argv[2]);
	opcion = atoi(argv[4]);
	np = atoi(argv[5]);
	
	//Dividimos la carga de trabajo
	int filasxbloque = filas/np;
	int resto = filas%np;
	int *vfilasxbloque = (int *) malloc (np*(sizeof(int)));
	for(int i=0 ; i<np-1 ; i++){
		vfilasxbloque[i] = filasxbloque;
	}
	vfilasxbloque[np-1] = filasxbloque+resto;
	// Guardamos los inicios en un vector
	int *vinicios= (int *) malloc (np*(sizeof(int)));
	for(int i=0 ; i<np ; i++){
		vinicios[i] = filasxbloque*i;
	}
	
	
	if(opcion==1){ //Filtrado por media 5x5
		// Abrimos el fichero de entrada
		f = fopen(argv[3], "rb");
		
		//Reservamos memoria para la matriz
		matriz = (unsigned char **) malloc ((filas+4)*sizeof(unsigned char*));
		for(int i=0; i<filas+4; i++){
			matriz[i] = (unsigned char*) malloc ((columnas+4)*sizeof(unsigned char));
		}
		
		// Guardamos los datos del archivo raw en la matriz
		for(int i=2; i<filas+2; i++){
			for(int j=2; j<columnas+2; j++){
				fread(&matriz[i][j], sizeof(unsigned char), 1, f);
			}
		}
		
		// Reservamos memoria para el resultado
		matrizFiltrada = (unsigned char **) malloc (filas*sizeof(unsigned char*));
		for(int i=0; i<filas; i++){
			matrizFiltrada[i] = (unsigned char*) malloc (columnas*sizeof(unsigned char));
		}
		
		// Realizamos la extension simetrica
		extensionSimetrica(matriz, filas, columnas);
		
		int i=0;
		inicio = clock();
		#pragma omp parallel num_threads(np) shared(matriz, matrizFiltrada, columnas, vinicios, vfilasxbloque) private(iam, i) default(none)
		{
			iam = omp_get_thread_num();
			media(matriz, matrizFiltrada, vfilasxbloque[iam], columnas, vinicios[iam]);
		}
		fin = clock();
		
		fout = fopen("ImagenFiltradaMedia5x5.raw", "wb");
		for(int i=0; i<filas ; i++){
			fwrite(matrizFiltrada[i], sizeof(unsigned char), columnas, fout);
		}
		
		printf("\nImagen filtrada por media con exito!\n\n");
	}
	else if(opcion==2){ //Filtrado por sobel
		
		// Abrimos el fichero de entrada
		f = fopen(argv[3], "rb");
		
		matriz = (unsigned char**) malloc (filas*sizeof(unsigned char*));
		for(int i=0; i<filas ; i++){
			matriz[i] = (unsigned char *) malloc (columnas*sizeof(unsigned char));
		}
		
		matrizFiltrada = (unsigned char**) malloc (filas*sizeof(unsigned char*));
		for(int i=0; i<filas ; i++){
			matrizFiltrada[i] = (unsigned char *) malloc (columnas*sizeof(unsigned char));
		}
		
		// Guardamos los datos del archivo raw en la matriz
		for(int i=0; i<filas; i++){
			for(int j=0; j<columnas; j++){
				fread(&matriz[i][j], sizeof(unsigned char), 1, f);
			}
		}
		inicio = clock();
		#pragma omp parallel num_threads(np) shared(matriz,matrizFiltrada,filas,columnas) private(iam) default(none)
		{
			iam = omp_get_thread_num();
			
			#pragma omp for schedule(static)
			for(int i=1; i<filas-1 ; i++)
			{
				sobel(matriz, matrizFiltrada, i, columnas);
			}
		}
		fin = clock();
		
		//Realizar extension simetrica
		exSimetricaSobel(matrizFiltrada, filas-2, columnas-2);
		
		fout = fopen("ImagenFiltradaSobel.raw", "wb");
		
		for(int i=0; i<filas ; i++){
			fwrite(matrizFiltrada[i], sizeof(unsigned char), columnas, fout);
		}

		printf("\nImagen filtrada por sobel con exito!\n\n");	
	}
	else if(opcion==3){ //Histograma
		// Reservamos memoria para la matriz
		matriz = (unsigned char**) malloc (filas*sizeof(unsigned char*));
		for(int i=0; i<filas ; i++){
			matriz[i] = (unsigned char*) malloc (columnas*sizeof(unsigned char));
		}
		// Guardamos la informacion en la matriz
		f = fopen(argv[3], "rb");
		for(int i=0 ; i<filas ; i++){
			for(int j=0; j<columnas ; j++){
				fread(&matriz[i][j], sizeof(unsigned char), 1, f);
			}
		}
		
		int histograma[256];
		for(int i=0; i<256; i++){
			histograma[i] = 0;
		}
		
		filasxbloque = 256/np;
		resto = 256%np;
		for(int i=0; i<np ; i++){
			vfilasxbloque[i] = filasxbloque;
		}
		for(int i=0; i<resto ; i++){
			vfilasxbloque[i]++;
		}
		
		for(int i=0; i<np ; i++){
			vinicios[i] = i*vfilasxbloque[i];
		}
		int maximo=0;
		int minimo=70000;
		int num_max = 0;
		int num_min = 0;
		
		
		inicio = clock();
		// Calculamos el histograma
		#pragma omp parallel num_threads(np) shared(filas, columnas, matriz, histograma) private(iam) default(none)
		{
			iam = omp_get_thread_num();
			
			#pragma omp for schedule(static)
			for(int i=0 ; i<256 ; i++){
				for(int j=0 ; j<filas ; j++){
					for(int k=0 ; k<columnas ; k++){
						if(matriz[j][k] == i){
							histograma[i]++;
						}
					}
				}
			}			
		}
		
		// Calculamos el max y el minimo
		#pragma omp parallel num_threads(np) shared(histograma, vinicios, vfilasxbloque) private(iam) reduction(max:maximo) reduction(min:minimo) default(none)
		{
			iam = omp_get_thread_num();
			for(int i=vinicios[iam] ; i<vfilasxbloque[iam]+vinicios[iam] ; i++){
				if(histograma[i]>maximo){
					maximo = histograma[i];
				}
				if(histograma[i]<minimo){
					minimo = histograma[i];
				}
			}
		}
		
		// Calculamos el numero de veces que se repite el numero minimo y maximo
		#pragma omp parallel num_threads(np) shared(histograma, maximo, minimo) reduction(+:num_max) reduction(+:num_min) default(none)
		{
			#pragma omp for schedule(static)
			for(int i=0; i<256 ; i++){
				if(histograma[i] == maximo){
					num_max++;
				}
				if(histograma[i] == minimo){
					num_min++;
				}
			}
		}
		fin = clock();
		
		//imprimimos en fichero
		fout = fopen("histograma.txt", "wb");
		for(int i=0 ; i<256 ; i++){
			fprintf(fout, "%d -> %d\n",i ,histograma[i]);
		}
		fprintf(fout, "\n\nEl numero maximo es %d\n",maximo);
		fprintf(fout, "El numero minimo es %d\n",minimo);
		fprintf(fout, "\n\nEl numero maximo se ha repetido %d veces\n",num_max);
		fprintf(fout, "El numero minimo se ha repetido %d veces\n",num_min);
		
		printf("\nEl histograma has sido calculado con exito\n");
	}
	else{
		printf("Error, opcion incorrecta\n");
		return 1;
	}
	
	fclose(f);
	fclose(fout);
	
	double tiempo;
	tiempo = (double) (fin-inicio)/CLOCKS_PER_SEC;
	
	printf("El proceso ha durado %f s\n", tiempo);
	return 0;
}