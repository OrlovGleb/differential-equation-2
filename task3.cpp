#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <math.h>
//#include <сmath>
#include <malloc.h>
const static double  pi=3.14159265358979323846;



void progonka(unsigned int n, double *Ab,double *Ac,double *At,double *BB, double *X){
    double y;
    double *a;
    double *B;
    unsigned m=n+1;
    a = static_cast<double *> (malloc((m)  * sizeof(double)));
    B = static_cast<double *> (malloc((m)  * sizeof(double)));
     // double a[n+1];
     // double B[n+1];
    y = Ac[0];
  a[0] = -At[1] / y;
  B[0] = BB[0] / y  ;
  for (unsigned i = 1; i < n; i++) {
    y = Ac[i] + Ab[i] * a[i - 1];
    a[i] = -At[i] / y;
    B[i] = (BB[i]- Ab[i] * B[i - 1]) / y;
  }

  X[n] = (BB[n] - Ab[n] * B[n - 1]) / (Ac[n] + Ab[n] * a[n - 1]);
  for (unsigned i = n; i >= 1; i--) {
    X[i - 1] = a[i - 1] * X[i] + B[i - 1];
  }
  free(a);
  free(B);
}


int main(int argc, char *argv[])
{
    
    unsigned long int t1;
    unsigned long int n1;
    n1 = (strtoul(argv[1], NULL, 0));
    t1 = (strtoul(argv[2], NULL, 0));
    unsigned  int t=static_cast<unsigned int>(t1);
    unsigned  int n=static_cast<unsigned int>(n1);
    if(argc<2){return 0;}
    
    
    
    
    unsigned j;
    unsigned i;
    unsigned k;

    unsigned i1;

    //unsigned int n=25;
    unsigned  m=n+1;
    //unsigned int t=320;
    unsigned  q=t+1;
    double err=0;
    double h;
    double tau;
    tau=pow(t,-1);
    h=pow(n,-1);
   // printf("%lf", h);
   // printf("\n");
    
    double *AA;
    double *F;
    F = static_cast<double *> (malloc((m)*(m)*(q)  * sizeof(double)));
    AA = static_cast<double *> (malloc((m)*(m)*(q)  * sizeof(double)));
    //зануляем матрицу
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
            for (k = 0; k <= t; k++){
              AA[k*(n+1)*(n+1)+i*(n+1)+j] = 0;
                
            }           
        }     
    }
    

    double *A1;
    A1 = static_cast<double *> (malloc((m)*(m)  * sizeof(double))); //матрица для слоя 1/2
    double *Ac;
    double *At;
    double *Ab;
    double *BB;
    double *X;
    Ac = static_cast<double *> (malloc((m)  * sizeof(double)));
    At = static_cast<double *> (malloc((m)  * sizeof(double))); 
    Ab = static_cast<double *> (malloc((m)  * sizeof(double)));
    BB = static_cast<double *> (malloc((m)  * sizeof(double)));
    X = static_cast<double *> (malloc((m)  * sizeof(double)));
    //double X1[m+1]; //матрица иксов
    //double Y1[m+1]; //матрица игреков
    //double T1[m+1]; //матрица т
    //scanf("%d", &a1);
    
    //заполняем матрицу иксов
    //for (i1 = 0; i1 < n+1; i1++) {
    //    X1[i1] = i1*h;
    //}
    //заполняем матрицу игреков
    //for (i1 = 0; i1 < n+1; i1++) {
    //    Y1[i1] = i1*h;
        
    //}
    //заполняем матрицу т
    //for (i1 = 0; i1 < t+1; i1++) {
    //    T1[i1] = i1*tau;
        
    //}
    // заполняем нулевой уровень
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
            
              AA[0*(n+1)*(n+1)+i*(n+1)+j] = (sin(j*h*pi))*(sin(i*h*pi));
              
                       
        }     
    }
    //заполняем F
    for (i = 0; i < n+1; i++) {
        for (j = 0; j <= n; j++) {
            for (k = 0; k <= t; k++){
               F[k*(n+1)*(n+1)+i*(n+1)+j] = -(exp(-pi*pi*k*tau))*(-pi*pi)*(sin(j*h*pi))*(sin(i*h*pi));
            }
        }     
    }



    for (k = 1; k <= t; k++){
    for (j = 0; j <= n; j++){

  //  for (i1 = 0; i1 < n+1; i1++) {
  //      for (j1 = 0; j1 < n+2; j1++) {
  //          A[i1*(n+2)+j1]=0;
  //      }     
  //  } 
    //задаем значения на границах
    Ac[0] = 1;//A[0]=1;
    BB[0] = 0;//A[n+1]=0;
    Ac[n] = 1;//A[n*(n+2)+n]=1;
    BB[n] = 0;//A[n*(n+2)+n+1]=0;
    //заполняем оставшуюся часть матрицы
    for (i=1; i<=n-1; i++){
        Ab[i]=tau;//  A[(n+2)*i+i-1]=tau;
        Ac[i]=-(h*h+2*tau);//A[(n+2)*i+i]=-(h*h+2*tau);
        At[i]=tau;//A[(n+2)*i+i+1]=tau;
        BB[i]=-h*h*AA[(k-1)*(n+1)*(n+1)+i*(n+1)+j]-(tau/2)*h*h*F[(k-1)*(n+1)*(n+1)+i*(n+1)+j];//A[(n+2)*i+n+1]=-h*h*AA[i][j][k-1]-(tau/2)*h*h*F[i][j]; 
    }
    
    //решаем слу
    //f(n+1, m+1, A, X);
    progonka(n,Ab,Ac,At,BB,X);
    for(i1=0;i1<=n;i1++){
        A1[(n+1)*j+i1] = X[i1];//заполняем матрицу промежуточного уровня
    } 
    }
    //for(i1=0;i1<n+1;i1++){
    //    A1[(n+1)*0+i1] = 0;//заполняем грaницу промежуточной матрицы
    //} 
    //for(i1=0;i1<n+1;i1++){
    //    A1[(n+1)*n+i1] = 0;//заполняем грaницу промежуточной матрицы
    //} 

    
    for (i = 0; i <= n; i++){

 //   for (i1 = 0; i1 < n+1; i1++) {
 //       for (j1 = 0; j1 < n+2; j1++) {
 //           A[i1*(n+2)+j1]=0;
 //       }     
 //  }
    //задаем значения на границах
    Ac[0] = 1;//A[0]=1;
    BB[0] = 0;//A[n+1]=0;
    Ac[n] = 1;//A[n*(n+2)+n]=1;
    BB[n] = 0;//A[n*(n+2)+n+1]=0;
    //заполняем оставшуюся часть матрицы
    for (i1=1; i1<=n-1; i1++){
        Ab[i1]=tau;//  A[(n+2)*i+i-1]=tau;
        Ac[i1]=-(h*h+2*tau);//A[(n+2)*i+i]=-(h*h+2*tau);
        At[i1]=tau;//A[(n+2)*i+i+1]=tau;
        BB[i1]=-h*h*A1[(n+1)*i1+i]-(tau/2)*h*h*F[k*(n+1)*(n+1)+i*(n+1)+i1];//A[(n+2)*i+n+1]=-h*h*AA[i][j][k-1]-(tau/2)*h*h*F[i][j]; 
        /*
        A[(n+2)*i1+i1-1]=tau;
        A[(n+2)*i1+i1]=-(h*h+2*tau);
        A[(n+2)*i1+i1+1]=tau;
        A[(n+2)*i1+n+1]=-h*h*A1[(n+1)*i1+i]-(tau/2)*h*h*F[i][i1]; */
    }

    //решаем слу
    //f(n+1, m+1, A, X);
    progonka(n,Ab,Ac,At,BB,X);
    for(i1=0;i1<=n;i1++){
        AA[k*(n+1)*(n+1)+i*(n+1)+i1] = X[i1];//заполняем матрицу следующего уровня
    }
    }
    // for(i1=0;i1<n+1;i1++){
    //    AA[0][i1][k] = 0;//заполняем границу матрицы следующего уровня
    //}
    // for(i1=0;i1<n+1;i1++){
    //    AA[n][i1][k] = 0;//заполняем границу матрицы следующего уровня
    //}
    //free(A1);
    //free(Ac);
    //free(At);
    //free(Ab);
    //free(BB);
    //free(X);

    }









    //находим ошибку
    //if (t==1){
    // for (i = 0; i < n+1; i++){
    //  if(fabs(X[i]-(a1*sin(pi*X1[i]/2)*sin(pi*X1[i]/2)))>err) 
    //  err=fabs(X[i]-(a1*sin(pi*X1[i]/2)*sin(pi*X1[i]/2))) ;
    //}
    //}

    //printf("%lf",err);

    //открываем файл и записываем в него ответ        
    FILE *f;
    //FILE *fx;
    //FILE *fy;
    //FILE *ft;
    //printf("\n");
    //printf("%.3lf", n*h);
    //printf("\n");
    f=fopen("x.txt", "w");
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
            for (k = 0; k <= t; k++){
                fprintf(f,"%lf", i*h);
                fprintf(f," ");
                fprintf(f,"%lf", j*h);
                fprintf(f," ");
                fprintf(f,"%lf", k*tau);
                fprintf(f," ");
                fprintf(f,"%lf", AA[k*(n+1)*(n+1)+i*(n+1)+j]);
                fprintf(f,"\n");
                if(err<fabs(AA[k*(n+1)*(n+1)+i*(n+1)+j] - (exp(-pi*pi*k*tau))*(sin(j*h*pi))*(sin(i*h*pi)) )){
                    err = fabs(AA[k*(n+1)*(n+1)+i*(n+1)+j] - (exp(-pi*pi*k*tau))*(sin(j*h*pi))*(sin(i*h*pi)) );
                }
            }           
        }     
    }
    FILE *fx;
    FILE *ft;
    fx=fopen("datax.txt", "w"); 
    for (i = 0; i <= n; i++){
        fprintf(fx,"%lf", i*h);
        fprintf(fx," ");
        fprintf(fx,"%lf", AA[(t/2)*(n+1)*(n+1)+i*(n+1)+n/2]);
        fprintf(fx,"\n");
    } 
    ft=fopen("datat.txt", "w"); 
    for (k = 0; k <= t; k++){
        fprintf(ft,"%lf", k*tau);
        fprintf(ft," ");
        fprintf(ft,"%lf", AA[k*(n+1)*(n+1)+(n/2)*(n+1)+n/2]);
        fprintf(ft,"\n");
    } 


    printf("%e", err);
    printf("\n");
    //теперь решаем систему для другого равенства, 
    //чтобы аппроксимировать результат, имея частное решение
    
    
    
    fclose(f);
    fclose(fx);
    fclose(ft);
    free(Ac);
    free(Ab);
    free(At);
    free(A1);
    free(BB);
    free(X);
    return 0;
}