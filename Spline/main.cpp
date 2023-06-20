#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <locale>
#include<fstream>
#include <string>
#include <cmath>
#define KOL 51//������������ ����� ����� �����

using namespace std;

const double PI = 3.141592653589793;
double Xravn[KOL],Xcheb[KOL],Yravn[KOL],Ycheb[KOL],Xcheb1[KOL],Ycheb1[KOL];
double a = 0.0, b = 1.0;//�������� ����������������

//��������������� �������
double f(double x)
{
    return x*x-asin(x-0.2);
}

//������������ ����������
double Maxf(int N, double Ax[],double Bx[])
{
    double maxb=0;
    for (int j=0; j<N; j++)
        {
        if (maxb<abs(Ax[j]-Bx[j])) maxb=abs(Ax[j]-Bx[j]);
        }
return maxb;
}

//����������� ����������� �����
void GetRavn(int n, double a, double b)
{
    double h = (b-a)/n;
    for(int i=0;i<=n;i++)
    {
        Xravn[i]=a+h*i;
        Yravn[i]=f(Xravn[i]);
    }
}

//����������� ����������� ����� (�� b �� a)
void GetCheb(int n, double a, double b)
{
    for(int i=0;i<=n;i++)
        {
            Xcheb[i]=0.5*((b-a)*cos((2.0*i+1.0)*PI/(2.0*(n+1.0)))+b+a);
            Ycheb[i]=f(Xcheb[i]);
        }
}

//����������� ����������� ����� (�� a �� b)
void GetCheb1(int n, double a, double b)
{
    for(int i=0;i<=n;i++)
        {
            Xcheb1[i]= Xcheb[n-i];
            Ycheb1[i]= Ycheb[n-i];
        }
}

//���������� ������ ������� n
double* my_S3(int n, double x[], double y[])
{
    double *h = new double[n+1];
    double *l = new double[n+1];
    double *delta = new double[n+1];
    double *lambda = new double[n+1];
    double *c = new double[n+1];
    double *d = new double[n+1];
    double *bb = new double[n+1];
    for(int k=1; k<=n; k++)
    {
       h[k] = x[k] - x[k-1];//��� �����, ����� ���� ����������
       l[k] = (y[k] - y[k-1])/h[k];//���������� ��������
    }
   //������ ������� ��������, delta � lambda - ����������� ������������
   delta[1] = - h[2]/(2*(h[1]+h[2]));
   lambda[1] = 1.5*(l[2] - l[1])/(h[1]+h[2]);
   for(int k=3; k<=n; k++)
    {
      delta[k-1] = - h[k]/(2*h[k-1] + 2*h[k] + h[k-1]*delta[k-2]);
      lambda[k-1] = (3*l[k] - 3*l[k-1] - h[k-1]*lambda[k-2])/(2*h[k-1] + 2*h[k] + h[k-1]*delta[k-2]);
    }
   c[0] = 0;
   c[n] = 0;
   //���������� ��������� ������������ �������
   for(int k=n; k>=2; k--){
      c[k-1] = delta[k-1]*c[k] + lambda[k-1];
   }
   for(int k=1; k<=n; k++){
      d[k] = (c[k] - c[k-1])/(3*h[k]);
      bb[k] = l[k] + (2*c[k]*h[k] + h[k]*c[k-1])/3;
   }

    int NN = 10*n;//����� ����������� �����
    double hh = (b-a)/(10*n);
    double *Xp=new double[NN+1];//����������� �����
    double *S3=new double[NN+1];//�������� ������� 3-�� �������
    //������� ��� ����������� ����� �������
    for (int i=0; i<=NN; i++)
    {
        Xp[i]=a+hh*i;
        int k;
        for(k=1; k<=n; k++)//������� k, ��� ������� Xp[i] ����� � [x_k-1; x_k]
            if(Xp[i]>=x[k-1] && Xp[i]<=x[k])  break;
        if (Xp[i]<x[0]) k=1;//��� ������ ������������� �����
        if (Xp[i]>x[n]) k=n;//��� ������ ������������� ������
        //k - ����� ������� �������
        S3[i] = y[k] + bb[k]*(Xp[i]-x[k]) + c[k]*pow(Xp[i]-x[k], 2) + d[k]*pow(Xp[i]-x[k], 3);
    }

    delete[] Xp;
    delete[] h;delete[] l;delete[] delta;delete[] lambda;delete[] c; delete[] d;delete[] bb;
    return S3;
    }


    //������� ������� ������������ ������� ������
    //��� ������� ������� ��� ����������� �������
    double* Gauss(int N, double **A,double B[])
{
    int n=N, i, j, k;
    double d, s,er[N];
    double *X = new double[n];//������� � ��������

    //��������� ����� ����� ������� ������� �� ��� ����� ����
    if (A[0][0]==0)
        {
         double temp = B[2]; B[2]=B[0]; B[0] = temp;
         for (i=0;i<N;i++)
            {
            temp = A[2][i];A[2][i]=A[0][i]; A[0][i]=temp;
            }
        }

    for (k = 0; k < n; k++) // ������ ���
    {
    for (j = k + 1; j < n; j++)
        {
        if (A[j][k]!=0)
            {
            d = A[j][k] / A[k][k];
            for (i = k; i < n; i++)
                {
                A[j][i] = A[j][i] - d * A[k][i];
                }
            B[j] = B[j] - d * B[k];
            }
        }
    }

for (k = n-1; k >= 0; k--) // �������� ���
    {
    d = 0;
    for (j = k + 1; j < n; j++)
        {
        s = A[k][j] * X[j];
        d = d + s;
        }
    X[k] = (B[k] - d) / A[k][k];
    }
   return X;
}

//������������ ������
double* my_S2(int n, double x[], double y[])
{
    int N =3*n;//����������� ����
    double *B = new double[N];//������� ��������� ������
    double **A = new double*[N];
    for (int i = 0; i < N; i++)
        A[i] = new double[N];
    double *coef = new double[N];//������������ �������

   //��������� ������� ������� ������
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++) A[i][j]=0;
    }

   //��������� ������� ������� �������� ��������
    for (int i=0; i<n; i++)
    {
        A[3*i][3*i]=x[i]*x[i];A[3*i][3*i+1]=x[i];A[3*i][3*i+2]=1;
        A[3*i+1][3*i]=x[i+1]*x[i+1];A[3*i+1][3*i+1]=x[i+1];A[3*i+1][3*i+2]=1;
        A[3*i+2][3*i]=2*x[i+1];A[3*i+2][3*i+1]=1;
        if (3*i+3<N)
        {
            A[3*i+2][3*i+3]=-2*x[i+1];A[3*i+2][3*i+4]=-1;
        }
        B[3*i]=y[i];B[3*i+1]=y[i+1];B[3*i+2]=0;
    }
    //��������� ������� � ����� ������� (����������� ��������� � ������������ ����������� �������)
    B[N-1]=(y[n]-y[n-1])/(x[n]-x[n-1]);

    //������ ������� ������� ������
    coef=Gauss(N,A,B);

    int NN = 10*n;
    double hh = (b-a)/(10*n);
    double *Xp=new double[NN+1];//����������� �����
    double *S2=new double[NN+1];//�������� ������� 2-�� �������

    for (int i=0; i<=NN; i++)
    {
        Xp[i]=a+hh*i;
        int k;
        for(k=1; k<=n; k++)//������� k, ��� ������� Xp[i] ����� � [x_k-1; x_k]
            if(Xp[i]>=x[k-1] && Xp[i]<=x[k])  break;
        if (Xp[i]<x[0]) k=1;
        if (Xp[i]>x[n]) k=n;
        k--;
        S2[i] = coef[3*k]*Xp[i]*Xp[i] + coef[3*k+1]*Xp[i]+coef[3*k+2];
    }

    delete [] A;delete [] B;delete [] coef;delete [] Xp;
    return S2;

}

//�������� ������
double* my_S1(int n, double x[], double y[])
{
    int NN = 10*n;
    double hh = (b-a)/(10*n);
    double *Xp=new double[NN+1];//����������� �����
    double *S1=new double[NN+1];//�������� ������� 1-�� �������

    for (int i=0; i<=NN; i++)
    {
        Xp[i]=a+hh*i;//�����
        int k;
        for(k=1; k<=n; k++)//������� k, ��� ������� Xp[i] ����� � [x_k-1; x_k]
            if(Xp[i]>=x[k-1] && Xp[i]<=x[k])  break;
        if (Xp[i]<x[0]) k=1;
        if (Xp[i]>x[n]) k=n;
        k--;//k - ��������� ���� �������
        //������ �� ���� ������
        S1[i] = y[k]+(y[k+1]-y[k])*(Xp[i]-x[k])/(x[k+1]-x[k]);
    }
    delete [] Xp;
    return S1;
}


//����������������
void Interpool(int n, char s)
{
    int NN = 10*n+1;//����� ����������� �����
    double h = (b-a)/(10*n);
    double *X=new double[NN];//����������� �����
    double *FX=new double[NN];//�������� ������� � ����������� ������
    double *S11=new double[NN];//�������� ������� 1-�� �������
    double *S1opt=new double[NN];//�������� ������� 1-�� ������� �� ����������� �����
    double *S2=new double[NN];//�������� ������� 2-�� �������
    double *S2opt=new double[NN];//�������� ������� 2-�� ������� �� ����������� �����
    double *S3=new double[NN];//�������� ������� 3-�� �������
    double *S3opt=new double[NN];//�������� ������� 3-�� ������� �� ����������� �����
    char s1[] = "outputX.txt";
    s1[6]=s;
    ofstream fout(s1);//����� � ����
    S11=my_S1(n,Xravn,Yravn);//�������� �� ����������� �����
    S1opt=my_S1(n,Xcheb1,Ycheb1);//�������� �� ����������� �����
    S2=my_S2(n,Xravn,Yravn);//���������� �� ����������� �����
    S2opt=my_S2(n,Xcheb1,Ycheb1);//���������� �� ����������� �����
    S3=my_S3(n,Xravn,Yravn);//���������� �� ����������� �����
    S3opt=my_S3(n,Xcheb1,Ycheb1);//���������� �� ����������� �����


    //����� � ����
    for (int i=0; i<NN; i++)
    {
        X[i]=a+h*i;
        FX[i]=f(X[i]);
        fout << X[i] << '\t' << FX[i] << '\t'<< S11[i] <<'\t' <<S1opt[i] << '\t'<< S2[i] <<'\t' <<S2opt[i]<< '\t'<< S3[i] <<'\t' <<S3opt[i]<< endl;
    }
    fout.close();
    //����� ������� �� �������
    cout << n << '\t' << NN << '\t'<< Maxf(NN,FX,S11) <<'\t' <<Maxf(NN,FX,S1opt)<< '\t'<< Maxf(NN,FX,S2) <<'\t' <<Maxf(NN,FX,S2opt)<<'\t'<< Maxf(NN,FX,S3) <<'\t' <<Maxf(NN,FX,S3opt)<<endl;
    //������� ������
    delete[] X;delete[] FX;delete[] S11;delete[] S1opt;delete[] S2;delete[] S2opt;delete[] S3;delete[] S3opt;
}


int main()
{
    setlocale(LC_ALL, "russian");
    int U[6] ={3, 10, 20, 30, 40, 50};//����� �������� ����������������
    int k = 6;
    char s[] = "123456";
    cout << "��������� �������� ��� ���������� ���������� ����� ����������������"<< endl;
    cout <<"\t\t\t"<< "��������" << "\t\t\t"<<"����������" << "\t\t\t" <<"����������"<<endl;
    for (int i=0; i<k; i++)
    {
        GetRavn(U[i],a,b);//�������� ����������� �����
        GetCheb(U[i],a,b);//�������� ����������� ����� (� �������� �������)
        GetCheb1(U[i],a,b);//�������� ����������� ����� (�� �������)
        Interpool(U[i],s[i]);//����� ������������
    }

    return 0;
}
