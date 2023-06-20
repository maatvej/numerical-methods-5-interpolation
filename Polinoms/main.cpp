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
double Xravn[KOL],Xcheb[KOL],Yravn[KOL],Ycheb[KOL],Xcheb1[KOL],Ycheb1[KOL],Xravn1[KOL],Yravn1[KOL];
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
void GetRavn1(int n, double a, double b)
{
    for(int i=0;i<=n;i++)
        {
            Xravn1[i]= Xravn[n-i];
            Yravn1[i]= Yravn[n-i];
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


//������� �������� (n-����� �����)
//x_arr[] � y_arr[] - ����� ����������������
double my_lagrange(double x, int n, double x_arr[], double y_arr[])
{

    double sum = 0;
    for (int i = 0; i < n; i++){

        double l = 1;//��������� �������� ������������ lj
        for (int j = 0; j < n; j++)
            if (j != i){
                l = l* (x - x_arr[j]) / (x_arr[i] - x_arr[j]);

            }

        sum += y_arr[i] * l;//����� lj*yj
    }
    return sum;
}

double my_newton(double x, int n, double x_arr[], double y_arr[])
{

    double sum = y_arr[0];
    for(int i = 1; i < n; ++i)//�� ���� �����
        {
        double F = 0;//i-�� ��������� �������� j ������������
        for(int j = 0; j <= i; ++j)
            {
            double den = 1;
            for(int k = 0; k <= i; ++k)
                if (k != j)
                //������������ �������� ��������� ����� ������� i ������
                    den *= (x_arr[j] - x_arr[k]);
            F += y_arr[j]/den;
        }
        //����� ������� ������������ (x-x[0])...(x-x[i-1])
        for(int k = 0; k < i; ++k)
            F *= (x - x_arr[k]);
        sum += F;//����� ��������� i-�� ���������
        }
    return sum;
}


//����������������
void Interpool(int n, char s)
{
    int N = 10*n+1;
    double h = (b-a)/(10*n);

    double *X=new double[N];//����������� �����
    double *FX=new double[N];//�������� ������� � ����������� ������
    double *LX=new double[N];//�������� �������� �������� � ����������� ������
    double *NX=new double[N];//�������� �������� ������� � ����������� ������
    double *LXopt=new double[N];//�������� �������� �������� � ����������� ������(�� ����������� ������)
    double *NXopt=new double[N];//�������� �������� ������� � ����������� ������(�� ����������� ������)
    char s1[] = "outputX.txt";
    s1[6]=s;
    ofstream fout(s1);//������� ������� ����������� ���� ��� ���������� ��������
    for (int i=0; i<N; i++)
    {
        X[i]=a+h*i;//����� �� x
        FX[i]=f(X[i]);//������ �������� �������
        LX[i]=my_lagrange(X[i],n+1,Xravn,Yravn);//������� �������� �� ����������� �����
        LXopt[i]=my_lagrange(X[i],n+1,Xcheb,Ycheb);//������� �������� �� ����������� �����
        //����� ���������� ����� ������� ������� ������������
        if (X[i]<(a+b)/2) //��� ������ �������� ������� 1-��
        {
            NX[i]=my_newton(X[i],n+1,Xravn,Yravn);//����������� �����
            NXopt[i]=my_newton(X[i],n+1,Xcheb1,Ycheb1);//����������� �����
        }
        else //��� ������ �������� ������� 2-��
        {
            NX[i]=my_newton(X[i],n+1,Xravn1,Yravn1);
            NXopt[i]=my_newton(X[i],n+1,Xcheb,Ycheb);
        }
        //����� � ����
        fout << X[i] << '\t' << FX[i] << '\t'<< LX[i] << '\t' <<LXopt[i] << '\t' << NX[i] << '\t' <<NXopt[i]<< endl;
    }
    fout.close();
    //����� ������� �� �������
    cout << n << '\t' << N << '\t'<< Maxf(N,FX,LX) << '\t' <<Maxf(N,FX,LXopt)<< "\t\t\t";
    cout << n << '\t' << N << '\t'<< Maxf(N,FX,NX) << '\t' <<Maxf(N,FX,NXopt)<< endl;
    //������������ ������
    delete[] X;delete[] FX;delete[] LX;delete[] NX;delete[] LXopt;delete[] NXopt;
}



int main()
{
    setlocale(LC_ALL, "russian");
    int U[6] ={3, 10, 20, 30, 40, 50};
    int k = 6;
    char s[] = "123456";
    cout << "��������� ����������������� �������� ��� ���������� ���������� ����� ����������������"<< endl;
    cout << "������� ��������" << "\t\t\t\t\t\t"<<"������� �������" <<endl;
    for (int i=0; i<k; i++)
    {
        GetRavn(U[i],a,b);
        GetRavn1(U[i],a,b);
        GetCheb(U[i],a,b);
        GetCheb1(U[i],a,b);
        Interpool(U[i],s[i]);
    }

    return 0;
}
