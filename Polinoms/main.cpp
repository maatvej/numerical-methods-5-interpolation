#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <locale>
#include<fstream>
#include <string>
#include <cmath>
#define KOL 51//максимальное число узлов сетки

using namespace std;

const double PI = 3.141592653589793;
double Xravn[KOL],Xcheb[KOL],Yravn[KOL],Ycheb[KOL],Xcheb1[KOL],Ycheb1[KOL],Xravn1[KOL],Yravn1[KOL];
double a = 0.0, b = 1.0;//интервал интерполирования

//интерполируемая функция
double f(double x)
{
    return x*x-asin(x-0.2);
}

//максимальное отклонение
double Maxf(int N, double Ax[],double Bx[])
{
    double maxb=0;
    for (int j=0; j<N; j++)
        {
        if (maxb<abs(Ax[j]-Bx[j])) maxb=abs(Ax[j]-Bx[j]);
        }
return maxb;
}

//определение равномерной сетки
void GetRavn(int n, double a, double b)
{
    double h = (b-a)/n;
    for(int i=0;i<=n;i++)
    {
        Xravn[i]=a+h*i;
        Yravn[i]=f(Xravn[i]);
    }
}

//определение равномерной сетки (от b до a)
void GetRavn1(int n, double a, double b)
{
    for(int i=0;i<=n;i++)
        {
            Xravn1[i]= Xravn[n-i];
            Yravn1[i]= Yravn[n-i];
        }
}

//определение оптимальной сетки (от b до a)
void GetCheb(int n, double a, double b)
{
    for(int i=0;i<=n;i++)
        {
            Xcheb[i]=0.5*((b-a)*cos((2.0*i+1.0)*PI/(2.0*(n+1.0)))+b+a);
            Ycheb[i]=f(Xcheb[i]);
        }
}

//определение оптимальной сетки (от a до b)
void GetCheb1(int n, double a, double b)
{
    for(int i=0;i<=n;i++)
        {
            Xcheb1[i]= Xcheb[n-i];
            Ycheb1[i]= Ycheb[n-i];
        }
}


//полином Лагранжа (n-число точек)
//x_arr[] и y_arr[] - сетка интерполирования
double my_lagrange(double x, int n, double x_arr[], double y_arr[])
{

    double sum = 0;
    for (int i = 0; i < n; i++){

        double l = 1;//начальное значение произведения lj
        for (int j = 0; j < n; j++)
            if (j != i){
                l = l* (x - x_arr[j]) / (x_arr[i] - x_arr[j]);

            }

        sum += y_arr[i] * l;//сумма lj*yj
    }
    return sum;
}

double my_newton(double x, int n, double x_arr[], double y_arr[])
{

    double sum = y_arr[0];
    for(int i = 1; i < n; ++i)//по всем узлам
        {
        double F = 0;//i-ое слагаемое содержит j произведений
        for(int j = 0; j <= i; ++j)
            {
            double den = 1;
            for(int k = 0; k <= i; ++k)
                if (k != j)
                //произведение конечных разностей между первыми i узлами
                    den *= (x_arr[j] - x_arr[k]);
            F += y_arr[j]/den;
        }
        //здесь считаем произведения (x-x[0])...(x-x[i-1])
        for(int k = 0; k < i; ++k)
            F *= (x - x_arr[k]);
        sum += F;//здесь посчитали i-ое слагаемое
        }
    return sum;
}


//интерполирование
void Interpool(int n, char s)
{
    int N = 10*n+1;
    double h = (b-a)/(10*n);

    double *X=new double[N];//проверочные точки
    double *FX=new double[N];//значение функции в проверочных точках
    double *LX=new double[N];//значение полинома Лагранжа в проверочных точках
    double *NX=new double[N];//значение полинома Ньютона в проверочных точках
    double *LXopt=new double[N];//значение полинома Лагранжа в проверочных точках(по оптимальной сетеке)
    double *NXopt=new double[N];//значение полинома Ньютона в проверочных точках(по оптимальной сетеке)
    char s1[] = "outputX.txt";
    s1[6]=s;
    ofstream fout(s1);//выводим таблицу результатов файл для построения графиков
    for (int i=0; i<N; i++)
    {
        X[i]=a+h*i;//сетка по x
        FX[i]=f(X[i]);//точное значение функции
        LX[i]=my_lagrange(X[i],n+1,Xravn,Yravn);//полином Лагранжа по равномерной сетке
        LXopt[i]=my_lagrange(X[i],n+1,Xcheb,Ycheb);//полином Лагранжа по оптимальной сетке
        //здесь определяем какую формулу Ньютона использовать
        if (X[i]<(a+b)/2) //для первой половины отрезка 1-ую
        {
            NX[i]=my_newton(X[i],n+1,Xravn,Yravn);//равномерная сетка
            NXopt[i]=my_newton(X[i],n+1,Xcheb1,Ycheb1);//оптимальная сетка
        }
        else //для второй половины отрезка 2-ую
        {
            NX[i]=my_newton(X[i],n+1,Xravn1,Yravn1);
            NXopt[i]=my_newton(X[i],n+1,Xcheb,Ycheb);
        }
        //вывод в файл
        fout << X[i] << '\t' << FX[i] << '\t'<< LX[i] << '\t' <<LXopt[i] << '\t' << NX[i] << '\t' <<NXopt[i]<< endl;
    }
    fout.close();
    //вывод таблицы по заданию
    cout << n << '\t' << N << '\t'<< Maxf(N,FX,LX) << '\t' <<Maxf(N,FX,LXopt)<< "\t\t\t";
    cout << n << '\t' << N << '\t'<< Maxf(N,FX,NX) << '\t' <<Maxf(N,FX,NXopt)<< endl;
    //освобождение памяти
    delete[] X;delete[] FX;delete[] LX;delete[] NX;delete[] LXopt;delete[] NXopt;
}



int main()
{
    setlocale(LC_ALL, "russian");
    int U[6] ={3, 10, 20, 30, 40, 50};
    int k = 6;
    char s[] = "123456";
    cout << "Поведение интерполяционного полинома при увеличении количества узлов интерполирования"<< endl;
    cout << "Полином Лагранжа" << "\t\t\t\t\t\t"<<"Полином Ньютона" <<endl;
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
