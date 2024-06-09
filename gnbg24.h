/// Author: Vladimir Stanovov (vladimirstanovov@yandex.ru)
/// Last edited: March 7th, 2024
/// C++ implementation of Generalized Numerical Benchmark Generator (GNBG) instance for GECCO 2024 competition
/// Includes implementation of simple Differential Evolution (DE) with rand/1 strategy and binomial crossover
/// Problem parameters are read from f#.txt files which should be prepared with python script convert.py from f#.mat
/// Competition page: https://competition-hub.github.io/GNBG-Competition/
/// Reference:
/// D. Yazdani, M. N. Omidvar, D. Yazdani, K. Deb, and A. H. Gandomi, "GNBG: A Generalized
///   and Configurable Benchmark Generator for Continuous Numerical Optimization," arXiv prepring	arXiv:2312.07083, 2023.
/// A. H. Gandomi, D. Yazdani, M. N. Omidvar, and K. Deb, "GNBG-Generated Test Suite for Box-Constrained Numerical Global
///   Optimization," arXiv preprint arXiv:2312.07034, 2023.
/// MATLAB and Python versions: https://github.com/Danial-Yazdani/GNBG-Instances
#include <fstream>
using namespace std;
class GNBG
{
public:
    int FEval;
    int MaxEvals;
    int Dimension;
    int CompNum;
    double MinCoordinate;
    double MaxCoordinate;
    double AcceptanceThreshold;
    double OptimumValue;
    double BestFoundResult;
    double fval;
    double AcceptanceReachPoint;
    double* CompSigma;
    double* Lambda;
    double* OptimumPosition;
    double* FEhistory;
    double* temp;
    double* a;
    double** CompMinPos;
    double** CompH;
    double** Mu;
    double** Omega;
    double*** RotationMatrix;
    GNBG(int func_num);
    ~GNBG();
    double Fitness(double* xvec);
};
GNBG::GNBG(int func_num)
{
    char buffer[15];
    sprintf(buffer,"f%d.txt",func_num);
    ifstream fin(buffer);
    AcceptanceReachPoint = -1;
    FEval = 0;
    fin>>MaxEvals;
    fin>>AcceptanceThreshold;
    fin>>Dimension;
    fin>>CompNum;
    fin>>MinCoordinate;
    fin>>MaxCoordinate;
    FEhistory = new double[MaxEvals];
    a = new double[Dimension];
    temp = new double[Dimension];
    OptimumPosition = new double[Dimension];
    CompSigma = new double[CompNum];
    Lambda = new double[CompNum];
    CompMinPos = new double*[CompNum];
    CompH = new double*[CompNum];
    Mu = new double*[CompNum];
    Omega = new double*[CompNum];
    RotationMatrix = new double**[CompNum];
    for(int i=0;i!=CompNum;i++)
    {
        CompMinPos[i] = new double[Dimension];
        CompH[i] = new double[Dimension];
        Mu[i] = new double[2];
        Omega[i] = new double[4];
        RotationMatrix[i] = new double*[Dimension];
        for(int j=0;j!=Dimension;j++)
            RotationMatrix[i][j] = new double[Dimension];
    }
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=Dimension;j++)
            fin>>CompMinPos[i][j];
    for(int i=0;i!=CompNum;i++)
        fin>>CompSigma[i];
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=Dimension;j++)
            fin>>CompH[i][j];
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=2;j++)
            fin>>Mu[i][j];
    for(int i=0;i!=CompNum;i++)
        for(int j=0;j!=4;j++)
            fin>>Omega[i][j];
    for(int i=0;i!=CompNum;i++)
        fin>>Lambda[i];
    for(int j=0;j!=Dimension;j++)
        for(int k=0;k!=Dimension;k++)
            for(int i=0;i!=CompNum;i++)
                fin>>RotationMatrix[i][j][k];
    fin>>OptimumValue;
    for(int i=0;i!=Dimension;i++)
        fin>>OptimumPosition[i];
}
double GNBG::Fitness(double* xvec)
{
    double res = 0;
    for(int i=0;i!=CompNum;i++)
    {
        for(int j=0;j!=Dimension;j++)
            a[j] = xvec[j] - CompMinPos[i][j];
        for(int j=0;j!=Dimension;j++)
        {
            temp[j] = 0;
            for(int k=0;k!=Dimension;k++)
                temp[j] += RotationMatrix[i][j][k]*a[k]; //matmul rotation matrix and (x - peak position)
        }
        for(int j=0;j!=Dimension;j++)
        {
            if(temp[j] > 0)
                a[j] = exp(log( temp[j])+Mu[i][0]*(sin(Omega[i][0]*log( temp[j]))+sin(Omega[i][1]*log( temp[j]))));
            else if(temp[j] < 0)
                a[j] =-exp(log(-temp[j])+Mu[i][1]*(sin(Omega[i][2]*log(-temp[j]))+sin(Omega[i][3]*log(-temp[j]))));
            else
                a[j] = 0;
        }
        fval = 0;
        for(int j=0;j!=Dimension;j++)
            fval += a[j]*a[j]*CompH[i][j];
        fval = CompSigma[i] + pow(fval,Lambda[i]);
        res = (i == 0)*fval + (i != 0)*min(res,fval);//if first iter then save fval, else take min
    }
    if(FEval > MaxEvals)
        return res;
    FEhistory[FEval] = res;
    BestFoundResult = (FEval == 0)*res + (FEval != 0)*min(res,BestFoundResult);
    if(FEhistory[FEval] - OptimumValue < AcceptanceThreshold && AcceptanceReachPoint == -1)
       AcceptanceReachPoint = FEval;
    FEval++;
    return res;
}
GNBG::~GNBG()
{
    delete a;
    delete temp;
    delete OptimumPosition;
    delete CompSigma;
    delete Lambda;
    for(int i=0;i!=CompNum;i++)
    {
        delete CompMinPos[i];
        delete CompH[i];
        delete Mu[i];
        delete Omega[i];
        for(int j=0;j!=Dimension;j++)
            delete RotationMatrix[i][j];
        delete RotationMatrix[i];
    }
    delete CompMinPos;
    delete CompH;
    delete Mu;
    delete Omega;
    delete RotationMatrix;
    delete FEhistory;
}
