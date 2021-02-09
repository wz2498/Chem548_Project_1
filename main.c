//
//  main.c
//  Jacobi
//
//  Created by 朱文清 on 2021/2/9.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 100

void Jacobi(double matrix[N][N], double epslion, int times, double Eigenvectors[N][N]);
void Sort(double matrix[N][N]); //将特征值从大到小排列；

int main()
{
    double mass;        /* Mass of the particle */
    double v[N];        /* Potential */
    double h[N][N];     /* Hamiltonian matrix */
    double evec[N][N];  /* Eigenvectors (wave functions) */
    double intval;      /* Interval between grid points = delta_x */
    double x;
    int i,j,k;
    
    mass = 100.0;
    //printf("mass = %f\n", mass);
    intval = 1.0/(N-1);
    
    /* particle in box potential */
    for (i = 0; i < N; i++)
    {
        x = (double) i / (double) (N - 1);
        v[i] = 0.0;
    }
    
    /* 构造哈密顿矩阵 */
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            h[i][j] = 0.0;
        }
        h[i][i] = v[i]+1.0/mass/intval/intval;
        if (i != 0)
        {
            h[i][i-1] = -0.5/mass/intval/intval;
        }
        
        if (i != N-1)
        {
            
            h[i][i+1] = -0.5/mass/intval/intval;
            
        }
    }
    
    printf("\nCalculating......please wait ^-^ \n");
    printf("\n**********************\n\n");
    
    Jacobi(h,1e-10f,10000000,evec); // 用Jacobi算法将Hamiltonina矩阵对角化
    Sort(h);
    
    printf("Successfully finish!!!\n\n");

    printf("analytical\n");
    for (i = 0; i < 5; i++)
    {
        x = i+1;
        printf("%d %15.10f\n",i+1,x*x/(8.0*mass)*pow(2.0*3.1416,2));
    }
    printf("numerical\n");
    for (i = 0; i < 5; i++)
    {
        printf("%d %15.10f\n",i+1,h[i][i]);
    }
    printf("\n**********************\n");
    
}


void Jacobi(double matrix[N][N], double epslion, int times, double Eigenvectors[N][N])
{
    int i,j;
    int count = 0; // 计数，当计数超过规定次数，停止迭代；
    int p = 0,q = 1;
    double max;
    double theta,costheta,sintheta;
    
    
    /* 初始化特征向量 */
    for(i = 0; i < N && count==0; i++)
    {
        for (j=0; j<N; j++)
        {
            Eigenvectors[i][j] = 0.0;
        }
        Eigenvectors[i][i] = 1.0;
    }
    
    while(1)
    {
        /* 遍历矩阵，得到最大非对角元素及其位置 */
        max = fabs(matrix[0][1]);
        p = 0; q = 1;
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                if (max < fabs(matrix[i][j]) && i!=j)
                {
                    max = fabs(matrix[i][j]);
                    p = i; q = j;
                }
            }
        }
        
        count++;
        if (max < epslion || count > times)
        {
            //printf("***%.10lf,%d***\n",max,count);
            break;
        }
        
# if(0)
        /* 给出计算进度 */
        if (count%10000 == 0)
        {
            printf("***%0.1f,%d***\n",max,count);
        }
        
# endif
        
        /* 计算旋转角度 */
        theta = 1.0/2.0 * atanf(2*matrix[p][q]/(matrix[p][p]-matrix[q][q]));
        costheta = cos(theta);
        sintheta = sin(theta);
        
        // 记录旋转前矩阵会改变的两行；
        double matrixip[N];
        double matrixiq[N];
        for (i = 0; i < N; i++)
        {
            matrixip [i] = matrix[i][p];
            matrixiq [i] = matrix[i][q];
        }
        
        /* 计算旋转后矩阵 */
        for (i = 0; i < N; i++)
        {
            matrix[i][p] = matrixip[i]*costheta + matrixiq[i]*sintheta;
            matrix[i][q] = -matrixip[i]*sintheta + matrixiq[i]*costheta;
            matrix[p][i] = matrix[i][p];
            matrix[q][i] = matrix[i][q];
        }
        matrix[p][p]=costheta*costheta*matrixip[p]+2.0*sintheta*costheta*matrixip[q]+sintheta*sintheta*matrixiq[q];
        matrix[q][q]=costheta*costheta*matrixiq[q]-2.0*sintheta*costheta*matrixip[q]+sintheta*sintheta*matrixip[p];
        matrix[q][p]=0.0;
        matrix[p][q]=0.0;
        
        /* 计算特征向量 */
        double Eigenvectorsip;
        for (i = 0; i < N; i++)
        {
            Eigenvectorsip=Eigenvectors[i][p];
            Eigenvectors[i][p]=Eigenvectors[i][p]*costheta+Eigenvectors[i][q]*sintheta;
            Eigenvectors[i][q]=-Eigenvectorsip*sintheta+Eigenvectors[i][q]*costheta;
        }
    }
}


void Sort(double matrix[N][N])
{
    int i,j;
    for (i = 0; i < N; i++ )
    {
        for (j = 0; j < N-2; j++)
        {
            double temp;
            if (matrix[j][j] > matrix[j+1][j+1])
            {
                temp = matrix[j][j];
                matrix[j][j] = matrix[j+1][j+1];
                matrix[j+1][j+1] = temp;
            }
        }
    }
}
