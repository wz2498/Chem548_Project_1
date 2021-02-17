//
//  Jacobi.c
//  Jacobi
//
//  Created by 朱文清 on 2021/2/17.
//

#include <stdio.h>
#include <math.h>
#include "FDM_Schrodinger.h"

#define N 100

void Jacobi(double matrix[N][N], double epslion, int times, double Eigenvectors[N][N])
{
    int i,j;
    int count = 0;
    int p = 0,q = 1;
    double max;
    double theta,costheta,sintheta;
    
    
    /* Initialize the eigenvector matrix */
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
        /* Traverse the matrix to get the largest non-diagonal element and its position */
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
            break;
        }
        
        /* Calculate the rotation angle \theta */
        theta = 1.0/2.0 * atanf(2*matrix[p][q]/(matrix[p][p]-matrix[q][q]));
        costheta = cos(theta);
        sintheta = sin(theta);
        
        /* Calculate the matrix after rotation transformation */
        double matrixip[N];
        double matrixiq[N];
        for (i = 0; i < N; i++)
        {
            matrixip [i] = matrix[i][p];
            matrixiq [i] = matrix[i][q];
        }
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
        
        /* Calculate the eigenvector matrix */
        double Eigenvectorsip;
        for (i = 0; i < N; i++)
        {
            Eigenvectorsip=Eigenvectors[i][p];
            Eigenvectors[i][p]=Eigenvectors[i][p]*costheta+Eigenvectors[i][q]*sintheta;
            Eigenvectors[i][q]=-Eigenvectorsip*sintheta+Eigenvectors[i][q]*costheta;
        }
    }
}
