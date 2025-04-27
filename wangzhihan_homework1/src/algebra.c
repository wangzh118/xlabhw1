#include "algebra.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>  // 正确声明 malloc

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    if (row > MAX_MATRIX_SIZE || col > MAX_MATRIX_SIZE)
    {
        printf("Error: The matrix size exceeds the maximum size.\n");
        return m;
    }
    int i, j;
   for( i = 0; i < row; i++)
   {
       for ( j = 0; j < col; j++)
       {
           m.data[i][j] = 0;
       }
   }
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, a.cols);
    int i, j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return c;

    return create_matrix(0, 0);
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    if(a.rows != b.rows || a.cols != b.cols)
    {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, a.cols);
    int i, j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return c;
    return create_matrix(0, 0);
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    if(a.cols != b.rows)
    {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix c = create_matrix(a.rows, b.cols);
    int i, j,k;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < b.cols; j++)
        {
            c.data[i][j] = 0;
            for ( k = 0; k < a.cols; k++)
            {
                c.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return c;
    return create_matrix(0, 0);
}

Matrix scale_matrix(Matrix a, double k)
{
    Matrix c = create_matrix(a.rows, a.cols);
    int i, j;

    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[i][j] = a.data[i][j] * k;
        }
    }
    return c;
    return create_matrix(0, 0);
}


Matrix transpose_matrix(Matrix a)
{
    Matrix c = create_matrix(a.cols, a.rows);
    int i, j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            c.data[j][i] = a.data[i][j];
        }
    }
    return c;
    return create_matrix(0, 0);
}
// 计算矩阵的余子式
// 该函数用于计算矩阵的余子式，输入参数为原始矩阵mat、临时矩阵temp、行索引p、列索引q和矩阵的大小n
void getCofactor(int n, double mat[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], double temp[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE], int p, int q) {
    int i = 0, j = 0;
    int row, col;
    for (row = 0; row < n; row++) {
        for (col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];
                if (j == n - 1) { j = 0; i++; }
            }
        }
    }
}

double det_matrix(Matrix a) {
    if (a.rows != a.cols) return 0;
    if (a.rows == 1) return a.data[0][0];
    if (a.rows == 2) return a.data[0][0] * a.data[1][1] - a.data[0][1] * a.data[1][0];

    int n = a.rows;
    double temp[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];

    double det = 0;
    int sign = 1;
    int i,j,f;
    for (f = 0; f < n; f++) {
        getCofactor(n, a.data, temp, 0, f);
        Matrix subMatrix = create_matrix(n - 1, n - 1);
        for ( i = 0; i < n - 1; i++)
            for ( j = 0; j < n - 1; j++)
                subMatrix.data[i][j] = temp[i][j];
        det += sign * a.data[0][f] * det_matrix(subMatrix);
        sign = -sign;
    }
    // 释放临时矩阵的内存
    // free(subMatrix.data);  // 注意：这里的内存释放需要在实际使用中进行管理
    // 这里的subMatrix.data是一个静态数组，不需要释放内存
    return det;
}


// 计算矩阵的逆矩阵
Matrix inv_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }

    if (a.rows == 1 && a.cols == 1) {
        Matrix inv = create_matrix(1, 1);
        inv.data[0][0] = 1 / a.data[0][0];
        return inv;
    }

    double det = det_matrix(a);
    if (det == 0) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }

    Matrix inv = create_matrix(a.rows, a.cols);
    int i, j;
    for ( i = 0; i < a.rows; i++) {
        for (j = 0; j < a.cols; j++) {
            Matrix subMatrix = create_matrix(a.rows - 1, a.cols - 1);
            getCofactor(a.rows, a.data, subMatrix.data, i, j);
            inv.data[j][i] = pow(-1, i + j) * det_matrix(subMatrix) / det;
        }
    }

    return inv;
}
int rank_matrix(Matrix a) {
    int rank = a.rows; // 初始秩为矩阵的行数
    int i, j, row, col;

    // 创建一个副本矩阵以避免修改原矩阵
    Matrix temp = create_matrix(a.rows, a.cols);
    for (i = 0; i < a.rows; i++) {
        for (j = 0; j < a.cols; j++) {
            temp.data[i][j] = a.data[i][j];
        }
    }

    // 高斯消元法
    for (row = 0, col = 0; row < rank && col < a.cols; col++) {
        // 找到当前列中绝对值最大的元素作为主元
        int maxRow = row;
        for (i = row + 1; i < rank; i++) {
            if (fabs(temp.data[i][col]) > fabs(temp.data[maxRow][col])) {
                maxRow = i;
            }
        }

        // 如果主元为 0，则跳过当前列
        if (fabs(temp.data[maxRow][col]) < 1e-9) {
            rank--; // 秩减少
            for (i = row; i < rank; i++) {
                for (j = 0; j < a.cols; j++) {
                    temp.data[i][j] = temp.data[i + 1][j];
                }
            }
            row--; // 保持当前行不变
            continue;
        }

        // 交换当前行和主元行
        if (maxRow != row) {
            for (j = 0; j < a.cols; j++) {
                double tempVal = temp.data[row][j];
                temp.data[row][j] = temp.data[maxRow][j];
                temp.data[maxRow][j] = tempVal;
            }
        }

        // 消去当前列的其他行
        for (i = row + 1; i < rank; i++) {
            double factor = temp.data[i][col] / temp.data[row][col];
            for (j = col; j < a.cols; j++) {
                temp.data[i][j] -= factor * temp.data[row][j];
            }
        }

        row++; // 处理下一行
    }


    return rank;
}
double trace_matrix(Matrix a)
{
    int i;
    for(i = 0; i < a.rows; i++)
    {
        if (a.rows != a.cols)
        {
            printf("Error: The matrix must be a square matrix.\n");
            return 0;
        }
    }
    double trace = 0;
    for(i = 0; i < a.rows; i++)
    {
        trace += a.data[i][i];
    }
    return trace;
        
    
    
    return 0;
}

void print_matrix(Matrix a)
{
    int i,j;
    for (i = 0; i < a.rows; i++)
    {
        for (j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}