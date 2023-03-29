#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

void printMatrix(unsigned int n, unsigned char matrix[n][n]) {
    // Print the entries of the given matrix: used for debugging.
    for (unsigned int i = 0; i < n; i++) {
        printf("[");
        for (unsigned int j = 0; j < n; j++) {
            printf("%d", matrix[i][j]);
            if (j < n - 1) {
                printf(", ");
            }
        }
        printf("]");
        if (i < n - 1) {
            printf("\n");
        }
    }
    printf("\n\n");
}

void transposeMatrix(unsigned int n, unsigned char matrix[n][n]) {
    // Transpose the given `matrix` in-place.
    unsigned int temp;
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i+1; j < n; j++) {
            temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

void writeToFile(unsigned int len, const unsigned int data[len], FILE *file) {
    // Write `data` as a line to the given `file`.
    for (int i = 0; i < len - 1; i++) {
        fprintf(file, "%d,", data[i]);
    }
    fprintf(file, "%d\n", data[len-1]);
}

void randomParityMatrix(unsigned int n, unsigned char matrix[n][n], unsigned int size) {
    // Start with the identity matrix.
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            matrix[i][j] = (i == j);
        }
    }

    // Iterate for the given number of CNOT gates to generate.
    for (unsigned int i = 0; i < size; i++) {
        unsigned int row1 = rand() % n;
        unsigned int row2 = rand() % n;

        // Make sure that the rows are not the same.
        while (row1 == row2) {
            row1 = rand() % n;
            row2 = rand() % n;
        }
        // Add row2 to row1 to simulate the action of a CNOT gate.
        // I.e., elementary matrix multiplication.
        addRows(n, matrix, row1, row2);
    }
}

void addRows(unsigned int n, unsigned char matrix[n][n], unsigned int i, unsigned int j) {
    // XOR the entries in row j with the entries in row i.
    for (unsigned int col = 0; col < n; col++) {
        matrix[j][col] ^= matrix[i][col];
    }
}

void copyMatrix(unsigned int n, unsigned char matrix1[n][n], unsigned char matrix2[n][n]) {
    // Copy the elements of matrix1 to matrix2.
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            matrix2[i][j] = matrix1[i][j];
        }
    }
}

unsigned int getSubrowPattern(unsigned int n, unsigned char matrix[n][n], unsigned int m, unsigned int section, int row) {
    // Map a sub-row to an integer value by computing the decimal value corresponding to the binary vector.
    unsigned int pattern = 0;
    for (unsigned int i = 0; i < m; i++) {
        // Make sure that the column is not beyond `n`.
        unsigned int col = (section - 1) * m + i;
        if (col >= n) {
            break;
        }
        pattern += matrix[row][col] * (1 << i);
    }
    return pattern;
}

unsigned int getMin(unsigned int n, const unsigned int data[n]) {
    // Get the minimum value in the given array.
    unsigned int minVal = UINT_MAX;
    for (unsigned int i = 0; i < n; i++) {
        if (data[i] < minVal) {
            minVal = data[i];
        }
    }
    return minVal;
}

unsigned int getMax(unsigned int n, const unsigned int data[n]) {
    // Get the maximum value in the given array.
    unsigned int maxVal = 0;
    for (unsigned int i = 0; i < n; i++) {
        if (data[i] > maxVal) {
            maxVal = data[i];
        }
    }
    return maxVal;
}

double getMean(unsigned int n, const unsigned int data[n]) {
    // Compute the mean of the given array.
    unsigned int sum = 0;
    for (unsigned int i = 0; i < n; i++) {
        sum += data[i];
    }
    return sum / (double) n;
}

double getSD(unsigned int n, const unsigned int data[n], double mean) {
    // Compute the standard deviation of the given array assuming the mean has already been computed.
    double sd = 0;
    for (unsigned int i = 0; i < n; i++) {
        sd += pow(data[i] - mean, 2);
    }
    return sqrt(sd / n);
}
