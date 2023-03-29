#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

unsigned int sythGauss(unsigned int n, unsigned char matrix[n][n]) {
    // Count for the number of CNOT gates.
    unsigned int count = 0;

    // First, eliminate ones from the lower half of the matrix.
    for (unsigned int col = 0; col < n; col++) {
        unsigned char diag = matrix[col][col];

        // Iterate over the rows below the diagonal.
        for (unsigned int row = col + 1; row < n; row++) {
            // If there is a one below the diagonal, it must be eliminated.
            if (matrix[row][col] == 1) {
                // Use the one to place a one on the diagonal, if necessary.
                if (diag == 0) {
                    addRows(n, matrix, row, col);
                    count++;
                    diag = 1;
                }
                // Eliminate the one below the diagonal.
                addRows(n, matrix, col, row);
                count++;
            }
        }
    }

    // Eliminate ones from the upper half of the matrix.
    for (int col = n - 1; col > 0; col--) {
        for (int row = col - 1; row >= 0; row--) {
            if (matrix[row][col] == 1) {
                addRows(n, matrix, col, row);
                count++;
            }
        }
    }

    return count;
}

unsigned int sythPMHLower(unsigned int n, unsigned char matrix[n][n], unsigned int m)  {
    unsigned int count = 0;

    // Iterate over each section of size `m`.
    unsigned int numSections = ceil(n / (double) m);
    for (int section = 1; section <= numSections; section++) {
        int sectionStart = (section - 1) * m;

        // Map sub-row patterns to the corresponding row in which they first appear.
        unsigned int patternToRow[1 << m];
        for (int row = n - 1; row >= sectionStart; row--) {
            unsigned int pattern = getSubrowPattern(n, matrix, m, section, row);
            patternToRow[pattern] = row;
        }

        // Iterate over the rows below the first diagonal for the section.
        for (int row = sectionStart + 1; row < n; row++) {
            // Get the row corresponding to the first occurrence of this sub-row's pattern.
            unsigned int pattern = getSubrowPattern(n, matrix, m, section, row);

            // Check if this sub-row contains a duplicate pattern (cannot be all zeros).
            unsigned int match = patternToRow[pattern];
            if (pattern != 0 && match != row) {
                addRows(n, matrix, match, row);
                count++;
            }
        }

        // Check for the case where `m` does not nicely divide `n`.
        unsigned int sectionEnd = sectionStart + m;
        if (sectionEnd > n) {
            sectionEnd = n;
        }

        // Iterate over the columns in the section.
        for (unsigned int col = sectionStart; col < sectionEnd; col++) {
            unsigned char diag = matrix[col][col];

            // Iterate over the rows below the diagonal for the column.
            for (unsigned int row = col + 1; row < n; row++) {
                // If there is a one below the diagonal, it must be eliminated.
                if (matrix[row][col] == 1) {
                    // Use the one to place a one on the diagonal if necessary.
                    if (diag == 0) {
                        addRows(n, matrix, row, col);
                        count++;
                        diag = 1;
                    }
                    // Eliminate the one below the diagonal.
                    addRows(n, matrix, col, row);
                    count++;
                }
            }
        }
    }

    return count;
}

unsigned int sythPMH(unsigned int n, unsigned char matrix[n][n], unsigned int m) {
    // Reduce the lower and upper triangles of the matrix.
    unsigned int count = sythPMHLower(n, matrix, m);
    transposeMatrix(n, matrix);
    count += sythPMHLower(n, matrix, m);
    return count;
}

void experiment1() {
    // Create CSV files to save results to.
    FILE *fileGauss = fopen("scalability_gauss.csv", "w+");
    FILE *filePHM = fopen("scalability_pmh.csv", "w+");
    if (fileGauss == NULL || filePHM == NULL) {
        exit(-1);
    }

    // Maximum number of qubits to consider, qubit step size and the number of runs per qubit amount.
    unsigned int nMax = 100;
    unsigned int nStep = 20;
    unsigned int numRuns = 100;

    // Two extra entries for the number of qubits and circuit size.
    unsigned int dataGauss[numRuns + 2];
    unsigned int dataPMH[numRuns + 2];

    // Iterate over the range of qubits to consider.
    for (unsigned int n = nStep; n <= nMax; n += nStep) {
        // Log progress.
        printf("%d/%d\n", n, nMax);

        // Get the section size and the number of random row operations to perform
        // when randomly generating parity matrices.
        unsigned int m = (unsigned int) round(log2(n) / 2.0);
        unsigned int size = (unsigned int) floor(pow(n, 2) / 4.0);

        // Record the number of qubits and the circuit size.
        dataGauss[0] = n;
        dataGauss[1] = size;
        dataPMH[0] = n;
        dataPMH[1] = size;

        // Allocate space on the heap for the parity matrix.
        // Make a copy so both algorithms can be run.
        unsigned char (*matrix1)[n] = calloc(n * sizeof(*matrix1), sizeof(unsigned char));
        unsigned char (*matrix2)[n] = calloc(n * sizeof(*matrix2), sizeof(unsigned char));

        // Iterate for the given number of runs.
        for (unsigned int i = 0; i < numRuns; i++) {
            // Generate a random parity matrix, overwriting any previous matrix entries.
            randomParityMatrix(n, matrix1, size);
            copyMatrix(n, matrix1, matrix2);

            // Run the Gaussian elimination and PMH algorithms and record the CNOT counts.
            dataGauss[i+2] = sythGauss(n, matrix1);
            dataPMH[i+2] = sythPMH(n, matrix2, m);
        }

        // Free the memory allocated for the matrices since the next run will use matrices of a different size.
        free(matrix1);
        free(matrix2);

        // Save the results to a local CSV file.
        writeToFile(numRuns + 2, dataGauss, fileGauss);
        writeToFile(numRuns + 2, dataPMH, filePHM);
    }

    // Close the CSV files.
    fclose(fileGauss);
    fclose(filePHM);
}

void experiment2() {
    // Create a CSV file to save results to.
    FILE *file = fopen("optimisation_raw.csv", "w+");
    if (file == NULL) {
        exit(-1);
    }

    // Define the range of qubits to consider.
    unsigned int nMin = 10;
    unsigned int nMax = 230;
    unsigned int nStep = 20;

    // Define the maximum section size and number of runs per value of n and m.
    unsigned int mMax = 8;
    unsigned int numRuns = 1000;

    // Two extra entries for the values of n and m.
    unsigned int data[numRuns + 2];

    // Iterate over the range of qubits to consider.
    for (unsigned int n = nMin; n <= nMax; n += nStep) {
        // Log progress.
        printf("%d/%d\n", n, nMax);

        // Allocate space on the heap for the parity matrix.
        unsigned char (*matrix)[n] = calloc(n * sizeof(*matrix), sizeof(unsigned char));

        // Iterate over the range of section sizes to consider.
        // Make sure that `m` does not go beyond the n // 2.
        for (unsigned int m = 2; m <= floor(n / 2.0) && m <= mMax; m++) {
            // Record the values of n and m for the synthesised circuit sizes.
            data[0] = n;
            data[1] = m;

            // Iterate for the given number of runs.
            for (unsigned int i = 0; i < numRuns; i++) {
                // Pick a random number of row operations to perform uniformly from [1, n^2]
                unsigned int size = (rand() % (unsigned int) pow(n, 2)) + 1;
                randomParityMatrix(n, matrix, size);

                // Apply the PMH algorithm and get the number of CNOTs from the synthesised circuit.
                data[i + 2] = sythPMH(n, matrix, m);
            }

            // Write the results to a local CSV file.
            writeToFile(numRuns + 2, data, file);
        }

        // Free the memory allocated to the matrix since the next run will use matrices of a different size.
        free(matrix);
    }

    // Close the CSV file.
    fclose(file);
}

int main() {
    // Use a different random seed each time the code is executed.
    time_t t;
    srand((unsigned) time(&t));

    // Run either experiment 1 or 2.
    experiment1();
    // experiment2();

    return 0;
}
