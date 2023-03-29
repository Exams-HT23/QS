// Headers for the utility functions in utils.c
void printMatrix(unsigned int n, unsigned char matrix[n][n]);
void transposeMatrix(unsigned int n, unsigned char matrix[n][n]);
void writeToFile(unsigned int len, const unsigned int data[len], FILE *fp);
void randomParityMatrix(unsigned int n, unsigned char matrix[n][n], unsigned int size);
void addRows(unsigned int n, unsigned char matrix[n][n], unsigned int i, unsigned int j);
void copyMatrix(unsigned int n, unsigned char matrix1[n][n], unsigned char matrix2[n][n]);
unsigned int getSubrowPattern(unsigned int n, unsigned char matrix[n][n], unsigned int m, unsigned int section, int row);

// Headers for functions that compute summary statistics: used for testing.
unsigned int getMin(unsigned int n, const unsigned int data[n]);
unsigned int getMax(unsigned int n, const unsigned int data[n]);
double getMean(unsigned int n, const unsigned int data[n]);
double getSD(unsigned int n, const unsigned int data[n], double mean);
