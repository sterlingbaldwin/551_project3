//see psudocode http://www.math.sjsu.edu/~foster/m143m/gaussian_elimination_algorithms_4.pdf
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <sstream>
using namespace std;

vector<double> gauss(vector<vector<double> > A){
  int n = A.size();

  for (int i=0; i<n; i++) {
      // Search for maximum in this column
      double maxEl = abs(A[i][i]);
      int maxRow = i;
      for (int k=i+1; k<n; k++) {
          //if the value is larger then max, set max to val
          if (abs(A[k][i]) > maxEl) {
              maxEl = abs(A[k][i]);
              maxRow = k;
          }
      }
  }

  for (int k=i; k<n+1;k++) {
      double tmp = A[maxRow][k];
      A[maxRow][k] = A[i][k];
      A[i][k] = tmp;
  }

  for (int k=i+1; k<n; k++) {
      double c = -A[k][i]/A[i][i];
      for (int j=i; j<n+1; j++) {
          if (i==j) {
              A[k][j] = 0;
          } else {
              A[k][j] += c * A[i][j];
          }
      }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  vector<double> x(n);
  for (int i=n-1; i>=0; i--) {
      x[i] = A[i][n]/A[i][i];
      for (int k=i-1;k>=0; k--) {
          A[k][n] -= A[k][i] * x[i];
      }
  }
  return x;
}

int main(int argc, char *argv[]) {
    int n;
    int threads;
    if(argc < 3){
      return -1;
    }
    stringstream ss(argv[1]);
    ss >> n;
    stringstream ss1(argv[2]);
    ss1 >> threads;
    cout << "theads should be:" << threads << "\n";

    vector<double> line(n+1,0);
    vector< vector<double> > A(n,line);

    srand48(time(NULL));
    #pragma omp parallel for
    for (int i=0; i<n; i++) {
        for (int j=0; j<=n; j++) {
            A[i][j] = drand48() * 2000000 - 1000000;
        }
    }

}
