//see psudocode http://www.math.sjsu.edu/~foster/m143m/gaussian_elimination_algorithms_4.pdf
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <sstream>
using namespace std;

vector<double> gauss(vector<vector<double> > A){
  //Do gauss
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
