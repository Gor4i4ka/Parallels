#include <cstdio>
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

// Nx, Ny, K1, K2

int Errormsg(int mode);
int argv_read(int argc, char **argv, int* Nx_point, int* Ny_point, int* K1_point, int* K2_point);
int readData(char **argv, int *Nx, int *Ny, int *K1, int *K2, bool mode);
// LEGACY: returns Graph object
bool ** fillGraph(int Nx, int Ny, int K1, int K2, int *size);
void printGraph(bool **graph, int Ny, int Nx);
// LEGACY: returns Portrait object
int ** fillPortrait(bool **graph, int size, int NumNodes);
void printPortrait(int **result);
int fillSystem(int N, int* IA, int* JA, double* A, double* b);
int printSystem(int N, int* IA, int* JA, double* A, double* b);

// STORING PROGRAMM MESSAGES BEGIN

int 
ErrorMsg(int mode) {
	const char * msg[6] = {
		"Read the task\n",
		"The path to the file is wrong or its data is tainted\n",
		"The data from program arguments is tainted\n",
		"Broken fill_system\n",
		"Broken print_system\n",
		"Call the program with no params. Vital info.\n"
	};
	cout << msg[mode];
	return mode;
}

// STORING PROGRAMM MESSAGES END
// READING DATA BEGIN

int
argv_read(int argc, char **argv, int* Nx_point, int* Ny_point, int* K1_point, int* K2_point) {
        if (argc == 1)
                return ErrorMsg(0);

        if (argc == 2)
        if (readData(argv, Nx_point, Ny_point, K1_point, K2_point, 0))
                return ErrorMsg(1);

        if (argc == 5)
        if (readData(argv, Nx_point, Ny_point, K1_point, K2_point, 1))
                return ErrorMsg(2);
        return 0;
}


int 
readData(char **argv, int *Nx, int *Ny, int *K1, int *K2, bool mode) {
	int nx, ny, k1, k2;
	if (mode) {
		// read from argv
		*Nx = atoi(argv[1]);
		*Ny = atoi(argv[2]);
		*K1 = atoi(argv[3]);
		*K2 = atoi(argv[4]);
		if (!*Nx or !*Ny or !*K1 or !*K2)
			return 1;
	}
	else {
		if (auto fhandler = fopen(argv[1], "r"))
		if (fscanf(fhandler, "%d,%d,%d,%d", Nx, Ny, K1, K2) == EOF)
			return 1;
		if (!*Nx or !*Ny or !*K1 or !*K2)
                        return 1;
		return 0;
	}
	return 0;

}

// READING DATA END
// GRAPH OBJECT BEGIN

bool **
fillGraph(int Nx, int Ny, int K1, int K2, int *size) {
	bool** graph = (bool**)malloc((Ny + 1)*(Nx + 1) * sizeof(bool*));
	for (auto i=0; i < (Nx + 1)*(Ny + 1); i++)
		graph[i] = (bool*)malloc((Ny + 1)*(Nx + 1)* sizeof(bool));

	for (int i=0; i < (Ny + 1); i++)
	for (int j=0; j < (Nx + 1); j++) 
		graph[i][j] = false;

	// i - row; j - column
	for (int i=0; i < (Ny + 1)*(Nx + 1); i++)
	{
	   	// diagonal
		graph[i][i] = true;
		(*size) += 1;
		// vertical line - up
		if (i > Nx) {
			graph[i][i - (Nx + 1)] = true;
			graph[i - (Nx + 1)][i] = true;
			*size += 2;
		}
		// horizontal line - right
		if ((i % (Nx + 1)) < Nx) {
			graph[i][i + 1] = true;
			//cout << "PEPEGA2.1\n";
			graph[i + 1][i] = true;
			//cout << "PEPEGA2.2\n";
			*size += 2;
		}
		// slash line - up
		if (i > Nx + 1)
		if ((i % (Nx + 1)) > 0) {
			int i1 = (i - (Nx + 1));
			int i2 = (i1 - i1 / (Nx + 1) - 1);
			int period_id = i2 % (K1 + K2);
			if (period_id > (K1 - 1)) {
				graph[i][i - Nx -2] = true;
				graph[i - Nx - 2][i] = true;
				*size += 2;
			}
		}
	}
			
	return graph;
}

void
printGraph(bool **graph, int Ny, int Nx) {
for (auto line = 0; line < (Ny + 1)*(Nx + 1); line++) {
        for (auto pil = 0; pil < (Ny + 1)*(Nx + 1); pil++) {
                if (pil == 0) {
                        if (line < 10 or line > 99)
                                printf("%d  ", line);
                        else 
                                printf("%d ", line);
		}

                printf("%d ", graph[line][pil]);
        }
        cout << "\n";
        }

}

// GRAPH OBJECT END
// PORTRAIT OBJECT BEGIN

int **
fillPortrait(bool **graph, int size, int NumNodes) {
	cout << "SIZE " << size << "\n";
	int ** result = (int **)malloc(3*sizeof(int*));
	// N
	result[0] = (int*)malloc(sizeof(int));
      	*result[0] = NumNodes;
	// Alloc IA
	result[1] = (int*)malloc((NumNodes + 1)*sizeof(int));
	result[1][NumNodes] = size;
	// Alloc JA
	result[2] = (int*)malloc(size*sizeof(int));
	
	// The cicle
	int JA_pos = 0;
	for (int i=0; i < (NumNodes); i++) {
		result[1][i] = JA_pos;
		for (int j=0; j < (NumNodes); j++) {
			if (graph[i][j]) {
				result[2][JA_pos] = j;
				JA_pos++;
			}
		}
	}
	return result;
}

void
printPortrait(int **result) {


        int NumNodes = *result[0];
        cout << "Amount of Nodes: " << NumNodes << "\n";

        int size = result[1][NumNodes];
        cout << "Amount of lines: " << size << "\n";

        for (int node_id = 0; node_id < NumNodes; node_id++) {
                int start, end;
                cout << node_id << ": ";
                if (node_id != NumNodes - 1) {
                        start = result[1][node_id];
                        end = result[1][node_id + 1];
                }
                else {
                        start = result[1][node_id];
                        end = size;
                }
                for (int column_id = start; column_id < end; column_id++)
                        cout << result[2][column_id] << " ";
                cout << "\n";

        }
        return;
}

// PORTRAIT OBJECT END
// SYSTEM BY DOCUMENTATION BEGIN

int
fillSystem(int N, int* IA, int* JA, double* A, double* b) {

        int size = IA[N];
        for (int i=0; i < N; i++) {
                int start = IA[i];
                int end = size;
                if (i < N -1)
                        end = IA[i+1];

                double diagonal_sum = 0;
                int diag_ind_ja = -1;

                // Non-diagonal A
                for (int ja_pos=start; ja_pos < end; ja_pos++) {
                        int j = JA[ja_pos];
                        if (i == j)
                                diag_ind_ja = ja_pos;
                        if (i != j) {
                                double val = cos(i*j + i + j);
                                //val = i;
                                //cout << val << "\n";
                                A[ja_pos] = val;
                                diagonal_sum += abs(val);
                        }
                }
                //cout << i << "\n";
                // diagonal A
                A[diag_ind_ja] = 1.234 * diagonal_sum;
                // vector b
                b[i] = sin(i);
        }
        return 0;
}

int
printSystem(int N, int* IA, int* JA, double* A, double* b) {
	
	int size = IA[N];
	for (int i=0; i < N; i++) {
		int start = IA[i];
		int end = size;
		if (i < N - 1)
			end = IA[i+1];
		cout.precision(3);
		cout << "IA " << i << ": ";
		for (int ja_pos = start; ja_pos < end; ja_pos++)
			cout << JA[ja_pos] << ": " << A[ja_pos] << ", ";
		cout << "b: " << b[i] << "\n";	
	}
	return 0;
}

// SYSTEM BY DOCUMENTATION END
// MAIN CALLING FUNCTION

int 
main (int argc, char **argv)
{
	int Nx, Ny, K1, K2;
	int size = 0;

	if (int error_id = argv_read(argc, argv, &Nx, &Ny, &K1, &K2))
		return ErrorMsg(error_id);

	bool **graph = fillGraph(Nx, Ny, K1, K2, &size);
	int **portrait = fillPortrait(graph, size, (Nx+1)*(Ny+1));
	
	//printGraph(graph, Ny, Nx);
	//printPortrait(portrait);

	int N = *portrait[0];
	int* IA = portrait[1];
	int* JA = portrait[2];

	double* A = (double*)malloc(size*sizeof(double));
	double* b = (double*)malloc(N*sizeof(double));

	if (fillSystem(N, IA, JA, A, b))
		return ErrorMsg(3);
	if (printSystem(N, IA, JA, A, b))
		return ErrorMsg(4);


	
	return 0;
}
