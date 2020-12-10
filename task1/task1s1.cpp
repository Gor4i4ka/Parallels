#include <cstdio>
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>

using namespace std;

// Nx, Ny, K1, K2
int 
ErrorMsg(int mode) {
	const char * msg[4] = {
		"Read the task\n",
		"The path to the file is wrong or its data is tainted\n",
		"The data from program arguments is tainted\n",
		"Call the program with no params. Vital info.\n"
	};
	cout << msg[mode];
	return mode;
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

bool **
fillGraph(int Nx, int Ny, int K1, int K2, int *NumOnes) {
	bool** graph = (bool**)malloc((Ny + 1)*(Nx + 1) * sizeof(bool*));
	for (auto i=0; i < (Nx + 1)*(Ny + 1); i++)
		graph[i] = (bool*)malloc((Ny + 1)*(Nx + 1)* sizeof(bool));

	for (auto i=0; i < (Ny + 1); i++)
	for (auto j=0; j < (Nx + 1); j++) 
		graph[i][j] = false;

	// i - row; j - column
	for (auto i=0; i < (Ny + 1)*(Nx + 1); i++)
	{
	   	// diagonal
		graph[i][i] = true;
		(*NumOnes) += 1;
		// vertical line - up
		if (i > Nx) {
			graph[i][i - (Nx + 1)] = true;
			graph[i - (Nx + 1)][i] = true;
			*NumOnes += 2;
		}
		// horizontal line - right
		if ((i % (Nx + 1)) < Nx) {
			graph[i][i + 1] = true;
			//cout << "PEPEGA2.1\n";
			graph[i + 1][i] = true;
			//cout << "PEPEGA2.2\n";
			*NumOnes += 2;
		}
		// slash line - up
		if (i > Nx)
		if ((i % (Nx + 1)) < Nx) {
			auto i1 = (i - (Nx + 1));
			auto i2 = i1 - (i1 / (Nx + 1));
			auto period_id = i2 % (K1 + K2);
			if (period_id > (K1 - 1)) {
				graph[i][i - Nx] = true;
				graph[i - Nx][i] = true;
				*NumOnes += 2;
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
void
printPortrait(int **result) {
	

	int NumNodes = *result[0];
	cout << "Amount of Nodes: " << NumNodes << "\n";

	int size = result[1][NumNodes];
	cout << "Amount of lines: " << size << "\n";

	for (auto node_id = 0; node_id < NumNodes; node_id++) {
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
		for (auto column_id = start; column_id < end; column_id++)
                	cout << result[2][column_id] << " ";
		cout << "\n";
	
	}
	return;
}

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

int 
main (int argc, char **argv)
{
	int Nx, Ny, K1, K2;
	int NumOnes = 0;

	if (argc == 1) 
		return ErrorMsg(0); 
	
	if (argc == 2) 
	if (readData(argv, &Nx, &Ny, &K1, &K2, 0))
		return ErrorMsg(1); 
	
	if (argc == 5) 
	if (readData(argv, &Nx, &Ny, &K1, &K2, 1))
		return ErrorMsg(2); 
	
	bool **graph = fillGraph(Nx, Ny, K1, K2, &NumOnes);
	int **portrait = fillPortrait(graph, NumOnes, (Nx+1)*(Ny+1));

	
	printGraph(graph, Ny, Nx);
	printPortrait(portrait);

	return 0;
}
