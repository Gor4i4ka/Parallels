#include <cstdio>
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

using namespace std;

// MPI VALUES FOR EACH PROCESS
int MPI_processes_amount;
int MPI_process_rank;
int MPI_process_i;
int MPI_process_j;

////////////////////////////////////////////PROTOTYPES BEGIN/////////////////////////////////////////

// READING PARAMETERS BEGIN
int Errormsg(int mode);
int argv_read(int argc, char **argv, int* Nx_point, int* Ny_point, int* K1_point, int* K2_point, int* Px_point, int* Py_point);
int readData(char **argv, int *Nx, int *Ny, int *K1, int *K2, int* Px, int* Py, bool mode);
// READING PARAMETERS END

// COMMON BEGIN
int G2L(int i_g, int* L2G, int N);
// COMMON END


// GENERATE BEGIN (including debug prints)
void Generate(int Nx, int Ny, int K1, int K2, int Px, int Py, int* N, int* N0, int** Part, int** L2G, int** IA, int** JA);
void printGraph(bool **graph, int size_p_galo, int N_g, int* L2G);
void printPortrait(int* IA, int* JA, int* L2G, int N0);
// GENERATE END


int fillSystem(int N, int* IA, int* JA, double* A, double* b);
int printSystem(int N, int* IA, int* JA, double* A, double* b);

/////////////////////////////////////////PROTOTYPES END/////////////////////////////////////////

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
argv_read(int argc, char **argv, int* Nx_point, int* Ny_point, int* K1_point, int* K2_point, int* Px_point, int* Py_point) {
        if (argc == 1)
                return ErrorMsg(0);

        if (argc == 2)
        if (readData(argv, Nx_point, Ny_point, K1_point, K2_point, Px_point, Py_point, 0))
                return ErrorMsg(1);

        if (argc == 7)
        if (readData(argv, Nx_point, Ny_point, K1_point, K2_point, Px_point, Py_point, 1))
                return ErrorMsg(2);
        return 0;
}


int 
readData(char** argv, int* Nx, int* Ny, int* K1, int* K2, int* Px, int* Py, bool mode) {
	int nx, ny, k1, k2;
	if (mode) {
		// read from argv
		*Nx = atoi(argv[1]);
		*Ny = atoi(argv[2]);
		*K1 = atoi(argv[3]);
		*K2 = atoi(argv[4]);
		*Px = atoi(argv[5]);
		*Py = atoi(argv[6]);
		if (!*Nx or !*Ny or !*K1 or !*K2 or !*Px or !*Py)
			return 1;
	}
	else {
		if (auto fhandler = fopen(argv[1], "r"))
		if (fscanf(fhandler, "%d,%d,%d,%d,%d,%d", Nx, Ny, K1, K2, Px, Py) == EOF)
			return 1;
		if (!*Nx or !*Ny or !*K1 or !*K2 or !*Px or !*Py)
                        return 1;
		return 0;
	}
	return 0;

}

// READING DATA END
// COMMON BEGIN
int 
G2L (int i_g, int* L2G, int N) {
	for (int i=0; i < N; i++)
		if (L2G[i] == i_g)
			return i;
	return -1;

}
// COMMON END

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////GENERATE BEGIN/////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Generate(int Nx, int Ny, int K1, int K2, int Px, int Py, int* N, int* N0, int** Part, int** L2G, int** IA, int** JA) {
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes_amount);
        MPI_Comm_rank(MPI_COMM_WORLD, &MPI_process_rank);
	MPI_process_i = MPI_process_rank / Py;
	MPI_process_j = MPI_process_rank % Px;

	/////////////////////////////////////////////////////////
	// GENERATE PHASE 1: FILLING GRAPH
	/////////////////////////////////////////////////////////


	// Getting the boundaries BEGIN

	int ib = ((Ny + 1) / Py) * MPI_process_i;
	int ie = 0;
	int jb = ((Nx + 1) / Px) * MPI_process_j;
	int je = 0;

	{
		int left_out_i = (Ny + 1) % Py;
		int left_out_j = (Nx + 1) % Px;

		int extended_dist_i = left_out_i - MPI_process_i;
		int extended_dist_j = left_out_j - MPI_process_j;

		if (extended_dist_i > 0) 
			ie++;
		
		if (extended_dist_j > 0 ) 
			je++;

		ie += ib + (Ny + 1) / Py;
		je += jb + (Nx + 1) / Px;
	}

	// Getting the boundaries END


	// Getting width, length and size (without GALO)
	int columns_p = je - jb;
	int lines_p = ie - ib;
	int size_p = columns_p * lines_p;

	// N0 can be already found and added to N
	
	*N0 = size_p;
	*N = size_p;
	
	// Taking GALO into consideration (with GALO)
	int columns_p_galo = columns_p + 2;
	int lines_p_galo = lines_p + 2;
	int size_p_galo = lines_p_galo * columns_p_galo;
	
	// Initializing Part & L2G
	*Part = (int*) malloc(size_p_galo * sizeof(int));
	*L2G = (int*) malloc(size_p_galo * sizeof(int));

	for (int i=0; i < size_p_galo; i++) {
		(*L2G)[i] = -1;
		(*Part)[i] = -1;
	}

	// Setting GALO

	// Getting L2G & Part from Allpoints
	for (int i_l=0; i_l < size_p_galo; i_l++)
	{

		// Global indicies
		int i_line_g;
		int i_column_g;
		int i_g;
		
		// For local points
		if (i_l < size_p) {
		 	
			// Getting the point's local indicies
                	int i_line_l = i_l / columns_p;
                	int i_column_l = i_l % columns_p;


			// Getting global values
			i_line_g = ib + i_line_l;
			i_column_g = jb + i_column_l;
			i_g = i_line_g * (Nx + 1) + i_column_g;

			(*Part)[i_l] = MPI_process_rank;
			(*L2G)[i_l] = i_g;
		}
		// For GALO points
		else {
			
			int i_galo = i_l - size_p;
			
			// Getting global values
			bool is_top = false;
			bool is_bottom = false;
			bool is_left = false;
			bool is_right = false;

			// Check default cases

			if (i_galo - columns_p_galo < 0)
				is_top = true;
			else
			if (i_galo - 2*columns_p_galo < 0)
				is_bottom = true;
			else
			if (i_galo - 2*columns_p_galo - lines_p_galo + 2 < 0)
				is_left = true;
			else
				is_right = true;
			
			if (is_top and (MPI_process_i > 0)) {
				
				i_line_g = ib - 1;
				i_column_g = (jb - 1) + i_galo;
				i_g = i_line_g * (Nx + 1) + i_column_g;

				// CHECK special cases
                                if (i_galo == 0) 
				if (MPI_process_j > 0){
					int i_period_pos = i_g % (K1 + K2 + 1);
                                	if (!(i_period_pos >= K1)) 
                                       		continue;
				}
				else
					continue;

				if (i_galo - columns_p_galo + 1 == 0)
					continue;


				(*Part)[i_l] = (MPI_process_i - 1) * Px + MPI_process_j;
				(*L2G)[i_l] = i_g;

			}
			
			else
			if (is_bottom and (MPI_process_i < (Py - 1))) {
				
				i_line_g = ie;
				i_column_g = (jb - 1) + (i_galo - columns_p_galo);
				i_g = i_line_g * (Nx + 1) + i_column_g;

				// CHECK special cases
				if (i_galo - 2 * columns_p_galo + 1 == 0) 
				if (MPI_process_j < (Px - 1)){

					int i_con = i_g - (Nx + 2);
                        		int i_con_period_pos = i_con % (K1 + K2 + 1);
                        		
					if (!(i_con_period_pos >= K1))
				       		continue;
				}
				else
					continue;

				if (i_galo - columns_p_galo == 0)
					continue;	

				(*Part)[i_l] = (MPI_process_i + 1) * Px + MPI_process_j;
				(*L2G)[i_l] = i_g;
			}

			else
			if (is_left and (MPI_process_j > 0)) {
				
				i_column_g = jb -1;
				i_line_g = ib + (i_galo - 2 * columns_p_galo);
				i_g = i_line_g * (Nx + 1) + i_column_g;

				(*Part)[i_l] = MPI_process_rank - 1;
				(*L2G)[i_l] = i_g;
			}

			else
			if (is_right and (MPI_process_j < (Px - 1))) {
				
				i_column_g = je;
				i_line_g = ib + (i_galo - 2 * columns_p_galo - lines_p_galo + 2);
				i_g = i_line_g * (Nx + 1) + i_column_g;

				(*Part)[i_l] = MPI_process_rank + 1;
				(*L2G)[i_l] = i_g;
			}
			else		
				continue;

			// The GALO point is valid and we can count it
			*N += 1;

		}
	}
	
	// compress L2G & Part
        {
            	int* L2G_buf = (int*) malloc(size_p_galo * sizeof(int));
       	        int* Part_buf = (int*) malloc(size_p_galo * sizeof(int));	
		memcpy(L2G_buf, *L2G, size_p_galo * sizeof(int));
		memcpy(Part_buf, *Part, size_p_galo * sizeof(int));

                free(*L2G);
                free(*Part);

                *L2G = (int*) malloc(*N * sizeof(int));
                *Part = (int*) malloc(*N * sizeof(int));

		memcpy(*L2G, L2G_buf, size_p * sizeof(int));
		memcpy(*Part, Part_buf, size_p * (sizeof(int)));

                int i_buf = *N0;
                for (int i=*N0; i < *N; i++) {	
      			while(i_buf < size_p_galo) {
      				if (L2G_buf[i_buf] == -1) {
					i_buf++;
                                 }
                                 else {
      	                                 (*L2G)[i] = L2G_buf[i_buf];
                                          (*Part)[i] = Part_buf[i_buf];
					  i_buf++;
                                          break;
                                 }
                      	}
               	}
      	}

        // Initializing graph
        bool** graph = (bool**)malloc(*N * sizeof(bool*));
        for (int i=0; i < *N0; i++)
                graph[i] = (bool*) calloc(*N, sizeof(bool));


	int connection_count = 0;

	// Filling the graph with local points data
	for (int i_l=0; i_l< size_p; i_l++) {

		// Global indicies
                int i_line_g;
                int i_column_g;
                int i_g;

                // For local points

                // Getting the point's local indicies
                int i_line_l = i_l / columns_p;
                int i_column_l = i_l % columns_p;

		// Getting global values
       	        i_line_g = ib + i_line_l;
                i_column_g = jb + i_column_l;
                i_g = i_line_g * (Nx + 1) + i_column_g;

		// COUNTING CONNECTIONS
		
	   	// diagonal
		graph[i_l][G2L(i_g, *L2G, *N)] = true;
		connection_count += 1;

		// vertical line - up
		if (i_line_g > 0) {
			graph[i_l][G2L(i_g - (Nx + 1), *L2G, *N)] = true;
			connection_count += 1;
		}

		// vertical line - down
		if (i_line_g < Ny) {
			graph[i_l][G2L(i_g + (Nx + 1), *L2G, *N)] = true;
			connection_count += 1;
		}
		
		// horizontal line - left
		if (i_column_g > 0) {
			graph[i_l][G2L(i_g - 1, *L2G, *N)] = true;
			connection_count += 1;
		}

		// horizontal line - right
		if (i_column_g < Nx) {
			graph[i_l][G2L(i_g + 1, *L2G, *N)] = true;
			connection_count += 1;
		}
		
		// slash line - up 
		if (i_column_g > 0)
		if (i_line_g > 0) {
			int i_con = i_g - (Nx + 2);
			int i_con_period_pos = i_con % (K1 + K2 + 1);
			if (i_con_period_pos >= K1) {
				graph[i_l][G2L(i_con, *L2G, *N)] = true;
				connection_count += 1;
			}

		}
		// slash line - down
		if (i_column_g < Nx)
		if (i_line_g < Ny) {
			int i_period_pos = i_g % (K1 + K2 + 1);
			int i_con = i_g + Nx + 2;
			if (i_period_pos >= K1) {
				graph[i_l][G2L(i_con, *L2G, *N)] = true;
				connection_count += 1;
			}

		}
		// FILLING VALUES END
	}
	
	// Filling the graph END

        /////////////////////////////////////////////////////////////
        // THE GRAPH IS NOW FILLED AND CAN BE PRINTED
        /////////////////////////////////////////////////////////////
 	
	//printGraph(graph, *N0, *N, *L2G);
	

	////////////////////////////////////////////////
        // GENERATE: PHASE 2
        ///////////////////////////////////////////////

        // Alloc IA
        *IA = (int*)malloc((*N0 + 1)*sizeof(int));
        (*IA)[*N0] = connection_count;
        // Alloc JA
        *JA = (int*)malloc(connection_count*sizeof(int));

        // The cicle
        int JA_pos = 0;
        for (int i=0; i < *N0; i++) {
                (*IA)[i] = JA_pos;
                for (int j=0; j < *N; j++) {
                        if (graph[i][j]) {
                                (*JA)[JA_pos] = j;
                                JA_pos++;
                        }
                }
        }
	
	/////////////////////////////////////////////////////////////
        // THE PORTRAIT IS NOW FILLED AND CAN BE PRINTED
        /////////////////////////////////////////////////////////////

	//printPortrait(*IA, *JA, *L2G, *N0);

	return;
}

void
printGraph(bool **graph, int N0, int N, int* L2G) {

for (auto i_l = 0; i_l < N0; i_l++) {
	int i_g = L2G[i_l];

        for (auto pil = 0; pil < N; pil++) {
		
                if (pil == 0) {
                        if ((i_g > -1 and i_g < 10) or i_g > 99)
                                printf("%d  ", i_g);
                        else 
                                printf("%d ", i_g);
		}

                printf("%d ", graph[i_l][pil]);
        }
        cout << "\n";
        }

}


void
printPortrait(int* IA, int* JA, int* L2G, int N0) {

        for (int i = 0; i < N0; i++) {
                int start, end;
                cout << L2G[i] << ": ";
                
               	start = IA[i];
              	end = IA[i + 1];
               

                for (int j = start; j < end; j++)
                        cout << L2G[JA[j]] << " ";
                cout << "\n";

        }
        return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////GENERATE END/////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



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

	// GLOBAL STUFF
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &MPI_processes_amount);
        MPI_Comm_rank(MPI_COMM_WORLD, &MPI_process_rank);

	// PROGRAMM PARAMETERS
	int Nx, Ny, K1, K2, Px, Py;

	// GENERATE VARIABLES
	int N;
	int N0;
	int* Part;
	int* L2G;
	int* IA;
	int* JA;

	if (int error_id = argv_read(argc, argv, &Nx, &Ny, &K1, &K2, &Px, &Py))
		return ErrorMsg(error_id);


	int buf = 0;
	if (MPI_process_rank > 0)
		MPI_Recv(&buf, 1, MPI_INT, MPI_process_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
			
	Generate(Nx, Ny, K1, K2, Px, Py,  &N, &N0, &Part, &L2G, &IA, &JA);
	printf("E:proc %d/%d\n", MPI_process_rank, MPI_processes_amount);
	if (MPI_process_rank < MPI_processes_amount - 1)
		MPI_Send(&buf, 1, MPI_INT, MPI_process_rank + 1, 0, MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
/*	//bool** graph;
	int size = 0;
	int **portrait = fillPortrait(graph, size, (Nx+1)*(Ny+1));
	
	//printGraph(graph, Ny, Nx);
	//printPortrait(portrait);

	int Q = *portrait[0];
	int* IA = portrait[1];
	int* JA = portrait[2];

	double* A = (double*)malloc(size*sizeof(double));
	double* b = (double*)malloc(N*sizeof(double));

	if (fillSystem(Q, IA, JA, A, b))
		return ErrorMsg(3);
	if (printSystem(Q, IA, JA, A, b))
		return ErrorMsg(4);


*/	//MPI_Finalize();
	return 0;
}
