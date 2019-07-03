#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/SparseLU>
#include <chrono>
#include <ctime>
#include <fstream>
#include <ratio>
#include <string>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/IterativeLinearSolvers>
#include <vector>
#include "bicgstabell.h"
#include <complex>
#include <type_traits>

/*
        Goal of script:

        Script to process a series of linear systems Ax=b for benchmarking purposes.
        A and b are read from file, the solution vector x is computed either iteratively
        via Eigen::BiCGStab, or via Eigen::BiCGStabL. Diagnostic information produced by
        the solver is witten to an output file.
*/

//#define MAJORITY Eigen::ColMajor
/*
        ColMajor is the default in Eigen
        Direct does not work with RowMajor
        BiCGStab is able to work in parallel for RowMajor, but not for ColMajor.
*/
//#define DEBUG false //true: print additional information for debugging purposes
#define SAVEX false //true: save computed result to file
//#define NUMTHREADS 2 //0: use default option

std::vector<std::string> read_inputfile(std::string filename)
{
	/*
	        Reads in an inputfile, returns the input data as a vector of strings.
	        inputfile is assumed to be one input on every new line.
	*/
	std::ifstream inputfile(filename);
	std::vector<std::string> input_list;

	if (inputfile.is_open())
	{
		std::string line;

		while (std::getline(inputfile, line))
		{
			input_list.push_back(line);

		}

		inputfile.close();
		return input_list;
	}
	else
	{
		std::cout << "Input file cannot be opened" << std::endl;
		exit(1);
	}
}
class Result
{
		/*
		        Keeps track of the solver performance for a given input.
		*/
	public:
		Result(std::string A_input_filename, std::string b_input_filename)
		{
			m_list_of_A = read_inputfile(A_input_filename);
			m_list_of_b = read_inputfile(b_input_filename);

			m_amount_of_results = m_list_of_A.size();
			m_time.resize(m_amount_of_results);
			m_residual.resize(m_amount_of_results);
			m_size.resize(m_amount_of_results);
			m_nnz.resize(m_amount_of_results);
			m_iterations.resize(m_amount_of_results);
			m_info.resize(m_amount_of_results);

			#if DEBUG==1
			{
				/*
				        Print out the input info to make sure the correct thing is read in.
				*/
				std::cout << "A_input_filename: " << A_input_filename << std::endl;

				for (int iter_A = 0; iter_A < m_list_of_A.size(); ++iter_A)
				{
					std::cout << m_list_of_A[iter_A] << std::endl;
				}

				std::cout << "b_input_filename: " << b_input_filename << std::endl;

				for (int iter_A = 0; iter_A < m_list_of_A.size(); ++iter_A)
				{
					std::cout << m_list_of_b[iter_A] << std::endl;
				}
			}
			#endif
		}
		void print() const
		{
			/*
			        Print the results to screen
			*/
			std::cout << "A\t"
				<< "b\t"
				<< "size\t"
				<< "nnz\t"
				<< "time\t";
			std::cout << "iterations\t";
			std::cout << "||Ax-b||_2/||b||_2" << "\t";
			std::cout << "info" << std::endl;

			for (unsigned iter_A = 0; iter_A < m_amount_of_results; ++iter_A)
			{
				std::cout << m_list_of_A[iter_A] << "\t"
					<< m_list_of_b[iter_A] << "\t"
					<< m_size[iter_A] << "\t"
					<< m_nnz[iter_A] << "\t"
					<< m_time[iter_A] << "\t";
				std::cout << m_iterations[iter_A] << "\t";
				std::cout << m_residual[iter_A] << "\t";
				std::cout << m_info[iter_A] << std::endl;
			}
		}
		void save_to_file(std::string result_output_filename) const
		{
			/*
			        Save the results to a file given by result_output_filename.
			*/
			std::ofstream output_file(result_output_filename);
			output_file << "A\t"
				<< "b\t"
				<< "size\t"
				<< "nnz\t"
				<< "time\t";
			output_file << "iterations\t";
			output_file << "||Ax-b||_2/||b||_2" << "\t";
			output_file << "info" << std::endl;

			if (output_file.is_open())
			{
				for (unsigned iter_A = 0; iter_A < m_amount_of_results; ++iter_A)
				{
					output_file << m_list_of_A[iter_A] << "\t"
						<< m_list_of_b[iter_A] << "\t"
						<< m_size[iter_A] << "\t"
						<< m_nnz[iter_A] << "\t"
						<< m_time[iter_A] << "\t";
					output_file << m_iterations[iter_A] << "\t";
					output_file << m_residual[iter_A] << "\t";
					output_file << m_info[iter_A] << std::endl;
				}

				output_file.close();
			}
			else
			{
				std::cout << "Output file could not be opened!" << std::endl;
			}
		}
		void get_current_A_b(size_t iter_A, std::string& current_A, std::string& current_b) const
		{
			current_A = m_list_of_A[iter_A];
			current_b = m_list_of_b[iter_A];
		}

		size_t m_amount_of_results;
		std::vector<double> m_time;
		std::vector<double> m_residual;
		std::vector<double> m_size;
		std::vector<double> m_nnz;
		std::vector<double> m_iterations;
		std::vector<double> m_info;
		Eigen::Matrix<double, Eigen::Dynamic, 1> m_real_previous;
		Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> m_complex_previous;
	private:
		std::vector<std::string> m_list_of_A;
		std::vector<std::string> m_list_of_b;

};
template<typename T>
void solve_system(Result& result, size_t iter_A, bool recycle_previous)
{
	/*
	        IN:
	        result:			reference to result object to store the performance details
	        iter_A			number of testcase to process
	        recycle_previous 	If the input is a sequence of matrices, this flag indicates if the previous solution vector is used as an initial guess to solve the current system
	*/
	typedef Eigen::SparseMatrix<T, MAJORITY> SparseMatrix;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> ColumnVector;
	bool solver_succes = true;
	SparseMatrix A;
	ColumnVector b, x, tmp;
	std::string current_A, current_b;
	result.get_current_A_b(iter_A, current_A, current_b);

	/*
		Load in A and b and do some basic consistency checks.
	*/
	std::cout << "Loading: " << current_A << std::endl;
	Eigen::loadMarket(A, current_A);

	if (A.rows() != A.cols())
	{
		std::cout << "Matrix is not square" << std::endl;
		solver_succes = false;
	}

	std::cout << "#Rows: " << A.rows() << " #Cols: " << A.cols() << std::endl;
	std::cout << "Loading: " << current_b << std::endl;
	Eigen::loadMarketVector(b, current_b);

	if (A.cols() != b.size())
	{
		std::cout << "Matrix and vector incompatible" << std::endl;
		solver_succes = false;
	}

	/*
		Setup of the iterative solver
	*/
	#if ALGORITHM==1
	Eigen::BiCGSTAB<SparseMatrix, PRECONDITIONER> solver;
	//Eigen::BiCGSTAB<SparseMatrix, Eigen::IdentityPreconditioner> solver;
	//Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<T>> solver;
	//Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<T>> solver;
	//solver.preconditioner().setFillfactor(2);
	std::cout << "Using BiCGSTAB" << std::endl;
	#elif ALGORITHM==2
	Eigen::BicgstabEll<SparseMatrix, PRECONDITIONER> solver;
	//Eigen::BicgstabEll<SparseMatrix, Eigen::IdentityPreconditioner> solver;
	//Eigen::BicgstabEll<SparseMatrix, Eigen::DiagonalPreconditioner<T>> solver;
	//Eigen::BicgstabEll<SparseMatrix, Eigen::IncompleteLUT<T>> solver;
	//solver.preconditioner().setFillfactor(10);
	const int L = ELL;
	std::cout << "Using BiCGStab(" << L << ")" << std::endl;
	#if DEBUG==1
	std::cout << "Printing debug information" << std::endl;
	#endif
	#endif
	solver.setTolerance(1e-12);
	solver.setMaxIterations(2000);

	#if DEBUG==1
	std::cout << "A:\n" << A << std::endl;
	std::cout << "b:\n" << b << std::endl;
	#endif

	std::cout << "Starting the computation" << std::endl;
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

	#if ALGORITHM==2
	solver.setL(L);
	#endif
	/*
		"compute" step, at this stage the preconditioner would be set up.
		If for some reason this fails, solver.info() will not be Eigen::Success
	*/
	solver.compute(A);

	/*
		Solve the system
	*/
	if (solver.info() == Eigen::Success)
	{
		if (recycle_previous == false)
		{
			std::cout << "Starting without previous solution" << std::endl;
			x = solver.solve(b);
		}
		else
		{
			std::cout << "Recycling previous solution" << std::endl;

			if constexpr (std::is_same_v<T, double>)
			{
				x = solver.solveWithGuess(b, result.m_real_previous);
			}
			else if constexpr (std::is_same_v<T, std::complex<double>>)
			{
				x = solver.solveWithGuess(b, result.m_complex_previous);
			}

		}
	}
	else
	{
		std::cout << "The solver failed " << std::endl;
		solver_succes = false;
	}

	#if DEBUG==1
	std::cout << "x:\n" << x << std::endl;
	#endif


	/*
		Stop the clock, print out performance information
	*/
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>
		(t2 - t1);
	std::cout << "Solver time: " << time_span.count() << " seconds" << std::endl;
	std::cout << "estimated error: " << solver.error() << std::endl;
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "info:" << solver.info() << std::endl;
	ColumnVector residual = b - A * x;
	double relative_residual = residual.norm() / b.norm();
	std::cout << "||Ax-b||_2/||b||_2 : " << relative_residual << std::endl;

	/*
		Save the results into the result object.
	*/
	if (solver_succes)
	{
		result.m_time[iter_A] = time_span.count();
		result.m_residual[iter_A] = relative_residual;
	}
	else
	{
		result.m_time[iter_A] = -1;
		result.m_residual[iter_A] = -1;
	}

	result.m_info[iter_A] = solver.info();
	result.m_size[iter_A] = A.cols();
	result.m_nnz[iter_A] = A.nonZeros();
	result.m_iterations[iter_A] = solver.iterations();
	#if SAVEX
	Eigen::saveMarket(x, "x_" + std::to_string(
			iter_A)); //The output precision is hardcoded in the Eigen library
	#endif

	/*
		This construction is to make sure the recycling of the previous solution uses the correct type.
	*/
	if constexpr (std::is_same_v<T, double>)
	{
		result.m_real_previous = x;
	}
	else if constexpr (std::is_same_v<T, std::complex<double>>)
	{
		result.m_complex_previous = x;
	}

}

int main(int argc, char* argv[])
{

	std::cout << "Number of input arguments: " << argc << std::endl;

	if (argc < 5 || argc > 6)
	{
		std::cout <<
			"Input is expected as: A_input(file) b_input(file) result_output(file) recycle(integer 0 or 1)";
		exit(1);
	}

	if (argc == 6)
	{
		int number_of_threads = atoi(argv[5]);
		Eigen::setNbThreads(number_of_threads);
		std::cout << "Number of Eigen's threads set to: " << number_of_threads << "" << std::endl;

		if (number_of_threads == 0)
		{
			std::cout << "number_of_threads is set to 0, this results in the default number of threads" <<
				std::endl;
		}
	}
	else
	{
		std::cout << "Using the default number of openMP threads: " << Eigen::nbThreads( ) << std::endl;
	}

	if (MAJORITY == Eigen::ColMajor)
	{
		std::cout << "Using the column major format" << std::endl;
	}
	else if (MAJORITY == Eigen::RowMajor)
	{
		std::cout << "Using the row major format " << std::endl;
	}


	std::string A_input_filename(argv[1]);
	std::string b_input_filename(argv[2]);
	std::string result_output_filename = argv[3];
	int recycle_x0 = atoi(argv[4]);
	/*
	        if recycle_x0 !=0, the computed x from the previous linear system will be used as
	        an initial guess for the next linear system.
	*/

	Result result(A_input_filename, b_input_filename);

	for (unsigned iter_A = 0; iter_A < result.m_amount_of_results; ++iter_A)
	{

		std::string current_A, current_b;
		result.get_current_A_b(iter_A, current_A, current_b);
		bool iscomplex, isvector;
		int sym;
		Eigen::getMarketHeader(current_A, sym, iscomplex, isvector);

		bool recycle_previous = (iter_A != 0 && recycle_x0 == 1) ? true : false;

		if (iscomplex)
		{
			std::cout << "Input matrix is complex" << std::endl;
			solve_system<std::complex<double>>(result, iter_A, recycle_previous);
		}
		else
		{
			std::cout << "input matrix is real" << std::endl;
			solve_system<double>(result, iter_A, recycle_previous);
		}

	}


	result.print();
	result.save_to_file(result_output_filename);

	return 0;
}
