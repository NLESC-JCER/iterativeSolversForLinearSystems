/*
	C.E.M. Schoutrop
	"The tactical nuke of iterative solvers"
	-working concept method on 26/08/2019

	Right-preconditioning is applied here;	Ax=b->AMinvu=b.
	instead of solving for x, we solve for MInv*u
	Then at the end x can be obtained via x=MInv*x
	See Algorithm 9.5 from Iterative methods for sparse linear systems for inspiration.

	Known issues:
	-None (for now)

	Possible optimizations (PO):
	-See //PO: notes in the code

	This implementation of IDR(S)Stab(L) is based on
		1. Aihara, K., Abe, K., & Ishiwata, E. (2014). A variant of IDRstab with reliable update strategies for solving sparse linear systems. Journal of Computational and Applied Mathematics, 259, 244-258. doi:10.1016/j.cam.2013.08.028
		2. Aihara, K., Abe, K., & Ishiwata, E. (2015). Preconditioned IDRStab Algorithms for Solving Nonsymmetric Linear Systems. International Journal of Applied Mathematics, 45(3).

	Special acknowledgement to Mischa Senders for his work on an initial reference implementation of this algorithm in MATLAB.

	Based on table 1 in ref 2, for L=2, S=4 left-preconditioning should be slightly less AXPY (71.5 compared to 67.5), however it is concluded in the paper that this difference is minimal in practice.
*/

#ifndef idrstab_h
#define idrstab_h

#include <eigen3/Eigen/QR>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>

#define IDRSTAB_DEBUG_INFO 1 	//Print info to console about the problem being solved.

#if IDRSTAB_DEBUG_INFO
#include <chrono>
#endif
namespace Eigen
{

	namespace internal
	{

		template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
		bool idrstab(const MatrixType& mat, const Rhs& rhs, Dest& x,
			const Preconditioner& precond, Index& iters,
			typename Dest::RealScalar& tol_error, Index L, Index S)
		{
			#if IDRSTAB_DEBUG_INFO
			//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
			std::cout << "Matrix size: " << mat.rows() << std::endl;
			#endif
			/*
				Setup
			*/
			typedef typename Dest::Scalar Scalar;
			typedef typename Dest::RealScalar RealScalar;
			typedef Matrix<Scalar, Dynamic, 1> VectorType;
			typedef Matrix<Scalar, Dynamic, Dynamic, ColMajor> DenseMatrixTypeCol;
			typedef Matrix<Scalar, Dynamic, Dynamic, RowMajor> DenseMatrixTypeRow;

			//If the matrix is small, or L,S are too large. Limit their size to size of the matrix.
			const Index N = x.rows();
			L = L < N ? L : N;
			S = S < N ? S : N;

			Index k = 0; //Iteration counter
			const Index maxIters = iters;

			//PO: sqrNorm() saves 1 sqrt calculation
			const RealScalar tol = tol_error;
			const RealScalar rhs_norm = rhs.norm();

			//Define maximum sizes to prevent any reallocation later on
			VectorType u(N * (L + 1));
			VectorType r(N * (L + 1));
			DenseMatrixTypeCol V(N * (L + 1), S);
			DenseMatrixTypeCol U(N * (L + 1), S);
			DenseMatrixTypeCol rHat(N, L + 1);
			VectorType alpha(S);
			VectorType gamma(L);
			VectorType update(N);

			/*
				IDR(S)Stab(L) algorithm
			*/
			//Set up the initial residual
			r.block(0, 0, N, 1) = rhs - mat * precond.solve(x);
			tol_error = r.head(N).norm();

			/*
				Select an initial (N x S) matrix R0
			*/
			//Generate random R0, orthonormalize the result
			//This results in R0, however to save memory and compute we only need the adjoint of R0
			//PO: To save on memory consumption identity can be sparse
			HouseholderQR<DenseMatrixTypeCol> qr(DenseMatrixTypeCol::Random(N, S));
			const DenseMatrixTypeRow R_T = (qr.householderQ() * DenseMatrixTypeCol::Identity(N, S)).adjoint();

			//Additionally, the matrix (mat.adjoint()*R_tilde).adjoint()=R_tilde.adjoint()*mat by the anti-distributivity property of the adjoint.
			//This results in AR_T, which is constant and can be precomputed.
			const DenseMatrixTypeRow AR_T = DenseMatrixTypeRow(R_T * mat);

			#if IDRSTAB_DEBUG_INFO
			std::cout << "Check orthonormality R_T\n" <<
				R_T* R_T.adjoint() << std::endl;
			#endif

			//Determine an initial U matrix of size N x S
			for (Index q = 0; q < S; ++q)
			{
				//By default S=4, q!=0 case is more likely, this ordering is better for branch prediction.
				//Arnoldi-like process to generate a set of orthogonal vectors spanning {u,A*u,A*A*u,...,A^(S-1)*u}.
				if (q != 0)
				{
					u.head(N) = mat * precond.solve(u.head(N));
					u.head(N) -= U.topLeftCorner(N, q) * (U.topLeftCorner(N, q).adjoint() * u.head(N));
				}
				else
				{
					u.head(N) = r.head(N);
				}

				u.head(N) /= u.head(N).norm();
				U.block(0, q, N, 1) = u.head(N);
			}

			#if IDRSTAB_DEBUG_INFO
			//Columns of U should be orthonormal
			std::cout << "Check orthonormality U\n" <<
				U.block(0, 0, N, S).adjoint()*U.block(0, 0, N, S) << std::endl;
			#endif

			//PO: Is there some way to avoid storing R_T and AR_T for left-preconditioning version of Ref. 2
			//, while still keeping the performance benefit of the precompute?

			//Pre-allocate sigma, this space will be recycled without additional allocations.
			//Also construct a LU-decomposition object beforehand.
			DenseMatrixTypeCol sigma(S, S);
			FullPivLU<DenseMatrixTypeCol> lu_sigma;

			//while (tol_error > tol * rhs_norm && k < maxIters)
			while (k < maxIters)
			{

				for (Index j = 1; j <= L; ++j)
				{
					//Cache some indexing variables that occur frequently and are constant
					const Index Nj = N * j;
					const Index Nj_plus_1 = N * (j + 1);
					const Index Nj_min_1 = N * (j - 1);

					/*
						The IDR Step
					*/

					//Construction of the sigma-matrix, and the LU decomposition of sigma.
					#pragma omp parallel for num_threads(S)
					for (Index i = 0; i < S; i++)
					{
						sigma.col(i).noalias() = AR_T * precond.solve(U.block(Nj_min_1, i, N, 1));
					}
					lu_sigma.compute(sigma);

					//Obtain the update coefficients alpha
					if (j != 1)
					{
						//alpha=inverse(sigma)*(AR_T*r_{j-2})
						alpha.noalias() = lu_sigma.solve(AR_T * precond.solve(r.block(N * (j - 2), 0, N, 1)));
					}
					else
					{
						//alpha=inverse(sigma)*(R_T*r_0)
						alpha.noalias() = lu_sigma.solve(R_T * r.head(N));
					}

					//Obtain new solution and residual from this update
					update.noalias() = U.topRows(N) * alpha;
					x += update;
					r.head(N) -=  mat * precond.solve(update);

					//It is possible to early-exit here, at the cost of computing L additional dot products per cycle.
					//However by continuing the residual is expected to be lowered even further, this gives a safety margin to counteract the residual gap.
					//Based on the matrices from Ref. 2 this early-exit was deemed not worth the extra cost.

					for (Index i = 1; i <= j - 2; ++i)
					{
						//This only affects the case L>2
						r.block(N * i, 0, N, 1) -= U.block(N * (i + 1), 0, N, S) * alpha;
					}
					if (j > 1)
					{
						//r=[r;A*r_{j-2}]
						r.block(Nj_min_1, 0, N, 1).noalias() = mat * precond.solve(r.block(N * (j - 2), 0, N, 1));
					}

					for (Index q = 1; q <= S; ++q)
					{
						if (q != 1)
						{
							//u=[u_1;u_2;...;u_j]
							u.head(Nj) = u.block(N, 0, Nj, 1);
						}
						else
						{
							//u = r;
							u.head(Nj_plus_1) = r.topRows(Nj_plus_1);
						}
						//Obtain the update coefficients beta implicitly
						//beta=lu_sigma.solve(AR_T * u.block(Nj_min_1, 0, N, 1)
						u.head(Nj) -=  U.topRows(Nj) * lu_sigma.solve(AR_T * precond.solve(u.block(Nj_min_1, 0, N, 1)));

						//u=[u;Au_{j-1}]
						u.block(Nj, 0, N, 1).noalias() = mat * precond.solve(u.block(Nj_min_1, 0, N, 1));

						//Orthonormalize u_j to the columns of V_j(:,1:q-1)
						if (q > 1)
						{
							//Gram-Schmidt-like procedure to make u orthogonal to the columns of V.
							//The vector mu from Ref. 1 is obtained implicitly:
							//mu=V.block(Nj, 0, N, q - 1).adjoint() * u.block(Nj, 0, N, 1).
							u.head(Nj_plus_1) -= V.topLeftCorner(Nj_plus_1, q - 1) * (V.block(Nj, 0, N, q - 1).adjoint() * u.block(Nj, 0, N, 1));
						}

						//Normalize u and assign to a column of V
						u.head(Nj_plus_1) /= u.block(Nj, 0, N, 1).norm();
						V.block(0, q - 1, Nj_plus_1, 1) = u.head(Nj_plus_1);
						//Since the segment u.head(Nj_plus_1) is not needed next q-iteration this can be combined into:
						//V.block(0, q - 1, Nj_plus_1, 1).noalias() = u.head(Nj_plus_1) / u.block(Nj, 0, N, 1).norm();

						#if IDRSTAB_DEBUG_INFO
						std::cout << "New u should be orthonormal to the columns of V" << std::endl;
						std::cout << V.block(Nj, 0, N, q).adjoint()*u.block(Nj, 0, N, 1) << std::endl; //OK
						#endif
					}

					#if IDRSTAB_DEBUG_INFO
					//This should be identity, since the columns of V are orthonormalized
					std::cout << "This should be identity matrix:" << std::endl;
					std::cout << V.block(Nj, 0, N, S).adjoint()* V.block(Nj, 0, N, S) << std::endl;
					#endif

					U = V;
				}

				//r=[r;mat*r_{L-1}]
				//Save this in rHat, the storage form for rHat is more suitable for the argmin step than the way r is stored.
				//In Eigen 3.4 this step can be compactly done via: rHat = r.reshaped(N, L + 1);
				r.block(N * L, 0, N, 1).noalias() = mat * precond.solve(r.block(N * (L - 1), 0, N, 1));
				for (Index i = 1; i <= L + 1; ++i)
				{
					rHat.col(i - 1) = r.block(N * (i - 1), 0, N, 1);
				}

				/*
					The polynomial step
				*/
				gamma.noalias() = rHat.rightCols(L).fullPivHouseholderQr().solve(r.head(N)); //Argmin step

				//Update solution and residual using the "minimized residual coefficients"
				update.noalias() = rHat.leftCols(L) * gamma;
				x += update;
				r.head(N) -= mat * precond.solve(update);

				//Update iteration info
				++k;
				tol_error = r.head(N).norm();

				#if IDRSTAB_DEBUG_INFO
				std::cout << "\nResidual: " << std::endl;
				std::cout << tol_error / rhs_norm << std::endl;
				std::cout << "True error:      " << (rhs - mat * precond.solve(x)).norm() / rhs.norm()<< std::endl;
				std::cout << "Iter" << std::endl;
				std::cout << k << std::endl;
				#endif

				if (tol_error < tol * rhs_norm)
				{
					//Slightly early exit by moving the criterion before the update of U,
					//after the main while loop the result of that calculation would not be needed.
					iters = k;
					x = precond.solve(x);
					return true;
				}

				/*
					U=U0-sum(gamma_j*U_j)
					Consider the first iteration. Then U only contains U0, so at the start of the while-loop
					U should be U0. Therefore only the first N rows of U have to be updated.
				*/
				//PO: Is there a way to do this without producing Umatrix?
				#pragma omp parallel for num_threads(L)
				for (Index i = 1; i <= L; ++i)
				{
					DenseMatrixTypeCol Umatrix(N, S);
					Umatrix.noalias() = U.block(N * i, 0, N, S) * gamma(i - 1);

					#pragma omp critical
					U.topRows(N) -= Umatrix;
				}

			}

			iters = k;
			x = precond.solve(x);
			return true;
		}

	}

	template< typename _MatrixType,
		typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
	class IDRStab;

	namespace internal
	{

		template< typename _MatrixType, typename _Preconditioner>
		struct traits<IDRStab<_MatrixType, _Preconditioner> >
		{
			typedef _MatrixType MatrixType;
			typedef _Preconditioner Preconditioner;
		};

	}

	template< typename _MatrixType, typename _Preconditioner>
	class IDRStab : public IterativeSolverBase<IDRStab<_MatrixType, _Preconditioner> >
	{
			typedef IterativeSolverBase<IDRStab> Base;
			using Base::matrix;
			using Base::m_error;
			using Base::m_iterations;
			using Base::m_info;
			using Base::m_isInitialized;
			Index m_L = 8;
			Index m_S = 2;
		public:
			typedef _MatrixType MatrixType;
			typedef typename MatrixType::Scalar Scalar;
			typedef typename MatrixType::RealScalar RealScalar;
			typedef _Preconditioner Preconditioner;

		public:

			/** Default constructor. */
			IDRStab() : Base() {}

			/**     Initialize the solver with matrix \a A for further \c Ax=b solving.

			        This constructor is a shortcut for the default constructor followed
			        by a call to compute().

			        \warning this class stores a reference to the matrix A as well as some
			        precomputed values that depend on it. Therefore, if \a A is changed
			        this class becomes invalid. Call compute() to update it with the new
			        matrix A, or modify a copy of A.
			*/
			template<typename MatrixDerived>
			explicit IDRStab(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}

			~IDRStab() {}

			/** \internal */
			/**     Loops over the number of columns of b and does the following:
			        1. Sets the tolerance and maxIterations
			        2. Calls the function that has the core solver routine
			*/
			#if BICGSTABL_IN_LIBRARY==0
			template<typename Rhs, typename Dest>
			void _solve_with_guess_impl(const Rhs& b, Dest& x) const
			{
				_solve_vector_with_guess_impl(b, x);
			}
			#endif

			template<typename Rhs, typename Dest>
			void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const
			{
				m_iterations = Base::maxIterations();
				m_error = Base::m_tolerance;

				bool ret = internal::idrstab(matrix(), b, x, Base::m_preconditioner, m_iterations, m_error,
						m_L, m_S);
				#if IDRSTAB_DEBUG_INFO
				std::cout << "ret: " << ret << std::endl;
				std::cout << "m_error: " << m_error << std::endl;
				std::cout << "Base::m_tolerance: " << Base::m_tolerance << std::endl;
				std::cout << "True error:      " << (matrix() * x - b).norm() / b.norm()<< std::endl;
				#endif

				// if (ret == false)
				// {
				// 	m_info = NumericalIssue;
				// 	x.setZero(); //x=nan does not pass Eigen's tests even if m_info=NumericalIssue :)
				// 	m_error = ((std::isfinite)(m_error) ? m_error : 1.0);
				// }
				// else
				// {
				// 	m_info = (m_error <= Base::m_tolerance) ? Success
				// 		: NoConvergence;
				// }

				m_info = (!ret) ? NumericalIssue
					: m_error <= Base::m_tolerance ? Success
					: NoConvergence;
				//m_info=NumericalIssue;
				#if IDRSTAB_DEBUG_INFO
				std::cout << "m_error_returned: " << m_error << std::endl;
				std::cout << "m_info: " << m_info << std::endl;
				#endif
				// m_info = (!ret) ? NumericalIssue
				// 	: m_error <= Base::m_tolerance ? Success
				// 	: NoConvergence;
				m_isInitialized = true;
			}

			/** \internal */
			/** Resizes the x vector to match the dimension of b and sets the elements to zero*/
			#if BICGSTABL_IN_LIBRARY==0
			using Base::_solve_impl;
			template<typename Rhs, typename Dest>
			void _solve_impl(const MatrixBase<Rhs>& b, Dest& x) const
			{

				x.resize(this->rows(), b.cols());
				x.setZero();

				_solve_with_guess_impl(b, x);
			}
			#endif
			void setL(Index L)
			{
				if (L < 1)
				{
					L = 2;
				}

				m_L = L;
			}
			void setS(Index S)
			{
				if (S < 1)
				{
					S = 4;
				}

				m_S = S;
			}

		protected:

	};

}

#endif /* idrstab_h */

