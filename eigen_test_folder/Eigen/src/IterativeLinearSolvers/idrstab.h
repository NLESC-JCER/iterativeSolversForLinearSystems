/*
	C.E.M. Schoutrop
	"The tactical nuke of iterative solvers"
	-working concept method on 26/08/2019
	-Version that consistently passes the Eigen unit-tests on 25/09/2019

	Right-preconditioning is applied here;	Ax=b->AMinvu=b.
	instead of solving for x, we solve for MInv*u
	Then at the end x can be obtained via x=MInv*x
	See Algorithm 9.5 from Iterative methods for sparse linear systems for inspiration.

	Known issues:
	-None (for now)

	Possible optimizations (PO):
	-See //PO: notes in the code

	This implementation of IDRStab is based on
		1. Aihara, K., Abe, K., & Ishiwata, E. (2014). A variant of IDRstab with reliable update strategies for solving sparse linear systems. Journal of Computational and Applied Mathematics, 259, 244-258. doi:10.1016/j.cam.2013.08.028
		2. Aihara, K., Abe, K., & Ishiwata, E. (2015). Preconditioned IDRStab Algorithms for Solving Nonsymmetric Linear Systems. International Journal of Applied Mathematics, 45(3).
		3. Saad, Y. (2003). Iterative Methods for Sparse Linear Systems: Second Edition. Philadelphia, PA: SIAM.

	Special acknowledgement to Mischa Senders for his work on an initial reference implementation of this algorithm in MATLAB and to Adithya Vijaykumar for providing the framework for this solver.

	Based on table 1 in ref 2, for L=2, S=4 left-preconditioning should be slightly less AXPY (71.5 compared to 67.5), however it is concluded in the paper that this difference is minimal in practice.
*/

#ifndef idrstab_h
#define idrstab_h

#include <Eigen/QR>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <cmath>

#define IDRSTAB_DEBUG_INFO 2 	//Print info to console about the problem being solved.

#if IDRSTAB_DEBUG_INFO>0
#include <chrono>
#include <Eigen/Eigenvalues>
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
			#if IDRSTAB_DEBUG_INFO >0
			std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
			std::cout << "Matrix size: " << mat.rows() << std::endl;
			#endif

			/*
				Setup and type definitions.
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
			const RealScalar rhs_norm = rhs.norm();
			const RealScalar tol2 = tol_error * rhs_norm;

			if (rhs_norm == 0)
			{
				/*
					If b==0, then the exact solution is x=0.
					rhs_norm is needed for other calculations anyways, this exit is a freebie.
				*/
				x.setZero();
				tol_error = 0.0;
				return true;
			}

			//Define maximum sizes to prevent any reallocation later on.
			VectorType u(N * (L + 1));
			VectorType r(N * (L + 1));
			DenseMatrixTypeCol V(N * (L + 1), S);
			DenseMatrixTypeCol U(N * (L + 1), S);
			DenseMatrixTypeCol rHat(N, L + 1);
			VectorType alpha(S);
			VectorType gamma(L);
			VectorType update(N);

			//Booleans to signal whether or not internal variables have become unfit to continue computing.
			bool valid_u = true;
			bool valid_r = true;
			bool valid_U = true;

			/*
				Main IDRStab algorithm
			*/
			//Set up the initial residual
			r.head(N) = rhs - mat * x;
			tol_error = r.head(N).norm();

			/*
				Select an initial (N x S) matrix R0.
				1. Generate random R0, orthonormalize the result.
				2. This results in R0, however to save memory and compute we only need the adjoint of R0. This is given by the matrix R_T.
			*/
			//PO: To save on memory consumption identity can be sparse
			HouseholderQR<DenseMatrixTypeCol> qr(DenseMatrixTypeCol::Random(N, S));
			DenseMatrixTypeRow R_T = (qr.householderQ() * DenseMatrixTypeCol::Identity(N, S)).adjoint();

			/*
				Additionally, the matrix (mat.adjoint()*R_tilde).adjoint()=R_tilde.adjoint()*mat by the anti-distributivity property of the adjoint.
				This results in AR_T, which is constant and can be precomputed.
			*/
			DenseMatrixTypeRow AR_T = DenseMatrixTypeRow(R_T * mat);

			#if IDRSTAB_DEBUG_INFO >1
			std::cout << "Check orthonormality R_T\n" <<
				R_T* R_T.adjoint() << std::endl;
			#endif

			DenseMatrixTypeCol h_FOM(S, S - 1);
			h_FOM.setZero();

			/*
				Determine an initial U matrix of size N x S
			*/
			for (Index q = 0; q < S; ++q)
			{
				//By default S=4, q!=0 case is more likely, this ordering is better for branch prediction.
				//Arnoldi-like process to generate a set of orthogonal vectors spanning {u,A*u,A*A*u,...,A^(S-1)*u}.
				if (q != 0)
				{
					/*
						Original Gram-Schmidt orthogonalization strategy from Ref. 1:
					*/
					//u.head(N) -= U.topLeftCorner(N, q) * (U.topLeftCorner(N, q).adjoint() * u.head(N));

					/*
					        Modified Gram-Schmidt strategy:
						Note that GS and MGS are mathematically equivalent, they are NOT numerically equivalent.

						Eventough h is a scalar, converting the dot product to Scalar is not supported:
						http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1610
					*/
					VectorType w = mat * precond.solve(u.head(N));
					for (Index i = 0; i < q; ++i)
					{
						//"Normalization factor" (v is normalized already)
						VectorType v = U.block(0, i, N, 1);

						//"How much do v and w have in common?"
						h_FOM(i, q - 1) = v.adjoint() * w;

						//"Subtract the part they have in common"
						w = w - h_FOM(i, q - 1) * v;
					}
					u.head(N) = w;
					h_FOM(q, q - 1) = u.head(N).norm();
					if (h_FOM(q, q - 1) - h_FOM(q, q - 1) != 0.0)
					{
						/*
							If this component of the Hessenberg matrix has become NaN, then the orthonormalization cannot continue.
							This implies that it is (numerically) not possible to construct a Krylov subspace of dimension S.
							Such cases occur if:
							1. The basis of dimension <S is sufficient to exactly solve the linear system.
								I.e. the current residual is in span{r,Ar,...A^{m-1}r}, where m<S.
							2. Two vectors vectors generated from r, Ar,... have are (practically) parallel.

							In case 1, the exact solution to the system can be obtained from the "Full Orthogonalization Method" (Algorithm 6.4 in the book of Saad).
						*/
						#if IDRSTAB_DEBUG_INFO >2
						std::cout << "FOM EXIT" << std::endl;
						#endif

						/*
							Complete the FOM-algorithm
						*/
						Scalar beta = r.head(N).norm(); //This is expected to be tol_error at this point!
						VectorType e1(q);
						e1(0) = beta;
						DenseMatrixTypeCol y = h_FOM.topLeftCorner(q, q).colPivHouseholderQr().solve(e1);
						x += U.topLeftCorner(N, q) * y;

						/*
							Exit
						*/
						iters = k;
						x = precond.solve(x);
						tol_error = (rhs - mat * x).norm() / rhs_norm;

						#if IDRSTAB_DEBUG_INFO >2
						std::cout << "tol_error: " << tol_error << std::endl;
						std::cout << "x:\n" << x << std::endl;
						std::cout << "h_FOM:\n" << h_FOM << std::endl;
						std::cout << "U:\n" << U << std::endl;
						std::cout << "u:\n" << u << std::endl;
						std::cout << "q:\n" << q << std::endl;
						#endif

						return true;
					}
					if (std::abs(h_FOM(q, q - 1)) != 0.0)
					{
						/*
							This only happens if u is NOT exactly zero. In case it is exactly zero
							it would imply that that this u has no component in the direction of the current residual.

							By then setting u to zero it will not contribute any further (as it should).
							Whereas attempting to normalize results in division by zero.

							Contrary to what one would suspect, the comparison with ==0.0 for floating-point types is intended here.
							Any arbritary non-zero u is fine to continue, however if it contains either NaN or Inf the algorithm will break down.
						*/
						u.head(N) /= h_FOM(q, q - 1);
					}

				}
				else
				{
					u.head(N) = r.head(N);
					u.head(N) /= u.head(N).norm();
				}

				U.block(0, q, N, 1) = u.head(N);
			}

			#if IDRSTAB_DEBUG_INFO >1
			//Columns of U should be orthonormal
			std::cout << "Check orthonormality U\n" <<
				U.block(0, 0, N, S).adjoint()*U.block(0, 0, N, S) << std::endl;
			//h_FOM should not contain any NaNs
			std::cout << "h_FOM:\n" << h_FOM << std::endl;
			#endif

			//Pre-allocate sigma, this space will be recycled without additional allocations.
			//Also construct an LU-decomposition object beforehand.
			DenseMatrixTypeCol sigma(S, S);
			FullPivLU<DenseMatrixTypeCol> lu_sigma;

			bool reset_while = false;

			while (k < maxIters)
			{

				for (Index j = 1; j <= L; ++j)
				{
					//Cache some indexing variables that occur frequently and are constant.
					const Index Nj = N * j;
					const Index Nj_plus_1 = N * (j + 1);
					const Index Nj_min_1 = N * (j - 1);

					/*
						The IDR Step
					*/

					//Construction of the sigma-matrix, and the LU decomposition of sigma.
					for (Index i = 0; i < S; ++i)
					{
						sigma.col(i).noalias() = AR_T * precond.solve(U.block(Nj_min_1, i, N, 1));
					}
					lu_sigma.compute(sigma);

					//Obtain the update coefficients alpha
					if (j != 1)
					{
						//alpha=inverse(sigma)*(AR_T*r_{j-2})
						alpha.noalias() = lu_sigma.solve(AR_T * precond.solve(r.segment(N * (j - 2), N)));
					}
					else
					{
						//alpha=inverse(sigma)*(R_T*r_0)
						alpha.noalias() = lu_sigma.solve(R_T * r.head(N));
					}

					//Obtain new solution and residual from this update
					update.noalias() = U.topRows(N) * alpha;
					r.head(N) -=  mat * precond.solve(update);
					x += update;


					for (Index i = 1; i <= j - 2; ++i)
					{
						//This only affects the case L>2
						r.segment(N * i, N) -= U.block(N * (i + 1), 0, N, S) * alpha;
					}
					if (j > 1)
					{
						//r=[r;A*r_{j-2}]
						r.segment(Nj_min_1, N).noalias() = mat * precond.solve(r.segment(N * (j - 2), N));
					}

					if (r.head(N).norm() < tol2)
					{
						//If at this point the algorithm has converged, the orthonormalization of U will fail.
						reset_while = true;
						break;
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
						u.head(Nj) -=  U.topRows(Nj) * lu_sigma.solve(AR_T * precond.solve(u.segment(Nj_min_1, N)));

						//u=[u;Au_{j-1}]
						u.segment(Nj, N).noalias() = mat * precond.solve(u.segment(Nj_min_1, N));

						//Orthonormalize u_j to the columns of V_j(:,1:q-1)
						if (q > 1)
						{
							/*
								Original Gram-Schmidt-like procedure to make u orthogonal to the columns of V from Ref. 1.

								The vector mu from Ref. 1 is obtained implicitly:
								mu=V.block(Nj, 0, N, q - 1).adjoint() * u.block(Nj, 0, N, 1).
							*/
							//u.head(Nj_plus_1) -= V.topLeftCorner(Nj_plus_1, q - 1) * (V.block(Nj, 0, N, q - 1).adjoint() * u.segment(Nj, N));

							/*
								The same, but using MGS instead of GS
							*/
							for (Index i = 0; i <= q - 2; ++i)
							{
								//"Normalization factor"
								DenseMatrixTypeCol h2 = V.block(Nj, i, N, 1).adjoint() * V.block(Nj, i, N, 1);

								//"How much do u and V have in common?"
								DenseMatrixTypeCol h = V.block(Nj, i, N, 1).adjoint() * u.segment(Nj, N) / h2(0, 0);

								//"Subtract the part they have in common"
								u.head(Nj_plus_1) = u.head(Nj_plus_1) - h(0, 0) * V.block(0, i, Nj_plus_1, 1);
							}
						}

						//Normalize u and assign to a column of V
						Scalar normalization_constant = u.block(Nj, 0, N, 1).norm();
						if (normalization_constant != 0.0)
						{
							/*
								If u is exactly zero, this will lead to a NaN. Small, non-zero u is fine. In the case of NaN the algorithm breaks down,
								eventhough it could have continued, since u zero implies that there is no further update in a given direction.
								In exact arithmetic this would imply the algorithm converged exactly.
							*/
							u.head(Nj_plus_1) /= normalization_constant;
						}
						V.block(0, q - 1, Nj_plus_1, 1) = u.head(Nj_plus_1);
						//Since the segment u.head(Nj_plus_1) is not needed next q-iteration this may be combined into one (Only works for GS method, not MGS):
						//V.block(0, q - 1, Nj_plus_1, 1).noalias() = u.head(Nj_plus_1) / u.segment(Nj, N).norm();

						#if IDRSTAB_DEBUG_INFO >1
						std::cout << "New u should be orthonormal to the columns of V" << std::endl;
						std::cout << V.block(Nj, 0, N, q).adjoint()*u.block(Nj, 0, N, 1) << std::endl; //OK
						#endif

					}

					#if IDRSTAB_DEBUG_INFO >1
					//This should be identity, since the columns of V are orthonormalized.
					std::cout << "Check if the columns of V are orthonormalized" << std::endl;
					std::cout << V.block(Nj, 0, N, S).adjoint()* V.block(Nj, 0, N, S) << std::endl;
					#endif

					U = V;

					//It is sufficient to check if the first component has become NaN or Inf.
					valid_u = (u(0, 0) - u(0, 0) == 0.0);
					valid_U = (U(0, 0) - U(0, 0) == 0.0);
					if ((valid_u && valid_U) == false)
					{
						/*
							One of the intermediate quantities has become invalid, either because
							the algorithm has converged, or an unfixable error has occured.

							Therefore take the most recent solution and return.
						*/
						iters = k;
						x = precond.solve(x);
						tol_error = (mat * x - rhs).norm() / rhs_norm;
						std::cout << "tol_error:" << tol_error << std::endl;
						return true;
					}

				}
				if (reset_while)
				{
					reset_while = false;
					tol_error = r.head(N).norm();
					if (tol_error < tol2)
					{
						/*
							Slightly early exit by moving the criterion before the update of U,
							after the main while loop the result of that calculation would not be needed.
						*/
						break;
					}
					continue;
				}
				valid_r = r(0, 0) - r(0, 0) == 0.0;
				if (valid_r == false)
				{
					/*
						The residual vector has become invalid, therefore the polynomial step cannot be completed.
						//TODO: Of course it can, if you just use the parts that are not NaN.
					*/
					iters = k;
					x = precond.solve(x);
					tol_error = (mat * x - rhs).norm() / rhs_norm;
					return true;
				}

				//r=[r;mat*r_{L-1}]
				//Save this in rHat, the storage form for rHat is more suitable for the argmin step than the way r is stored.
				//In Eigen 3.4 this step can be compactly done via: rHat = r.reshaped(N, L + 1);
				r.segment(N * L, N).noalias() = mat * precond.solve(r.segment(N * (L - 1), N));
				for (Index i = 1; i <= L + 1; ++i)
				{
					rHat.col(i - 1) = r.segment(N * (i - 1), N);
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
				if (tol_error < tol2)
				{
					//Slightly early exit by moving the criterion before the update of U,
					//after the main while loop the result of that calculation would not be needed.
					break;
				}

				valid_U = (U(0, 0) - U(0, 0) == 0.0);
				if (valid_U == false)
				{
					/*
					        The polynomial step could be completed, however the matrix U is not valid. Therefore the next iteration cannot take place.
					*/
					iters = k;
					x = precond.solve(x);
					tol_error = (mat * x - rhs).norm() / rhs_norm;
					return true;
				}

				/*
					U=U0-sum(gamma_j*U_j)
					Consider the first iteration. Then U only contains U0, so at the start of the while-loop
					U should be U0. Therefore only the first N rows of U have to be updated.
				*/
				for (Index i = 1; i <= L; ++i)
				{
					U.topRows(N) -= U.block(N * i, 0, N, S) * gamma(i - 1);
				}

			}

			/*
				Exit after the while loop terminated.
			*/
			iters = k;
			x = precond.solve(x);
			tol_error = tol_error / rhs_norm;
			#if IDRSTAB_DEBUG_INFO >0
			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>
				(t2 - t1);
			std::cout << "Solver time: " << time_span.count() << " seconds" << std::endl;
			std::cout << "#iterations:     " << k << std::endl;
			std::cout << "Estimated relative residual: " << tol_error << std::endl;
			std::cout << "True relative residual:      " << (mat * x - rhs).norm() / rhs.norm() << std::endl;
			#endif

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
			Index m_L = 2;
			Index m_S = 4;
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
			// template<typename Rhs, typename Dest>
			// void _solve_with_guess_impl(const Rhs& b, Dest& x) const
			// {
			// 	_solve_vector_with_guess_impl(b, x);
			// }
			// using Base::_solve_impl;
			// template<typename Rhs, typename Dest>
			// void _solve_impl(const MatrixBase<Rhs>& b, Dest& x) const
			// {

			// 	x.resize(this->rows(), b.cols());
			// 	x.setZero();

			// 	_solve_with_guess_impl(b, x);
			// }

			template<typename Rhs, typename Dest>
			void _solve_vector_with_guess_impl(const Rhs& b, Dest& x) const
			{
				m_iterations = Base::maxIterations();
				m_error = Base::m_tolerance;
				bool ret = internal::idrstab(matrix(), b, x, Base::m_preconditioner, m_iterations, m_error,
						m_L, m_S);

				m_info = (!ret) ? NumericalIssue
					: m_error <= Base::m_tolerance ? Success
					: NoConvergence;

				#if IDRSTAB_DEBUG_INFO >0
				std::cout << "ret: " << ret << std::endl;
				std::cout << "m_error: " << m_error << std::endl;
				std::cout << "Base::m_tolerance: " << Base::m_tolerance << std::endl;
				std::cout << "m_info: " << m_info << std::endl;
				#endif


			}

			/** \internal */
			/** Resizes the x vector to match the dimension of b and sets the elements to zero*/

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

