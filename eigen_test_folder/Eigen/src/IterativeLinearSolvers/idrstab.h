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

#if IDRSTAB_DEBUG_INFO
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
			const RealScalar rhs_norm = rhs.norm();
			const RealScalar tol2 = tol_error * rhs_norm;

			if (rhs_norm == 0)
			{
				x.setZero();
				return true;
			}

			//Define maximum sizes to prevent any reallocation later on
			VectorType u(N * (L + 1));
			VectorType r(N * (L + 1));
			DenseMatrixTypeCol V(N * (L + 1), S);
			DenseMatrixTypeCol U(N * (L + 1), S);
			DenseMatrixTypeCol rHat(N, L + 1);
			VectorType alpha(S);
			VectorType gamma(L);
			VectorType update(N);

			//VectorType h(1);
			/*
				IDR(S)Stab(L) algorithm
			*/
			//Set up the initial residual
			r.head(N) = rhs - mat * x;
			tol_error = r.head(N).norm();

			/*
				Select an initial (N x S) matrix R0
			*/
			//Generate random R0, orthonormalize the result
			//This results in R0, however to save memory and compute we only need the adjoint of R0
			//PO: To save on memory consumption identity can be sparse
			HouseholderQR<DenseMatrixTypeCol> qr(DenseMatrixTypeCol::Random(N, S));
			DenseMatrixTypeRow R_T = (qr.householderQ() * DenseMatrixTypeCol::Identity(N, S)).adjoint();

			//Additionally, the matrix (mat.adjoint()*R_tilde).adjoint()=R_tilde.adjoint()*mat by the anti-distributivity property of the adjoint.
			//This results in AR_T, which is constant and can be precomputed.
			DenseMatrixTypeRow AR_T = DenseMatrixTypeRow(R_T * mat);

			#if IDRSTAB_DEBUG_INFO >1
			std::cout << "Check orthonormality R_T\n" <<
				R_T* R_T.adjoint() << std::endl;
			#endif

			bool valid_u=true;
			bool valid_r=true;
			bool valid_V=true;
			bool valid_U=true;

			DenseMatrixTypeCol h_FOM(S,S-1);
			h_FOM.setZero();
			U.setZero(); //PO: Not needed, only for testing
			//Determine an initial U matrix of size N x S
			for (Index q = 0; q < S; ++q)
			{
				//By default S=4, q!=0 case is more likely, this ordering is better for branch prediction.
				//Arnoldi-like process to generate a set of orthogonal vectors spanning {u,A*u,A*A*u,...,A^(S-1)*u}.
				if (q != 0)
				{
					/*
					        Gram-Schmidt orthogonalization:
					*/
					//u.head(N) -= U.topLeftCorner(N, q) * (U.topLeftCorner(N, q).adjoint() * u.head(N));
					/*
					        Modified Gram-Schmidt method:
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
						//DenseMatrixTypeCol h2 = v.adjoint() * w;
						h_FOM(i,q-1)=v.adjoint() * w;

						//"Subtract the part they have in common"
						//w = w - h2(0, 0) * v;
						w = w - h_FOM(i, q-1) * v;
					}
					u.head(N) = w;
					h_FOM(q+1-1,q-1)=u.head(N).norm();
					//if(h_FOM(q+1-1,q-1)-h_FOM(q+1-1,q-1) !=0.0 || std::abs(h_FOM(q+1-1,q-1))<1e-16)
					//if(h_FOM(q+1-1,q-1)-h_FOM(q+1-1,q-1) !=0.0 || std::abs(h_FOM(q+1-1,q-1))==0.0)
					if(h_FOM(q+1-1,q-1)-h_FOM(q+1-1,q-1) !=0.0)
					//TODO: SET TOL TO SQRT EPS,
					//If the orthogonalization cannot continue, exit and use the FOM algorithm to obtain a best estimate.
					{
						std::cout << "FOM EXIT"<<std::endl;
						//Apply the FOM algorithm and exit
						//Not happy with this criterion
						Scalar beta=r.head(N).norm(); //This is expected to be tol_error at this point!
						VectorType e1(q);
						e1(0)=1.0;
						DenseMatrixTypeCol y=h_FOM.block(0,0,q,q).colPivHouseholderQr().solve(beta*e1);
						x=x+U.block(0, 0, N, q)*y;
						iters = k;
						x = precond.solve(x);
						tol_error = (mat * x - rhs).norm() / rhs_norm;
						std::cout << "tol_error: " <<tol_error<< std::endl;
						std::cout << "x:\n" <<x<< std::endl;
						std::cout<< "h_FOM:\n"<< h_FOM<<std::endl;
						std::cout<< "U:\n"<< U<<std::endl;
						std::cout<< "u:\n"<< u<<std::endl;
						std::cout<< "q:\n"<< q<<std::endl;
						return true;
					}
					if(std::abs(h_FOM(q+1-1,q-1))!=0.0){
						//This only happens if u is exactly zero.
						u.head(N) /= h_FOM(q+1-1,q-1);
					}

				}
				else
				{
					u.head(N) = r.head(N);
					u.head(N)/=u.head(N).norm();
				}

				U.block(0, q, N, 1) = u.head(N);
			}
/*
TODO: For some matrices it is not possible to span a basis with dimension S
In that case U.block(0, 0, N, S).adjoint()*U.block(0, 0, N, S) will not be identity
but partially NaN.
The dimension that should be spannable should also indicate the amount of iterations needed to some extent. If the solution exists in a Krylov subspace consisting of r and Ar, then there is no need
to also check A^2r!

Het is als de Krylov subspace een lineair afhankelijke vectoren krijgt, waardoor het een dimensie lager dan S opspant

This seems to be the final remaining issue :D
*/
			#if IDRSTAB_DEBUG_INFO >1
			//Columns of U should be orthonormal
			std::cout << "Check orthonormality U\n" <<
				U.block(0, 0, N, S).adjoint()*U.block(0, 0, N, S) << std::endl;
			std::cout<< "h_FOM:\n"<< h_FOM<<std::endl;
			#endif
			/*
			DenseMatrixTypeRow titanic = U.block(0, 0, N, S).adjoint()*U.block(0, 0, N, S);
			//Scalar tol_titanic=1e5;
			if(std::abs(titanic.block(0,0,1,S).norm())>1+1e-2 || titanic.block(0,0,1,S).norm()-titanic.block(0,0,1,S).norm()!=0.0)
			{
				std::cout<<"Titanic 4"<<std::endl;

				//	Perform argmin, return exact result
				//	eth3.m part with U.block(0, 0, N, argmin_L) as p


				//Figure out which part is bad
				Index argmin_L=0;
				for(Index t=0;t<S;++t){
					if(std::abs(titanic(t,t))>1e-2){
						argmin_L++;
					}
				}
				DenseMatrixTypeCol B(N, argmin_L);
				for (Index i = 0; i < argmin_L; ++i)
				{
					B.col(i).noalias() = mat * precond.solve(U.block(0, i, N, 1));
				}

				VectorType delta=B.fullPivHouseholderQr().solve(r.head(N));
				x=x+U.block(0, 0, N, argmin_L)*delta;
				iters = k;
				x = precond.solve(x);
				tol_error = (mat * x - rhs).norm() / rhs_norm;
				return true;
			}
			*/

/*
men kan het ook via householder rank revealing ofzo
*/

			//Pre-allocate sigma, this space will be recycled without additional allocations.
			//Also construct a LU-decomposition object beforehand.
			DenseMatrixTypeCol sigma(S, S);
			FullPivLU<DenseMatrixTypeCol> lu_sigma;

			bool reset_while = false;

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


					// for(auto &a: alpha){
					// 	if(a-a!=0.0){
					// 		a=0;
					// 	}
					// }
					if (alpha.norm()-alpha.norm()!=0.0)
					{
						/*
							Everything will fail, take the best estimate for x and abandon ship.
						*/
						std::cout << "TITANIC 305" << std::endl;

						//Obtain new solution and residual from this update
						for(auto &a: alpha){
							if(a-a!=0.0){
								a=0;
							}
						}
						update.noalias() = U.topRows(N) * alpha;
						r.head(N) -=  mat * precond.solve(update);
						x += update;


						iters = k;
						x = precond.solve(x);
						tol_error = (mat * x - rhs).norm() / rhs_norm;
						std::cout << "tol_error:" << tol_error<< std::endl;
						return true;
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
					//It is possible to early-exit here, at the cost of computing L additional dot products per cycle.
					//However by continuing the residual is expected to be lowered even further, this gives a safety margin to counteract the residual gap.
					//Based on the matrices from Ref. 2 this early-exit was deemed not worth the extra cost.
					//one has to early exist here, else u becomes NaN
					if (r.head(N).norm() < tol2)
					{
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
								Gram-Schmidt-like procedure to make u orthogonal to the columns of V.
								The vector mu from Ref. 1 is obtained implicitly:
								mu=V.block(Nj, 0, N, q - 1).adjoint() * u.block(Nj, 0, N, 1).
							*/
							//u.head(Nj_plus_1) -= V.topLeftCorner(Nj_plus_1, q - 1) * (V.block(Nj, 0, N, q - 1).adjoint() * u.segment(Nj, N));

							/*
								The same, but using MGS instead of GS!
								This results in NaN if uhead is zero. This only happens if something went wrong in u already
								I.e. when the residual is zero!
							*/
							u.head(Nj_plus_1) /= u.head(Nj_plus_1).norm(); //This is not strictly necessary
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
						u.head(Nj_plus_1) /= u.block(Nj, 0, N, 1).norm();
						V.block(0, q - 1, Nj_plus_1, 1) = u.head(Nj_plus_1);
						//Since the segment u.head(Nj_plus_1) is not needed next q-iteration this may be combined into (Only works for GS method, not MGS):
						//V.block(0, q - 1, Nj_plus_1, 1).noalias() = u.head(Nj_plus_1) / u.segment(Nj, N).norm();

						// #if IDRSTAB_DEBUG_INFO >1
						// std::cout << "New u should be orthonormal to the columns of V" << std::endl;
						// std::cout << V.block(Nj, 0, N, q).adjoint()*u.block(Nj, 0, N, 1) << std::endl; //OK
						// DenseMatrixTypeCol h=V.block(Nj, 0, N, q).adjoint()*u.block(Nj, 0, N, 1);
						// if(h(0,0)-h(0,0)!=0.0)
						// {
						// 	std::cout << "mat:\n" << mat << std::endl;
						// 	std::cout << "rhs:\n" << rhs << std::endl;
						// 	std::cout << "u:\n" << u << std::endl;
						// 	std::cout << "V:\n" << V << std::endl;
						// 	std::cout << "r:\n" << r << std::endl;
						// 	std::cout << "q:" << q << std::endl;
						// 	std::cout << "j:" << j << std::endl;
						// 	std::cout << "k:" << k << std::endl;
						// 	DenseMatrixTypeCol mat_dense=DenseMatrixTypeCol(mat);
						// 	ComplexEigenSolver<DenseMatrixTypeCol> es(mat_dense);
						// 	std::cout << "The eigenvalues of A are:" << std::endl << es.eigenvalues() << std::endl;
						// 	std::cout << "The matrix of eigenvectors, V, is:" << std::endl << es.eigenvectors() << std::endl << std::endl;
						// 	std::cout << "x:\n" << x << std::endl;
						// 	x = precond.solve(x);
						// 	std::cout << "x:\n" << x << std::endl;
						// 	std::cout << "True relative residual:      " << (mat * x - rhs).norm() / rhs.norm() << std::endl;
						// 	*(int*)0 = 0;
						// }

						// #endif
					}

					#if IDRSTAB_DEBUG_INFO >1
					//This should be identity, since the columns of V are orthonormalized
					std::cout << "This should be identity matrix:" << std::endl;
					std::cout << V.block(Nj, 0, N, S).adjoint()* V.block(Nj, 0, N, S) << std::endl;
					#endif

					U = V;

					valid_u=u(0,0)-u(0,0)==0.0;
					valid_r=r(0,0)-r(0,0)==0.0;
					valid_V=V(0,0)-V(0,0)==0.0;
					valid_U=U(0,0)-U(0,0)==0.0;
					if ((valid_u && valid_r && valid_V && valid_U)==false)
					{
						/*
							Everything will fail, take the best estimate for x and abandon ship.
						*/
						std::cout << "TITANIC 410" << std::endl;
						iters = k;
						x = precond.solve(x);
						tol_error = (mat * x - rhs).norm() / rhs_norm;
						std::cout << "tol_error:" << tol_error<< std::endl;
						return true;
					}

				}
				if (reset_while)
				{
					tol_error = r.head(N).norm();
					if (tol_error < tol2)
					{
						//Slightly early exit by moving the criterion before the update of U,
						//after the main while loop the result of that calculation would not be needed.
						break;
					}
					continue;
				}
				//bool valid_u=u(0,0)-u(0,0)==0.0;
				valid_r=r(0,0)-r(0,0)==0.0;
				//bool valid_U=U(0,0)-U(0,0)==0.0;
				if (valid_r==false)
				{
					/*
						Everything will fail, take the best estimate for x and abandon ship.
					*/
					std::cout << "TITANIC 427" << std::endl;
					iters = k;
					x = precond.solve(x);
					tol_error = (mat * x - rhs).norm() / rhs_norm;
					std::cout << "tol_error:" << tol_error<< std::endl;
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

				if(gamma.norm()-gamma.norm()!=0.0){
					for(auto &g:gamma){
						if(g-g!=0.0){
							g=0.0;
						}
					}
				}

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

				valid_u=u(0,0)-u(0,0)==0.0;
				valid_r=r(0,0)-r(0,0)==0.0;
				valid_V=V(0,0)-V(0,0)==0.0;
				valid_U=U(0,0)-U(0,0)==0.0;
				if ((valid_u && valid_r && valid_V && valid_U)==false)
				{
					/*
						Everything will fail, take the best estimate for x and abandon ship.
					*/
					std::cout << "TITANIC 524" << std::endl;
					iters = k;
					x = precond.solve(x);
					tol_error = (mat * x - rhs).norm() / rhs_norm;
					std::cout << "tol_error:" << tol_error<< std::endl;
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
			iters = k;
			x = precond.solve(x);
			tol_error = tol_error / rhs_norm;
			#if IDRSTAB_DEBUG_INFO >0
			//Print experimental info
			std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>
				(t2 - t1);
			std::cout << "Solver time: " << time_span.count() << " seconds" << std::endl;
			std::cout << "#iterations:     " << k << std::endl;
			std::cout << "Estimated relative residual: " << tol_error << std::endl;
			std::cout << "True relative residual:      " << (mat * x - rhs).norm() / rhs.norm() << std::endl;
			// tol_error = (mat * x - rhs).norm() / rhs.norm();
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
				//std::cout << "ret: " << ret << std::endl;
				//std::cout << "m_error: " << m_error << std::endl;
				//std::cout << "Base::m_tolerance: " << Base::m_tolerance << std::endl;
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

