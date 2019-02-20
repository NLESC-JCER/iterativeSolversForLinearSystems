//
//  bicgstabell.h
//  conjugateGradient
//
//  Created by Adithya Vijaykumar on 19/02/2019.
//  Copyright © 2019 Adithya Vijaykumar. All rights reserved.
//

#ifndef bicgstabell_h
#define bicgstabell_h

namespace Eigen {
    
    namespace internal {
        
        template<typename MatrixType, typename Rhs, typename Dest, typename Preconditioner>
        bool bicgstabell(const MatrixType& mat, const Rhs& rhs, Dest& xguess,
                                const Preconditioner& precond, Index& iters,
                         typename Dest::RealScalar& tol_error, typename Dest::RealScalar& l)
        {
            using std::sqrt;
            using std::abs;
            typedef typename Dest::RealScalar RealScalar;
            typedef typename Dest::Scalar Scalar;
            
            //start with k=-l or k+l=0
            RealScalar L = l;
            RealScalar tol = tol_error;
            Index maxIters = iters;
            
            RealScalar k = -L;
            
            typedef Matrix<Scalar,Dynamic,1> VectorType;
           
            //We start with an initial guess x_0 and let us set the shadow residual as r_0 (residual calculated from x_0)
            VectorType x0 = xguess;
            VectorType r0  = rhs - mat * x0; //r_0
            VectorType u0;
            u0.setZero();
            VectorType rShadow = r0;
            
            VectorType x = x0;
            VectorType r = r0;
            VectorType u = u0;
            
            RealScalar deltaNew = r0.squaredNorm();//r.r
            RealScalar delta0 = deltaNew;
            RealScalar rhs_sqnorm = rhs.squaredNorm();
            
            //Other vectors and scalars initialisation
            
            RealScalar rho0 = 1;
            RealScalar alpha = 0;
            RealScalar omega = 0;
            
            if(rhs_sqnorm == 0)
            {
                x.setZero();
                return true;
            }
            
            
            RealScalar tol2 = tol*tol*delta0;
            Index i = 0;
            
            while ( r.squaredNorm() > tol2 && i<iters )
            {
                k = k+L;
                rho0 = -omega*rho0;
                
                VectorType rJloop = r;
                for(Index j=0;j<=L-1;++j)
                {
                    RealScalar rho1 = rJloop.dot(rShadow);
                    RealScalar beta = alpha * (rho1/rho0);
                    rho0 = rho1;
                    VectorType rIloop = r;
                    for(Index i=0; i<=j; ++i)
                    {
                        
                    }
                    
                }
                
                
            }
            
            tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);
            iters = i;
            return true;
        }
        
    }
    
    template< typename _MatrixType,
    typename _Preconditioner = DiagonalPreconditioner<typename _MatrixType::Scalar> >
    class BicgstabEll;
    
    namespace internal {
        
        template< typename _MatrixType, typename _Preconditioner>
        struct traits<BicgstabEll<_MatrixType,_Preconditioner> >
        {
            typedef _MatrixType MatrixType;
            typedef _Preconditioner Preconditioner;
        };
        
    }
    
    template< typename _MatrixType, typename _Preconditioner>
    class BicgstabEll : public IterativeSolverBase<BicgstabEll<_MatrixType,_Preconditioner> >
    {
        typedef IterativeSolverBase<BicgstabEll> Base;
        using Base::matrix;
        using Base::m_error;
        using Base::m_iterations;
        using Base::m_info;
        using Base::m_isInitialized;
    public:
        typedef _MatrixType MatrixType;
        typedef typename MatrixType::Scalar Scalar;
        typedef typename MatrixType::RealScalar RealScalar;
        typedef _Preconditioner Preconditioner;
        
    public:
        
        /** Default constructor. */
        BicgstabEll() : Base() {}
        
        /** Initialize the solver with matrix \a A for further \c Ax=b solving.
         *
         * This constructor is a shortcut for the default constructor followed
         * by a call to compute().
         *
         * \warning this class stores a reference to the matrix A as well as some
         * precomputed values that depend on it. Therefore, if \a A is changed
         * this class becomes invalid. Call compute() to update it with the new
         * matrix A, or modify a copy of A.
         */
        template<typename MatrixDerived>
        explicit BicgstabEll(const EigenBase<MatrixDerived>& A) : Base(A.derived()) {}
        
        ~BicgstabEll() {}
        
        /** \internal */
        /** Loops over the number of columns of b and does the following:
         1. sets the tolerence and maxIterations
         2. Calls the function that has the core solver routine
         */
        template<typename Rhs,typename Dest>
        void _solve_with_guess_impl(const Rhs& b, Dest& x) const
        {
            bool failed = false;
            for(Index j=0; j<b.cols(); ++j)
            {
                m_iterations = Base::maxIterations();
                //******************MANUALLY SET NUM ITERATIONS
                //m_iterations = 30;
                m_error = Base::m_tolerance;
                
                typename Dest::ColXpr xj(x,j);
                if(!internal::bicgstabell(matrix(), b.col(j), xj, Base::m_preconditioner, m_iterations, m_error))
                    failed = true;
            }
            m_info = failed ? NumericalIssue
            : m_error <= Base::m_tolerance ? Success
            : NoConvergence;
            m_isInitialized = true;
        }
        
        /** \internal */
        /** Resizes the x vector to match the dimenstion of b and sets the elements to zero*/
        using Base::_solve_impl;
        template<typename Rhs,typename Dest>
        void _solve_impl(const MatrixBase<Rhs>& b, Dest& x) const
        {
            x.resize(this->rows(),b.cols());
            x.setZero();
            _solve_with_guess_impl(b,x);
        }
        
    protected:
        
    };
    
}

#endif /* bicgstabell_h */

