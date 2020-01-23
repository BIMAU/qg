#include "QG.hpp"

#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <iostream>

namespace QG
{

    QG::QG(int m, int n, int perio)
        :
        m_(m),
        n_(n),
        nun_(2),
        ndim_(nun_*m*n),
        adim_(25*ndim_),
        hdim_(6.0e+02),
        f0dim_(1.0e-04),
        beta0dim_(1.6e-11),
        Lxdim_(1.0e+06),
        Lydim_(1.0e+06),
        udim_(1.6e-02),
        gdim_(2.0e-02),
        rhodim_(1.0e+03),
        taudim_(1.5e-01),
        Ahdim_(1.0e+03),
        bfdim_(0.0e+00),
        xmin_(0.0),
        xmax_(1.0),
        ymin_(0.0),
        ymax_(Lydim_/Lxdim_),
        periodic_(perio),
        x_(n),
        y_(m),
        tx_(n,m),
        ty_(n,m),
        par_(21),
        Llzz_(n,m,10),Llzp_(n,m,10),
        Nlzz_(n,m,10),Nlzp_(n,m,10),
        Nlpp_(n,m,10),Nlpz_(n,m,10),
        Tlzz_(n,m,10),Tlzp_(n,m,10),
        Tlpz_(n,m,10),Tlpp_(n,m,10),
        Llpz_(n,m,10),Llpp_(n,m,10),
        Alzz_(n,m,10),Alzp_(n,m,10),
        Alpz_(n,m,10),Alpp_(n,m,10),
        A_(ndim_, adim_), B_(ndim_, adim_),
        lu_(ndim_, ndim_ * 800),
        diag_(ndim_), scaling_(ndim_)
    {
        dx_ = (xmax_ - xmin_) / (n_ - 1);
        dy_ = (ymax_ - ymin_) / (m_ - 1);

        for (int i = 0; i < n_; i++)
            x_[i] = i * dx_ + xmin_;
        for (int i = 0; i < m_; i++)
            y_[i] = i * dy_ + ymin_;

        // if not periodic we only use interior points for many computations
        jmin_ = (periodic_) ? 0 : 1;
        jmax_ = (periodic_) ? m_ : m_-1;
        imin_ = (periodic_) ? 0 : 1;
        imax_ = (periodic_) ? n_ : n_-1;

        /******************************************************************************
         *     1:   alpha_tau
         *     2:   beta parameter
         *     3:   bottom friction r
         *     4:   aspect ratio (B/L)
         *     5:   Reynolds number
         *     6:   No (free) slip parameter EW  - 1.0 (0.0)
         *     7:   No (free) slip parameter NS  - 1.0 (0.0)
         *     10:  F
         ******************************************************************************/
        // par_(1) = taudim*Lxdim/(rhodim*hdim*udim**2)     // par 1&2 omgedraaid
        par_(1) = 1.0e+03;
        par_(2) = beta0dim_*Lxdim_*Lxdim_/udim_;
        // par_(3) = bfdim*Lxdim/udim;
        par_(3) = 0.0; // no bottom friction
        par_(4) = Lydim_/Lxdim_;
        par_(5) = udim_*Lxdim_/Ahdim_;
        // par_(10) = (f0dim*Lxdim)**2/(gdim*hdim);
        par_(10) = 0.0;
        par_(11) = 0.0;

        par_(6) = 1.0;
        par_(7) = 0.0;

        par_(19) = 0.0;

        compute_linear();
    }

    QG::~QG()
    {
    }

    void QG::compute_linear()
    {
        compute_forcing(); // definition windstress pattern (sin 2 pi y) on grid

        lin(); // linear part which is independent of the state
    }

    // fully 0-based (including XX) find row function
    int QG::findRow(int i, int j, int XX)
    {
        return nun_*(n_*j+i)+XX;
    }

    // shift: give coordinates of neighbour at location loc. Indexing:
    //        2 5     
    //        1 4 7
    //          3 6
    void QG::shift(int i, int j, int &i2, int &j2, int loc)
    {
        j2 = j + (loc%3)-1;
        i2 = i + (loc/3)-1;

        // resolve periodic boundaries
        if (periodic_)
        {
            if (i2 == n_) i2 = 0;
            if (i2 == -1) i2 = n_-1;
            if (j2 == m_) j2 = 0;
            if (j2 == -1) j2 = m_-1;
        }
    }

    void QG::assemble(std::vector<std::vector<Vector3D> > const &stencil, Matrix &A)
    {
        for (int i = 0; i < adim_; i++)
            A.co[i] = 0.0;
        
        int ZZ = 0; // state component zeta index
        int PP = 1; // state component psi index

        int np = 7; // number of possible neighbouring positions

        int i2,j2;
        double val;
        int col;

        int elm_ctr = 0;
        int beg_ctr = 0;
        
        for (int j = 0; j < m_; ++j)
            for (int i = 0; i < n_; ++i)
                for (int a = ZZ; a <= PP; ++a)
                {
                    A.beg[beg_ctr] = elm_ctr;
                    ++beg_ctr;

                    for (int loc = 1; loc <= np; ++loc)
                    {
                        shift(i, j, i2, j2, loc);
                        
                        // now I know the neighbour location i2,j2 and
                        // I can get the dependency for Z and P
                        for (int b = ZZ; b <= PP; ++b)
                        {
                            val = stencil[a][b](i, j, loc);
                            if (std::abs(val) > 1e-12)
                            {
                                A.co[elm_ctr] = val;
                                col = findRow(i2, j2, b);
                                A.jco[elm_ctr] = col;
                                ++elm_ctr;
                            }
                        }
                    }
                }
        A.beg[beg_ctr] = elm_ctr;        
    }
    
    void QG::assembleA()
    {        
        // Let's make a 5D array out of 4 3D arrays
        std::vector<std::vector<Vector3D> > Al = {
            {Alzz_, Alzp_},
            {Alpz_, Alpp_} };

        // Assemble the 5D stencil into the matrix A_;
        assemble(Al, A_);

    }

    void QG::assembleB()
    {
        // Let's make a 5D array out of 4 3D arrays
        std::vector<std::vector<Vector3D> > Bl = {
            {Tlzz_, Tlzp_},
            {Tlpz_, Tlpp_} };

        // Assemble the 5D stencil into the matrix B_;
        assemble(Bl, B_);        
    }

    void QG::jacob(double const *un, double sig)
    {
        Vector2D om(n_, m_);
        Vector2D ps(n_, m_);

        u_to_psi(un, om, ps);
        nlin_jac(om, ps);
        timedep();

        for (int i = 0; i < Alzz_.size(); i++)
        {
            Alzz_[i] = Llzz_[i] + Nlzz_[i] - sig * Tlzz_[i];
            Alzp_[i] = Llzp_[i] + Nlzp_[i] - sig * Tlzp_[i];
            Alpz_[i] = Llpz_[i];
            Alpp_[i] = Llpp_[i];
        }

        boundaries();
        assembleA();
        Asort();
    }

    void QG::jacobian(double const *un, double sig, int *beg, int *jco, double *co )
    {
        Vector2D om(n_, m_);
        Vector2D ps(n_, m_);

        u_to_psi(un, om, ps);
        nlin_jac(om, ps);
        timedep();

        for (int i = 0; i < Alzz_.size(); i++)
        {
            Alzz_[i] = Llzz_[i] + Nlzz_[i] - sig * Tlzz_[i];
            Alzp_[i] = Llzp_[i] + Nlzp_[i] - sig * Tlzp_[i];
            Alpz_[i] = Llpz_[i];
            Alpp_[i] = Llpp_[i];
        }

        boundaries();
        assembleA();
        Asort();

        for (int i = 0; i < ndim_+1; i++ )
            beg[i]=A_.beg[i];
        for (int i = 0; i < beg[ndim_]; i++ )
        {
            co[i]=A_.co[i] ;
            jco[i]=A_.jco[i] ;
        }
    }

  
    void QG::writeA(char const *name, double const *un, double sig)
    {
        jacob(un, sig);

        std::ofstream f(name);
        f << "%%MatrixMarket matrix coordinate real general\n";
        f << ndim_ << ' ' << ndim_ << ' ' << A_.beg[ndim_] << std::endl;
        for (int i = 0; i < ndim_; i++)
        {
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                f << i + 1 << " " << A_.jco[v] + 1 << " " <<  A_.co[v] << '\n';
        }
        f << std::flush;
        f.close();
    }

    void QG::writeM(char const *name)
    {
        timedep();
        assembleB();

        std::ofstream f(name);
        f << "%%MatrixMarket matrix coordinate real general\n";
        f << ndim_ << ' ' << ndim_ << ' ' << B_.beg[ndim_] << std::endl;
        for (int i = 0; i < ndim_; i++)
        {
            for (int v = B_.beg[i]; v < B_.beg[i+1]; v++)
                f << i + 1 << " " << B_.jco[v] + 1 << " " <<  B_.co[v] << '\n';
        }
        f << std::flush;
        f.close();
    }

    void QG::rhs(double const *un, double *b)
    {
        Vector2D om(n_, m_);
        Vector2D ps(n_, m_);
        Vector1D Frc(ndim_);

        double alpha_tau = par_(11) * par_(1);

        for (int i = 0; i < ndim_; i++)
            Frc[i] = 0.0;

        u_to_psi(un, om, ps);
        nlin_rhs(om, ps);

        for (int j = jmin_; j < jmax_; j++)
        {
            for (int i = imin_; i < imax_; i++)
            {
                int row = 2 * (n_ * j + i);
                Frc(row) = alpha_tau * tx_(i,j);
                Frc(row+1) = 0.0;
            }
        }

        for (int i = 0; i < Alzz_.size(); i++)
        {
            Alzz_[i] = Llzz_[i] + Nlzz_[i];
            Alzp_[i] = Llzp_[i] + Nlzp_[i];
            Alpz_[i] = Llpz_[i];
            Alpp_[i] = Llpp_[i];
        }

        boundaries();
        assembleA();
        Asort();

        // This is the reverse of the fortran version!!!
        // B = Au - Frc
        for (int i = 0; i < ndim_; i++)
        {
            b[i] = -Frc(i);
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                b[i] += A_.co[v] * un[A_.jco[v]];
        }
    }

    void QG::bilin(double const *un,double const *vn, double *b)
    {
        Vector2D om(n_, m_);
        Vector2D ps(n_, m_);
    
        u_to_psi(un, om, ps);
        nlin_rhs(om, ps);

        for (int i = 0; i < Alzz_.size(); i++)
        {
            Alzz_[i] = Nlzz_[i];
            Alzp_[i] = Nlzp_[i];
            Alpz_[i] = 0.0;
            Alpp_[i] = 0.0;
        }

        // boundariesNl();
        assembleA();
        Asort();

        // This is the reverse of the fortran version!!!
        // B = Au 
        for (int i = 0; i < ndim_; i++)
        {
            b[i] = 0;
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                b[i] -= A_.co[v] * vn[A_.jco[v]];
        }
    }

    void QG::apply(double const *x, double *y)
    {
        for (int i = 0; i < ndim_; i++)
        {
            y[i] = 0.0;
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                y[i] += A_.co[v] * x[A_.jco[v]];
        }
    }

    void QG::mass(double *M)
    {
        timedep();
        assembleB();

        for (int i = 0; i < ndim_; i++)
        {
            M[i] = 0.0;
            for (int j = B_.beg[i]; j < B_.beg[i+1]; j++)
            {
                if (B_.jco[j] != i)
                {
                    if (B_.co[j] != 0.0)
                    {
                        std::cerr << "Mass matrix has wrong format" << std::endl;
                        std::cerr << "i = " << i << ", j = " << B_.jco[j];
                        std::cerr << ", val = " << B_.co[j] << std::endl;
                        exit(1);
                    }
                }
                else
                    M[i] = B_.co[j];
            }
        }
    }

    double QG::windFun(double x, double y)
    {
        double asym = par_(19);
        double y2   = (y-ymin_) / (ymax_-ymin_);
        return - (1. - asym) * sin(2 * M_PI * y2)
            - asym * sin(M_PI * y2);
    }

    // wind forcing according to [Edeling, 2019]
    double QG::windFun2(double x, double y)
    {
        double asym = par_(19);
        double y2   = (y-ymin_) / (ymax_-ymin_);
        double x2   = (x-xmin_) / (xmax_-xmin_);
        return -sin(5*2*M_PI*x2+M_PI/2)*sin(5*2*M_PI*y2+M_PI/2);
    }

    void QG::lin()
    {
/*     Quasi-geostrophic equations, barotropic
 *     Produce local element matrices for linear operators
 */
        Vector3D z(n_,m_,10);
        Vector3D dxx(n_,m_,10);
        Vector3D dyy(n_,m_,10);
        Vector3D cor(n_,m_,10);

        double Re   = par_(5);
        double Beta = par_(2);
        double rbf  = par_(3);

        deriv(1, z);   // ALLEEN BIJ GEBRUIK BODEMFRICTIE NIET NUL
        deriv(2, dxx); // ALLEEN GETALLEN: VOOR BEIDE VGL TE GEBRUIKEN
        deriv(3, dyy); 
        deriv(4, cor);

        for (int i = 0; i < Llzz_.size(); i++)
            Llzz_[i] = -(dxx[i] + dyy[i]) / Re + rbf * z[i];
        for (int i = 0; i < Llzp_.size(); i++)
            Llzp_[i] = Beta * cor[i];

        deriv(1, Llpz_);

        for (int i = 0; i < Llpp_.size(); i++)
            Llpp_[i] = -(dxx[i] + dyy[i]);
    }

    void QG::nlin_rhs(Vector2D const &om, Vector2D const &ps)
    {
        // Produce local matrices for nonlinear operators for calc of Rhs
        Vector3D Udx(n_,m_,10), Vdy(n_,m_,10);

        double F = par_(10);

        nonlin(1,Udx,om,ps);
        nonlin(2,Vdy,om,ps);

        for (int i = 0; i < Nlzz_.size(); i++)
            Nlzz_[i] = Udx[i] + Vdy[i];
        for (int i = 0; i < Nlzp_.size(); i++)
            Nlzp_[i] = - F*(Udx[i] + Vdy[i]);
        // for (int i = 0; i < Nlpp_.size(); i++)
        //     Nlpp_[i] = 0.0;
        // for (int i = 0; i < Nlpz_.size(); i++)
        //     Nlpz_[i] = 0.0;
    }

    void QG::nlin_jac(Vector2D const &om, Vector2D const &ps)
    {
        // Produce local matrices for nonlinear operators for calc of Jacobian

        Vector3D Udx(n_,m_,10), Vdy(n_,m_,10);
        Vector3D udZx(n_,m_,10),vdZy(n_,m_,10);
        // Vector3D udPx(n_,m_,10),vdPy(n_,m_,10);

        double F = par_(10);

        nonlin(1,Udx, om,ps);
        nonlin(2,Vdy, om,ps);
        nonlin(3,udZx,om,ps);
        nonlin(4,vdZy,om,ps);
        // nonlin(5,udPx,om,ps);
        // nonlin(6,vdPy,om,ps);

        for (int i = 0; i < Nlzz_.size(); i++)
            Nlzz_[i] = Udx[i] + Vdy[i];
        for (int i = 0; i < Nlzp_.size(); i++)
            // Nlzp_[i] = udZx[i] + vdZy[i] - F*(Udx[i] + Vdy[i]) - F*(udPx[i] + vdPy[i]);
            Nlzp_[i] = udZx[i] + vdZy[i] - F*(Udx[i] + Vdy[i]);
    }

    void QG::timedep()
    {
        // Produce local matrices for time-dependent operators
        for (int i = 0; i < Tlzz_.size(); i++)
        {
            Tlzz_[i] = 0.0;
            Tlzp_[i] = 0.0;
            Tlpz_[i] = 0.0; // not used
            Tlpz_[i] = 0.0; // not used
        }

        // Residual formulation: Bdx/dt + F(x) = 0        
        double F = par_(10);
        for (int j = jmin_; j < jmax_; j++)
        {
            for (int i = imin_; i < imax_; i++)
            {
                Tlzz_(i,j,4) = -1.0;
                Tlzp_(i,j,4) = F;
            }
        }
    }

    void QG::u_to_psi(double const *un, Vector2D &om, Vector2D &ps)
    {
        for (int j = 0; j < m_; j++)
        {
            for (int i = 0; i < n_; i++)
            {
                int row = 2 * (n_ * j + i);
                om(i, j) = un[row];
                ps(i, j) = un[row+1];
            }
        }
    }

    void QG::boundaries()
    {
// insert conditions at the 'real' boundaries of the domain
/*
 *     2 5
 *     1 4 7
 *       3 6
 */
        if (periodic_) return;
        
        double oml2 =    par_(6) / 2.0;
        double oml3 = -3*par_(6) / (dx_*dx_); // for equidistant grids no differences
        double omr2 =    par_(6) / 2.0;
        double omr3 = -3*par_(6) / (dx_*dx_);
        double omb2 =    par_(7) / 2.0;
        double omb3 = -3*par_(7) / (dy_*dy_);
        double omt2 =    par_(7) / 2.0;
        double omt3 = -3*par_(7) / (dy_*dy_);

        for (int i = 0; i < n_; i++)
        {
            Alzz_(i,0,4) = 1.0; // SOUTH:i=(1,N), j=1
            Alzz_(i,0,5) = omb2;
            Alzp_(i,0,5) = omb3;
            Alpp_(i,0,4) = 1.0;

            Alzz_(i,m_-1,4) = 1.0; // NORTH:i=(1,N), j=M
            Alzz_(i,m_-1,3) = omt2;
            Alzp_(i,m_-1,3) = omt3;
            Alpp_(i,m_-1,4) = 1.0;
        }

        for (int j = 1; j < m_-1; j++)
        {
            Alzz_(0,j,4) = 1.0; // WEST: i=1      , j=(2,M-1)
            Alzz_(0,j,7) = oml2;
            Alzp_(0,j,7) = oml3;
            Alpp_(0,j,4) = 1.0;

            Alzz_(n_-1,j,4) = 1.0; //EAST: i=N      , j=(2,M-1)
            Alzz_(n_-1,j,1) = omr2;
            Alzp_(n_-1,j,1) = omr3;
            Alpp_(n_-1,j,4) = 1.0;
        }
    }
    
    void QG::compute_forcing()
    {
        for (int i = 0; i < tx_.size(); i++)
        {
            tx_[i] = 0.0;
            ty_[i] = 0.0;
        }

        for (int j = jmin_; j < jmax_; j++)
        {
            for (int i = imin_; i < imax_; i++)
            {
                if (periodic_)
                    tx_(i, j) = windFun2(x_(i),y_(j));
                else
                    tx_(i, j) = windFun(x_(i),y_(j));
                
                ty_(i, j) = 0.0;
            }
        }
    }

/*
 *  dx = coef(i,j,pos,unkn)
 *
 *  position: 4 central(i,j)
 *                    1 west   (i-1)           2  5
 *                    7 east   (i+1)           1  4  7
 *                    3 south  (j-1)              3  6
 *                    5 north  (j+1)
 */
    void QG::nonlin(int type, Vector3D &atom, Vector2D const &om, Vector2D const &ps)
    {
        /*
         *     nonlinear terms for the zeta and psi equation
         *     1:  Udx     (udx)
         *     2:  Vdy     (vdy)
         *     3:  udZx
         *     4:  vdZy
         *     5:  udPx
         *     6:  vdPy
         */

        for (int i = 0; i < atom.size(); i++)
            atom[i] = 0.0;

        double r4dxdy = 1.0/(4*dy_*dx_);

        switch (type)
        {
        case 1: // Udx
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,1) = +r4dxdy*(ps(i,j+1)-ps(i,j-1));
                    atom(i,j,7) = -r4dxdy*(ps(i,j+1)-ps(i,j-1));
                }
            }
        }
        break;
        case 2: // Vdy
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,3) = -r4dxdy*(ps(i+1,j)-ps(i-1,j));
                    atom(i,j,5) =  r4dxdy*(ps(i+1,j)-ps(i-1,j));
                }
            }
        }
        break;
        case 3: // udZx
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,3) = + r4dxdy*(om(i+1,j)-om(i-1,j));
                    atom(i,j,5) = - r4dxdy*(om(i+1,j)-om(i-1,j));
                }
            }
        }
        break;
        case 4: // vdZy
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,1) = - r4dxdy*(om(i,j+1)-om(i,j-1));
                    atom(i,j,7) = + r4dxdy*(om(i,j+1)-om(i,j-1));
                }
            }
        }
        break;
        case 5: // udPx
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,3) = + r4dxdy*(ps(i+1,j)-ps(i-1,j));
                    atom(i,j,5) = - r4dxdy*(ps(i+1,j)-ps(i-1,j));
                }
            }
        }
        break;
        case 6: // vdPy
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,1) = - r4dxdy*(ps(i,j+1)-ps(i,j-1));
                    atom(i,j,7) = + r4dxdy*(ps(i,j+1)-ps(i,j-1));
                }
            }
        }
        break;
        }
    }

    void QG::deriv(int type, Vector3D &atom)
    {
        /*
         *     1:  z, p
         *     2:  zxx, pxx
         *     3:  zyy, pyy
         *     4:  cor
         */

        for (int i = 0; i < atom.size(); i++)
            atom[i] = 0.0;

        switch (type)
        {
        case 1:
        {
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,4) = 1.0;
                }
            }
        }
        break;
        case 2:
        {
            double rdx2i = (1.0/dx_) * (1.0/dx_);
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,1) =    rdx2i;
                    atom(i,j,4) = -2*rdx2i;
                    atom(i,j,7) =    rdx2i;
                }
            }
        }
        break;
        case 3:
        {
            double rdy2i = (1.0/dy_) * (1.0/dy_);
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,3) =    rdy2i;
                    atom(i,j,4) = -2*rdy2i;
                    atom(i,j,5) =    rdy2i;
                }
            }
        }
        break;
        case 4:
        {
            double r2dx = 1.0 / (2*dx_);
            for (int j = jmin_; j < jmax_; j++)
            {
                for (int i = imin_; i < imax_; i++)
                {
                    atom(i,j,1) = -r2dx;
                    atom(i,j,7) =  r2dx;
                }
            }
        }
        break;
        }
    }

    void QG::readfort(int irs, double *u)
    {
        std::ifstream f("fort.4");
        while (f.good())
        {
            int lbl, contpar, evecs, ndim, nskip;
            std::string dummys;
            f >> lbl;
            f >> contpar;
            f >> evecs;
            f >> ndim;
            f >> nskip;

            if (ndim != ndim_)
            {
                std::cerr << "Wrong dimension " << ndim << " != " << ndim_ << std::endl;
                exit(1);
            }

            if (lbl != irs)
            {
                for (int i = 0; i < nskip + 1; i++)
                    std::getline(f, dummys);
                continue;
            }

            double dummy;
            for (int i = 0; i < evecs + 2; i++)
                std::getline(f, dummys);

            for (int i = 0; i < ndim; i++)
            {
                for (int j = 0; j < evecs + 2; j++)
                {
                    if (j == 0)
                        f >> u[i];
                    else
                        f >> dummy;
                }
            }

            for (int i = 0; i < 20; i++)
            {
                f >> par_[i+1];
            }
        }

        std::cout << "Read parameters: ";
        for (int i = 0; i < 20; i++)
            std::cout << par_[i+1] << " ";
        std::cout << std::endl;

        f.close();

        compute_linear();
    }

    double QG::get_par(int par)
    {
        return par_(par);
    }

    void QG::set_par(int par, double val)
    {
        par_(par) = val;

        compute_linear();
    }

    int QG::solve(double *x)
    {
        apply_scaling(x);

        int ierr = cgstab(x, 1e-6);

        remove_scaling();

        return ierr;
    }

    int QG::compute_precon()
    {
        double eps = 1e-4;

        compute_scaling();

        apply_scaling();
 
        int *volg = new int[ndim_];
        for (int i = 0; i < ndim_; i++)
            volg[i] = -1;

        double *kop = new double[ndim_];
        for (int i = 0; i < ndim_; i++)
            kop[i] = 0.0;

        A_.jco[A_.beg[ndim_]] = -1;
        lu_.beg[0] = 0;
        int tellu = -1;

        for (int i = 0; i < ndim_; i++)
        {
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
            {
                volg[A_.jco[v]] = A_.jco[v+1];
                kop[A_.jco[v]] = A_.co[v];
            }
            volg[A_.jco[A_.beg[i+1]-1]] = ndim_;
            int j = A_.jco[A_.beg[i]];
            double ff = 0.0;

            while (j < i)
            {
                if (std::abs(kop[j]) > eps * diag_[i])
                {
                    tellu++;
                    lu_.co[tellu] = kop[j] * lu_.co[lu_.di[j]];
                    lu_.jco[tellu] = j;
                    for (int v = lu_.di[j]+1; v < lu_.beg[j+1]; v++)
                    {
                        int k = lu_.jco[v];
                        if (volg[k] > -1)
                            kop[k] -= lu_.co[tellu] * lu_.co[v];
                        else
                        {
                            int zoek = lu_.jco[v-1];
                            while (volg[zoek] <= k)
                                zoek = volg[zoek];
                            volg[k] = volg[zoek];
                            volg[zoek] = k;
                            kop[k] = -lu_.co[tellu] * lu_.co[v];
                        }
                    }
                }
                else
                    ff += kop[j];
                int jn = volg[j];
                volg[j] = -1;
                j = jn;
            }

            tellu++;
            lu_.jco[tellu] = i;
            lu_.co[tellu] = kop[i];
            lu_.di[i] = tellu;
            j = volg[i];
            volg[i] = -1;

            while (j < ndim_)
            {
                if (std::abs(kop[j]) > eps * diag_[i])
                {
                    tellu++;
                    lu_.co[tellu] = kop[j];
                    lu_.jco[tellu] = j;
                }
                else
                    ff += kop[j];
                int jn = volg[j];
                volg[j] = -1;
                j = jn;
            }

            ff += kop[i];
            if (std::abs(ff) < 1e-8)
                ff = 1.0;
            ff = 1.0 / ff;
            lu_.co[lu_.di[i]] = ff;
            lu_.beg[i+1] = tellu+1;
        }

        delete[] kop;
        delete[] volg;

        remove_scaling();

        return 0;
    }

    int QG::lusolve(double *x)
    {
        for (int i = 1; i < ndim_; i++)
        {
            double tmp = 0.0;
            for (int j = lu_.beg[i]; j < lu_.di[i]; j++)
                tmp += lu_.co[j]*x[lu_.jco[j]];
            x[i] -= tmp;
        }
        x[ndim_-1] = x[ndim_-1] * lu_.co[lu_.di[ndim_-1]];
        for (int i = ndim_-2; i > -1; i--)
        {
            double tmp = 0.0;
            for (int j = lu_.di[i]+1; j < lu_.beg[i+1]; j++)
                tmp += lu_.co[j]*x[lu_.jco[j]];
            x[i] = (x[i]-tmp) * lu_.co[lu_.di[i]];
        }
        return 0;
    }

    double QG::dot(double *x, double *y)
    {
        double out = 0.0;
        for (int i = 0; i < ndim_; i++)
            out += x[i] * y[i];
        return out;
    }

    int QG::cgstab(double *b, double tol)
    {
        if (dot(b, b) < 1e-20)
            return 1;
 
        double *p = new double[ndim_];
        double *q = new double[ndim_];
        double *r = new double[ndim_];
        double *s = new double[ndim_];
        double *t = new double[ndim_];
        double *v = new double[ndim_];

        bool nyet = true;
        double maxb = 0;
        for (int i = 0; i < ndim_; i++)
            maxb = std::max(maxb, std::abs(b[i]));
        tol *= maxb;

        lusolve(b);
        memcpy(r, b, ndim_ * sizeof(double));
        memset(v, 0.0, ndim_ * sizeof(double));
        memset(q, 0.0, ndim_ * sizeof(double));
        memcpy(p, r, ndim_ * sizeof(double));
        memset(b, 0.0, ndim_ * sizeof(double));

        int ret = -1;
        double wdak = 1.0;
        double bet = 1.0;
        double alf = 1.0;

        for (int tel = 0; tel < 100; tel++)
        {
            double betdak = dot(r, p);
            double w = (betdak / bet) * (wdak / alf);
            bet = betdak;
            for (int i = 0; i < ndim_; i++)
                q[i] = r[i] + w*(q[i]-alf*v[i]);

            apply(q, v);
            lusolve(v);

            wdak = dot(p, v);
            wdak = bet / wdak;
            for (int i = 0; i < ndim_; i++)
                s[i] = r[i] - wdak * v[i];

            apply(s, t);
            lusolve(t);

            double ts = dot(t, s);
            double tt = dot(t, t);

            if (tt < 1e-30)
            {
                alf = 0.0;
                nyet = false;
            }
            else
                alf = ts / tt;

            double rmax = 0.0;
            for (int i = 0; i < ndim_; i++)
            {
                b[i] += wdak * q[i] + alf * s[i];
                r[i] = s[i] - alf * t[i];
                rmax = std::max(rmax, std::abs(r[i]));
            }
            if (rmax < tol || !nyet)
            {
                if (rmax > tol * 1e3)
                {
                    ret = -2;
                    break;
                }
                ret = 0;
                break;
            }
        }

        delete[] p;
        delete[] q;
        delete[] r;
        delete[] s;
        delete[] t;
        delete[] v;

        return ret;
    }

    void QG::compute_scaling()
    {
        // Row-wise scaling of A
        for (int i = 0; i < ndim_; i++)
        {
            diag_[i] = 0.0;
            scaling_[i] = 0.0;
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
            {
                if (A_.jco[v] == i)
                    diag_[i] = A_.co[v];
                scaling_[i] = std::max(scaling_[i], std::abs(A_.co[v]));
            }
            if (diag_[i] < 0.0)
                scaling_[i] = -scaling_[i];
            diag_[i] /= scaling_[i];
        }
    }

    void QG::apply_scaling()
    {
        for (int i = 0; i < ndim_; i++)
        {
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                A_.co[v] /= scaling_[i];
        }
    }

    void QG::apply_scaling(double *x)
    {
        for (int i = 0; i < ndim_; i++)
        {
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                A_.co[v] /= scaling_[i];
            x[i] /= scaling_[i];
        }
    }

    void QG::remove_scaling()
    {
        for (int i = 0; i < ndim_; i++)
        {
            for (int v = A_.beg[i]; v < A_.beg[i+1]; v++)
                A_.co[v] *= scaling_[i];
        }
    }

    void QG::Asort()
    {
        for (int i = 0; i < ndim_; i++)
        {
            for (int v = A_.beg[i]; v < A_.beg[i+1]-1; v++)
            {
                int k = v;
                for (int vv = k+1; vv < A_.beg[i+1]; vv++)
                    if (A_.jco[vv] < A_.jco[k])
                        k = vv;
                double dum = A_.co[k];
                int jdum = A_.jco[k];
                A_.co[k] = A_.co[v];
                A_.jco[k] = A_.jco[v];
                A_.jco[v] = jdum;
                A_.co[v] = dum;
            }
        }
    }
}
