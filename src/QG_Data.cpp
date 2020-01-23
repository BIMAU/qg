#include "QG_Data.hpp"

#include <cstring>
#include <cmath>

namespace QG
{

    Matrix::Matrix(int ndim, int adim)
        :
        ndim_(ndim),
        adim_(adim)
    {
        beg = new int[ndim+1];
        jco = new int[ndim*800];
        co  = new double[ndim*800];
        di  = new int[ndim];
    }

    Matrix::Matrix(Matrix const &other)
        :
        Matrix(other.ndim_, other.adim_)
    {
        memcpy(beg, other.beg, (ndim_ + 1) * sizeof(int));
        memcpy(jco, other.jco, adim_ * sizeof(int));
        memcpy(co, other.co, adim_ * sizeof(double));
        memcpy(di, other.di, ndim_ * sizeof(int));
    }

    Matrix::~Matrix()
    {
        delete[] beg;
        delete[] jco;
        delete[] co;
        delete[] di;
    }

    void Matrix::pack()
    {
        // Remove entries smaller than 1e-12
        int vv = 0;
        for (int i = 0; i < ndim_; i++)
        {
            int begin = vv;
            for (int v = beg[i]; v < beg[i+1]; v++)
            {
                if (std::abs(co[v]) > 1e-12)
                {
                    co[vv] = co[v];
                    jco[vv] = jco[v];
                    vv++;
                }
            }
            beg[i] = begin;
        }
        beg[ndim_] = vv;
    }

    Vector1D::Vector1D(int n)
        :
        n_(n)
    {
        data_ = std::vector<double>(n, 0.0);
    }

    Vector1D::Vector1D(Vector1D const &other)
        :
        Vector1D(other.n_)
    {
        data_ = other.data_;
    }

    Vector1D::~Vector1D()
    {}


    double &Vector1D::operator[](int i)
    {
        return data_[i];
    }

    double const &Vector1D::operator[](int i) const
    {
        return data_[i];
    }

    double &Vector1D::operator()(int i)
    {
        return data_[i];
    }

    double const &Vector1D::operator()(int i) const
    {
        return data_[i];
    }

    int Vector1D::size()
    {
        return n_;
    }

    Vector2D::Vector2D(int m, int n)
        :
        m_(m),
        n_(n)
    {
        data_ = std::vector<double>(m * n, 0.0);
    }

    Vector2D::Vector2D(Vector2D const &other)
        :
        Vector2D(other.m_, other.n_)
    {
        data_ = other.data_;
    }

    Vector2D::~Vector2D()
    {}
    
    double &Vector2D::operator[](int i)
    {
        return data_[i];
    }

    double const &Vector2D::operator[](int i) const
    {
        return data_[i];
    }

    double &Vector2D::operator()(int i, int j)
    {
        i = (n_ + i) % n_;
        j = (m_ + j) % m_;
        return data_[i + j*m_];
    }

    double const &Vector2D::operator()(int i, int j) const
    {
        i = (n_ + i) % n_;
        j = (m_ + j) % m_;
        return data_[i + j*m_];
    }

    int Vector2D::size()
    {
        return m_ * n_;
    }

    Vector3D::Vector3D(int l, int m, int n)
        :
        l_(l),
        m_(m),
        n_(n)
    {
        data_ = std::vector<double>(l * m * n, 0.0);
    }

    Vector3D::Vector3D(Vector3D const &other)
        :
        Vector3D(other.l_, other.m_, other.n_)
    {
        data_ = other.data_;
    }
    
    Vector3D::~Vector3D()
    {}

    double &Vector3D::operator[](int i)
    {
        return data_[i];
    }

    double const &Vector3D::operator[](int i) const
    {
        return data_[i];
    }

    double &Vector3D::operator()(int i, int j, int k)
    {
        return data_[i + j*l_ + k*l_*m_];
    }

    double const &Vector3D::operator()(int i, int j, int k) const
    {
        return data_[i + j*l_ + k*l_*m_];
    }

    int Vector3D::size()
    {
        return l_ * m_ * n_;
    }

}
