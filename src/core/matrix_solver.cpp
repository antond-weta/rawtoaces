// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#ifdef ENABLE_EIGEN
#    include <Eigen/Core>
#    include <Eigen/Dense>
#endif // ENABLE_EIGEN

#ifdef ENABLE_CERES
#    include <ceres/ceres.h>
#endif // ENABLE_CERES

#include "matrix_solver.hpp"

namespace rta
{
namespace core
{

bool use_eigen = true;
bool use_ceres = true;

inline void check_dependencies()
{
#ifndef ENABLE_EIGEN
    if ( use_eigen )
    {
        std::cerr << "The library was built without Eigen support. "
                  << "The build-in solver will be used instead. " << std::endl;
        use_eigen = false;
    }
#endif // ENABLE_EIGEN

#ifndef ENABLE_CERES
    if ( use_ceres )
    {
        std::cerr << "The library was built without Ceres-solver support. "
                  << "The build-in solver will be used instead." << std::endl;
        use_ceres = false;
    }
#endif // ENABLE_EIGEN
}

#ifdef ENABLE_EIGEN

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> vec_from_std( const std::vector<T> &vec )
{
    size_t n = vec.size();
    assert( n > 0 );

    Eigen::Matrix<T, Eigen::Dynamic, 1> result( vec.size() );
    for ( size_t i = 0; i < n; i++ )
        result( i ) = vec[i];
    return result;
}

template <typename T>
std::vector<T> vec_to_std( const Eigen::Matrix<T, Eigen::Dynamic, 1> &vec )
{
    size_t n = vec.rows();
    assert( n > 0 );

    std::vector<T> result( n );

    for ( size_t i = 0; i < n; i++ )
        result[i] = vec( i );
    return result;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
mat_from_std( const std::vector<std::vector<T>> &mat )
{
    size_t n = mat.size();
    assert( n > 0 );

    size_t m = mat[0].size();
    assert( m > 0 );

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result( n, m );

    for ( size_t i = 0; i < n; i++ )
        for ( size_t j = 0; j < m; j++ )
            result( i, j ) = mat[i][j];
    return result;
}

template <typename T>
std::vector<std::vector<T>>
mat_to_std( const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &mat )
{
    size_t n = mat.rows();
    assert( n > 0 );

    size_t m = mat.cols();
    assert( m > 0 );

    std::vector<std::vector<T>> result( n, std::vector<T>( m ) );

    for ( size_t i = 0; i < n; i++ )
        for ( size_t j = 0; j < m; j++ )
            result[i][j] = mat( i, j );
    return result;
}

#endif // ENABLE_EIGEN

template <typename T_dst, typename T_src>
std::vector<T_dst> transform_vec( const std::vector<T_src> &vec )
{
    std::vector<T_dst> result( vec.size() );
    std::transform(
        vec.begin(), vec.end(), result.begin(), []( const T_src &v ) {
            return T_dst( v );
        } );
    return result;
}

template <typename T_dst, typename T_src>
std::vector<std::vector<T_dst>>
transform_mat( const std::vector<std::vector<T_src>> &mat )
{
    std::vector<std::vector<T_dst>> result( mat.size() );
    std::transform(
        mat.begin(),
        mat.end(),
        result.begin(),
        []( const std::vector<T_src> &v ) {
            return transform_vec<T_dst>( v );
        } );
    return result;
}

template <> std::vector<double> transform_vec( const std::vector<double> &vec )
{
    return vec;
}

template <>
std::vector<std::vector<double>>
transform_mat( const std::vector<std::vector<double>> &mat )
{
    return mat;
}

template <typename T>
std::vector<std::vector<T>> transposed( std::vector<std::vector<T>> mat )
{
    size_t n = mat.size();
    assert( n > 0 );

    size_t m = mat[0].size();
    assert( m > 0 );

    if ( use_eigen )
    {
#ifdef ENABLE_EIGEN
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat2 =
            mat_from_std( mat );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat3 =
            mat2.transpose();
        return mat_to_std( mat3 );
#else  // ENABLE_EIGEN
        std::cerr << "The library was built without Eigen support."
                  << std::endl;
        exit( -1 );
#endif // ENABLE_EIGEN
    }
    else
    {
        std::vector<std::vector<T>> result( m, std::vector<T>( n ) );

        for ( size_t i = 0; i < n; i++ )
            for ( size_t j = 0; j < m; j++ )
                result[j][i] = mat[i][j];

        return result;
    }
}

template <typename T>
std::vector<std::vector<T>>
product( std::vector<std::vector<T>> mat1, std::vector<std::vector<T>> mat2 )
{
    size_t n1 = mat1.size();
    assert( n1 > 0 );

    size_t m1 = mat1[0].size();
    assert( m1 > 0 );

    size_t n2 = mat2.size();
    assert( n2 > 0 );

    size_t m2 = mat2[0].size();
    assert( m2 > 0 );

    assert( m1 == n2 );

    std::vector<std::vector<T>> result( n1, std::vector<T>( m2 ) );

    if ( use_eigen )
    {
#ifdef ENABLE_EIGEN
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat3 =
            mat_from_std( mat1 );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat4 =
            mat_from_std( mat2 );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat5 = mat3 * mat4;
        return mat_to_std( mat5 );
#else  // ENABLE_EIGEN
        std::cerr << "The library was built without Eigen support."
                  << std::endl;
        exit( -1 );
#endif // ENABLE_EIGEN
    }
    else
    {
        std::vector<std::vector<T>> result( n1, std::vector<T>( m2 ) );

        for ( size_t i = 0; i < n1; i++ )
            for ( size_t j = 0; j < m2; j++ )
                for ( size_t k = 0; k < m1; k++ )
                    result[i][j] += mat1[i][k] * mat2[k][j];

        return result;
    }
}

template <typename T>
std::vector<T> product( std::vector<std::vector<T>> mat, std::vector<T> vec )
{
    size_t n1 = mat.size();
    assert( n1 > 0 );

    size_t m1 = mat[0].size();
    assert( m1 > 0 );

    size_t n2 = vec.size();
    assert( n2 > 0 );

    assert( m1 == n2 );

    if ( use_eigen )
    {
#ifdef ENABLE_EIGEN
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat2 =
            mat_from_std( mat );
        Eigen::Matrix<T, Eigen::Dynamic, 1> vec2 = vec_from_std( vec );
        Eigen::Matrix<T, Eigen::Dynamic, 1> vec3 = mat2 * vec2;
        return vec_to_std( vec3 );
#else  // ENABLE_EIGEN
        std::cerr << "The library was built without Eigen support."
                  << std::endl;
        exit( -1 );
#endif // ENABLE_EIGEN
    }
    else
    {
        std::vector<T> result( n1 );

        for ( size_t i = 0; i < n1; i++ )
            for ( size_t j = 0; j < m1; j++ )
                result[i] += mat[i][j] * vec[j];

        return result;
    }
}

template <typename T>
bool inverse(
    const std::vector<std::vector<T>> &in, std::vector<std::vector<T>> &out )
{
    if ( use_eigen )
    {
#ifdef ENABLE_EIGEN
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat2 =
            mat_from_std( in );
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat3 = mat2.inverse();
        out = mat_to_std( mat3 );
        return true;
#else  // ENABLE_EIGEN
        std::cerr << "The library was built without Eigen support."
                  << std::endl;
        exit( -1 );
#endif // ENABLE_EIGEN
    }
    else
    {
        assert( in.size() == 3 );
        assert( in[0].size() == 3 );

        T a = in[0][0];
        T b = in[0][1];
        T c = in[0][2];
        T d = in[1][0];
        T e = in[1][1];
        T f = in[1][2];
        T g = in[2][0];
        T h = in[2][1];
        T i = in[2][2];

        T A = e * i - f * h;
        T B = f * g - d * i;
        T C = d * h - e * g;

        T det = A * a + B * b + C * c;
        if ( det == 0 )
            return false;

        T scale = 1.0 / det;

        T D = c * h - b * i;
        T E = a * i - c * g;
        T F = b * g - a * h;
        T G = b * f - c * e;
        T H = c * d - a * f;
        T I = a * e - b * d;

        out.resize( 3 );

        out[0].resize( 3 );
        out[0][0] = A * scale;
        out[0][1] = D * scale;
        out[0][2] = G * scale;

        out[1].resize( 3 );
        out[1][0] = B * scale;
        out[1][1] = E * scale;
        out[1][2] = H * scale;

        out[2].resize( 3 );
        out[2][0] = C * scale;
        out[2][1] = F * scale;
        out[2][2] = I * scale;

        return true;
    }
}

bool solve_linear( std::vector<std::vector<double>> &a, std::vector<double> &b )
{
    size_t n = a.size();
    assert( n > 0 );
    assert( a[0].size() == n );
    assert( b.size() == n );

    for ( size_t i = 0; i < n; i++ )
    {
        double max_val = a[i][i];
        size_t pivot   = i;

        for ( size_t j = i; j < n; j++ )
        {
            double v = std::abs( a[i][j] );
            if ( v > max_val )
            {
                max_val = v;
                pivot   = j;
            }
        }

        if ( max_val == 0 )
            return false;

        if ( pivot != i )
        {
            std::swap( a[i], a[pivot] );
            std::swap( b[i], b[pivot] );
        }

        double scale = 1.0 / a[i][i];
        for ( size_t k = i + 1; k < n; k++ )
        {
            a[i][k] *= scale;
        }
        a[i][i] = 1.0;
        b[i] *= scale;

        for ( size_t j = i + 1; j < n; j++ )
        {
            if ( a[j][i] != 0 )
            {
                double scale = -a[j][i];
                for ( size_t k = i; k < n; k++ )
                {
                    a[j][k] += a[i][k] * scale;
                }
                b[j] += b[i] * scale;
            }
        }
    }

    for ( size_t i = n - 1; i > 0; i-- )
    {
        for ( size_t j = 0; j < i; j++ )
        {
            double scale = -a[j][i];
            //            a[j][i] += a[i][i] * scale;
            b[j] += b[i] * scale;
        }
    }

    return true;
}

constexpr double v_6_29_3 = 6.0 * 6.0 * 6.0 / 29.0 / 29.0 / 29.0;
constexpr double v_29_6_3 = 29.0 * 29.0 / 6.0 / 6.0;

template <typename T> T lab_F( T val )
{
    if ( val > v_6_29_3 )
    {
#ifdef ENABLE_CERES
        if ( use_ceres )
        {
            return ceres::pow( val, 1.0 / 3.0 );
        }
        else
#endif //ENABLE_CERES
        {
            return pow( val, 1.0 / 3.0 );
        }
    }
    else
        return val * v_29_6_3 / 3.0;
}

template <typename T> T lab_F_( T val, T val_prime )
{
    if ( val > v_6_29_3 )
#ifdef ENABLE_CERES
        if ( use_ceres )
        {
            return T( 1.0 / 3.0 ) * ceres::pow( val, -2.0 / 3.0 ) * val_prime;
        }
        else
#endif //ENABLE_CERES
        {
            return T( 1.0 / 3.0 ) * pow( val, -2.0 / 3.0 ) * val_prime;
        }
    else
        return val_prime * v_29_6_3 / 3.0;
}

template <typename T> void xyz_to_lab( const T ( &in )[3], T ( &out )[3] )
{
    out[0] = T( 116.0 ) * lab_F( in[1] ) - T( 16.0 );
    out[1] = T( 500.0 ) * ( lab_F( in[0] ) - lab_F( in[1] ) );
    out[2] = T( 200.0 ) * ( lab_F( in[1] ) - lab_F( in[2] ) );
}

template <typename T>
void xyz_to_lab_( const T ( &in )[3], const T ( &in_ )[3], T ( &out )[3] )
{
    out[0] = T( 116.0 ) * lab_F_( in[1], in_[1] );
    out[1] = T( 500.0 ) * ( lab_F_( in[0], in_[0] ) - lab_F_( in[1], in_[1] ) );
    out[2] = T( 200.0 ) * ( lab_F_( in[1], in_[1] ) - lab_F_( in[2], in_[2] ) );
}

bool solve_nonlinear(
    double ( &x )[],
    double ( *objective )(
        const double ( &x )[],
        std::vector<std::vector<double>> &J,
        std::vector<double>              &F,
        void                             *context ),
    void  *context,
    double error_threshold,
    size_t max_steps,
    size_t verbosity )
{
    size_t step     = 0;
    double old_cost = std::numeric_limits<double>::max();

    std::vector<std::vector<double>> J;
    std::vector<double>              F;

    while ( true )
    {
        if ( step == max_steps )
        {
            if ( verbosity > 0 )
                std::cerr << "Reached the maximum iterations number."
                          << std::endl;
            break;
        }

        double new_cost = objective( x, J, F, context );

        auto J_t = transposed( J );
        auto A   = product( J_t, J );
        auto b   = product( J_t, F );

        if ( use_eigen )
        {
#ifdef ENABLE_EIGEN
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> AA =
                mat_from_std( A );
            Eigen::Matrix<double, Eigen::Dynamic, 1> bb = vec_from_std( b );
            Eigen::Matrix<double, Eigen::Dynamic, 1> x =
                AA.colPivHouseholderQr().solve( bb );
            b = vec_to_std( x );
#endif // ENABLE_EIGEN
        }
        else
        {
            if ( !solve_linear( A, b ) )
            {
                return false;
            }
        }

        for ( size_t i = 0; i < 6; i++ )
        {
            x[i] -= b[i];
        }

        if ( verbosity > 1 )
        {
            std::cerr << "Step " << std::setw( 3 ) << step << std::setw( 0 )
                      << ": cost = " << std::scientific << new_cost / 2.0
                      << std::fixed << std::endl;
        }

        if ( old_cost - new_cost < old_cost * error_threshold )
        {
            if ( verbosity > 0 )
                std::cerr
                    << "Reached the objective function tolerance threshold."
                    << std::endl;
            break;
        }
        old_cost = new_cost;

        step++;
    }

    return true;
}

template <typename Objective>
bool solve_nonlinear(
    double ( &x )[6],
    const Objective *objective,
    double           error_threshold,
    size_t           max_steps,
    size_t           verbosity )
{
    size_t n = objective->camera_points.size();
    assert( n > 0 );

    size_t step     = 0;
    double old_cost = std::numeric_limits<double>::max();

    std::vector<std::vector<double>> J( n * 3, std::vector<double>( 6, 0 ) );
    std::vector<double>              F( n * 3, 0 );

    while ( true )
    {
        if ( step == max_steps )
        {
            if ( verbosity > 0 )
                std::cerr << "Reached the maximum iterations number."
                          << std::endl;
            break;
        }

        double *B   = x;
        bool    res = ( *objective )( B, F.data(), &J );

        double new_cost = 0; //objective(x, J, F, context);
        for ( auto &i: F )
        {
            new_cost += i * i;
        }

        auto J_t = transposed( J );
        auto A   = product( J_t, J );
        auto b   = product( J_t, F );

        if ( use_eigen )
        {
#ifdef ENABLE_EIGEN
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> AA =
                mat_from_std( A );
            Eigen::Matrix<double, Eigen::Dynamic, 1> bb = vec_from_std( b );
            Eigen::Matrix<double, Eigen::Dynamic, 1> x =
                AA.colPivHouseholderQr().solve( bb );
            b = vec_to_std( x );
#endif // ENABLE_EIGEN
        }
        else
        {
            if ( !solve_linear( A, b ) )
            {
                return false;
            }
        }

        for ( size_t i = 0; i < 6; i++ )
        {
            x[i] -= b[i];
        }

        if ( verbosity > 1 )
        {
            std::cerr << "Step " << std::setw( 3 ) << step << std::setw( 0 )
                      << ": cost = " << std::scientific << new_cost / 2.0
                      << std::fixed << std::endl;
        }

        if ( old_cost - new_cost < old_cost * error_threshold )
        {
            if ( verbosity > 0 )
                std::cerr
                    << "Reached the objective function tolerance threshold."
                    << std::endl;
            break;
        }
        old_cost = new_cost;

        step++;
    }

    return true;
}

struct IDTSolverContext
{
    const std::vector<std::vector<double>> &camera_points;
    const std::vector<std::vector<double>> &xyz_points;
    const std::vector<double>              &white_rgb;
    const std::vector<std::vector<double>> &out_to_xyz_mat;

    IDTSolverContext(
        const std::vector<std::vector<double>> &rgb,
        const std::vector<std::vector<double>> &xyz,
        const std::vector<double>              &white_rgb_out,
        const std::vector<std::vector<double>> &out_to_xyz_matrix )
        : camera_points( rgb )
        , xyz_points( xyz )
        , white_rgb( white_rgb_out )
        , out_to_xyz_mat( out_to_xyz_matrix )
    {}
};

class Objective_LAB : public IDTSolverContext
{
public:
    using IDTSolverContext::IDTSolverContext;

    // Calculates the residuals (F) and their Jacobian (J). The residuals are the
    // colour distance in the CIELAB space between the point (L,A,B) calculated
    // via the matrix being solved, and the true value (Lt,At,Bt) calculated using
    // the standard observer curves.
    //
    // F:                       J:
    //      L(p0) - Lt(p0)          dL(p0)/dx0, dL(p0)/dx1, ... dL(p0)/dx5
    //      A(p0) - At(p0)          dA(p0)/dx0, dA(p0)/dx1, ... dA(p0)/dx5
    //      B(p0) - Bt(p0)          dB(p0)/dx0, dB(p0)/dx1, ... dB(p0)/dx5
    //      L(p1) - Lt(p1)          dL(p1)/dx0, dL(p1)/dx1, ... dL(p1)/dx5
    //      A(p1) - At(p1)          dA(p1)/dx0, dA(p1)/dx1, ... dA(p1)/dx5
    //      B(p1) - Bt(p1)          dB(p1)/dx0, dB(p1)/dx1, ... dB(p1)/dx5
    //      ...                     ...
    // Returns the cost value, which as the sum of squares of the vector F.
    template <typename T>
    bool operator()(
        const T                     *B,
        T                           *residuals,
        std::vector<std::vector<T>> *jacobian = nullptr ) const
    {
        size_t n = camera_points.size();
        assert( n > 0 );

        size_t m = camera_points[0].size();
        assert( m > 0 );

        std::vector<std::vector<T>> rgb =
            transform_mat<T, double>( camera_points );
        std::vector<std::vector<T>> xyz =
            transform_mat<T, double>( xyz_points );

        std::vector<std::vector<T>> mat_idt = {
            { B[0], B[1], T( white_rgb[0] ) - B[0] - B[1] },
            { B[2], B[3], T( white_rgb[1] ) - B[2] - B[3] },
            { B[4], B[5], T( white_rgb[2] ) - B[4] - B[5] },
        };

        std::vector<std::vector<T>> mat_idt_transposed = transposed( mat_idt );

        std::vector<std::vector<T>> mat_out_to_xyz_transposed(
            3, std::vector<T>( 3 ) );
        for ( size_t i = 0; i < 3; i++ )
        {
            for ( size_t j = 0; j < 3; j++ )
            {
                mat_out_to_xyz_transposed[i][j] = T( out_to_xyz_mat[j][i] );
            }
        }

        std::vector<std::vector<T>> rgb_W = {
            { T( white_rgb[0] ), T( white_rgb[1] ), T( white_rgb[2] ) }
        };

        std::vector<std::vector<T>> out_rgb =
            product( rgb, mat_idt_transposed );
        std::vector<std::vector<T>> out_xyz =
            product( out_rgb, mat_out_to_xyz_transposed );
        std::vector<std::vector<T>> xyz_W =
            product( rgb_W, mat_out_to_xyz_transposed );

        for ( size_t i = 0; i < n; i++ )
        {
            T xyz1[3] = { out_xyz[i][0], out_xyz[i][1], out_xyz[i][2] };
            T xyz2[3] = { xyz[i][0], xyz[i][1], xyz[i][2] };

            for ( size_t j = 0; j < 3; j++ )
            {
                xyz1[j] /= xyz_W[0][j];
                xyz2[j] /= xyz_W[0][j];
            }

            T lab1[3];
            T lab2[3];

            xyz_to_lab( xyz1, lab1 );
            xyz_to_lab( xyz2, lab2 );

            residuals[i * 3 + 0] = lab1[0] - lab2[0];
            residuals[i * 3 + 1] = lab1[1] - lab2[1];
            residuals[i * 3 + 2] = lab1[2] - lab2[2];

            if ( jacobian != nullptr )
            {
                // Partial derivatives of the out_rgb with regards to the IDT
                // matrix variables (x)

                // Ro = x0*(Rs-Bs) + x1*(Gs-Bs) + Rw*Bs
                // Go = x2*(Rs-Bs) + x3*(Gs-Bs) + Gw*Bs
                // Bo = x4*(Rs-Bs) + x5*(Gs-Bs) + Bw*Bs

                // Ro_ = (Rs-Bs) when j == 0, 0 otherwise
                // Ro_ = (Gs-Bs) when j == 1, 0 otherwise
                // Go_ = (Rs-Bs) when j == 2, 0 otherwise
                // Go_ = (Gs-Bs) when j == 3, 0 otherwise
                // Bo_ = (Rs-Bs) when j == 4, 0 otherwise
                // Bo_ = (Gs-Bs) when j == 5, 0 otherwise

                T &Rs = rgb[i][0];
                T &Gs = rgb[i][1];
                T &Bs = rgb[i][2];

                std::vector<std::vector<T>> out_rgb_ = {
                    { Rs - Bs, T( 0 ), T( 0 ) }, { Gs - Bs, T( 0 ), T( 0 ) },
                    { T( 0 ), Rs - Bs, T( 0 ) }, { T( 0 ), Gs - Bs, T( 0 ) },
                    { T( 0 ), T( 0 ), Rs - Bs }, { T( 0 ), T( 0 ), Gs - Bs },
                };

                // Partial derivatives of the out_xyz with regards to the IDT
                // matrix variables (x), which is just the out_to_xyz matrix
                // applied to out_rgb_.
                std::vector<std::vector<T>> out_xyz_ =
                    product( out_rgb_, mat_out_to_xyz_transposed );

                for ( size_t j = 0; j < 6; j++ )
                {
                    T xyz_[3] = { out_xyz_[j][0],
                                  out_xyz_[j][1],
                                  out_xyz_[j][2] };

                    for ( size_t c = 0; c < 3; c++ )
                    {
                        xyz_[j] /= xyz_W[0][c];
                    }

                    T lab_[3];
                    xyz_to_lab_( xyz1, xyz_, lab_ );

                    for ( size_t c = 0; c < 3; c++ )
                        ( *jacobian )[i * 3 + c][j] = T( lab_[c] );
                }
            }
        }

        return true;
    }
};

#ifdef ENABLE_CERES

template <typename Objective>
bool solve_nonlinear_ceres(
    double ( &x )[6],
    Objective *objective,
    double     error_threshold,
    size_t     max_steps,
    size_t     verbosity )
{
    const auto &points          = objective->camera_points;
    int         residuals_count = (int)( points.size() * points[0].size() );
    ceres::CostFunction *cost_function =
        new ceres::AutoDiffCostFunction<Objective, ceres::DYNAMIC, 6>(
            objective, residuals_count );

    ceres::Problem problem;
    problem.AddResidualBlock( cost_function, NULL, x );

    ceres::Solver::Options options;
    options.linear_solver_type        = ceres::DENSE_QR;
    options.parameter_tolerance       = error_threshold;
    options.function_tolerance        = error_threshold;
    options.min_line_search_step_size = error_threshold;
    options.max_num_iterations        = (int)max_steps;

    if ( verbosity > 2 )
        options.minimizer_progress_to_stdout = true;

    ceres::Solver::Summary summary;
    ceres::Solve( options, &problem, &summary );

    if ( verbosity >= 2 )
        std::cout << summary.FullReport() << std::endl;
    else if ( verbosity > 1 )
        std::cout << summary.BriefReport() << std::endl;

    return true;
}
#endif // ENABLE_CERES

bool solve_IDT_LAB(
    const std::vector<std::vector<double>> &src_points,
    const std::vector<std::vector<double>> &xyz_points,
    const std::vector<double>              &white_rgb,
    const std::vector<std::vector<double>> &out_to_xyz_mat,
    std::vector<std::vector<double>>       &out,
    int                                     verbosity )
{
    check_dependencies();

    if ( verbosity > 0 )
        std::cerr << "MatrixSolver: starting..." << std::endl;

    // We only need to solve for two unknowns per row, calculating the third
    // one so the lines add up to the destination white point, usually (1,1,1).
    // The initial guess is an identity matrix.
    double x[6] = { 1, 0, 0, 1, 0, 0 };

    bool result;

    // Ceres-solver takes ownership of this object. No need to delete.
    Objective_LAB *objective =
        new Objective_LAB( src_points, xyz_points, white_rgb, out_to_xyz_mat );

#ifdef ENABLE_CERES
    if ( use_ceres )
    {
        result = solve_nonlinear_ceres<Objective_LAB>(
            x, objective, 1e-20, 100, verbosity );
    }
    else
#endif // ENABLE_CERES
    {
        result = solve_nonlinear<Objective_LAB>(
            x, objective, 1e-20, 100, verbosity );
        delete objective;
    }

    if ( !result )
        return false;

    out.resize( 3 );

    out[0].resize( 3 );
    out[0][0] = x[0];
    out[0][1] = x[1];
    out[0][2] = white_rgb[0] - x[0] - x[1];

    out[1].resize( 3 );
    out[1][0] = x[2];
    out[1][1] = x[3];
    out[1][2] = white_rgb[1] - x[2] - x[3];

    out[2].resize( 3 );
    out[2][0] = x[4];
    out[2][1] = x[5];
    out[2][2] = white_rgb[2] - x[4] - x[5];

    if ( verbosity > 0 )
        std::cerr << "MatrixSolver: done!" << std::endl;

    return true;
}

std::vector<std::vector<double>> calculate_CAT(
    const std::vector<double> &src_white_XYZ,
    const std::vector<double> &dst_white_XYZ )
{
    assert( src_white_XYZ.size() == 3 );
    assert( dst_white_XYZ.size() == 3 );

    std::vector<double> src_white_LMS = product( CAT02, src_white_XYZ );
    std::vector<double> dst_white_LMS = product( CAT02, dst_white_XYZ );

    std::vector<std::vector<double>> mat( 3, std::vector<double>( 3, 0 ) );
    for ( size_t i = 0; i < 3; i++ )
        mat[i][i] = dst_white_LMS[i] / src_white_LMS[i];

    return product( CAT02_inv, product( mat, CAT02 ) );
}

// Instantiate some templates to avoid linking errors.
void __dummy()
{
    std::vector<std::vector<double>> mat1, mat2;
    inverse( mat1, mat2 );

    std::vector<std::vector<float>> mat3;
    transposed( mat3 );
}

} // namespace core
} // namespace rta
