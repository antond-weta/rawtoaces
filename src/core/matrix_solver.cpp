// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#include <assert.h>

#include <rawtoaces/matrix_solver.hpp>
#include <rawtoaces/mathOps.h>

namespace rta
{
namespace core
{

constexpr double v_6_29_3 = 6.0 * 6.0 * 6.0 / 29.0 / 29.0 / 29.0;
constexpr double v_29_6_3 = 29.0 * 29.0 / 6.0 / 6.0;

double F( double val )
{
    if ( val > v_6_29_3 )
        return pow( val, 1.0 / 3.0 );
    else
        return val * v_29_6_3 / 3.0;
}

double F_( double val, double val_prime )
{
    if ( val > v_6_29_3 )
        return 1.0 / 3.0 * pow( val, -2.0 / 3.0 ) * val_prime;
    else
        return val_prime * v_29_6_3 / 3.0;
}

void lab( const double ( &in )[3], double ( &out )[3] )
{
    out[0] = 116.0 * F( in[1] ) - 16.0;
    out[1] = 500.0 * ( F( in[0] ) - F( in[1] ) );
    out[2] = 200.0 * ( F( in[1] ) - F( in[2] ) );
}

void lab_(
    const double ( &in )[3], const double ( &in_ )[3], double ( &out )[3] )
{
    out[0] = 116.0 * F_( in[1], in_[1] );
    out[1] = 500.0 * ( F_( in[0], in_[0] ) - F_( in[1], in_[1] ) );
    out[2] = 200.0 * ( F_( in[1], in_[1] ) - F_( in[2], in_[2] ) );
}

// Calculates the residuals of the i-th equation when j == -1,
// or the partial derivative of the residual with regards to the j-th variable
// if j >= 0. The residuals are the colour distance in the CIELAB space between
// the point (L,A,B) calculated via the matrix being solved, and the true value
// (Lt,At,Bt) calculated using the standard observer curves.
//    (L-Lt, A-At, B-Bt) when j == -1
//    (L'  , A'  , B'  ) when j >=  0
void MatrixSolver::objective( size_t i, int j, double ( &result )[3] )
{
    // First calculate the values which are required in all cases.

    auto  &src = _camera_points[i];
    double Rs  = src[0];
    double Gs  = src[1];
    double Bs  = src[2];

    auto  &dst = _xyz_points[i];
    double Xd  = dst[0];
    double Yd  = dst[1];
    double Zd  = dst[2];

    double Rw = _white_rgb[0];
    double Gw = _white_rgb[1];
    double Bw = _white_rgb[2];

    // Calculate the values in the output space by applying the colour matrix
    // which we are solving for to the values in the camera space.
    // The matrix is
    //     | x0 x1 (Rw-x0-x1) |
    //     | x2 x3 (Gw-x2-x3) |
    //     | x4 x5 (Bw-x4-x5) |

    // Ro = x0*Rs + x1*Gs + (Rw-x0-x1)*Bs
    // Go = x2*Rs + x3*Gs + (Gw-x2-x3)*Bs
    // Bo = x4*Rs + x5*Gs + (Zw-x4-x5)*Bs

    // Ro = x0*(Rs-Bs) + x1*(Gs-Bs) + Rw*Bs
    // Go = x2*(Rs-Bs) + x3*(Gs-Bs) + Gw*Bs
    // Bo = x4*(Rs-Bs) + x5*(Gs-Bs) + Bw*Bs

    double Ro = _x[0] * ( Rs - Bs ) + _x[1] * ( Gs - Bs ) + Rw * Bs;
    double Go = _x[2] * ( Rs - Bs ) + _x[3] * ( Gs - Bs ) + Gw * Bs;
    double Bo = _x[4] * ( Rs - Bs ) + _x[5] * ( Gs - Bs ) + Bw * Bs;

    // Convert the values from the output colour space to CIEXYZ using the
    // matrix provided by the user, that is usually the ACES_to_XYZ, when
    // solving for an ACES IDT, or an identity matrix when solving for a
    // camera RGB to XYZ matrix.

    // X = m00*Rd + m01*Gd + m02*Bd
    // Y = m10*Rd + m11*Gd + m12*Bd
    // Z = m20*Rd + m21*Gd + m22*Bd

    const double( &m )[3][3] = _out_to_xyz_mat;

    double XYZ[3];
    double XYZw[3];
    double XYZww[3];
    double XYZdww[3];
    for ( size_t i = 0; i < 3; i++ )
    {
        XYZ[i]    = m[i][0] * Ro + m[i][1] * Go + m[i][2] * Bo;
        XYZw[i]   = m[i][0] * Rw + m[i][1] * Gw + m[i][2] * Bw;
        XYZww[i]  = XYZ[i] / XYZw[i];
        XYZdww[i] = dst[i] / XYZw[i];
    }

    double LAB[3];
    lab( XYZww, LAB );

    double LABd[3];
    lab( XYZdww, LABd );

    if ( j == -1 )
    {
        for ( size_t i = 0; i < 3; i++ )
            result[i] = LAB[i] - LABd[i];
        return;
    }

    // Ro = x0*(Rs-Bs) + x1*(Gs-Bs) + Rw*Bs
    // Go = x2*(Rs-Bs) + x3*(Gs-Bs) + Gw*Bs
    // Bo = x4*(Rs-Bs) + x5*(Gs-Bs) + Bw*Bs

    // Ro_ = (Rs-Bs) when j == 0, 0 otherwise
    // Ro_ = (Gs-Bs) when j == 1, 0 otherwise
    // Go_ = (Rs-Bs) when j == 2, 0 otherwise
    // Go_ = (Gs-Bs) when j == 3, 0 otherwise
    // Bo_ = (Rs-Bs) when j == 4, 0 otherwise
    // Bo_ = (Gs-Bs) when j == 5, 0 otherwise

    double Ro_ = 0;
    double Go_ = 0;
    double Bo_ = 0;

    switch ( j )
    {
        case 0: Ro_ = Rs - Bs; break;
        case 1: Ro_ = Gs - Bs; break;
        case 2: Go_ = Rs - Bs; break;
        case 3: Go_ = Gs - Bs; break;
        case 4: Bo_ = Rs - Bs; break;
        case 5: Bo_ = Gs - Bs; break;
    }

    double XYZ_[3];
    double XYZ_ww[3];

    for ( size_t i = 0; i < 3; i++ )
    {
        XYZ_[i]   = m[i][0] * Ro_ + m[i][1] * Go_ + m[i][2] * Bo_;
        XYZ_ww[i] = XYZ_[i] / XYZw[i];
    }

    lab_( XYZww, XYZ_ww, result );
}

std::vector<std::vector<double>> MatrixSolver::solve()
{
    // We only need to solve for two unknowns per row, calculating the third
    // one so the lines add up to the destination white point, usually (1,1,1).
    // The initial guess is an identity matrix.
    _x[0] = 1;
    _x[1] = 0; /* 0 */
    _x[2] = 0;
    _x[3] = 1; /* 0 */
    _x[4] = 0;
    _x[5] = 0; /* 1 */

    // Pre-calculate the values which don't change between the steps.

    double Xw = _white_xyz[0];
    double Yw = _white_xyz[1];
    double Zw = _white_xyz[2];

    _lab_points.resize( _xyz_points.size() );
    for ( size_t i = 0; i < _xyz_points.size(); i++ )
    {
        const auto &xyz = _xyz_points[i];
        double      XYZww[3];
        for ( size_t j = 0; j < 3; j++ )
            XYZww[j] = xyz[j] / _white_xyz[j];

        double LAB[3];
        lab( XYZww, LAB );

        _lab_points[i] = std::vector<double>( LAB, LAB + 3 );
    }

    size_t                           step  = 0;
    double                           error = std::numeric_limits<double>::max();
    std::vector<std::vector<double>> jacobian(
        _camera_points.size() * 3, std::vector<double>( 6, 0 ) );
    std::vector<double> residual( _camera_points.size() * 3 );

    while ( step < max_steps && error > error_threshold )
    {

        double curr_error = 0;

        for ( size_t i = 0; i < _camera_points.size(); i++ )
        {
            for ( int j = -1; j < 6; j++ )
            {
                double val[3];
                objective( i, j, val );

                for ( size_t k = 0; k < 3; k++ )
                {
                    double v = val[k];
                    if ( j == -1 )
                    {
                        residual[i * 3 + k] = v;
                        curr_error += v * v;
                    }
                    else
                    {
                        jacobian[i * 3 + k][j] = v;
                    }
                }
            }
        }

        error = curr_error;

        auto jacobian_t = transposeVec( jacobian );
        auto A          = mulVector( jacobian_t, jacobian_t );
        auto B          = mulVector( jacobian_t, residual );

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> AA;
        AA.resize( A.size(), A[0].size() );
        FORIJ( AA.rows(), AA.cols() )
        AA( i, j ) = A[i][j];

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> BB;
        BB.resize( B.size(), 1 );
        FORI( BB.rows() )
        BB( i, 0 ) = B[i];

        Eigen::VectorXd S = AA.colPivHouseholderQr().solve( BB );

        for ( size_t i = 0; i < 6; i++ )
        {
            _x[i] -= S( i );
        }

        //        std::cerr << "Iteration: " << std::setw(3) << step << std::setw(0);
        //        std::cerr << std::scientific;
        //        std::cerr << "  Error: " << error * 10000 / 2;
        //        std::cerr << "  Curr error: " << curr_error * 10000 / 2;
        //        std::cerr << "  mult: " << mult << std::endl;
        //        std::cerr << std::fixed;

        step++;
    }

    std::vector<std::vector<double>> result( 3, std::vector<double>( 3, 0 ) );
    result[0][0] = _x[0];
    result[0][1] = _x[1];
    result[0][2] = _white_rgb[0] - _x[0] - _x[1];
    result[1][0] = _x[2];
    result[1][1] = _x[3];
    result[1][2] = _white_rgb[1] - _x[2] - _x[3];
    result[2][0] = _x[4];
    result[2][1] = _x[5];
    result[2][2] = _white_rgb[2] - _x[4] - _x[5];

    return result;
}

} // namespace core
} // namespace rta
