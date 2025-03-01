// Copyright Contributors to the rawtoaces project.
// SPDX-License-Identifier: Apache-2.0
// https://github.com/AcademySoftwareFoundation/rawtoaces

#pragma once

#include <compare>
#include <vector>
#include <assert.h>

namespace rta
{
namespace core
{

struct Spectrum
{
    struct Shape
    {
        float first = 0, last = 0, step = 0;
        auto  operator<=>( const Shape  & ) const = default;
    } shape;

    inline static Shape ReferenceShape = { 380, 780, 5 };

    std::vector<double> values;

    Spectrum( double value = 0, const Shape &reference_shape = ReferenceShape )
        : shape( reference_shape )
    {
        if ( shape.step > 0 )
            values.resize(
                ( shape.last - shape.first + shape.step ) / shape.step, value );
    }

    template <typename Val_or_Ref, typename F>
    friend Val_or_Ref op( Val_or_Ref lhs, const Spectrum &rhs, F func )
    {
        assert( lhs.shape == rhs.shape );
        assert( lhs.values.size() == rhs.values.size() );

        auto l = lhs.values.begin();
        auto r = rhs.values.begin();
        while ( l != lhs.values.end() )
            *l++ = func( *l, *r++ );
        return lhs;
    }

    friend Spectrum operator+( Spectrum lhs, const Spectrum &rhs )
    {
        return op<Spectrum>( lhs, rhs, std::plus<double>() );
    }

    friend Spectrum operator-( Spectrum lhs, const Spectrum &rhs )
    {
        return op<Spectrum>( lhs, rhs, std::minus<double>() );
    }

    friend Spectrum operator*( Spectrum lhs, const Spectrum &rhs )
    {
        return op<Spectrum>( lhs, rhs, std::multiplies<double>() );
    }

    friend Spectrum operator/( Spectrum lhs, const Spectrum &rhs )
    {
        return op<Spectrum>( lhs, rhs, std::divides<double>() );
    }

    Spectrum &operator+=( const Spectrum &rhs )
    {
        return op<Spectrum &>( *this, rhs, std::plus<double>() );
    }

    Spectrum &operator-=( const Spectrum &rhs )
    {
        return op<Spectrum &>( *this, rhs, std::minus<double>() );
    }

    Spectrum &operator*=( const Spectrum &rhs )
    {
        return op<Spectrum &>( *this, rhs, std::multiplies<double>() );
    }

    Spectrum &operator/=( const Spectrum &rhs )
    {
        return op<Spectrum &>( *this, rhs, std::divides<double>() );
    }

    void reshape()
    {
        if ( shape == ReferenceShape )
            return;

        std::vector<double> temp;
        size_t              src    = 0;
        float               wl_src = shape.first;
        float               wl_dst = ReferenceShape.first;

        while ( wl_dst <= ReferenceShape.last )
        {
            if ( wl_src < wl_dst )
            {
                if ( src < values.size() - 1 )
                {
                    float next_wl_src = shape.first + shape.step * ( src + 1 );
                    if ( next_wl_src <= wl_dst )
                    {
                        // The next source wavelength is still not big enough,
                        // advancing.
                        src++;
                        wl_src = next_wl_src;
                    }
                    else
                    {
                        // The target wavelength is between two source samples,
                        // linearly interpolating.
                        double ratio =
                            ( wl_dst - wl_src ) / ( next_wl_src - wl_src );
                        double vv = values[src] * ( 1.0 - ratio ) +
                                    values[src + 1] * ratio;
                        temp.push_back( vv );
                        wl_dst = ReferenceShape.first +
                                 ReferenceShape.step * temp.size();
                    }
                }
                else
                {
                    // We have passed all available source samples,
                    // copying the last sample.
                    temp.push_back( values[src] );
                    wl_dst = ReferenceShape.first +
                             ReferenceShape.step * temp.size();
                }
            }
            else if ( wl_src == wl_dst )
            {
                // Found an exact match, just copy it over.
                temp.push_back( values[src] );
                wl_dst =
                    ReferenceShape.first + ReferenceShape.step * temp.size();
            }
            else
            {
                // Haven't reached the available source range yet, advancing.
                temp.push_back( values[src] );
                wl_dst =
                    ReferenceShape.first + ReferenceShape.step * temp.size();
            }
        }

        values = temp;
        shape  = ReferenceShape;
    }

    double integrate()
    {
        double result = 0;
        for ( auto &v: values )
            result += v;
        return result;
    }

    double max() const
    {
        double result = std::numeric_limits<float>::lowest();
        for ( auto &v: values )
            result = std::max( result, v );
        return result;
    }
};

} // namespace core
} // namespace rta
