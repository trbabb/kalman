/* 
 * File:   lintypes.h
 * Author: tbabb
 *
 * Created on June 16, 2014, 11:20 PM
 */

#ifndef LINTYPES_H
#define	LINTYPES_H

#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Ray.h>
#include <geomc/linalg/Matrix.h>
#include <geomc/linalg/Quaternion.h>
#include <geomc/function/Raster.h>
#include <geomc/random/RandomTools.h>

using namespace geom;

typedef double real_t;

typedef Vec<real_t,  2>  vec2;
typedef Vec<real_t,  3>  vec3;
typedef Vec<real_t,  4>  vec4;
typedef Vec<index_t, 2> ivec2;
typedef Vec<index_t, 3> ivec3;
typedef Vec<index_t, 4> ivec4;

typedef Ray<real_t, 2> ray2;
typedef Ray<real_t, 3> ray3;
typedef Hit<real_t, 3> rayhit;

typedef AffineTransform<real_t,3> transform_t;
typedef SimpleMatrix<real_t,2,2>  matrix2;
typedef SimpleMatrix<real_t,3,3>  matrix3;
typedef SimpleMatrix<real_t,4,4>  matrix4;

typedef Plane<real_t,3>   plane3;
typedef Sphere<real_t,3> sphere3;
typedef Rect<real_t,2>     rect2;
typedef Rect<real_t,3>     rect3;
typedef Rect<index_t,2>   irect2;
typedef Rect<index_t,3>   irect3;

typedef Quat<real_t> quat;

typedef Raster<real_t,real_t,2,1> img1ch;
typedef Raster<real_t,real_t,2,2> img2ch;
typedef Raster<real_t,real_t,2,3> img3ch;
typedef Raster<real_t,real_t,2,4> img4ch;

typedef Sampler<real_t> sampler;

#endif	/* LINTYPES_H */

