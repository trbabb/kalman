/* 
 * File:   Predictor.h
 * Author: tbabb
 *
 * Created on February 4, 2015, 1:18 PM
 */

#ifndef PREDICTOR_H
#define	PREDICTOR_H

#include <geomc/function/Dual.h>
#include <geomc/linalg/Orthogonal.h>

template <typename T>
class Predictor {
  public:
    // todo: handle control input
    
    virtual void covariance(T *cov, T *state, T t, T dt) = 0;
    
    virtual void predict(Dual<T> *out, const Dual<T> *state, T t, T dt) = 0;
};


template <typename T>
class KinematicPredictor : public Predictor<T> {
  public:
    
    T accel_variance;
    T omega_variance;
    
    
    KinematicPredictor():
            accel_variance(1),
            omega_variance(1) {}
    KinematicPredictor(T accel_variance, T omega_variance):
            accel_variance(accel_variance),
            omega_variance(omega_variance) {}
    
    
    void predict(Dual<T> *out, const Dual<T> *state, T t, T dt) {
        typedef KinematicState< Dual<T> > kinstate;
        typedef Quat< Dual<T> > quatx;
        
        kinstate &i = *(kinstate*)state;
        kinstate &o = *(kinstate*)out;
        
        o.x += 0.5 * i.a * dt * dt + i.v * dt;
        o.v +=       i.a * dt;
        o.orient = std::exp(quatx(dt * i.omega / 2, 0)) * i.orient.unit();
    }
    
    
    void covariance(T *cov, T *state, T t, T dt) {
        const index_t n = KINSTATE_SIZE;
        
        WrapperMatrix<T,0,0> cov_m(cov, n, n);
        cov_m.setZero();
        
        for (index_t i = 0; i < 3; i++) {
            // we assume that accel is goverend by gaussian random forces with a 
            // particular fixed variance. Velocity will be the integral of these
            // forces and is thus a gaussian random walk, and will be gaussian
            // distributed with variance(v) = t * variance(a). We can apply this
            // logic again to the distribution of v to obtain:
            //    variance(x) = t^2 * variance(a).
            index_t r = KINSTATE_IDX_A + i;
            cov_m[r][r] = accel_variance;
            r = KINSTATE_IDX_V + i;
            cov_m[r][r] = accel_variance * dt;
            r = KINSTATE_IDX_X + i;
            cov_m[r][r] = accel_variance * dt * dt;
            r = KINSTATE_IDX_OMEGA + i;
            cov_m[r][r] = omega_variance;
            // todo: use the formula that covariance = var[x]*var[y] * cov(x,y)
            //       and a frequency cutoff for acceleration so that 
            //       the variance of a attenuates on short time scales.
            //       we can then assume that v is linearly related to a if a is 
            //       smooth over that scale; and also that x is similarly related to v.
        }
        
        // we construct a gaussian blob about the current orientation quaternion.
        // this blob is essentially flattened to lie tangent to the surface of
        // the unit 3-sphere near q. The following matrix takes a "white" gaussian
        // blob with variance 1 and re-orients it so that the q_0 basis is aligned
        // with q. We'll use this in a moment.
        Vec<T,4> bases[4];
        bases[0] = Vec<T,4>(state + KINSTATE_IDX_ORIENT).unit();
        nullspace(bases, 1, bases + 1);
        orthonormalize(bases, 4);
        
        T sigm = dt * omega_variance;
        
        // omega acts tangent to the orientation, so omega's variance will
        // cause a gaussian distribution to spread with linearly-increasing variance
        // along those tangent directions. Basis 0 is the "normal", so all the
        // other bases get scaled by dt * omega_variance.
        for (index_t i = 1; i < 4; i++) { bases[i] *= sigm; }
        
        // for short timescales, the advancement of q is strictly along the 
        // tangent plane. However, for longer scales, q creeps "over the horizon";
        // this introduces some variance along the normal direction. The more
        // q advances, the less we know about where it could be, so the gaussian
        // blob becomes more and more round, and more uniform over the unit 3-ball. 
        // we construct things so that thickness(blob) / radius(blob) -> 1
        // as time progresses; from the above we see that radius increases 
        // linearly. why erf? it's roughly the right shape and I pulled it outta
        // mah butt. the actual math is too hard for wolframalpha.
        bases[0] *= sigm * std::erf(sigm);
        
        SimpleMatrix<T,4,4> m(bases[0].begin());
        
        // If A is a matrix acting on a random vector v, then the covariance
        // matrix of Av is A * cov[v] * A^T. Thus the covariance matrix for a
        // transformed "unit" distribution:
        m = transpose(m) * m;
        
        for (index_t r = 0; r < 4; r++) {
            // copy this row of m to a (section of a) row in the covariance mtx.
            T *section_start = cov + n * (KINSTATE_IDX_ORIENT + n) + KINSTATE_IDX_ORIENT;
            std::copy(m.row(r), m.row(r) + 4, section_start);  // xxx debug
        }
        
        const T scale = 0.01;
        for (index_t i = 0; i < n*n; i++) cov[i] *= scale;
    }
    
    
};


#endif	/* PREDICTOR_H */

