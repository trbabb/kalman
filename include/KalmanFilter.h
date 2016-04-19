/* 
 * File:   KalmanFilter.h
 * Author: tbabb
 *
 * Created on February 3, 2015, 6:50 PM
 */

#ifndef KALMANFILTER_H
#define	KALMANFILTER_H

#include <geomc/function/Dual.h>
#include <geomc/linalg/Matrix.h>

using namespace geom;

#include "Sensor.h"
#include "Predictor.h"
#include "KalmanBuffer.h"

// todo: Control input not passed to Predictor
// todo: Should move to square root filter formulation for better stability
// todo: Consider adding an Unscented Kalman Filter (like kalman + particle filter)
// todo: can we get the process noise by looking at y~?
// todo: process noise can be a vector of different size.
//       it is an input to f().

template <typename T>
index_t measurement_size(Sensor<T>* s, index_t n) {
    index_t k = 0;
    for (index_t i = 0; i < n; i++) {
        k += s->reading_size();
    }
    return k;
}

template <typename T>
index_t measurement_maxwidth(Sensor<T>* s, index_t n) {
    index_t k = 0;
    for (index_t i = 0; i < n; i++) {
        k = std::max(k, s->reading_size());
    }
    return k;
}

/**
 * Estimates the true state of a possibly-multivariate system based on
 * gaussian-uncertain measurements of that state, using an extended Kalman 
 * filter (EKF).
 * 
 * One advantage of this class is that dynamically-varying numbers and types
 * of sensors are tolerated. The state can be updated as data becomes available,
 * and sensors with no new data may be omitted. Similarly, new sensors may be 
 * added to the estimation process as they come online.
 * 
 * This class works with Predictor and Sensor classes to automatically linearize
 * the state about the current estimate, in order to get accurate covariance
 * matrices.
 * 
 * If multiple sensors have noise with mutual covariance, they should be modeled 
 * as a single Sensor object.
 */
template <typename T>
class KalmanFilter {
    
  public:
    
    index_t  n;              // number of state fields
    T*       x;              // most recent state estimate
    SimpleMatrix<T,0,0> P;   // most recent estimate covariance
    Predictor<T> *predictor; // estimator of "next" state (pass to predict() instead?)
    KalmanBuffer<T> pool;    // buffer allocator
        
    KalmanFilter(index_t n=1, index_t m=1):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(NULL),
            pool(n,m,m) {
        std::fill(x, x+n, 0);
    }
    
    KalmanFilter(index_t n, index_t m, Predictor<T> *predictor):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(predictor),
            pool(n,m,m) {
        std::fill(x, x+n, 0);
    }
    
    KalmanFilter(index_t n, Sensor<T>* sensors, index_t n_sensors, Predictor<T> *predictor):
            n(n),
            x(new T[n]),
            P(n,n),
            predictor(predictor),
            pool(n, measurement_size(sensors, n_sensors), measurement_maxwidth(sensors, n_sensors) ) {
        std::fill(x, x+n, 0);
    }
    
    ~KalmanFilter() { delete [] x; }
    
    
    // todo: pass the state to/from sensors as a WrapperMatrix<T,0,1>
    // todo: handle control input to prediction
    // todo: some formulations factor the P matrix to maintain positive definite 
    //       symmetry. Can we just manually condition it? abs(max(elem, elem^T))?
    
    void predict(T t, T dt) {
        Dual<T>* x0 = pool.getX0(); // dual version of x
        Dual<T>* x1 = pool.getX1(); // predicted x
        WrapperMatrix<T,0,0>      F = pool.getPredictionJacobian();     // Jacobian of predictor fn
        WrapperMatrix<T,0,0> tmp_nn = pool.getPredictionMatrixBuffer(); // a temp
        
        std::copy(x, x + n, x0);
        // x_hat  <- f(x_{k-1}, u_{k-1})  (predict based on prev state and ctrl input)
        // F      <- J[f] dx              (compute jacobian of prediction with respect to state)
        
        // compute prediction and its covariance.
        for (index_t i = 0; i < n; i++) {
            // compute derivative in the ith direction
            x0[i].dx = 1;
            
            std::copy(x0, x0 + n, x1); // predictor might only update some variables.
            predictor->predict(x1, x0, t, dt);
            // copy derivatives to jacobian
            for (index_t j = 0; j < n; j++) {
                F[j][i] = x1[j].dx;
            }
            
            x0[i].dx = 0;
        }
         
        // P_k <- F * P_{k-1} * F^T + Q_k  (compute prediction covariance)
        
        mul(&tmp_nn, F, P);  // tmp <- F * P_{k-1}
        F.transpose();
        mul(&P, tmp_nn, F);  // P_k <- tmp * F^T
        
        // add covariance. We use F's memory as a temp.
        std::fill(F.begin(), F.end(), 0);
        predictor->covariance(F.begin(), x, t, dt);
        
        for (index_t i = 0; i < n * n; i++) P.begin()[i] += F.begin()[i];
        
        // xhat no longer needs to be dual
        for (index_t i = 0; i < n; i++) x[i] = x1[i].x;
    }
    
    
    void update(const Measurement<T>* observations, index_t n_obs) {
        // how many readings are we dealing with?
        index_t m = 0;
        for (index_t i = 0; i < n_obs; i++) {
            const Measurement<T>& z_i = observations[i];
            index_t m_i = z_i.sensor->reading_size();
            m += m_i;
        }
        
        if (m == 0) return; // no readings to use.
        
        // compute predicted sensor readings and covariance
        Dual<T>* z = pool.getZ(); // for holding derivatives                 (length m)
        WrapperMatrix<T,0,0>    H_k = pool.getHk(m);     // sensor jacobian      (m x n)
        WrapperMatrix<T,0,0>   H_kT = pool.getHkT(m);    // transpose of H_k     (n x m)
        WrapperMatrix<T,0,1>      y = pool.getY(m);      // measurement residual (m x 1)
        WrapperMatrix<T,0,0>    S_k = pool.getSk(m);     // residual covariance  (m x m)
        WrapperMatrix<T,0,0>    K_k = pool.getKk(m);     // Kalman gain          (n x m)
        WrapperMatrix<T,0,0> tmp_nn = pool.getTempNn();  // temporary mtx        (n x n)
        WrapperMatrix<T,0,1>  x_hat(x, n, 1);            // predicted state      (n x 1) (wraps x)
        
        // compute the sensor reading vector, and its derivative
        // with respect to the jth state variable. here we are linearizing
        // the readings about the measurement.
        Dual<T>* x0 = pool.getX0(); // state as a dual vector
        std::copy(x, x+n, x0);
        for (index_t j = 0; j < n; j++) {
            index_t i = 0;
            x0[j].dx  = 1; // query derivative in jth direction
            for (index_t k = 0; k < n_obs; k++) {
                const Measurement<T>& z_i = observations[k];
                index_t m_i = z_i.sensor->reading_size();
                z_i.sensor->measure(z + i, x0);
                i += m_i;
            }
            // copy derivatives to jacobian
            for (index_t a = 0; a < m; a++) { H_k[a][j] = z[a].dx; }
            x0[j].dx = 0;
        }

        // y <- z_k - h(x_hat)  (compute measurement residual)
        //                      (h was already invoked above and h(x_hat) remains in z)
        
        // copy measurments into y.
        index_t a = 0;
        for (index_t i = 0; i < n_obs; i++) {
            const Measurement<T>& z_i = observations[i];
            index_t m_i = z_i.sensor->reading_size();
            std::copy(z_i.data, z_i.data + m_i , y.begin() + a);
            a += m_i;
        }
        
        // subtract h(z).
        for (index_t i = 0; i < m; i++) {
            y[i][0] -= z[i].x;
        }
        
        // S_k <- H_k * P_k * H_k^T + R_k  (compute residual covariance)
        
        WrapperMatrix<T,0,0> tmp_nm = pool.getTempPHkT(m);
        transpose(&H_kT, H_k);
        mul(&tmp_nm, P, H_kT);    // tmp <- P_k * H_k^T
        mul(&S_k, H_k, tmp_nm);   // S_k <- H_k * tmp
        
        // add sensor noise covariance.
        // if sensors don't interfere with each other (an assumption),
        // the covariance matrix is block-diagonal.
        a = 0;
        for (index_t i = 0; i < n_obs; i++) {
            const Measurement<T>& z_i = observations[i];
            index_t m_i = z_i.sensor->reading_size();
            WrapperMatrix<T,0,0> mtx_tmp = pool.getTempSensorNoiseCovariance(m_i);
            
            z_i.sensor->covariance(mtx_tmp, z_i.data);
            
            index_t j = 0;
            auto    k = S_k.region_begin(MatrixRegion(a, a + m_i));
            auto  end = k.end();
            for (; k != end; k++) {
                *k += mtx_tmp.begin()[j++];
            }
            a += m_i;
        }

        // K_k <- P_k * H_k^T * S_k^-1
        
        WrapperMatrix<T,0,0> HkTSk_tmp = pool.getTempHkTSk(m); // (n x m)
        inv(&S_k, S_k);
        mul(&HkTSk_tmp, H_kT, S_k); // tmp <- H_k^T * S_k^-1 
        mul(&K_k, P, HkTSk_tmp);    // K_k <- P_k * tmp

        // x` <- x_hat + K_k * y  (compute new state estimate)
        
        WrapperMatrix<T,0,1> x_prime = pool.getTempX();
        mul(&x_prime, K_k, y);
        add(&x_hat, x_hat, x_prime);

        // P <- (I - K_k * H_k) * P_k  (compute new state covariance)

        WrapperMatrix<T,0,0> tmp_kk = pool.getTempKk();
        mul(&tmp_nn, K_k, H_k);
        for (index_t i = 0; i < n*n; i++) tmp_nn.begin()[i] *= -1;
        for (index_t i = 0; i < n; i++) {
            tmp_nn[i][i] += 1;
        }
        mul(&tmp_kk, tmp_nn, P);
        mtxcopy(&P, tmp_kk);
    }
    
    
    void advance(const Measurement<T>* observations, index_t n_obs, T t, T dt, index_t iters=5) {
        if (predictor) predict(t, dt);
        
        if (n_obs > 0) {
            T *xx = new T[n];
            for (int i = 0; i < iters; i++) {
                // tah-dah, we are an iterated extended kalman filter.
                std::copy(x, x + n, xx);
                update(observations, n_obs);
            }
            delete [] xx;
        }
    }
    
};


#endif	/* KALMANFILTER_H */

