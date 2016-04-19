#include <geomc/linalg/Vec.h>
#include <geomc/linalg/Matrix.h>

/**
 * A memory pool class for the internals of the Kalman filter class.
 * Computes the exact quantity of memory needed and pre-allocates it, and
 * will alias storage for certain temporaries that don't have overlapping usage.
 * 
 * This construct, while rather ugly and cumbersome, avoids expensive and
 * redundant memory reallocation on every frame update. The speed difference 
 * versus alloc'ing is a little over 2x.
 * 
 * We must tell the pool ahead of time what the maximum number of sensor readings
 * `m` could be. Since we allow incremental updates with subsets of sensors,
 * per-frame queries may ask for smaller sensor matrices.
 */

template <typename T>
class KalmanBuffer {
    index_t  m;
    index_t  n;
    index_t  mmax;
    T*       buf_real;
    Dual<T>* buf_dual;
    index_t  n_real;
    index_t  n_dual;
    T*       blk_1;
    T*       blk_2;
    
public:
    
    // todo: move sensor-list initialization here.
    
    KalmanBuffer(index_t n, index_t m, index_t mmax):
            m(m), n(n), mmax(mmax) {
        n_real = 2 * n * n  +  3 * m * n  +  m * m  +  m  +  mmax * mmax;
        n_dual = n + std::max(m, n);
        buf_real = new T[n_real];
        buf_dual = new Dual<T>[n_dual];
        blk_1 = buf_real + 2 * n * n;
        blk_2 = blk_1 + 3 * m * n;
    }
    
    ~KalmanBuffer() {
        delete [] buf_real;
        delete [] buf_dual;
        buf_real = NULL;
        buf_dual = NULL;
    }
    
    // dual layout:
    // [ n (x0) ] [ n (x1) ]
    //            [ m (z)  ]
    
    
    inline Dual<T>* getX0() {
        return buf_dual;
    }
    
    inline Dual<T>* getX1() {
        return buf_dual + n;
    }
    
    inline Dual<T>* getZ() {
        return buf_dual + n;
    }
    
    
    // prediction layout:
    // [ n x n (jacobian) ] [ n x n (buffer) ]
    
    inline WrapperMatrix<T,0,0> getPredictionJacobian() {
        return WrapperMatrix<T,0,0>(buf_real, n, n);
    }
    
    inline WrapperMatrix<T,0,0> getPredictionMatrixBuffer() {
        return WrapperMatrix<T,0,0>(buf_real + n * n, n, n);
    }
    
    
    // update layout:
    // [       2 (n x n) blk_0           ] [         3 (m x n) blk 1                            ] [       oddball blk 2                        ]
    // [ n x n (tmp_nn) ] [ n x n K_k_tmp] [ m x n (H_k) ] [ n x m (H_kT) ] [ n x m (K_k)       ] [ m x m (S_k) ] [ m (y) ] [ mmax^2 (mtx_tmp) ]
    // [ n     (x_tmp)  ]                                  [ PK_k_tmp     ] [ n x m (HkTSk tmp) ]
    //                                                                      [ n x m (PH_kT tmp) ]
    
    inline WrapperMatrix<T,0,0> getTempNn() {
        return WrapperMatrix<T,0,0>(buf_real, n, n);
    }
    
    inline WrapperMatrix<T,0,0> getTempKk() {
        return WrapperMatrix<T,0,0>(buf_real + n * n, n, n);
    }
    
    inline WrapperMatrix<T,0,1> getTempX() {
        return WrapperMatrix<T,0,1>(buf_real, n);
    }
    
    inline WrapperMatrix<T,0,0> getHk(index_t m0) {
        return WrapperMatrix<T,0,0>(blk_1, m0, n);
    }
    
    inline WrapperMatrix<T,0,0> getHkT(index_t m0) {
        return WrapperMatrix<T,0,0>(blk_1 + m * n, n, m0);
    }
    
    inline WrapperMatrix<T,0,0> getKk(index_t m0) {
        return WrapperMatrix<T,0,0>(blk_1 + 2 * m * n, n, m0);
    }
    
    inline WrapperMatrix<T,0,0> getTempPKk(index_t m0) {
        return getHkT(m0);
    }
    
    inline WrapperMatrix<T,0,0> getTempHkTSk(index_t m0) {
        return getKk(m0);
    }
    
    inline WrapperMatrix<T,0,0> getTempPHkT(index_t m0) {
        return getKk(m0);
    }
    
    inline WrapperMatrix<T,0,0> getSk(index_t m0) {
        return WrapperMatrix<T,0,0>(blk_2, m0, m0);
    }
    
    inline WrapperMatrix<T,0,1> getY(index_t m0) {
        return WrapperMatrix<T,0,1>(blk_2 + m * m, m0, 1);
    }
    
    inline WrapperMatrix<T,0,0> getTempSensorNoiseCovariance(index_t sensor_width) {
        return WrapperMatrix<T,0,0>(blk_2 + m * m + m, sensor_width, sensor_width);
    }
    
};
