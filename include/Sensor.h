/* 
 * File:   Sensor.h
 * Author: tbabb
 *
 * Created on February 3, 2015, 1:38 PM
 */

#ifndef SENSOR_H
#define	SENSOR_H

#include <algorithm>
#include <random>

#include <geomc/function/Dual.h>
#include <geomc/linalg/Vec.h>
#include <geomc/linalg/AffineTransform.h>

#include "lintypes.h"
#include "KinematicState.h"
 
 // todo: store covariance as a mtx square root?

/***************************************
 * Random number drawing               *
 ***************************************/


typedef std::mt19937_64                        rng_t;
typedef std::normal_distribution<real_t>       d_normal_t; 
typedef std::uniform_real_distribution<real_t> d_uniform_t;


// sample the value of a normal distribution pdf at x
real_t normal(real_t x, real_t mean, real_t variance) {
    variance = std::abs(variance);
    real_t p = x - mean;
    real_t q = (p * p) / (2 * variance); 
    return std::exp(-q) / std::sqrt(variance * 2 * M_PI);
}


vec3 draw_random_vector(const vec3 &ctr, const vec3 &variance, rng_t *rng) {
    vec3 o;
    for (index_t i = 0; i < 3; i++) {
        d_normal_t N = d_normal_t(ctr[i], std::sqrt(std::abs(variance[i])));
        o[i] = N(*rng);
    }
    return o;
}

quat draw_quat(rng_t *rng) {
    quat q;
    for (index_t i = 0; i < 4; i++) {
        q[i] = d_normal_t(0, 1)(*rng);
    }
    return q.unit();
}

/***************************************
 * Sensor base class                   *
 ***************************************/

/**
 * A Sensor is an object which maps an underlying (unobservable) state to a 
 * reading.
 * 
 * Characteristics of this mapping, as well as sensor noise, may be used by a 
 * state estimator to optimally solve for the underlying state.
 */
template <typename T>
class Sensor {
  public:
    
    index_t _state_size;
    index_t _measurement_size;
    index_t _id;
    
    /// Construct a new Sensor object
    Sensor():
            _state_size(1),
            _measurement_size(1),
            _id(0) {}
    /// Construct a new Sensor object
    Sensor(index_t state_size, index_t measurement_size):
            _state_size(state_size),
            _measurement_size(measurement_size),
            _id(0) {}
    /// Construct a new Sensor object
    Sensor(index_t state_size, index_t measurement_size, index_t id) :
            _state_size(state_size),
            _measurement_size(measurement_size),
            _id(id) {}
            
    virtual ~Sensor() {}
    
    /// Number of expected elements in a state.
    inline index_t state_size()   { return _state_size; }
    /// Number of elements in a reading.
    inline index_t reading_size() { return _measurement_size; }
    /// Sensor ID
    inline index_t id() { return _id; }
    
    /** 
     * Return an exact measurement of the given state. 
     * 
     * The Duals provide a mechanism for obtaining a linearization of the change
     * of measurement with respect to state about the given point (this makes
     * Sensor objects amenable to the Extended Kalman Filter). All arithmetic 
     * involving elements of the state should therefore be performed on Dual objects down to 
     * the level of the basic arithmetic operators and elementary functions in 
     * order to preserve a sensible derivative.
     */
    virtual void measure(Dual<T> *measurement, const Dual<T> *state) = 0;
    
    /**
     * Return the sensor noise covariance near a particular measurement as an 
     * M x M matrix. The covariance is taken to be in reading-space.
     * The current measurement is passed as a convenience, so that sensors
     * with nonlinear covariance may be modeled.
     */
    virtual void covariance(WrapperMatrix<T,0,0> cov, const T* measurement) const = 0;
    /// Return the likelihood that a particular state and measurement coexist.
    virtual    T likelihood(const T* state, const T* measurement) = 0;
    /// Return a stochastic measurement of the given state, with noise added according to the characteristics of the sensor.
    virtual void simulate(T* measurement, const T* state, rng_t* rng) = 0;
    
};


/***************************************
 * Oriented Sensor                     *
 ***************************************/

/**
 * An oriented sensor is one that is associated with a rigid body. Its alignment
 * with world-space quantities changes based on the orientation described by the state.
 * In addition to rigid body orientation, the sensor may have a fixed orientation 
 * relative to the body itself. This fixed orientation, along with sensor biases 
 * and gains, is encapsulated by the state2reading matrix.
 */
template <typename T>
class SensorOriented : public Sensor<T> {
  public:
    /// Matrix encapsulating fixed body-relative orientation, and sensor biases and gains.
    AffineTransform<T,3> state2reading;
    /// Covariance matrix, in reading-space. 
    SimpleMatrix<T,3,3>  covariance_mtx;
    
    /// A dual 3-vector.
    typedef Vec<Dual<T>,3> diff_measure_t;
    
    /// Construct a new OrientedSensor
    SensorOriented()          :Sensor<T>(KINSTATE_SIZE, 3)     {}
    /// Construct a new OrientedSensor
    SensorOriented(index_t id):Sensor<T>(KINSTATE_SIZE, 3, id) {}
    
    /**
     * Return the measured quanitity in the co-rotating coordinate frame.
     * This is called by measure() to return a quantity in measurement space.
     */
    virtual diff_measure_t body_space_state(const KinematicState< Dual<T> > &state) = 0;
    
    
    Vec<T,3> body_space_state(const KinematicState<T> &state) {
        KinematicState< Dual<T> > buf;
        std::copy(state.begin(), state.end(), buf.begin());
        diff_measure_t m = body_space_state(buf);
        
        return (Vec<T,3>)m;
    }
    
    /// Returns a dual quaternion orienting the body in world space.
    template <typename U>
    inline Quat<U> body2world(const KinematicState<U> &state) {
        return state.orient.unit();
    }
    
    void measure(Dual<T> *measurement, const Dual<T> *state) {
        // reinterpret `state` as a kinematic state object
        const KinematicState< Dual<T> > &ks = *((const KinematicState< Dual<T> >*)state);
        // measure it in body space
        Vec<Dual<T>,4> body_qty = Vec<Dual<T>,4>(body_space_state(ks), 1);
        // transform the reading to sensor space
        diff_measure_t sens_qty = (state2reading.mat * body_qty).template resized<3>(); // this works; T becomes the elem type of the vector
        // save the result to output
        std::copy(sens_qty.begin(), sens_qty.end(), measurement);
    }
    
    void covariance(WrapperMatrix<T,0,0> cov, const T *measurement) const {
        mtxcopy(&cov, covariance_mtx);
    }
    
    void set_variance(Vec<T,3> v) {
        covariance_mtx[0][0] = v[0];
        covariance_mtx[1][1] = v[1];
        covariance_mtx[2][2] = v[2];
    }
    
    T likelihood(const T *state, const T *measurement) {
        // todo: this.
        return 1;
    }
    
    void simulate(T *measurement, const T* state, rng_t *rng) {
        // XXX TODO. the full covariance matrix is not used. to do so
        //           we need the mtx sqrt.
        const KinematicState<T> &ks = *((const KinematicState<T>*)state);
        Vec<T,3> reading = state2reading * this->body_space_state(ks);
        Vec<T,3> v = draw_random_vector(reading, 
                                    Vec<T,3>(covariance_mtx[0][0],
                                             covariance_mtx[1][1],
                                             covariance_mtx[2][2]),
                                    rng);
        
        // when we have the matrix sqrt, we will draw with variance 1
        // and then transform the sample by said matrix.
        
        measurement[0] = v[0];
        measurement[1] = v[1];
        measurement[2] = v[2];
    }
    
};


/***************************************
 * Inertial Sensor                     *
 ***************************************/


template <typename T>
class SensorInertial : public Sensor<T> {
   public:
    SimpleMatrix<T,3,3> covariance_mtx;
    
    typedef Vec<Dual<T>,3> diff_measure_t;
    
    SensorInertial()           : Sensor<T>(KINSTATE_SIZE, 3) {}
    SensorInertial(index_t id) : Sensor<T>(KINSTATE_SIZE, 3, id) {}
    
    void set_variance(const Vec<T,3> &var) {
        covariance_mtx[0][0] = var[0];
        covariance_mtx[1][1] = var[1];
        covariance_mtx[2][2] = var[2];
    }
    
    virtual diff_measure_t inertial_measurement(const KinematicState< Dual<T> > &state) = 0;
    
    Vec<T,3> inertial_measurement(const KinematicState<T> &state) {
        KinematicState< Dual<T> > buf;
        std::copy(state.begin(), state.end(), buf.begin());
        diff_measure_t m = inertial_measurement(buf);
        
        return (Vec<T,3>)m;
    }
    
    void measure(Dual<T> *measurement, const Dual<T> *state) {
        const KinematicState< Dual<T> > &s = *((const KinematicState< Dual<T> >*)state);
        std::copy(s.v.begin(), s.v.end(), measurement);
    }
    
    void covariance(WrapperMatrix<T,0,0> cov, const T* measurement) const {
        mtxcopy(&cov, covariance_mtx);
    }
    
    T likelihood(const T *state, const T *measurement) {
        // xxx todo
        return 1;
    }
    
    void simulate(T *measurement, const T *state, rng_t *rng) {
        const KinematicState<T> &s = *reinterpret_cast<const KinematicState<T>*>(state);
        Vec<T,3> v = draw_random_vector(inertial_measurement(s), 
                                    Vec<T,3>(covariance_mtx[0][0],
                                             covariance_mtx[1][1],
                                             covariance_mtx[2][2]),
                                    rng);
        std::copy(v.begin(), v.end(), measurement);
    }
};


/***************************************
 * Kinematic Sensors                   *
 ***************************************/

/// Accelerometer sensor
template <typename T>
class SensorAccelerometer : public SensorOriented<T> {
  public:
    Vec<T,3> c_e;
    
    using diff_measure_t = typename SensorOriented<T>::diff_measure_t;
    
    SensorAccelerometer():
        c_e(0,0,-6378100) {}
    SensorAccelerometer(index_t id):
        SensorOriented<T>(id),
        c_e(0,0,-6378100) {}
    
    diff_measure_t body_space_state(const KinematicState< Dual<T> > &state) {// todo: model 1/r^2
        // todo: gravity model
        // could later depend on ECEF coordinates, i.e. WGS84 or EGM2008
        // http://earth-info.nga.mil/GandG/wgs84/gravitymod/
        // todo: model offset from cm?
        diff_measure_t g = ((diff_measure_t)c_e - state.x).unit() * 9.80;
        return (diff_measure_t(state.a /*0.0*/) - g) * this->body2world(state);
    }
};


/// Magnetometer/compass sensor
template <typename T>
class SensorMagnetometer : public SensorOriented<T> {
  public:
    
    using diff_measure_t = typename SensorOriented<T>::diff_measure_t;
    
    SensorMagnetometer() {}
    SensorMagnetometer(index_t id):SensorOriented<T>(id) {}
    
    template <typename U>
    Vec<U,3> field_inertial(const KinematicState<U> &state) {
        // todo: compute this
        // http://www.ngdc.noaa.gov/geomag/geomag.shtml
        typedef Vec<U,3> V; 
        // 61.8  = degrees dip in san francisco
        // 0.486 = field strength, gauss.
        return V(1,0,0).rotate(V(0,1,0), 61.8 * M_PI / 180.0) * 0.486; 
    }
    
    diff_measure_t body_space_state(const KinematicState< Dual<T> > &state) {
        return field_inertial(state) * this->body2world(state);
    }
};


/// Rate gyro sensor
template <typename T>
class SensorRateGyro : public SensorOriented<T> {
  public:
    
    using diff_measure_t = typename SensorOriented<T>::diff_measure_t;
    
    SensorRateGyro() {}
    SensorRateGyro(index_t id):SensorOriented<T>(id) {}
    
    diff_measure_t body_space_state(const KinematicState< Dual<T> > &state) {
        return state.omega * this->body2world(state);
    }
    
};

// GPS velocity sensor.
template <typename T>
class SensorGPSVelocity : public SensorInertial<T> {
   public:
    
       using diff_measure_t = typename SensorInertial<T>::diff_measure_t;
       
       SensorGPSVelocity() {}
       SensorGPSVelocity(index_t id):SensorInertial<T>(id) {}
       
       diff_measure_t inertial_measurement(const KinematicState< Dual<T> > &state) {
           return state.v;
       }
};

// GPS location.
template <typename T>
class SensorGPSLocation : public SensorInertial<T> {
   public:
    
       using diff_measure_t = typename SensorInertial<T>::diff_measure_t;
       
       SensorGPSLocation() {}
       SensorGPSLocation(index_t id):SensorInertial<T>(id) {}
       
       diff_measure_t inertial_measurement(const KinematicState< Dual<T> > &state) {
           return state.x;
       }
};


/***************************************
 * Measurement                         *
 ***************************************/

template <typename T>
struct Measurement {
    Sensor<T> *sensor;
    T         *data;
    
    Measurement(Sensor<T> *s, T *data) :
            sensor(s),
            data(data) {}
    Measurement() : 
            sensor(NULL), 
            data(NULL) {}
};

// Quat< Dual<real_t> > q(state + KINSTATE_IDX_ORIENT);

#endif	/* SENSOR_H */

