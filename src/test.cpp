
#include <iostream>
#include <algorithm>
#include <cmath>
#include <random>
#include <time.h>

#include <geomc/function/Dual.h>

#define SID_ACC  1
#define SID_VEL  2
#define SID_GYR  3
#define SID_MAG  4

#include "KalmanFilter.h"
#include "Solver.h"

#define OMEGA_PROCESS_VARIANCE   0.25
#define ACCEL_PROCESS_VARIANCE   0.025

#define MAX_READINGS 12


class SensorSimulator;


/********** simulation stuff **********/

vec3 cmp_accel(real_t t, rng_t *rng) {
    // white noise shall buffet our object
    return draw_random_vector(0.0, ACCEL_PROCESS_VARIANCE, rng);
}

vec3 cmp_omega(real_t t, vec3 base_omega, rng_t *rng) {
    // buffet with white noise
    return draw_random_vector(base_omega, OMEGA_PROCESS_VARIANCE, rng);
}

/********** simulation stuff **********/

struct delta_data {
    rng_t *rng;
    real_t dt;
    vec3 base_omega;
};

// accumulation function for kinematic states
//   s_1 = s_0 + d_ds * dt
// (operates on many states simultaneously)
void macc(KinematicState<real_t> *o, real_t k, const KinematicState<real_t> *d_ds, const KinematicState<real_t> *s, size_t n) {
    KinematicState<real_t> *end = o + n;
    for (; o != end; o++, d_ds++, s++) {
        o->x = s->x + k * d_ds->x;
        o->v = s->v + k * d_ds->v;
        o->a = s->a + k * d_ds->a;
        
        o->orient = std::exp(k * d_ds->orient / 2) * s->orient;
        o->omega  = s->omega + k * d_ds->omega;
    }
}


// delta function for kinematic states
//   d_ds = f(s0, t, data)
// A very simple model: inertially propagate current state, plus some noise on angular velocity.
void delta(KinematicState<real_t> *d_dt, const KinematicState<real_t> *s0, real_t t, size_t n, void *data) {
    delta_data *dat = (delta_data*) data;
    for (size_t i = 0; i < n; i++) {
        d_dt[i].x = s0[i].v;
        d_dt[i].v = s0[i].a;
        d_dt[i].a = (cmp_accel(t, dat->rng) - s0[i].a) / dat->dt;
        
        d_dt[i].orient = quat(s0[i].omega, 0);
        d_dt[i].omega  = (cmp_omega(t, dat->base_omega, dat->rng) - s0[i].omega) / dat->dt;
    }
}


/********** simulator **********/


class SensorSimulator {
    /*
     * This class keeps a ground-truth state and advances it according to physics
     * plus some random buffeting. The ground truth is used to generate simulated
     * sensor readings, which can then be fed to a KinematicSolver for state estimation.
     */
public:
    rng_t *rng;
    KalmanFilter<real_t> *filter;
    KinematicState<real_t> truth;
    
    SensorMagnetometer<real_t>   s_mag;
    SensorAccelerometer<real_t>  s_acc;
    SensorGPSVelocity<real_t>    s_vel;
    SensorRateGyro<real_t>       s_gyr;
    real_t t;
    vec3 base_omega;
    
    real_t measure_buf[3 * 5];
    Measurement<real_t> latest_obs;
    bool measure_available;
    
    SensorSimulator():
            // random number gen for fake sensor noise
            rng(new rng_t(11937294775LL)),
            // new kalman filter
            filter(new KalmanFilter<real_t>(KINSTATE_SIZE, 
                                            MAX_READINGS,
                                            // Predictor for advancing state:
                                            new KinematicPredictor<real_t>(
                                                    ACCEL_PROCESS_VARIANCE, 
                                                    OMEGA_PROCESS_VARIANCE))), 
            t(0) {
        init_sensors();
        init_state();
    }
    
    SensorSimulator(rng_t *rng):
            rng(rng),
            filter(new KalmanFilter<real_t>(KINSTATE_SIZE, 
                                            MAX_READINGS,
                                            new KinematicPredictor<real_t>(
                                                    ACCEL_PROCESS_VARIANCE, 
                                                    OMEGA_PROCESS_VARIANCE))), 
            t(0) {
        init_sensors();
        init_state();
    }
    
    void init_state() {
        measure_available = false;
        truth.x = 0.;
        truth.v = 0.;
        truth.a = 0.;
        truth.orient = quat(0,0,0,1);
        truth.omega  = base_omega = draw_random_vector(0., 0.25, rng);
        std::copy(truth.begin(), truth.end(), filter->x);
    }
    
    void init_sensors() {
        
        real_t sens_var = pow(0.025, 2);
        
        // store IDs in all the sensors, so we can identify them later
        s_acc._id = SID_ACC;
        s_vel._id = SID_VEL;
        s_gyr._id = SID_GYR;
        s_mag._id = SID_MAG;
        
        // set the gravitational center of earth
        s_acc.c_e = vec3(0,0,-6378100);
        
        // make up some sensor variances
        s_acc.set_variance(sens_var);
        s_vel.set_variance(1);
        s_gyr.set_variance(sens_var);
        s_mag.set_variance(sens_var);
    }
    
    real_t advance_dt(real_t dt) {
        int n_subsims = 8;
        KinematicState<real_t> s0 = truth;
        KinematicState<real_t> s1;
        KinematicState<real_t> b0, b1;
        
        delta_data dat = {rng, dt, base_omega};

        // advance the ground truth position/orientation
        real_t dt_i = dt / n_subsims;
        for (int i = 0; i < n_subsims; i++) {
            real_t t_i = t + i / (real_t)n_subsims;
            rk4_advance(&s1, 1, t_i, dt_i, &s0, &delta, &macc, &dat, &b0, &b1);
            s0 = s1;
        }
        truth = s0;

        index_t n_samples_per_frame = 3;
        // only sample velocity occasionally.
        if (std::uniform_int_distribution<int>(0,120)(*rng) == 1) n_samples_per_frame++;
        Sensor<real_t> *sens;
        std::vector< Measurement<real_t> > obs_list(n_samples_per_frame);
        
        // make a few random (fake) observations
        // we may not have data from all sensors, but that's OK
        for (int i = 0; i < n_samples_per_frame; i++) {
            real_t *b = measure_buf + (i * 3);
            switch (i) { //(std::uniform_int_distribution<int>(0,3)(*rng)) {
                case 0:
                    sens = &s_mag;
                    break;
                case 1:
                    sens = &s_acc;
                    break;
                case 2:
                    sens = &s_gyr;
                    break;
                case 3:
                    sens = &s_vel;
                    break;
            }
            Measurement<real_t> obs;
            obs.data = b;
            obs.sensor = sens;

            sens->simulate(b, truth.begin(), rng);

            latest_obs = obs_list[i] = obs;
            measure_available = true;
        }
        
        // update the filter's estimate, given some jittered sensor readings.
        clock_t start = clock();
        filter->advance(obs_list.data(), obs_list.size(), t, dt);
        clock_t end = clock();
        
        // xxx debug
        KinematicState<real_t> &fs = *(KinematicState<real_t>*)filter->x;
        fs.orient = fs.orient.unit();

        t += dt;
        
        return (end-start) / (real_t)CLOCKS_PER_SEC;
    }
    
};


/********** program state **********/


KinematicState<real_t> get_estimate(real_t *x) {
    KinematicState<real_t> s;
    std::copy(x, x + KINSTATE_SIZE, s.begin());
    return s;
}


void profile_kalman() {
    SensorSimulator ssim;
    real_t t_acc = 0;
    
    clock_t start = clock();
    for (index_t i = 0; i < 10000; i++) {
        t_acc += ssim.advance_dt(1/60.);
    }
    clock_t end = clock();
    real_t t_delt = (end-start) / (real_t)CLOCKS_PER_SEC;
    std::cout << "time: " << t_delt << " acc: " << t_acc << std::endl;
}


int main(int argc, char** argv) {
    profile_kalman();
    return 0;
}
