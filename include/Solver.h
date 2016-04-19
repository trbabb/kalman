/* 
 * File:   Solver.h
 * Author: tbabb
 *
 * Created on October 29, 2013, 8:44 PM
 */

#ifndef SOLVER_H
#define	SOLVER_H

#include <cstddef>

/**
 * Multiply / Accumulate.
 * Function for computing `m * x + b` over arrays `x` and `b`.
 * @param out Output array of ``State`s.
 * @param m Gain. 
 * @param x Array of `State`s to be scaled.
 * @param b Bias array of `State`s.
 * @param n Number of items.
 */
template <typename T, typename State>
void macc(State *out, T m, const State *x, const State *b, size_t n) {
    T *end = out + n;
    for (; out != end; out++, x++, b++) {
        *out = m * (*x) + (*b);
    }
}

// note that d_dt(state) is not necessarily the same shape as state.
// also, some of state may not be integrated. it would be better if dSdt had
// its own type, and we allowed for direct computation of intermediate and
// state variables.



/**
 * Advance a state with first-order euler integration.
 * 
 * @param s_final Output state.
 * @param n Number of parallel states.
 * @param t0 Time from which to advance.
 * @param dt Length of time to advance.
 * @param s0 State at time `t0`.
 * @param delta Function accepting a state `s`, a time `t`, an item count, and 
 *        a pointer to arbitrary external data, generating a rate of change for 
 *        each state variable at `t`.
 * @param accum Function for accumulating `n` (parameter 5) new states (parameter 3) 
 *        into existing states (parameter 4), with the new states having a 
 *        scaling factor of `k` (parameter 2) applied.
 * @param buf0 Optional pre-allocated buffer for intermediate work, of length `n`.
 */
template <typename T, typename State>
void euler_advance(State *s_final, size_t n,
                   T t0, T dt,
                   const State *s0, 
                   void (*delta)(State*, const State*, T, size_t, void*),
                   void (*accum)(State*, T, const State*, const State*, size_t),
                   void *external_data=NULL,
                   State *buf=NULL) {
    bool alloc_buf = external_data == NULL;
    if (alloc_buf) {
        buf = new State[n];
    }
    
    delta(buf, s0, t0, n, external_data);  //     buf <- dS/dt 
    accum(s_final, dt, buf, s0, n);        // s_final <- dt * dS/dt + s0
    
    if (alloc_buf) delete [] buf;
}

/**
 * Advance a state with 4th-order Runge-Kutta.
 * 
 * @param s_final Output state.
 * @param n Number of parallel states.
 * @param t0 Time from which to advance.
 * @param dt Length of time to advance.
 * @param s0 State at time `t0`.
 * @param delta Function accepting a state `s`, a time `t`, an item count, and 
 *        a pointer to arbitrary external data, generating a rate of change for 
 *        each state variable at `t`.
 * @param buf0 Optional pre-allocated buffer for intermediate work, of length `n`.
 * @param buf1 Optional pre-allocated buffer, distinct from `buf0`, of length `n`.
 */
template <typename T, typename State>
inline void rk4_advance(State *s_final, size_t n, 
                 T t0, T dt, 
                 const State *s0, 
                 void (*delta)(State*, const State*, T, size_t, void *),
                 void *external_data=NULL, 
                 State *buf0=NULL, 
                 State *buf1=NULL) {
    rk4_advance<T,State>(s_final, n, t0, dt, s0, delta, macc, external_data, buf0, buf1);
}


/**
 * Advance a state with 4th-order Runge-Kutta.
 * 
 * @param s_final Output state.
 * @param n Number of parallel states.
 * @param t0 Time from which to advance.
 * @param dt Length of time to advance.
 * @param s0 State at time `t0`.
 * @param delta Function accepting a state `s`, a time `t`, an item count, and 
 *        a pointer to arbitrary external data, generating a rate of change for 
 *        each state variable at `t`.
 * @param accum Function for accumulating `n` (parameter 5) new states (parameter 3) 
 *        into existing states (parameter 4), with the new states having a 
 *        scaling factor of `k` (parameter 2) applied.
 * @param buf0 Optional pre-allocated buffer for intermediate work, of length `n`.
 * @param buf1 Optional pre-allocated buffer, distinct from `buf0`, of length `n`.
 */
template <typename T, typename State>
void rk4_advance(State *s_final, size_t n, 
                 T t0, T dt, 
                 const State *s0, 
                 void (*delta)(State*, const State*, T, size_t, void *),
                 void (*accum)(State*, T, const State*, const State*, size_t),
                 void *external_data=NULL, 
                 State *buf0=NULL, 
                 State *buf1=NULL) {
    bool alloc_buf = (buf0 == NULL or buf1 == NULL);
    if (alloc_buf) {
        buf0 = new State[n*2];
        buf1 = buf0 + n;
    }
    
    T t1 = t0 + dt / 2.0;
    T t2 = t0 + dt;
    
    // k1 <- dSdt(t0, s0):
    delta(buf0, s0, t0, n, external_data);    // buf0 = k1
    
    // k2 <- dSdt(t0 + dt/2, s0 + k1 * dt/2):
    accum(s_final, dt/6, buf0, s0, n);        // s_final <- dt/6 * buf0 + s0
    accum(buf0, dt/2, buf0, s0, n);           // buf0    <- dt/2 * buf0 + s0
    delta(buf1, buf0, t1, n, external_data);  // buf1 = k2
    
    // k3 <- dSdt(t0 + dt/2, s0 + k2 * dt / 2):
    accum(s_final, dt/3, buf1, s_final, n);   // s_final <- dt/3 * buf1 + s_final
    accum(buf1, dt/2, buf1, s0, n);           // buf1    <- dt/2 * buf1 + s0
    delta(buf0, buf1, t1, n, external_data);  // buf0 = k3
    
    // k4 <- dSdt(t0 + dt, s0 + k3 * dt):
    accum(s_final, dt/3, buf0, s_final, n);   // s_final <- dt/3 * buf0 + s_final
    accum(buf0, dt, buf0, s0, n);             // buf0    <- dt   * buf0 + s0
    delta(buf1, buf0, t2, n, external_data);  // buf1 = k4
    
    accum(s_final, dt/6, buf1, s_final, n);   // s_final <- dt/6 * buf1 + s_final
    
    if (alloc_buf) delete [] buf0;
}


#endif	/* SOLVER_H */

