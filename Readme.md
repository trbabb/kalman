Kalman
======

Implements a generalized extended Kalman filter.

Makes heavy use of [automatic differentiation](http://alexey.radul.name/ideas/2013/introduction-to-automatic-differentiation/) to automatically linearize your arbitrary, possibly-nonlinear system model and sensor mappings. 

Automatic differentiation is very fast, and computes derivatives to machine precision while avoiding finite differencing or symbolic manipulations. As long as all intermidiate computations are carried out using the `Dual` class, a sensible derivative will be computed.

The pieces are:

- A `Sensor` takes a state (made of `Dual`s) and converts it to a reading, using arbitrary math.
- A `Predictor` takes a state (made of `Dual`s) and a time delta, and produces a new state for t + dt, using arbitrary math.
- A `KalmanFilter` takes a `Predictor` and a list of `Sensors`, and constructs an EKF model for the system they describe.
- In an update loop, the `KalmanFilter` is fed `Measurement`s (i.e. a sensor object and its reading) and `dt`, and updates itself to the new optimal estimate of state.
