#pragma once

static float kalman_x0 = 0; // predicted state
static float kalman_p0 = 0.1f; // predicted covariance

static float kalman_f = 1; // factor of real value to previous real value
static float kalman_q = 2; // environment noise
static float kalman_h = 1; // factor of measured value to real value
static float kalman_r = 15; // measurement noise

static float kalman_state = 192;
static float kalman_covariance = 0.1f;


void kalman_init(float init_state, float init_covariance, float measurement_noise, float environment_noise)
{
	kalman_state = init_state;
	kalman_covariance = init_covariance;
	kalman_q = environment_noise;
	kalman_r = measurement_noise;
}
float kalman_filter(float value)
{
    //time update - prediction
    kalman_x0 = kalman_f * kalman_state;
    kalman_p0 = kalman_f * kalman_covariance * kalman_f + kalman_q;

    //measurement update - correction
    float k = kalman_h * kalman_p0 / (kalman_h * kalman_p0 * kalman_h + kalman_r);
    kalman_state = kalman_x0 + k * (value - kalman_h * kalman_x0);
    kalman_covariance = (1 - k * kalman_h) * kalman_p0;

    return kalman_state;
}