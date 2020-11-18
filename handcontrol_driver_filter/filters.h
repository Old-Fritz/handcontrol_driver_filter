#pragma once

// kalman filter
void kalman_init(float init_state, float init_covariance, float measurement_noise, float environment_noise);
float kalman_filter(float value);

void lowpass_init(float init_output, float cut_off_frequency, float delta_time);
float lowpass_filter(float value);
