#pragma once

#define FILTERS_COUNT 3

// kalman filter
void kalman_init(float init_state, float init_covariance, float measurement_noise, float environment_noise);
float kalman_filter(float value);

// lowpass filter
void lowpass_init(float init_output, float cut_off_frequency, float delta_time);
float lowpass_filter(float value);

// butter low pass filter
void butter_init(float init_output, float order, float cut_off_frequency, float sample_rate);
float butter_filter(float value);

// adjusting filter
void adjust_update_func(float(*filter_func)(float));
void adjust_init(float(*filter_func)(float), float min_lerp_diff, float max_lerp_diff, float min_speed, float max_speed);
float adjust_filter(float filter);
