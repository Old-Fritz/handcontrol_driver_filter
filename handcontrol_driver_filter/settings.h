#pragma once

int settings_load_from_file(char* filename);

void settings_init_default();

void setting_get_input_filename(const char** filename);

void settings_get_lowpass_params(float* cut_off_frequency, float* delta_time, const char** output_filename);

void settings_get_kalman_params(float* init_covariance, float* measurement_noise, float* environment_noise, const char** output_filename);

void settings_get_butter_params(float* order, float* cut_off_frequency, float* sample_rate, const char** output_filename);

void settings_get_adjust_params(float* min_lerp_diff, float* max_lerp_diff, float* min_speed, float* max_speed);

