#pragma once

#include "settings.h"
#include "logger.h"

#include <stdio.h>
#include <string.h>

#define BUFFER_SIZE 512


static char  input_filename[BUFFER_SIZE] = "../full.csv";

static float kalman_init_covariance = 0.1f;
static float kalman_measurement_noise = 15.0f;
static float kalman_environment_noise = 2.0f;
static char  kalman_output_filename[BUFFER_SIZE] = "../full_filtered_kalman";

static float lowpass_cut_off_frequency = 200.0f;
static float lowpass_delta_time = 0.001f;
static char  lowpass_output_filename[BUFFER_SIZE] = "../full_filtered_lowpass";

static float butter_order = 6.0f;
static float butter_cut_off_frequency = 50.0f;
static float butter_sample_rate = 1000.0f;
static char  butter_output_filename[BUFFER_SIZE] = "../full_filtered_butter";

static float adjust_min_lerp_diff = 30.0f;
static float adjust_max_lerp_diff = 100.0f;
static float adjust_min_speed = 0.02f;
static float adjust_max_speed = 0.2f;

#define FAIL_LOAD(read_args, param_name) if (!read_args) { return FAIL; }

static int load_param(FILE * file, const char* param_name, void* param, char printf_type)
{
	char format[BUFFER_SIZE];
	
	sprintf_s(format, BUFFER_SIZE, "%s = %%%c\n", param_name, printf_type);

	int read_args = fscanf_s(file, format, (char*)param, BUFFER_SIZE);
	if (!read_args)
	{
		log("Failed to load param: %s", param_name);
		return FAIL;
	}
	
	sprintf_s(format, BUFFER_SIZE, "%s loaded from settings: %%%c", param_name, printf_type);

	switch (printf_type)
	{
	case 's':
		log(format, (char*)param);
		break;
	case 'f':
		log(format, *(float*)param);
	}

	return SUCCESS;
}

#define LOAD_PARAM(file, param_name, printf_type) if( load_param(file, #param_name, &param_name, printf_type) == FAIL ) { return FAIL;}
#define LOAD_FLOAT(file, param_name) LOAD_PARAM(file, param_name, 'f')
#define LOAD_STR(file, param_name) LOAD_PARAM(file, param_name, 's')


int settings_load_from_file(char* filename)
{
	FILE* file;
	int result = fopen_s(&file, filename, "r");
	if (!file || result)
	{
		log("Can't open settings file %s", filename);
		return FAIL;
	}

	// input filename
	LOAD_STR(file, input_filename);

	// lowpass init params
	LOAD_FLOAT(file, lowpass_cut_off_frequency);
	LOAD_FLOAT(file, lowpass_delta_time);
	LOAD_STR(file, lowpass_output_filename);

	// kalman init params
	LOAD_FLOAT(file, kalman_init_covariance);
	LOAD_FLOAT(file, kalman_measurement_noise);
	LOAD_FLOAT(file, kalman_environment_noise);
	LOAD_STR(file, kalman_output_filename);

	// butter init params
	LOAD_FLOAT(file, butter_order);
	LOAD_FLOAT(file, butter_cut_off_frequency);
	LOAD_FLOAT(file, butter_sample_rate);
	LOAD_STR(file, butter_output_filename);

	// adjust init params
	LOAD_FLOAT(file, adjust_min_lerp_diff);
	LOAD_FLOAT(file, adjust_max_lerp_diff);
	LOAD_FLOAT(file, adjust_min_speed);
	LOAD_FLOAT(file, adjust_max_speed);

	fclose(file);

	return SUCCESS;
}

void settings_init_default()
{
	strcpy_s(input_filename, BUFFER_SIZE, "../full.csv");

	kalman_init_covariance = 0.1f;
	kalman_measurement_noise = 15.0f;
	kalman_environment_noise = 2.0f;
	strcpy_s(kalman_output_filename, BUFFER_SIZE, "../full_filtered_kalman");

	lowpass_cut_off_frequency = 20.0f;
	lowpass_delta_time = 0.001f;
	strcpy_s(lowpass_output_filename, BUFFER_SIZE, "../full_filtered_lowpass");

	butter_order = 6.0f;
	butter_cut_off_frequency = 50.0f;
	butter_sample_rate = 1000.0f;
	strcpy_s(butter_output_filename, BUFFER_SIZE, "../full_filtered_butter");

	adjust_min_lerp_diff = 30.0f;
	adjust_max_lerp_diff = 100.0f;
	adjust_min_speed = 0.02f;
	adjust_max_speed = 0.2f;

	log("Use default input_filename: %s", input_filename);

	log("Use default lowpass_cut_off_frequency: %f", lowpass_cut_off_frequency);
	log("Use default lowpass_delta_time: %f", lowpass_delta_time);
	log("Use default lowpass_output_filename: %s", lowpass_output_filename);

	log("Use default kalman_init_covariance: %f", kalman_init_covariance);
	log("Use default kalman_measurement_noise: %f", kalman_measurement_noise);
	log("Use default kalman_environment_noise: %f", kalman_environment_noise);
	log("Use default kalman_output_filename: %s", kalman_output_filename);

	log("Use default butter_order: %f", butter_order);
	log("Use default butter_cut_off_frequency: %f", butter_cut_off_frequency);
	log("Use default butter_sample_rate: %f", butter_sample_rate);
	log("Use default butter_output_filename: %s", butter_output_filename);

	log("Use default adjust_min_lerp_diff: %f", adjust_min_lerp_diff);
	log("Use default adjust_max_lerp_diff: %f", adjust_max_lerp_diff);
	log("Use default adjust_min_speed: %f", adjust_min_speed);
	log("Use default adjust_max_speed: %f", adjust_max_speed);
}

void setting_get_input_filename(const char** filename)
{
	*filename = input_filename;
}

void settings_get_lowpass_params(float* cut_off_frequency, float* delta_time, const char** output_filename)
{
	*cut_off_frequency = lowpass_cut_off_frequency;
	*delta_time = lowpass_delta_time;
	*output_filename = lowpass_output_filename;
}

void settings_get_kalman_params(float* init_covariance, float* measurement_noise, float* environment_noise, const char** output_filename)
{
	*init_covariance = kalman_init_covariance;
	*measurement_noise = kalman_measurement_noise;
	*environment_noise = kalman_environment_noise;
	*output_filename = kalman_output_filename;
}

void settings_get_butter_params(float* order, float* cut_off_frequency, float* sample_rate, const char** output_filename)
{
	*order = butter_order;
	*cut_off_frequency = butter_cut_off_frequency;
	*sample_rate = butter_sample_rate;
	*output_filename = butter_output_filename;
}

void settings_get_adjust_params(float* min_lerp_diff, float* max_lerp_diff, float* min_speed, float* max_speed)
{
	*min_lerp_diff = adjust_min_lerp_diff;
	*max_lerp_diff = adjust_max_lerp_diff;
	*min_speed = adjust_min_speed;
	*max_speed = adjust_max_speed;
}

