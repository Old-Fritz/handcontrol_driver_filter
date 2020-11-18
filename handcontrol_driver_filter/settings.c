#pragma once

#include "settings.h"
#include "logger.h"

#include <stdio.h>
#include <string.h>

#define BUFFER_SIZE 512


static char  input_filename[BUFFER_SIZE] = "../full.csv";

static float kalman_init_covariance = 0.1f;
static float kalman_measurement_noise = 15;
static float kalman_environment_noise = 2;
static char  kalman_output_filename[BUFFER_SIZE] = "../full_filtered_kalman.csv";

static float lowpass_cut_off_frequency = 200;
static float lowpass_delta_time = 0.001f;
static char  lowpass_output_filename[BUFFER_SIZE] = "../full_filtered_lowpass.csv";


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

	fclose(file);

	return SUCCESS;
}

void settings_init_default()
{
	strcpy_s(input_filename, BUFFER_SIZE, "../full.csv");

	kalman_init_covariance = 0.1f;
	kalman_measurement_noise = 15;
	kalman_environment_noise = 2;
	strcpy_s(kalman_output_filename, BUFFER_SIZE, "../full_filtered_kalman.csv");

	lowpass_cut_off_frequency = 20;
	lowpass_delta_time = 0.001f;
	strcpy_s(lowpass_output_filename, BUFFER_SIZE, "../full_filtered_lowpass.csv");

	log("Use default input_filename: %s", input_filename);

	log("Use default lowpass_cut_off_frequency: %f", lowpass_cut_off_frequency);
	log("Use default lowpass_delta_time: %f", lowpass_delta_time);
	log("Use default lowpass_output_filename: %s", lowpass_output_filename);

	log("Use default kalman_init_covariance: %f", kalman_init_covariance);
	log("Use default kalman_measurement_noise: %f", kalman_measurement_noise);
	log("Use default kalman_environment_noise: %f", kalman_environment_noise);
	log("Use default kalman_output_filename: %s", kalman_output_filename);
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

