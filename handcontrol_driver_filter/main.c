#pragma once

#include "csv.h"
#include "filters.h"
#include "logger.h"
#include "settings.h"

#include <stdlib.h>


int main(int argc, char* argv[])
{
	int dataSize;
	ADCEntry* data;

	// try load settings
	char* settings_filename;
	if (argc > 1)
	{
		settings_filename = argv[0];
		log("Using settings filename from arguments: %s", settings_filename);
	}
	else
	{
		settings_filename = "../settings.txt";
		log("Using default settings filename: %s", settings_filename);
	}
	int result = settings_load_from_file(settings_filename);
	if (result == FAIL)
	{
		log("Failed to load settings.", settings_filename);
		settings_init_default();
	}

	// try read input data
	const char* input_filename;
	setting_get_input_filename(&input_filename);
	if (read_csv(input_filename, &dataSize, &data) == FAIL)
	{
		log("Failed to read input csv file: %s", input_filename);
		return 1;
	}

	// init filtered data
	ADCEntry* data_kalman = malloc(sizeof(ADCEntry) * dataSize);
	ADCEntry* data_low_pass = malloc(sizeof(ADCEntry) * dataSize);

	data_kalman[0] = data_low_pass[0] = data[0]; // init first element

	// get filter params
	float kalman_init_covariance, kalman_measurement_noise, kalman_environment_noise;
	const char* kalman_output_filename;
	settings_get_kalman_params(&kalman_init_covariance, &kalman_measurement_noise, &kalman_environment_noise, &kalman_output_filename);
	float lowpass_cut_off_frequency, lowpass_delta_time;
	const char* lowpass_output_filename;
	settings_get_lowpass_params(&lowpass_cut_off_frequency, &lowpass_delta_time, &lowpass_output_filename);

	// init filters
	kalman_init(data_kalman[0].value, kalman_init_covariance, kalman_measurement_noise, kalman_environment_noise);
	lowpass_init(data_low_pass[0].value, lowpass_cut_off_frequency, lowpass_delta_time);

	// perform filtrarion
	for (int i = 1; i < dataSize; i++)
	{
		data_kalman[i].tick = data_low_pass[i].tick = data[i].tick;

		data_kalman[i].value = kalman_filter(data[i].value);
		data_low_pass[i].value = lowpass_filter(data[i].value);
	}
	
	// save data
	if (write_csv(kalman_output_filename, dataSize, data_kalman) == FAIL)
	{
		log("Failed to save in file: %s", kalman_output_filename);
	}
	if(write_csv(lowpass_output_filename, dataSize, data_low_pass) == FAIL)
	{
		log("Failed to save in file: %s", lowpass_output_filename);
	}

	// free resources
	free(data_kalman);
	free(data_low_pass);
	free(data);
}