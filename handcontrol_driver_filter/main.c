#pragma once

#include "csv.h"
#include "filters.h"
#include "logger.h"
#include "settings.h"

#include <stdlib.h>
#include <string.h>

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

	const char* filters_output_names[FILTERS_COUNT];
	float(*filter_funcs[FILTERS_COUNT])(float) = {kalman_filter, lowpass_filter, butter_filter};

	// get filter params
	float kalman_init_covariance, kalman_measurement_noise, kalman_environment_noise;
	settings_get_kalman_params(&kalman_init_covariance, &kalman_measurement_noise, &kalman_environment_noise, &filters_output_names[0]);
	float lowpass_cut_off_frequency, lowpass_delta_time;
	settings_get_lowpass_params(&lowpass_cut_off_frequency, &lowpass_delta_time, &filters_output_names[1]);
	float butter_order, butter_cut_off_frequency, butter_sample_rate;
	settings_get_butter_params(&butter_order, &butter_cut_off_frequency, &butter_sample_rate, &filters_output_names[2]);
	float adjust_min_lerp_diff, adjust_max_lerp_diff, adjust_min_speed, adjust_max_speed;
	settings_get_adjust_params(&adjust_min_lerp_diff, &adjust_max_lerp_diff, &adjust_min_speed, &adjust_max_speed);


	// init filters
	kalman_init(data[0].value, kalman_init_covariance, kalman_measurement_noise, kalman_environment_noise);
	lowpass_init(data[0].value, lowpass_cut_off_frequency, lowpass_delta_time);
	butter_init(data[0].value, butter_order, butter_cut_off_frequency, butter_sample_rate);//6, 100, 1000);
	adjust_init(butter_filter, adjust_min_lerp_diff, adjust_max_lerp_diff, adjust_min_speed, adjust_max_speed);//30, 100, 0.02f, 0.2f);

	// pass all filtres
	for (int i = 0; i < FILTERS_COUNT; i++)
	{
		ADCEntry* data_filtered = malloc(sizeof(ADCEntry) * dataSize);
		ADCEntry* data_diffs = malloc(sizeof(ADCEntry) * dataSize);
		
		// select next filter
		adjust_update_func(filter_funcs[i]);

		// perform filtrarion
		for (int j = 0; j < dataSize; j++)
		{
			data_filtered[j].tick = data_diffs[j].tick = data[j].tick;
			data_filtered[j].value = adjust_filter(data[j].value);
			data_diffs[j].value = data_filtered[j].value - data[j].value;
		}

		char filename[512];

		// save filtered data
		strcpy_s(filename, 512, filters_output_names[i]);
		strcat_s(filename, 512, ".csv");
		if (write_csv(filename, dataSize, data_filtered) == FAIL)
		{
			log("Failed to save in file: %s", filename);
		}
		// save differences data
		strcpy_s(filename, 512, filters_output_names[i]);
		strcat_s(filename, 512, "_diffs.csv");
		if (write_csv(filename, dataSize, data_diffs) == FAIL)
		{
			log("Failed to save in file: %s", filename);
		}

		free(data_filtered);
		free(data_diffs);
	}
	
	free(data);
}