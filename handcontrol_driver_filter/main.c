#pragma once

#include "csv.h"
#include "filters.h"
#include "logger.h"
#include "settings.h"

#include <stdlib.h>
#include <string.h>

#define TEST_DATA_SIZE 2000
#define TEST_OFFSET 100

void add_sin(ADCEntry* data, int data_size, int period, float amplitude, float offset_amplitude)
{
	for (int i = 0; i < data_size; i++)
	{
		data[i].value += sin((float)i / (float)period) * amplitude + offset_amplitude;
	}
}

void add_mush(ADCEntry* data, int data_size, float min_value, float max_value)
{
	for (int i = 0; i < data_size; i++)
	{
		data[i].value += ((rand() / (float)RAND_MAX) * (max_value - min_value) + min_value);
	}
}

void gen_test_data()
{
	const int test_data_size = TEST_DATA_SIZE;
	ADCEntry data[TEST_DATA_SIZE];
	ADCEntry mushed_data[TEST_DATA_SIZE];
	ADCEntry filtred_data[TEST_DATA_SIZE];
	ADCEntry filtred_data_diffs[TEST_DATA_SIZE];

	// gen data
	for (int i = 0; i < test_data_size; i++)
	{
		data[i].tick = i;
		data[i].value = 0;
		mushed_data[i].tick = i;
		mushed_data[i].value = 0;
	}

	add_sin(data, test_data_size, 100, 500, 500);
	add_sin(data, test_data_size, 5, 20, 0);
	add_sin(mushed_data, test_data_size, 100, 500, 500);
	add_sin(mushed_data, test_data_size, 5, 20, 0);
	add_mush(mushed_data, test_data_size, -10, 10);

	// filter data
	adjust_update_func(butter_filter);
	for (int i = 0; i < test_data_size; i++)
	{
		filtred_data[i].tick = filtred_data_diffs[i].tick = mushed_data[i].tick;
		filtred_data[i].value = adjust_filter(mushed_data[i].value);
		filtred_data_diffs[i].value = filtred_data[i].value - data[i].value;
	}

	write_csv("../csv/test_data_filtered_diffs.csv", test_data_size - TEST_OFFSET, filtred_data_diffs + TEST_OFFSET);
	write_csv("../csv/test_data_filtered.csv", test_data_size - TEST_OFFSET, filtred_data + TEST_OFFSET);
	write_csv("../csv/test_data.csv", test_data_size - TEST_OFFSET, data + TEST_OFFSET);
	write_csv("../csv/test_data_mushed.csv", test_data_size - TEST_OFFSET, mushed_data + TEST_OFFSET);
}

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

	gen_test_data();
}