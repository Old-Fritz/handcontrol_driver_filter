#pragma once

#include "csv.h"
#include "filters.h"

int main()
{
	int dataSize;
	ADCEntry* data;
	if (read_csv("../full.csv", &dataSize, &data) == FAIL)
	{
		exit();
	}

	ADCEntry* data_kalman = malloc(sizeof(ADCEntry) * dataSize);
	ADCEntry* data_low_pass = malloc(sizeof(ADCEntry) * dataSize);

	data_kalman[0] = data_low_pass[0] = data[0]; // init first element

	// not ideal params
	kalman_init(data[0].value, 0.1f, 20, 1);
	lowpass_init(data[0].value, 10, 0.001f);

	for (int i = 1; i < dataSize; i++)
	{
		data_kalman[i].tick = data_low_pass[i].tick = data[i].tick;

		data_kalman[i].value = kalman_filter(data[i].value);
		data_low_pass[i].value = lowpass_filter(data[i].value);
	}
	
	write_csv("../full_filtered_kalman.csv", dataSize, data_kalman);
	write_csv("../full_filtered_low_pass.csv", dataSize, data_low_pass);

	free(data_kalman);
	free(data_low_pass);
	free(data);
}