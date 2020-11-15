#pragma once

#include "csv.h"
#include "filters.h"

int main()
{
	int dataSize;
	ADCEntry* data;
	if (read_csv("full.csv", &dataSize, &data) == FAIL)
	{
		exit();
	}

	kalman_init(data[0].value, 0.1f, 20, 1); // not ideal params
	for (int i = 1; i < dataSize; i++)
	{
		data[i].value = kalman_filter(data[i].value);
	}
	
	write_csv("full_filtered.csv", dataSize, data);

	free(data);
}