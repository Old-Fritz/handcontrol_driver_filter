#pragma once

#define FAIL 0
#define SUCCESS 1
#define CSV_SEPARATOR ';'

typedef struct ADCEntry
{
	float tick;
	float value;
} ADCEntry;

int read_csv(const char* filename, int* outDataSize, ADCEntry** outData);
int write_csv(const char* filename, int dataSize, ADCEntry* data);