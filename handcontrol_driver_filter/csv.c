#pragma once

#include "csv.h"
#include "logger.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static int get_lines_count(const char* filename)
{
	FILE* file;
	if (fopen_s(&file, filename, "r"))
	{
		log("Failed to open file: %s", filename);
		return -1;
	}

	int lines = 0;
	while (!feof(file))
	{
		char ch = fgetc(file);
		if (ch == '\n')
		{
			lines++;
		}
	}

	fclose(file);

	return lines;
}

static void read_csv_line(FILE* file, float* left_value, float* right_value)
{
	char left_valueBuffer[256], right_valueBuffer[256];
	int left_value_buffer_size = 0;
	int right_value_buffer_size = 0;

	while (!feof(file) && left_value_buffer_size < 255)
	{
		char ch = fgetc(file);
		if (ch == CSV_SEPARATOR)
		{
			break;
		}
		left_valueBuffer[left_value_buffer_size++] = ch;
	}
	left_valueBuffer[left_value_buffer_size] = 0;


	while (!feof(file) && right_value_buffer_size < 255)
	{
		char ch = fgetc(file);
		if (ch == '\n')
		{
			break;
		}
		right_valueBuffer[right_value_buffer_size++] = ch;
	}
	right_valueBuffer[right_value_buffer_size] = 0;

	if (left_value)
	{
		*left_value = atof(left_valueBuffer);
	}
	if (right_value)
	{
		*right_value = atof(right_valueBuffer);
	}
}

static void write_csv_line(FILE* file, float left_value, float right_value)
{
	char buffer[256];

	size_t string_length;
	sprintf_s(buffer, 256, "%f", left_value);
	string_length = strlen(buffer);

	buffer[string_length++] = CSV_SEPARATOR;
	buffer[string_length] = 0;

	sprintf_s(buffer + string_length, 256 - string_length, "%f", right_value);
	string_length = strlen(buffer);

	buffer[string_length++] = '\n';
	buffer[string_length] = 0;
	fwrite(buffer, sizeof(char), string_length, file);
}

int read_csv(const char* filename, int* outDataSize, ADCEntry** outData)
{
	int lines_count = get_lines_count(filename);
	if (lines_count < 0)
	{
		log("Failed to count lines in file: %s", filename);
		return FAIL;
	}
	int dataSize = lines_count - 1;   // first line contains column names
	ADCEntry* data = malloc(dataSize * sizeof(ADCEntry));

	FILE* file;
	if (fopen_s(&file, filename, "r"))
	{
		log("Failed to open file: %s", filename);
		return FAIL;
	}

	read_csv_line(file, NULL, NULL);      // first line contains column names
	for (int i = 0; i < dataSize; i++)
	{
		read_csv_line(file, &data[i].tick, &data[i].value);
	}

	fclose(file);

	*outDataSize = dataSize;
	*outData = data;

	return SUCCESS;
}

int write_csv(const char* filename, int dataSize, ADCEntry* data)
{
	FILE* file;
	if (fopen_s(&file, filename, "w"))
	{
		log("Failed to open file: %s", filename);
		return FAIL;
	}


	const char column_names[] = "tick;adc\n";
	fwrite(column_names, sizeof(char), strlen(column_names), file);
	for (int i = 0; i < dataSize; i++)
	{
		write_csv_line(file, data[i].tick, data[i].value);
	}

	fclose(file);

	return SUCCESS;
}
