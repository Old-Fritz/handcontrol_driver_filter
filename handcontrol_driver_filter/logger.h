#pragma once

#include <stdio.h>

#define FAIL 0
#define SUCCESS 1

#define log(format, ...) printf_s(format, __VA_ARGS__); printf_s("\n");