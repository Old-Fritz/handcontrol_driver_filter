#include "filters.h"
#define _USE_MATH_DEFINES
#include <math.h>

static float m_output = 0;
static float m_e_pow = 0.118085f;

void lowpass_init(float init_output, float cut_off_frequency, float delta_time)
{
	m_output = init_output;
	m_e_pow = 1 - exp(-delta_time * 2 * M_PI * cut_off_frequency);
}
float lowpass_filter(float value)
{
	m_output += (value - m_output) * m_e_pow;
	return m_output;
}