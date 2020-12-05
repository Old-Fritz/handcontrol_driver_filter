#pragma once
#include "filters.h"
#include <math.h>


static float(*m_filter_func)(float);
static float m_min_lerp_diff = 0.f;
static float m_max_lerp_diff = 0.f;

static float m_min_speed = 0.f;
static float m_max_speed = 0.f;

static float m_current_diff = 0.f;

void adjust_update_func(float(*filter_func)(float))
{
	m_filter_func = filter_func;
}

void adjust_init(float(*filter_func)(float), float min_lerp_diff, float max_lerp_diff, float min_speed, float max_speed)
{
	m_filter_func = filter_func;
	m_min_lerp_diff = min_lerp_diff; 
	m_max_lerp_diff = max_lerp_diff;

	m_min_speed = min_speed;
	m_max_speed = max_speed;

	m_current_diff = 0.0f;
}
float adjust_filter(float filter)
{
	float orig_value = filter;
	filter = m_filter_func(orig_value);
	filter += m_current_diff;
	float diff = orig_value - filter;
	float speed = 0.0f;
	float abs_diff = fabs(diff);

	// linear interpolation between areas
	if (abs_diff > m_max_lerp_diff)
	{
		speed = m_max_speed;
	}
	else if (abs_diff > m_min_lerp_diff)
	{
		speed = m_min_speed + (abs_diff - m_min_lerp_diff) / (m_max_lerp_diff - m_min_lerp_diff) * (m_max_speed - m_min_speed);
	}
	else
	{
		speed = abs_diff / m_min_lerp_diff * m_min_speed;
	}

	// adjust 
	m_current_diff += diff * speed;

	return filter;
}
