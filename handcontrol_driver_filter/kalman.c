#pragma once
#include "filters.h"


static float m_x0 = 0; // predicted state
static float m_p0 = 0.1f; // predicted covariance

static float m_f = 1; // factor of real value to previous real value
static float m_q = 2; // environment noise
static float m_h = 1; // factor of measured value to real value
static float m_r = 15; // measurement noise

static float m_state = 192;
static float m_covariance = 0.1f;


void kalman_init(float init_state, float init_covariance, float measurement_noise, float environment_noise)
{
	m_state = init_state;
	m_covariance = init_covariance;
	m_q = environment_noise;
	m_r = measurement_noise;
}
float kalman_filter(float value)
{
    //time update - prediction
    m_x0 = m_f * m_state;
    m_p0 = m_f * m_covariance * m_f + m_q;

    //measurement update - correction
    float k = m_h * m_p0 / (m_h * m_p0 * m_h + m_r);
    m_state = m_x0 + k * (value - m_h * m_x0);
    m_covariance = (1 - k * m_h) * m_p0;

    return m_state;
}