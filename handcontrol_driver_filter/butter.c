#pragma once
#include "filters.h"
#define _USE_MATH_DEFINES
#include <math.h>


#define MAX_ORDER 12
#define MAX_STAGES MAX_ORDER / 2 + 1

// Filter types
typedef float Sample;

typedef struct Complex
{
    Sample real;
    Sample imag;
} Complex;

typedef struct ComplexPair
{
    Complex first;
    Complex second;
} ComplexPair;

typedef struct PoleZeroPair
{
    ComplexPair poles;
    ComplexPair zeros;
} PoleZeroPair;

typedef struct LayoutBase
{
    int num_poles;
    PoleZeroPair pairs[MAX_STAGES];
    Sample normal_w;
    Sample normal_gain;
} LayoutBase;

typedef struct Biquad
{
    Sample a0;
    Sample a1;
    Sample a2;
    Sample b1;
    Sample b2;
    Sample b0;
} Biquad;

typedef struct State
{
    Sample v1;
    Sample v2;
} State;

typedef struct Cascade
{
    int num_stages;
    Biquad stages[MAX_STAGES];
    State states[MAX_STAGES];
} Cascade;

// constants
static const Complex c_infinity = { INFINITY, INFINITY };
static const Complex c_minus_one = { -1, 0 };

// params
static Cascade m_cascade;

// complex numbers operations
static Complex complex_polar(Sample rho, Sample theta)
{
    Complex out;
    out.real = rho * cos(theta);
    out.imag = rho * sin(theta);
    return out;
}
static Complex complex_conj(Complex value)
{
    Complex out;
    out.real = value.real;
    out.imag = -value.imag;
    return out;
}
static Complex complex_zero()
{
    Complex out = { 0,0 };
    return out;
}
static Complex complex_exp(Complex value)
{
    Complex out;
    Sample exp_part = exp(value.real);
    out.real = exp_part * cos(value.imag);
    out.imag = exp_part * sin(value.imag);
    return out;
}
static Sample complex_norm(Complex value)
{
    return value.real * value.real + value.imag * value.imag;
}
static Sample complex_abs(Complex value)
{
    return sqrt(complex_norm(value));
}
static Complex complex_addmul(Complex arg1, Sample v, Complex arg2)
{
    Complex out;
    out.real = arg1.real + v * arg2.real;
    out.imag = arg1.imag + v * arg2.imag;
    return out;
}
static Complex complex_div(Complex num, Complex den)
{
    Complex output;

    Sample den_coef = den.real * den.real + den.imag * den.imag;
    output.real = (num.real * den.real + num.imag * den.imag) / den_coef;
    output.imag = (num.imag * den.real - num.real * den.imag) / den_coef;

    return output;
}
static Complex complex_mul(Complex arg1, Complex arg2)
{
    Complex output;

    output.real = arg1.real * arg2.real - arg1.imag * arg2.imag;
    output.imag = (arg1.real * arg2.imag + arg1.imag * arg2.real);

    return output;
}
static Complex complex_transform(Complex value, Sample f)
{
    if (value.real == c_infinity.real && value.imag == c_infinity.imag)
        return c_minus_one;

    // frequency transform
    value.real = f * value.real;
    value.imag = f * value.imag;

    // bilinear low pass transform
    Complex num = { 1 + value.real, value.imag };
    Complex den = { 1 - value.real, value.imag };

    return complex_div(num, den);
}

// layout operations
static void layout_add_pz_conj(LayoutBase* layout, Complex pole, Complex zero)
{
    layout->pairs[layout->num_poles / 2].poles.first = pole;
    layout->pairs[layout->num_poles / 2].poles.second = complex_conj(pole);
    layout->pairs[layout->num_poles / 2].zeros.first = zero;
    layout->pairs[layout->num_poles / 2].zeros.second = complex_conj(zero);
    layout->num_poles += 2;
}
static void layout_add_pz(LayoutBase* layout, Complex pole, Complex zero)
{
    layout->pairs[layout->num_poles / 2].poles.first = pole;
    layout->pairs[layout->num_poles / 2].poles.second = complex_zero();
    layout->pairs[layout->num_poles / 2].zeros.first = zero;
    layout->pairs[layout->num_poles / 2].zeros.second = complex_zero();
    layout->num_poles += 2;
}
static void layout_design(LayoutBase* layout, int num_poles)
{
    layout->num_poles = 0;
    layout->normal_w = 0;
    layout->normal_gain = 1;

    const Sample n2 = 2 * num_poles;
    const int pairs_count = num_poles / 2;
    for (int i = 0; i < pairs_count; ++i)
    {
        Complex c = complex_polar(1, M_PI_2 + (2 * i + 1) * M_PI / n2);
        layout_add_pz_conj(layout, c, c_infinity);
    }
    if (num_poles & 1)
        layout_add_pz(layout, c_minus_one, c_infinity);
}
static void layout_transform(Sample fc, LayoutBase* digital, LayoutBase* analog)
{
    // reset digital
    digital->num_poles = 0;
    digital->normal_w = 0;
    digital->normal_gain = 1;

    Sample f = tan(M_PI * fc);
    const int num_poles = analog->num_poles;
    const int pairs = num_poles / 2;

    for (int i = 0; i < pairs; ++i)
    {
        PoleZeroPair pair = analog->pairs[i];
        layout_add_pz_conj(digital,
            complex_transform(pair.poles.first, f),
            complex_transform(pair.zeros.first, f));
    }

    if (num_poles & 1)
    {
        PoleZeroPair pair = analog->pairs[pairs];
        layout_add_pz(digital,
            complex_transform(pair.poles.first, f),
            complex_transform(pair.zeros.first, f));
    }

    digital->normal_w = analog->normal_w;
    digital->normal_gain = analog->normal_gain;
}

// biquad operations
static void bq_set_coefficients(Biquad* bq, Sample a0, Sample a1, Sample a2,
    Sample b0, Sample b1, Sample b2)
{
    bq->a0 = a0;
    bq->a1 = a1 / a0;
    bq->a2 = a2 / a0;
    bq->b0 = b0 / a0;
    bq->b1 = b1 / a0;
    bq->b2 = b2 / a0;
}
static void bq_init(Biquad* bq)
{
    bq_set_coefficients(bq, 1, 0, 0, 1, 0, 0);
}
static void bq_set_one_pole(Biquad* bq, Complex pole, Complex zero)
{
    const Sample a0 = 1;
    const Sample a1 = -pole.real;
    const Sample a2 = 0;
    const Sample b0 = -zero.real;
    const Sample b1 = 1;
    const Sample b2 = 0;

    bq_set_coefficients(bq, 0, a1, a2, b0, b1, b2);
}
static void bq_set_two_pole(Biquad* bq, Complex pole1, Complex zero1,
    Complex pole2, Complex zero2)
{
    const Sample a0 = 1;
    Sample a1;
    Sample a2;

    if (pole1.imag != 0)
    {
        a1 = -2 * pole1.real;
        a2 = complex_norm(pole1);
    }
    else
    {
        a1 = -(pole1.real + pole2.real);
        a2 = pole1.real * pole2.real;
    }

    const Sample b0 = 1;
    Sample b1;
    Sample b2;

    if (zero1.imag != 0)
    {
        b1 = -2 * zero1.real;
        b2 = complex_norm(zero1);
    }
    else
    {
        b1 = -(zero1.real + zero2.real);
        b2 = zero1.real * zero2.real;
    }

    bq_set_coefficients(bq, a0, a1, a2, b0, b1, b2);
}
static void bq_set_pole_zero_pair(Biquad* bq, PoleZeroPair pair)
{
    if (pair.poles.second.real == 0 &&
        pair.poles.second.imag == 0 &&
        pair.zeros.second.real == 0 &&
        pair.zeros.second.imag == 0)
    {
        bq_set_one_pole(bq, pair.poles.first, pair.zeros.first);
    }
    else
    {
        bq_set_two_pole(bq, pair.poles.first, pair.zeros.first,
            pair.poles.second, pair.zeros.second);
    }
}

// cascade operations
static Complex cascade_response(Cascade* cascade, Sample norm_freq)
{
    Sample w = 2 * M_PI * norm_freq;
    const Complex czn1 = complex_polar(1., -w);
    const Complex czn2 = complex_polar(1., -2 * w);
    Complex ch = { 1, 0 };
    Complex cbot = { 1, 0 };

    const Biquad* stage = cascade->stages;
    for (int i = cascade->num_stages; --i >= 0; ++stage)
    {
        Complex cb = { 1, 0 };
        Complex ct = { stage->b0 / stage->a0, 0 };
        ct = complex_addmul(ct, stage->b1 / stage->a0, czn1);
        ct = complex_addmul(ct, stage->b2 / stage->a0, czn2);
        cb = complex_addmul(cb, stage->a1 / stage->a0, czn1);
        cb = complex_addmul(cb, stage->a2 / stage->a0, czn2);
        ch = complex_mul(ch, ct);
        cbot = complex_mul(cbot, cb);
    }

    return complex_div(ch, cbot);
}
static void cascade_set_layout(Cascade* cascade, LayoutBase* layout)
{
    const int num_poles = layout->num_poles;
    cascade->num_stages = (num_poles + 1) / 2;

    Biquad* stage = cascade->stages;
    State* state = cascade->states;
    for (int i = 0; i < MAX_STAGES; ++i, ++stage, ++state)
    {
        bq_init(stage);
        state->v1 = 0;
        state->v2 = 0;
    }

    stage = cascade->stages;
    for (int i = 0; i < cascade->num_stages; ++i, ++stage)
        bq_set_pole_zero_pair(stage, layout->pairs[i]);

    Sample scale = layout->normal_gain /
        complex_abs(cascade_response(cascade, layout->normal_w / (2 * M_PI)));

    cascade->stages->b0 *= scale;
    cascade->stages->b1 *= scale;
    cascade->stages->b2 *= scale;
}
static Sample state_filter(State* state, const Sample value, const Biquad* stage)
{
    Sample w = value - stage->a1 * state->v1 - stage->a2 * state->v2;
    Sample out = stage->b0 * w + stage->b1 * state->v1 + stage->b2 * state->v2;

    state->v2 = state->v1;
    state->v1 = w;

    return out;
}

// butter filter specific
static void butter_init_output(float init_value, int deep)
{
    for (int i = 0; i < deep; i++)
    {
        butter_filter(init_value);
    }
}

float butter_filter(float value)
{
    Sample out = value;
    State* state = m_cascade.states;
    Biquad const* stage = m_cascade.stages;
    for (int i = m_cascade.num_stages; --i >= 0; ++state, ++stage)
    {
        out = state_filter(state, out, stage);
    }
    return out;
}


void butter_init(float init_output, float order, float cut_off_frequency, float sample_rate)
{
    LayoutBase analog_layout, digital_layout;

    layout_design(&analog_layout, order);
    layout_transform(cut_off_frequency / sample_rate * 2, &digital_layout, &analog_layout);
    cascade_set_layout(&m_cascade, &digital_layout);

    // init state to match first input
    butter_init_output(init_output, 100);
}


/*
typedef struct zpk 
{
    float z, p, k;
} zpk;

typedef struct ZerosPoles
{

} ZerosPoles;

float relative_degree(z, p)
{
    //Return relative degree of transfer function from zeros and poles
    float degree = len(p) - len(z);
    if degree < 0 :
        raise ValueError("Improper transfer function. "
            "Must have at least as many poles as zeros.")
    else:

    return degree;
}

zpk lp2hp_zpk(float zeros, float poles, float k, float wo)
{
    zpk output;

    z = atleast_1d(z);
    p = atleast_1d(p);

    float degree = relative_degree(z, p);

        // Invert positions radially about unit circle to convert LPF to HPF
        // Scale all points radially from origin to shift cutoff frequency
    output.z = wo / z;
    output.p = wo / p;

    //# If lowpass had zeros at infinity, inverting moves them to origin.
    output.z = append(output.z, zeros(degree));

    //# Cancel out gain change caused by inversion
    output.k = k * real(prod(-z) / prod(-p));

    return output;
}
*/

/*
#define FILTER_SECTIONS   2

typedef struct  {
    unsigned int length;       // size of filter 
    float history[2 * FILTER_SECTIONS];    // history in filter 
    float coef [4 * FILTER_SECTIONS];               // pointer to coefficients of filter
} filter;

typedef struct {
    double a0, a1, a2;       // numerator coefficients 
    double b0, b1, b2;       // denominator coefficients 
} buquad;

buquad m_proto_coef[FILTER_SECTIONS];      // Filter prototype coefficients, 1 for each filter section


filter   m_iir;


void prewarp(double* a0, double* a1, double* a2, double fc, double fs);
void bilinear(
    double a0, double a1, double a2,
    double b0, double b1, double b2,
    double* k,
    double fs,
    float* coef);


void prewarp(
    double* a0, double* a1, double* a2,
    double fc, double fs)
{
    double wp, pi;

    pi = 4.0 * atan(1.0);
    wp = 2.0 * fs * tan(pi * fc / fs);

    *a2 = (*a2) / (wp * wp);
    *a1 = (*a1) / wp;
}


void bilinear(
    double a0, double a1, double a2,
    double b0, double b1, double b2,
    double* k,
    double fs,
    float* coef
)
{
    double ad, bd;

    ad = 4. * a2 * fs * fs + 2. * a1 * fs + a0;
    bd = 4. * b2 * fs * fs + 2. * b1 * fs + b0;

    *k *= ad / bd;

    *coef++ = (float)((2. * b0 - 8. * b2 * fs * fs) / bd); // beta1
    *coef++ = (float)((4. * b2 * fs * fs - 2. * b1 * fs + b0) / bd); // beta2 

    *coef++ = (float)((2. * a0 - 8. * a2 * fs * fs) / ad);			 // alpha1 
    *coef = (float)((4. * a2 * fs * fs - 2. * a1 * fs + a0) / ad); // alpha2 
}


void szxform(
    double* a0, double* a1, double* a2, // numerator coefficients 
    double* b0, double* b1, double* b2, // denominator coefficients 
    double fc,         // Filter cutoff frequency 
    double fs,         // sampling rate 
    double* k,         // overall gain factor 
    float* coef)         // pointer to 4 iir coefficients 
{
    // Calculate a1 and a2 and overwrite the original values 
    prewarp(a0, a1, a2, fc, fs);
    prewarp(b0, b1, b2, fc, fs);
    bilinear(*a0, *a1, *a2, *b0, *b1, *b2, k, fs, coef);
}



float butter_filter(float input)
{
    unsigned int i;
    float* hist1_ptr, * hist2_ptr, * coef_ptr;
    float output, new_hist, history1, history2;
    static float dc = (float)1E-25;
    input += dc;
    dc = -dc;

    coef_ptr = m_iir.coef;               

    hist1_ptr = m_iir.history;           
    hist2_ptr = hist1_ptr + 1;           

    output = input * (*coef_ptr++);

    for (i = 0; i < m_iir.length; i++)
    {
        history1 = *hist1_ptr;           
        history2 = *hist2_ptr;

        output = output - history1 * coef_ptr[0];
        new_hist = output - history2 * coef_ptr[1]; 

        output = new_hist + history1 * coef_ptr[2];
        output = output + history2 * coef_ptr[3];   

        coef_ptr += 4;
        *hist2_ptr++ = *hist1_ptr;
        *hist1_ptr++ = new_hist;
        hist1_ptr++;
        hist2_ptr++;
    }

    return(output);
}

static void init_coeffs(float init_value, int deep)
{
    for (int i = 0; i < deep; i++)
    {
        butter_filter(init_value);
    }
}

void update_filter(float resonance, float cutoff, int samplerate)
{
    unsigned nInd;
    double   a0, a1, a2, b0, b1, b2;
    double   fs;    // Sampling frequency, cutoff frequency 
    double   k;          // overall gain factor 
    float* coef;

    k = 1.0;				// Set overall filter gain 
    coef = m_iir.coef + 1;		// Skip k, or gain 
    fs = (double)samplerate;    // Sampling frequency (Hz) 


    for (nInd = 0; nInd < m_iir.length; nInd++)
    {
        a0 = m_proto_coef[nInd].a0;
        a1 = m_proto_coef[nInd].a1;
        a2 = m_proto_coef[nInd].a2;

        b0 = m_proto_coef[nInd].b0;
        b1 = m_proto_coef[nInd].b1 / resonance;      // Divide by resonance or Q
        b2 = m_proto_coef[nInd].b2;
        szxform(&a0, &a1, &a2, &b0, &b1, &b2, cutoff, fs, &k, coef);
        coef += 4;                       // Point to next filter section 
    }

    m_iir.coef[0] = (float)k;
}


void butter_init(float init_output, float order, float cut_off_frequency, float sample_rate)
{
    m_proto_coef[0].a0 = 1.0;
    m_proto_coef[0].a1 = 0;
    m_proto_coef[0].a2 = 0;
    m_proto_coef[0].b0 = 1.0;
    m_proto_coef[0].b1 = 0.765367;
    m_proto_coef[0].b2 = 1.0;

    m_proto_coef[1].a0 = 1.0;
    m_proto_coef[1].a1 = 0;
    m_proto_coef[1].a2 = 0;
    m_proto_coef[1].b0 = 1.0;
    m_proto_coef[1].b1 = 1.847759;
    m_proto_coef[1].b2 = 1.0;

    m_iir.length = FILTER_SECTIONS;         // Number of filter sections 

    if (!m_iir.coef)
    {
        return;
    }

    update_filter(1, 100, 1000);

    init_coeffs(init_output, 100000);
    return 1;
}
*/

/*

//http://www.exstrom.com/journal/sigproc/

#define MAX_ORDER 12

int m_order = 6;

double m_A[MAX_ORDER];
double m_d1[MAX_ORDER];
double m_d2[MAX_ORDER];
double m_w0[MAX_ORDER];
double m_w1[MAX_ORDER];
double m_w2[MAX_ORDER];
double m_ep;



float butter_filter(float value)
{
    for (int i = 0; i < m_order; ++i)
    {
        m_w0[i] = m_d1[i] * m_w1[i] + m_d2[i] * m_w2[i] + value;
        value = m_A[i] * (m_w0[i] + 2.0 * m_w1[i] );
        m_w2[i] = m_w1[i];
        m_w1[i] = m_w0[i];
    }
    return m_ep * value;
}

static void init_coeffs(float init_value, int deep)
{
    for (int i = 0; i < deep; i++)
    {
        butter_filter(init_value);
    }
}

static void butter(float init_output, float order, float cut_off_frequency, float sample_rate)
{
    m_order = 6;
    float s = 1000.0f;
    float f = 50.0f;
    float a = tan(M_PI * f / s);
    float a2 = a * a;

    for (int i = 0; i < m_order; ++i) {
        float r = sin(M_PI * (2.0 * i + 1.0) / (4.0 * m_order));
        s = a2 + 2.0 * a * r + 1.0;
        m_A[i] = a2 / s;
        m_d1[i] = 2.0 * (1 - a2) / s;
        m_d2[i] = -(a2 - 2.0 * a * r + 1.0) / s;
    }

    init_coeffs(init_output, 300);
}

static void cheb(float init_output, float order, float cut_off_frequency, float sample_rate)
{
    m_order = 3;
    double ep = 0.6;
    double s = 1000.0f;
    double f = 20.0f;
    double a = tan(M_PI * f / s);
    double a2 = a * a;
    double u = log((1.0 + sqrt(1.0 + ep * ep)) / ep);
    double su = sinh(u / (double)m_order * 2.0);
    double cu = cosh(u / (double)m_order * 2.0);
    double b, c;

    for (int i = 0; i < m_order; ++i) {
        b = sin(M_PI * (2.0 * i + 1.0) / (4.0 * m_order)) * su;
        c = cos(M_PI * (2.0 * i + 1.0) / (4.0 * m_order)) * cu;
        c = b * b + c * c;
        s = a2 * c + 2.0 * a * b + 1.0;
        m_A[i] = a2 / (4.0 * s); // 4.0
        m_d1[i] = 2.0 * (1 - a2 * c) / s;
        m_d2[i] = -(a2 * c - 2.0 * a * b + 1.0) / s;
    }

    m_ep = 2.0 / ep;

    init_coeffs(init_output, 300);
}

void butter_init(float init_output, float order, float cut_off_frequency, float sample_rate)
{
    cheb(init_output, order, cut_off_frequency, sample_rate);
}



*/


/*

//https://github.com/nxsEdson/Butterworth-Filter/blob/master/butterworth.cpp
float coeff_a[MAX_ORDER * 2 + 1];
float coeff_b[MAX_ORDER * 2 + 1];
float zi[MAX_ORDER * 2 + 1];
int a_count;
int b_count;

void trinomial_multiply(int order, double* b, double* c, double* output)
{
    int i, j;

    output[2] = c[0];
    output[3] = c[1];
    output[0] = b[0];
    output[1] = b[1];

    for (i = 1; i < order; ++i)
    {
        output[2 * (2 * i + 1)] += c[2 * i] * output[2 * (2 * i - 1)] - c[2 * i + 1] * output[2 * (2 * i - 1) + 1];
        output[2 * (2 * i + 1) + 1] += c[2 * i] * output[2 * (2 * i - 1) + 1] + c[2 * i + 1] * output[2 * (2 * i - 1)];

        for (j = 2 * i; j > 1; --j)
        {
            output[2 * j] += b[2 * i] * output[2 * (j - 1)] - b[2 * i + 1] * output[2 * (j - 1) + 1] +
                c[2 * i] * output[2 * (j - 2)] - c[2 * i + 1] * output[2 * (j - 2) + 1];
            output[2 * j + 1] += b[2 * i] * output[2 * (j - 1) + 1] + b[2 * i + 1] * output[2 * (j - 1)] +
                c[2 * i] * output[2 * (j - 2) + 1] + c[2 * i + 1] * output[2 * (j - 2)];
        }

        output[2] += b[2 * i] * output[0] - b[2 * i + 1] * output[1] + c[2 * i];
        output[3] += b[2 * i] * output[1] + b[2 * i + 1] * output[0] + c[2 * i + 1];
        output[0] += b[2 * i];
        output[1] += b[2 * i + 1];
    }
}

void compute_den_coeffs(int filter_order, float cut_off_freq_low, float cut_off_freq_up)
{
    int k;            // loop variables
    double theta;     // PI * (Ucutoff - Lcutoff) / 2.0
    double cp;        // cosine of phi
    double st;        // sine of theta
    double ct;        // cosine of theta
    double s2t;       // sine of 2*theta
    double c2t;       // cosine 0f 2*theta
    double r_coeffs[MAX_ORDER * 2];  // z^-2 coefficients 
    double t_coeffs[MAX_ORDER * 2];  // z^-1 coefficients
    double pole_angle;      // pole angle
    double sin_pole_angle;     // sine of pole angle
    double cos_pole_angle;     // cosine of pole angle
    double a;         // workspace variables

    cp = cos(M_PI * (cut_off_freq_up + cut_off_freq_low) / 2.0);
    theta = M_PI * (cut_off_freq_up - cut_off_freq_low) / 2.0;
    st = sin(theta);
    ct = cos(theta);
    s2t = 2.0 * st * ct;        // sine of 2*theta
    c2t = 2.0 * ct * ct - 1.0;  // cosine of 2*theta

    for (k = 0; k < filter_order; ++k)
    {
        pole_angle = M_PI * (double)(2 * k + 1) / (double)(2 * filter_order);
        sin_pole_angle = sin(pole_angle);
        cos_pole_angle = cos(pole_angle);
        a = 1.0 + s2t * sin_pole_angle;
        r_coeffs[2 * k] = c2t / a;
        r_coeffs[2 * k + 1] = s2t * cos_pole_angle / a;
        t_coeffs[2 * k] = -2.0 * cp * (ct + st * sin_pole_angle) / a;
        t_coeffs[2 * k + 1] = -2.0 * cp * st * cos_pole_angle / a;
    }

    double mul_result[MAX_ORDER * 4];
    for (int i = 0; i < MAX_ORDER * 4; i++)
    {
        mul_result[i] = 0;
    }
    trinomial_multiply(filter_order, t_coeffs, r_coeffs, mul_result);

    mul_result[1] = mul_result[0];
    mul_result[0] = 1.0;
    for (k = 0; k < 3; ++k)
    {
        coeff_a[k] = mul_result[k];
    }
    for (k = 3; k <= 2 * filter_order; ++k)
    {
        coeff_a[k] = mul_result[2 * k - 2];
    }

}

// output size - max_order + 1
void compute_lp(int filter_order, double* output)
{
    int m;
    int i;

    output[0] = 1;
    output[1] = filter_order;
    m = filter_order / 2;
    for (i = 2; i <= m; ++i)
    {
        output[i] = (double)(filter_order - i + 1) * output[i - 1] / i;
        output[filter_order - i] = output[i];
    }
    output[filter_order - 1] = filter_order;
    output[filter_order] = 1;
}

// output size - max_order + 1
void compute_hp(int filter_order, double* output)
{
    int i;

    compute_lp(filter_order, output);

    for (i = 0; i <= filter_order; ++i)
        if (i % 2) output[i] = -output[i];

}

void compute_num_coeffs(int filter_order, float cut_off_freq_low, float cut_off_freq_up)
{
    double t_coeffs[MAX_ORDER + 1];
    double num_coeffs[MAX_ORDER * 2 + 2];
    Complex normalized_kernel[MAX_ORDER * 2 + 2];

    float numbers[MAX_ORDER * 2 + 2];
    for (int i = 0; i < filter_order * 2 + 1; i++)
        numbers[i] = (double)i;

    compute_hp(filter_order, t_coeffs);

    for (int i = 0; i < filter_order; ++i)
    {
        num_coeffs[2 * i] = t_coeffs[i];
        num_coeffs[2 * i + 1] = 0.0;
    }
    num_coeffs[2 * filter_order] = t_coeffs[filter_order];

    double cp[2];
    double bw, wn;
    cp[0] = 2 * 2.0 * tan(M_PI * cut_off_freq_low / 2.0);
    cp[1] = 2 * 2.0 * tan(M_PI * cut_off_freq_up / 2.0);

    bw = cp[1] - cp[0];
    //center frequency
    wn = sqrt(cp[0] * cp[1]);
    wn = 2 * atan2(wn, 4);
    double kern;

    for (int k = 0; k < filter_order * 2 + 1; k++)
    {
        const Complex result = { 0, wn * numbers[k] };
        normalized_kernel[k] = complex_exp(result);
    }
    double b = 0;
    double den = 0;
    for (int d = 0; d < filter_order * 2 + 1; d++)
    {
        b += normalized_kernel[d].real * num_coeffs[d];
        den += normalized_kernel[d].real * coeff_a[d];
    }
    for (int c = 0; c < filter_order * 2 + 1; c++)
    {
        coeff_b[c] = (num_coeffs[c] * den) / b;
    }
}



void butter_init(float init_output, float order, float cut_off_frequency, float sample_rate)
{
    double fps = 15;
    cut_off_frequency = 1.5 / fps * 2;
    sample_rate = 2.5 / fps * 2;
    compute_den_coeffs(order, cut_off_frequency, sample_rate);
    compute_num_coeffs(order, cut_off_frequency, sample_rate);
    a_count = b_count = order * 2 + 1;
    //a_count = 1;
    for (int i = 0; i < order + 1; i++)
    {
        zi[i] = 0;
    }
}

float butter_filter(float value)
{
    float result;

    if (a_count == 1)
    {
        result = coeff_b[0] * value + zi[0];
        for (int i = 1; i < b_count; i++)
        {
            zi[i - 1] = coeff_b[i] * value + zi[i];//-coeff_a[i]*filter_x[m];
        }
    }
    else
    {
        result = coeff_b[0] * value + zi[0];
        for (int i = 1; i < b_count; i++)
        {
            zi[i - 1] = coeff_b[i] * value + zi[i] - coeff_a[i] * result;
        }
    }

    return result;
}

*/