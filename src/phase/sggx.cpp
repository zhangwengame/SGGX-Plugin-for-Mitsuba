/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/frame.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>

/// Generate a few statistics related to the implementation?
// #define SGGX_STATISTICS 1


/// The following file implements the micro-flake distribution
/// for rough fibers
#include "sggx_fiber.h"

MTS_NAMESPACE_BEGIN
//#define SGGXDEBUGEVAL 1
//#define SGGXDEBUGSAMPLE 1
//#define SGGXDEBUG 1
#if defined(SGGX_STATISTICS)
static StatsCounter avgSampleIterations("Sggx model",
		"Average rejection sampling iterations", EAverage);
#endif

/*!\plugin{microflake}{Micro-flake phase function}
 * \parameters{
 *     \parameter{stddev}{\Float}{
 *       Standard deviation of the micro-flake normals. This
 *       specifies the roughness of the fibers in the medium.
 *     }
 * }
 *
 * \renderings{
 *     \rendering{\code{stddev}=0.2}{phase_microflakes_02}
 *     \rendering{\code{stddev}=0.05}{phase_microflakes_005}
 * }
 *
 * This plugin implements the anisotropic micro-flake phase function
 * described in
 * ``A radiative transfer framework for rendering materials with
 * anisotropic structure'' by Wenzel Jakob, Adam Arbree,
 * Jonathan T. Moon, Kavita Bala, and Steve Marschner
 * \cite{Jakob2010Radiative}.
 *
 * The implementation in this plugin is specific to rough fibers
 * and uses a Gaussian-type flake distribution. It is much faster
 * than the spherical harmonics approach proposed in the original
 * paper. This distribution, as well as the implemented sampling
 * method, are described in the paper
 * ``Building Volumetric Appearance Models of Fabric using
 * Micro CT Imaging'' by Shuang Zhao, Wenzel Jakob, Steve Marschner,
 * and Kavita Bala \cite{Zhao2011Building}.
 *
 * Note: this phase function must be used with a medium that specifies
 * the local fiber orientation at different points in space. Please
 * refer to \pluginref{heterogeneous} for details.
 *
 * \vspace{1cm}
 *
 * \renderings{
 *   \fbox{\includegraphics[width=\textwidth]{images/medium_knitpatterns}}
 *   \caption{A range of different knit patterns, rendered using the
 *   \pluginref{heterogeneous} and \pluginref{microflake} plugins. Courtesy
 *   of Yuksel et al. \cite{Yuksel2012Stitch}.}
 * }
 */
void buildOrthonormalBasis(vec3& omega_1, vec3& omega_2, const vec3& omega_3)
{
	if (omega_3.z < -0.9999999f)
	{
		omega_1 = vec3(0.0f, -1.0f, 0.0f);
		omega_2 = vec3(-1.0f, 0.0f, 0.0f);
	}
	else {
		const float a = 1.0f / (1.0f + omega_3.z);
		const float b = -omega_3.x*omega_3.y*a;
		omega_1 = vec3(1.0f - omega_3.x*omega_3.x*a, b, -omega_3.x);
		omega_2 = vec3(b, 1.0f - omega_3.y*omega_3.y*a, -omega_3.y);
	}
}
float D(vec3 wm,
	float S_xx, float S_yy, float S_zz,
	float S_xy, float S_xz, float S_yz)
{
	const float detS =
		S_xx*S_yy*S_zz - S_xx*S_yz*S_yz - S_yy*S_xz*S_xz - S_zz*S_xy*S_xy + 2.0f*S_xy*S_xz*S_yz;
	const float den = wm.x*wm.x*(S_yy*S_zz - S_yz*S_yz) + wm.y*wm.y*(S_xx*S_zz - S_xz*S_xz) + wm.z*wm.z*(S_xx*S_yy - S_xy*S_xy)
		+ 2.0f*(wm.x*wm.y*(S_xz*S_yz - S_zz*S_xy) + wm.x*wm.z*(S_xy*S_yz - S_yy*S_xz) + wm.y*wm.z*(S_xy*S_xz - S_xx*S_yz));
	const float D = powf(fabsf(detS), 1.50f) / (M_PI*den*den);
	return D;
}
vec3 sample_VNDF(const vec3 wi,
	const float S_xx, const float S_yy, const float S_zz,
	const float S_xy, const float S_xz, const float S_yz,
	const float U1, const float U2)
{
	// generate sample (u, v, w)
	const float r = sqrtf(U1);
	const float phi = 2.0f*M_PI*U2;
	const float u = r*cosf(phi);
	const float v = r*sinf(phi);
	const float w = sqrtf(1.0f - u*u - v*v);
	// build orthonormal basis
	vec3 wk, wj;
	buildOrthonormalBasis(wk, wj, wi);
	// project S in this basis
	const float S_kk = wk.x*wk.x*S_xx + wk.y*wk.y*S_yy + wk.z*wk.z*S_zz
		+ 2.0f * (wk.x*wk.y*S_xy + wk.x*wk.z*S_xz + wk.y*wk.z*S_yz);
	const float S_jj = wj.x*wj.x*S_xx + wj.y*wj.y*S_yy + wj.z*wj.z*S_zz
		+ 2.0f * (wj.x*wj.y*S_xy + wj.x*wj.z*S_xz + wj.y*wj.z*S_yz);
	const float S_ii = wi.x*wi.x*S_xx + wi.y*wi.y*S_yy + wi.z*wi.z*S_zz
		+ 2.0f * (wi.x*wi.y*S_xy + wi.x*wi.z*S_xz + wi.y*wi.z*S_yz);
	const float S_kj = wk.x*wj.x*S_xx + wk.y*wj.y*S_yy + wk.z*wj.z*S_zz
		+ (wk.x*wj.y + wk.y*wj.x)*S_xy
		+ (wk.x*wj.z + wk.z*wj.x)*S_xz
		+ (wk.y*wj.z + wk.z*wj.y)*S_yz;
	const float S_ki = wk.x*wi.x*S_xx + wk.y*wi.y*S_yy + wk.z*wi.z*S_zz
		+ (wk.x*wi.y + wk.y*wi.x)*S_xy + (wk.x*wi.z + wk.z*wi.x)*S_xz + (wk.y*wi.z + wk.z*wi.y)*S_yz;
	const float S_ji = wj.x*wi.x*S_xx + wj.y*wi.y*S_yy + wj.z*wi.z*S_zz
		+ (wj.x*wi.y + wj.y*wi.x)*S_xy
		+ (wj.x*wi.z + wj.z*wi.x)*S_xz
		+ (wj.y*wi.z + wj.z*wi.y)*S_yz;
	// compute normal
	float sqrtDetSkji = sqrtf(fabsf(S_kk*S_jj*S_ii - S_kj*S_kj*S_ii - S_ki*S_ki*S_jj - S_ji*S_ji*S_kk + 2.0f*S_kj*S_ki*S_ji));
	float inv_sqrtS_ii = 1.0f / sqrtf(S_ii);
	float tmp = sqrtf(S_jj*S_ii - S_ji*S_ji);
	vec3 Mk(sqrtDetSkji / tmp, 0.0f, 0.0f);
	vec3 Mj(-inv_sqrtS_ii*(S_ki*S_ji - S_kj*S_ii) / tmp, inv_sqrtS_ii*tmp, 0);
	vec3 Mi(inv_sqrtS_ii*S_ki, inv_sqrtS_ii*S_ji, inv_sqrtS_ii*S_ii);
	vec3 wm_kji = normalize(u*Mk + v*Mj + w*Mi);
	// rotate back to world basis
	return wm_kji.x * wk + wm_kji.y * wj + wm_kji.z * wi;
}
float sigma(vec3 wi,
	float S_xx, float S_yy, float S_zz,
	float S_xy, float S_xz, float S_yz)
{
	const float sigma_squared = wi.x*wi.x*S_xx + wi.y*wi.y*S_yy + wi.z*wi.z*S_zz
		+ 2.0f * (wi.x*wi.y*S_xy + wi.x*wi.z*S_xz + wi.y*wi.z*S_yz);
	return (sigma_squared > 0.0f) ? sqrtf(sigma_squared) : 0.0f; // conditional to avoid numerical errors
}
float eval_diffuse(vec3 wi, vec3 wo,
	const float S_xx, const float S_yy, const float S_zz,
	const float S_xy, const float S_xz, const float S_yz,
	const float U1, const float U2)
{
	// sample VNDF
	const vec3 wm = sample_VNDF(wi, S_xx, S_yy, S_zz, S_xy, S_xz, S_yz, U1, U2);
	// eval diffuse
	return 1.0f / M_PI * std::max(0.0f, dot(wo, wm));
}
float eval_specular(vec3 wi, vec3 wo,
	const float S_xx, const float S_yy, const float S_zz,
	const float S_xy, const float S_xz, const float S_yz)
{
	vec3 wh = normalize(wi + wo);
	return 0.25f * D(wh, S_xx, S_yy, S_zz, S_xy, S_xz, S_yz) / sigma(wi, S_xx, S_yy, S_zz, S_xy, S_xz, S_yz);
}
vec3 sample_specular(const vec3 wi,
	const float S_xx, const float S_yy, const float S_zz,
	const float S_xy, const float S_xz, const float S_yz,
	const float U1, const float U2)
{
	// sample VNDF
	const vec3 wm = sample_VNDF(wi, S_xx, S_yy, S_zz, S_xy, S_xz, S_yz, U1, U2);
	// specular reflection
	const vec3 wo = -1.0*wi + 2.0f * dot(wm, wi)* wm;
	return wo;
}
class SggxPhaseFunction : public PhaseFunction {
public:
	SggxPhaseFunction(const Properties &props) : PhaseFunction(props) {
		/// Standard deviation of the flake distribution
		m_fiberDistr = GaussianFiberDistribution(props.getFloat("stddev"));
		sigma3 = OsigmaDir(std::cos(0));
		sigma1 = sigma2 = OsigmaDir(0.0);
		Float rate = sigma1 / 1.0;
		sigma1 /= rate;
		sigma2 /= rate;
		sigma3 /= rate;
#if defined(SGGXDEBUG)
		if (count == 0){
			FILE *f = fopen("log_phase.txt", "w");
			fprintf_s(f, "original sigma %f %f %f\n", sigma1*rate, sigma2*rate, sigma3*rate);
			fprintf_s(f, "scaled sigma %f %f %f rate %f\n", sigma1, sigma2, sigma3,rate);
			count++;
			Float Sxx, Sxy, Sxz, Syy, Syz, Szz;
			vec3 omega1, omega2, omega3(0.0,0.0,1.0);
			Sxx = sigma3*sigma3*omega3.x*omega3.x + omega3.y*omega3.y + omega3.z*omega3.z;
			Sxy = sigma3*sigma3*omega3.x*omega3.y - omega3.x*omega3.y;
			Sxz = sigma3*sigma3*omega3.x*omega3.z - omega3.x*omega3.z;
			Syy = sigma3*sigma3*omega3.y*omega3.y + omega3.x*omega3.x + omega3.z*omega3.z;
			Syz = sigma3*sigma3*omega3.y*omega3.z - omega3.y*omega3.z;
			Szz = sigma3*sigma3*omega3.z*omega3.z + omega3.x*omega3.x + omega3.y*omega3.y;
			for (int j = 0; j <= 360; j++)
			{
				for (int i = 0; i <= 180; i++)
				{
					Float radi = 2 * M_PI*(1.0*i / 360);
					Float radj = 2 * M_PI*(1.0*i / 360);
					fprintf_s(f, "%f ", D(vec3(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(radi)), Sxx, Syy, Szz, Sxy, Sxz, Syz));
					//fprintf_s(f, "%d %f;\n", i, D(vec3(sinf(radi)*sinf(radj), sinf(radi)*cosf(radj), cosf(rad)), Sxx, Syy, Szz, Sxy, Sxz, Syz));
				}
				fprintf_s(f, ";\n");
			}
				
		//	omega3 = vec3(0.0, 0.0, 1.0);
		//	buildOrthonormalBasis(omega1, omega2, omega3);
			
			
			fclose(f);

		}
#endif
	}

	SggxPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		m_fiberDistr = GaussianFiberDistribution(stream->readFloat());
		sigma1 = stream->readFloat();
		sigma2 = stream->readFloat();
		sigma3 = stream->readFloat();
		configure();
	}

	virtual ~SggxPhaseFunction() { }

	void configure() {
		PhaseFunction::configure();
		m_type = EAnisotropic | ENonSymmetric;
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		stream->writeFloat(m_fiberDistr.getStdDev());
		stream->writeFloat(sigma1);
		stream->writeFloat(sigma2);
		stream->writeFloat(sigma3);
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		if (pRec.mRec.orientation.isZero()) {
			/* What to do when the local orientation is undefined */
			#if 0
				return 1.0f / (4 * M_PI);
			#else
				return 0.0f;
			#endif
		}
		Float Sxx, Sxy, Sxz, Syy, Syz, Szz;
		vec3 omega3(pRec.mRec.orientation.x, pRec.mRec.orientation.y, pRec.mRec.orientation.z);
		omega3 = normalize(omega3);
		Sxx = sigma3*sigma3*omega3.x*omega3.x + omega3.y*omega3.y + omega3.z*omega3.z;
		Sxy = sigma3*sigma3*omega3.x*omega3.y - omega3.x*omega3.y;
		Sxz = sigma3*sigma3*omega3.x*omega3.z - omega3.x*omega3.z;
		Syy = sigma3*sigma3*omega3.y*omega3.y + omega3.x*omega3.x + omega3.z*omega3.z;
		Syz = sigma3*sigma3*omega3.y*omega3.z - omega3.y*omega3.z;
		Szz = sigma3*sigma3*omega3.z*omega3.z + omega3.x*omega3.x + omega3.y*omega3.y;
		vec3 wi(pRec.wi.x, pRec.wi.y, pRec.wi.z);
		vec3 wo(pRec.wo.x, pRec.wo.y, pRec.wo.z);
		wi = normalize(wi);
		wo = normalize(wo);
		//if (fabs(length) < 0.0001) return 0.0f;
		float value = eval_specular(wi, wo, Sxx, Syy, Szz, Sxy, Sxz, Syz);

	//	Frame frame(pRec.mRec.orientation);
	//	Vector wi = frame.toLocal(pRec.wi);
	//	Vector wo = frame.toLocal(pRec.wo);
	//	Vector H = wi + wo;
	//	Float length = H.length();

	//	if (length == 0)
	//		return 0.0f;
		#if defined(SGGXDEBUGEVAL)
		if (value > 0.0000)
		{
			FILE *f = fopen("log_phase.txt", "a+");
		//	fprintf(f, "ori %f %f %f\n", omega3.x, omega3.y, omega3.z);
		//	fprintf(f, "wi %f %f %f\n", wi.x, wi.y, wi.z);
		//	fprintf(f, "wo %f %f %f\n", wo.x, wo.y, wo.z);
			fprintf(f, "%f\n", value);
			fclose(f);
		}
		#endif
		return value;
	}

	inline Float sample(PhaseFunctionSamplingRecord &pRec, Sampler *sampler) const {
		if (pRec.mRec.orientation.isZero()) {
			/* What to do when the local orientation is undefined */
			#if 0
				pRec.wo = warp::squareToUniformSphere(sampler->next2D());
				return 1.0f;
			#else
				return 0.0f;
			#endif
		}
		Float Sxx, Sxy, Sxz, Syy, Syz, Szz;
		vec3 omega3(pRec.mRec.orientation.x, pRec.mRec.orientation.y, pRec.mRec.orientation.z);
		omega3 = normalize(omega3);
		Sxx = sigma3*sigma3*omega3.x*omega3.x + omega3.y*omega3.y + omega3.z*omega3.z;
		Sxy = sigma3*sigma3*omega3.x*omega3.y - omega3.x*omega3.y;
		Sxz = sigma3*sigma3*omega3.x*omega3.z - omega3.x*omega3.z;
		Syy = sigma3*sigma3*omega3.y*omega3.y + omega3.x*omega3.x + omega3.z*omega3.z;
		Syz = sigma3*sigma3*omega3.y*omega3.z - omega3.y*omega3.z;
		Szz = sigma3*sigma3*omega3.z*omega3.z + omega3.x*omega3.x + omega3.y*omega3.y;
		vec3 wi(pRec.wi.x, pRec.wi.y, pRec.wi.z);		

		#if defined(SGGX_STATISTICS)
			avgSampleIterations.incrementBase();
		#endif
		wi=normalize(wi);
		Point2 point = sampler->next2D();
		vec3 wo=sample_specular(wi, Sxx, Syy, Szz, Sxy, Sxz, Syz, point.x, point.y);
		wo = normalize(wo);
		pRec.wo = Vector(wo.x,wo.y,wo.z);
#if defined(SGGXDEBUGSAMPLE)
		float length = sqrtf(dot(wo, wo));
		if (fabs(length - 1.0) > 0.0001)
		{
			FILE *f = fopen("log_phase_sample.txt", "a+");
			fprintf(f, "wo %f %f %f\n", wo.x, wo.y, wo.z);
			fprintf(f, "%f\n",length );
			fclose(f);
		}
#endif
		/*int iterations = 0, maxIterations = 1000;
		while (true) {
			Vector H = m_fiberDistr.sample(sampler->next2D());
			#if defined(SGGX_STATISTICS)
				++avgSampleIterations;
			#endif
			++iterations;

			if (sampler->next1D() < absDot(wi, H)) {
				Vector wo = H*(2*dot(wi, H)) - wi;
				pRec.wo = frame.toWorld(wo);
				break;
			}

			if (iterations >= maxIterations) {
				Log(EWarn, "Sample generation unsuccessful after %i iterations"
					" (dp=%f, fiberOrientation=%s, wi=%s)", iterations,
					absDot(pRec.wi, pRec.mRec.orientation),
					pRec.mRec.orientation.toString().c_str(),
					pRec.wi.toString().c_str());
				return 0.0f;
			}
		}*/

		return 1.0f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Float &pdf, Sampler *sampler) const {
		if (fabs(sample(pRec, sampler)) <1e-6) {
			pdf = 0; return 0.0f;
		}
		pdf = eval(pRec);
		return 1.0f;
	}

	bool needsDirectionallyVaryingCoefficients() const { return true; }

	Float sigmaDir(Float cosTheta) const {
		// Scaled such that replacing an isotropic phase function with an
		// isotropic microflake distribution does not cause changes
		Float Sxx, Sxy, Sxz, Syy, Syz, Szz;
		vec3 omega3(0.0, 0.0, 1.0);
		Sxx = sigma3*sigma3*omega3.x*omega3.x + omega3.y*omega3.y + omega3.z*omega3.z;
		Sxy = sigma3*sigma3*omega3.x*omega3.y - omega3.x*omega3.y;
		Sxz = sigma3*sigma3*omega3.x*omega3.z - omega3.x*omega3.z;
		Syy = sigma3*sigma3*omega3.y*omega3.y + omega3.x*omega3.x + omega3.z*omega3.z;
		Syz = sigma3*sigma3*omega3.y*omega3.z - omega3.y*omega3.z;
		Szz = sigma3*sigma3*omega3.z*omega3.z + omega3.x*omega3.x + omega3.y*omega3.y;
		vec3 wi = vec3(sqrtf(1 - (cosTheta*cosTheta)), 0, cosTheta);
		float ret = sigma(wi,Sxx,Syy,Szz,Sxy,Sxz,Syz);
		return ret;
		//return 2 * m_fiberDistr.sigmaT(cosTheta);
	}

	Float sigmaDirMax() const {
//#if defined(SGGXDEBUG)
//		FILE *f = fopen("log_phase_max.txt", "a+");
//		fprintf(f, "sigmadirmax %f\n", sigmaDir(0));
//		fprintf(f, "osigmadirmax %f\n", OsigmaDir(0));
//		fclose(f);
//#endif
		return sigmaDir(0);
	}
	Float OsigmaDir(Float cosTheta) const {
		// Scaled such that replacing an isotropic phase function with an
		// isotropic microflake distribution does not cause changes
		return 2 * m_fiberDistr.sigmaT(cosTheta);
	}
	std::string toString() const {
		std::ostringstream oss;
		oss << "SggxPhaseFunction[" << endl
			<< "   fiberDistr = " << indent(m_fiberDistr.toString()) << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	GaussianFiberDistribution m_fiberDistr;
	Float sigma1, sigma2, sigma3;
	static int count;
};
int SggxPhaseFunction::count = 0;
MTS_IMPLEMENT_CLASS_S(SggxPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(SggxPhaseFunction, "Sggx phase function");
MTS_NAMESPACE_END
