
#ifndef CORRESPOND_CU_
#define CORRESPOND_CU_

#include "correspond.hpp"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/transform.h>

namespace pt {

using namespace std;

class Part {
public:
	float X[2];
	int id;
};

class DistanceCalc : public thrust::binary_function<Part&, Part&, float>
{
public:
	__host__ __device__
		float operator()(Part& A, Part& B) { 
			float dx = A.X[0]-B.X[0];
			float dy = A.X[1]-B.X[1];
			return sqrtf(dx*dx + dy*dy); 
	}
};


class EpipolarLineCalc : public thrust::binary_function<Part&, float**, float3>
{
public:
//__host__ __device__
	float3 operator()(Part& A, float** F) { 
		float x = A.X[0];
		float y = A.X[1];
		float z = 1.0;
		float3 line;
		line.x = F[0][0] * x + F[1][0] * y + F[2][0] * z;
		line.y = F[0][1] * x + F[1][1] * y + F[2][1] * z;
		line.z = F[0][2] * x + F[1][2] * y + F[2][2] * z;
		float factor = line.x * line.x + line.y * line.y;
		factor = factor ? 1./sqrt(factor) : 1.;
		line.x *= factor;
		line.y *= factor;
		line.z *= factor;
		return line;
	}
};

class EpipolarResidualCalc : public thrust::binary_function<Part&, float3&, float>
{
public:
	__host__ __device__
		float operator()(Part& A, float3& line) { 
			return fabs( A.X[0] * line.x + A.X[1] * line.y + line.z );
	}
};

void calcDistancesHOST(thrust::host_vector<Part>& A, thrust::host_vector<Part>& B, thrust::host_vector<float>& distances);
void calcDistancesDEVICE(thrust::device_vector<Part>& A, thrust::device_vector<Part>& B, thrust::device_vector<float>& distances);


} /* NAMESPACE_PT */

#endif /* CORRESPOND_CU_ */
