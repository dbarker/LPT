#include <correspondcuda.h>

using namespace std;

__constant__ lpt::CameraPairCUDA pairs_k[100];

__constant__ int num_matches_k;

__global__ void calcEpipolarResidualAllInOneStreams_kernel(int p, float match_threshold, lpt::KernelArray<float> particles_x, lpt::KernelArray<float> particles_y, lpt::KernelArray<int> num_particles, lpt::KernelArray<lpt::MatchIDs> matches2way, lpt::KernelArray<int> num_matches )
{	
	p += blockIdx.x;
	int b = blockIdx.y * blockDim.x + threadIdx.x;
	int id_b = b;
	if (pairs_k[p].cam_b_id != 0)
		id_b += num_particles.data[pairs_k[p].cam_b_id - 1];
	
	if (id_b < num_particles.data[pairs_k[p].cam_b_id] ) {	
		float line[3];
		float x = particles_x.data[id_b];
		float y = particles_y.data[id_b];
		
		line[0] = pairs_k[p].F[0][0] * x + pairs_k[p].F[1][0] * y + pairs_k[p].F[2][0] * 1.f;
		line[1] = pairs_k[p].F[0][1] * x + pairs_k[p].F[1][1] * y + pairs_k[p].F[2][1] * 1.f;
		line[2] = pairs_k[p].F[0][2] * x + pairs_k[p].F[1][2] * y + pairs_k[p].F[2][2] * 1.f;
		float factor = line[0] * line[0] + line[1] * line[1];
		factor = factor ? 1.f/sqrtf(factor) : 1.f;
		line[0] *= factor;
		line[1] *= factor;
		line[2] *= factor;
     
		int id_a_start = 0;
		if (pairs_k[p].cam_a_id != 0) 
			id_a_start = num_particles.data[ pairs_k[p].cam_a_id - 1];
		int id_a_end = num_particles.data[ pairs_k[p].cam_a_id ];
		
		int match_id = b;
		if (p !=0 )
			match_id += num_matches.data[p - 1];
		int match_count = 0; 
		for (int id_a = id_a_start; id_a < id_a_end; ++id_a) {		
			bool matched = static_cast<int>( floor( match_threshold / fabs( particles_x.data[id_a] * line[0]  + particles_y.data[id_a] * line[1] + line[2] ) ) );
			if ( matched && match_count < num_matches_k ) { 
				matches2way.data[match_id].ids[match_count] = id_a;
				++match_count;
			}
		}
	}
}

__global__ void calcEpipolarResidualAllInOne_kernel( float match_threshold, lpt::KernelArray<float> particles_x, lpt::KernelArray<float> particles_y, lpt::KernelArray<int> num_particles, lpt::KernelArray<lpt::MatchIDs> matches2way, lpt::KernelArray<int> num_matches )
{	
	int p = blockIdx.x; 
	int b = blockIdx.y * blockDim.x + threadIdx.x;
	int id_b = b;
	if (pairs_k[p].cam_b_id != 0)
		id_b = num_particles.data[pairs_k[p].cam_b_id - 1] + b;
	
	if (id_b < num_particles.data[pairs_k[p].cam_b_id] ) {	
		float line[3];
		float x = particles_x.data[id_b]; 
		float y = particles_y.data[id_b]; 
		
		line[0] = pairs_k[p].F[0][0] * x + pairs_k[p].F[1][0] * y + pairs_k[p].F[2][0] * 1.f;
		line[1] = pairs_k[p].F[0][1] * x + pairs_k[p].F[1][1] * y + pairs_k[p].F[2][1] * 1.f;
		line[2] = pairs_k[p].F[0][2] * x + pairs_k[p].F[1][2] * y + pairs_k[p].F[2][2] * 1.f;
		float factor = line[0] * line[0] + line[1] * line[1];
		factor = factor ? 1.f/sqrtf(factor) : 1.f;
		line[0] *= factor;
		line[1] *= factor;
		line[2] *= factor;
     
		int id_a_start = 0;
		if (pairs_k[p].cam_a_id != 0) 
			id_a_start = num_particles.data[ pairs_k[p].cam_a_id - 1];
		int id_a_end = num_particles.data[ pairs_k[p].cam_a_id ];
		int match_id = b;
		if (p !=0 )
			match_id += num_matches.data[p - 1];
		int match_count = 0; 
		for (int id_a = id_a_start; id_a < id_a_end; ++id_a) {		
			bool matched = static_cast<int>( floor( match_threshold / fabs( particles_x.data[id_a] * line[0]  + particles_y.data[id_a] * line[1] + line[2] ) ) );
			if ( matched && match_count < num_matches_k ) { 
				matches2way.data[match_id].ids[match_count] = id_a;
				++match_count;
			}
		}
	}
}

__global__ void calcEpipolarLines_kernel(lpt::KernelArray<float> particles_x, lpt::KernelArray<float> particles_y, lpt::KernelArray<int> num_particles, lpt::KernelArray<float> lines_x, lpt::KernelArray<float> lines_y, lpt::KernelArray<float> lines_z, lpt::KernelArray<int> num_lines )
{	
	int p = blockIdx.x; 
	int b = blockIdx.y * blockDim.x + threadIdx.x;
	int cam_b_id = pairs_k[p].cam_b_id;
	int id_b = b;
	if (cam_b_id != 0)
		id_b = num_particles.data[cam_b_id - 1] + b;

	if (id_b < num_particles.data[cam_b_id] ) {

		float x = particles_x.data[id_b];
		float y = particles_y.data[id_b];
		int line_id = b;
		if (p !=0 )
			line_id += num_lines.data[p -1];

		lines_x.data[line_id] = pairs_k[p].F[0][0] * x + pairs_k[p].F[1][0] * y + pairs_k[p].F[2][0] * 1.f;
		lines_y.data[line_id] = pairs_k[p].F[0][1] * x + pairs_k[p].F[1][1] * y + pairs_k[p].F[2][1] * 1.f;
		lines_z.data[line_id] = pairs_k[p].F[0][2] * x + pairs_k[p].F[1][2] * y + pairs_k[p].F[2][2] * 1.f;
		float factor = lines_x.data[line_id] * lines_x.data[line_id] + lines_y.data[line_id] * lines_y.data[line_id];
		factor = factor ? 1.f/sqrtf(factor) : 1.f;
		lines_x.data[line_id] *= factor;
		lines_y.data[line_id] *= factor;
		lines_z.data[line_id] *= factor; 
	}
}

__global__ void calcEpipolarResiduals_kernel(float match_threshold, lpt::KernelArray<float> particles_x, lpt::KernelArray<float> particles_y, lpt::KernelArray<int> num_particles, lpt::KernelArray<float> lines_x, lpt::KernelArray<float> lines_y, lpt::KernelArray<float> lines_z, lpt::KernelArray<int> num_lines, lpt::KernelArray<lpt::MatchIDs> matches2way, lpt::KernelArray<int> num_matches)
{
	////__shared__ float line[3];
	//
	//int line_id = blockIdx.x;
	////float line[] = {lines_x.data[line_id], lines_y.data[line_id], lines_z.data[line_id]};
	//int cam_id;
	//int r_id = residuals.size;
	//for (int pair_id = 0; pair_id < num_lines.size; ++pair_id) {
	//	if( line_id < num_lines.data[pair_id] ) {
	//		cam_id = pairs_k[pair_id].cam_a_id;
	//		if (pair_id !=0 )
	//			r_id = num_residuals.data[pair_id - 1] + (line_id - num_lines.data[pair_id -1]) * (num_particles.data[cam_id] - num_particles.data[cam_id - 1]) + threadIdx.x;
	//		else
	//			r_id = (line_id - num_lines.data[pair_id -1]) * (num_particles.data[cam_id] - num_particles.data[cam_id - 1]) + threadIdx.x;
	//		break;
	//	}
	//}

	////if (threadIdx.x == 0) {
	////	line[0] = lines_x.data[line_id];
	////	line[1] = lines_y.data[line_id];
	////	line[2] = lines_z.data[line_id];
	////}	
	////__syncthreads(); //FIXME: put lines in constant memory if possible

	//int id_a = num_particles.data[cam_id - 1] + blockIdx.y * blockDim.x + threadIdx.x;
	//
	//if ( id_a < num_particles.data[cam_id] && r_id < residuals.size ) {
	//	residuals.data[r_id] = static_cast<int>( match_threshold / fabs( particles_x.data[id_a] * lines_x.data[line_id]  + particles_y.data[id_a] * lines_y.data[line_id] + lines_z.data[line_id] ) );
	//}

}

namespace lpt {

PointMatcherCUDA::PointMatcherCUDA() {
	cout << "Epipolor Point matcher created (CUDA Enabled)" << endl;
	int devcount = 0;
	cudaGetDeviceCount(&devcount);
	for (int i = 0; i < devcount; ++i) {
		cudaDeviceProp device_prop;
		cudaGetDeviceProperties(&device_prop, i);
		if (! device_prop.kernelExecTimeoutEnabled ) {
			cout << "Device " << i << ":  " <<  device_prop.name << "  added to available queue" <<endl;
		} else {
			cout << "Device " << i << ":  " << device_prop.name << " added to available queue (Kernel run time limited)" << endl;
		}
		compute_devices_available.push(i);
	}
}

int PointMatcherCUDA::getNextComputeDeviceID() {
	boost::mutex::scoped_lock(this->mutex);
	int id = this->compute_devices_available.front();
	this->compute_devices_available.pop();
	return id;
}

void PointMatcherCUDA::initializeEpipolarMatchThread(int thread_id) {
	auto& cameras = shared_objects->cameras;
	auto& camera_pairs = shared_objects->camera_pairs;
	int id = getNextComputeDeviceID();
	cudaSetDevice( id );
	cout << "PointMatcherCUDA Thread " << thread_id << " setting device " << id << endl;

	particles_x_h.resize(this->initial_max_particles_per_image * cameras.size(), 0.f );
	particles_x_d = particles_x_h;
	
	particles_y_h.resize(this->initial_max_particles_per_image * cameras.size(), 0.f );
	particles_y_d = particles_y_h;
	
	num_particles_h.resize(cameras.size(), 0 );
	num_particles_d = num_particles_h;

	camera_pairs_h.resize( camera_pairs.size() );
	
	num_matches_h.resize( camera_pairs.size(), 0 );
	num_matches_d = num_matches_h;

	for (int i = 0; i < NUM_MATCHES; ++i)
		this->match_initializer.ids[i] = -1;

	matches2way_h.resize( camera_pairs.size() * this->initial_max_particles_per_image, this->match_initializer);
	matches2way_d = matches2way_h;

	for (int i = 0; i < camera_pairs.size(); ++i) {
		for (int n = 0; n < 3; ++n)
			for (int m = 0; m < 3; ++m)
				camera_pairs_h[i].F[n][m] = static_cast<float>(camera_pairs[i].F[n][m]);
		camera_pairs_h[i].cam_a_id = camera_pairs[i].cam_A.id;
		camera_pairs_h[i].cam_b_id = camera_pairs[i].cam_B.id;
	}
			
	streams.clear();
	for (int f = 2; f < camera_pairs_h.size(); ++f) {
		//if ( camera_pairs_h.size() % f == 0 ) {
			streams.resize(f);
			//break;
		//}
	}
	cout << "Streams size = " << streams.size() << endl;
	for(int i = 0; i < streams.size(); ++i) 
        cudaStreamCreate(&(streams[i]));

	int num = NUM_MATCHES;

	cudaMemcpyToSymbol( "num_matches_k", &num,  sizeof(int) );
	cudaMemcpyToSymbol( "pairs_k", thrust::raw_pointer_cast(&camera_pairs_h[0]),  sizeof(CameraPairCUDA) * camera_pairs_h.size() );
}

void PointMatcherCUDA::initialize() {
	this->initializeMatchMap();
}

void PointMatcherCUDA::addControls() {
	void* matcher_void_ptr = static_cast<void*> ( this );
	cv::createTrackbar("Match Thresh", string() , &params.match_thresh_level, 100, callbackMatchThreshcuda, matcher_void_ptr);
}

void PointMatcherCUDA::findEpipolarMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap) {
	
	thrust::fill(matches2way_d.begin(), matches2way_d.end(), match_initializer);
	
	num_particles_h[0] = frame_group[0].particles.size();
	for(int p = 0; p < frame_group[0].particles.size(); ++p) {
			particles_x_h[p] = static_cast<float>(frame_group[0].particles[p]->x);
			particles_y_h[p] = static_cast<float>(frame_group[0].particles[p]->y);
	}

	int max_particles = num_particles_h[0];
	for(int i = 1; i < frame_group.size(); ++i) {
		num_particles_h[i] = frame_group[i].particles.size() + num_particles_h[i-1];
		for(int p = 0; p < frame_group[i].particles.size(); ++p) {
			particles_x_h[ num_particles_h[i-1] + p] = static_cast<float>(frame_group[i].particles[p]->x);
			particles_y_h[ num_particles_h[i-1] + p] = static_cast<float>(frame_group[i].particles[p]->y);
		}
		if (frame_group[i].particles.size() > max_particles) {
			max_particles = frame_group[i].particles.size();
			if ( max_particles > this->initial_max_particles_per_image) 
				cout << "WARNING IMAGE FRAME HAS EXCEEDED MAX PARTICLES:  correspondcuda.cu" << endl;
		}
	}

	cudaMemcpyAsync(thrust::raw_pointer_cast(&num_particles_d[0]), thrust::raw_pointer_cast(&num_particles_h[0]), num_particles_h.size() * sizeof(int), cudaMemcpyHostToDevice, streams[0]);
	cudaMemcpyAsync(thrust::raw_pointer_cast(&particles_x_d[0]), thrust::raw_pointer_cast(&particles_x_h[0]), *num_particles_h.rbegin() * sizeof(float), cudaMemcpyHostToDevice, streams[0]);
	cudaMemcpyAsync(thrust::raw_pointer_cast(&particles_y_d[0]), thrust::raw_pointer_cast(&particles_y_h[0]), *num_particles_h.rbegin() * sizeof(float), cudaMemcpyHostToDevice, streams[0]);

	num_matches_h[0] = frame_group[camera_pairs_h[0].cam_b_id].particles.size(); 
	for(int i = 1; i < this->camera_pairs_h.size(); ++i ) {
		num_matches_h[i] = frame_group[camera_pairs_h[i].cam_b_id].particles.size() + num_matches_h[i-1]; 
	}

	cudaMemcpyAsync(thrust::raw_pointer_cast(&num_matches_d[0]), thrust::raw_pointer_cast(&num_matches_h[0]), num_matches_h.size() * sizeof(int), cudaMemcpyHostToDevice, streams[0]);
	
	int num_pairs = camera_pairs_h.size();

	dim3 dimblock(128,1,1);
	dim3 dimgrid( static_cast<unsigned int>(num_pairs), (static_cast<unsigned int>(max_particles) / dimblock.x ) + 1 );
	
	calcEpipolarResidualAllInOne_kernel <<< dimgrid, dimblock, 0, streams[0] >>> (params.match_threshold, particles_x_d, particles_y_d, num_particles_d, matches2way_d, num_matches_d);
	
	cudaMemcpyAsync(thrust::raw_pointer_cast(&matches2way_h[0]), thrust::raw_pointer_cast(&matches2way_d[0]), *num_matches_h.rbegin() * sizeof(MatchIDs), cudaMemcpyDeviceToHost, streams[0]);
	
	cudaStreamSynchronize(streams[0]);
	
	int match_overload = 0;
	for (int p = 0; p < camera_pairs_h.size(); ++p) {
		int match_id = (p == 0 ? 0 : num_matches_h[p-1]);
		int cam_b = camera_pairs_h[p].cam_b_id; 
		int cam_a = camera_pairs_h[p].cam_a_id;
		int b_end = num_particles_h[cam_b];
		int b_start = (cam_b !=0 ? num_particles_h[cam_b-1] : 0);
		int a_start = (cam_a !=0 ? num_particles_h[cam_a-1] : 0);
		for (int b_id = b_start; b_id < b_end; ++b_id, ++match_id) {
			for (int m = 0; m < NUM_MATCHES; ++m) {
				int a_id = matches2way_h[match_id].ids[m];
				if (a_id >= 0) {
					//matchmap[b_id][cam_a][m] = a_id - a_start;
					//matchmap[a_id][cam_b][m] = b_id - b_start; 
					auto itb = std::find(matchmap[b_id][cam_a].begin(), matchmap[b_id][cam_a].end(), -1);
					auto ita = std::find(matchmap[a_id][cam_b].begin(), matchmap[a_id][cam_b].end(), -1);
					if (itb != matchmap[b_id][cam_a].end() && ita != matchmap[a_id][cam_b].end() ) {
						*itb = a_id - static_cast<int>(a_start);
						*ita = static_cast<int>(b_id - b_start);
					} else
						match_overload++;
				}
				else				
					break;
			}				
		}
	}
	if (match_overload > 0)
		;//cout << "WARNING: MORE MATCHES THAN ARRAY SIZE NUM_MATCHES: total overload = " << match_overload << endl;
}

void PointMatcherCUDA::findUniqueMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches) {
	vector<int> num_particles(frame_group.size());
	num_particles[0] = frame_group[0].particles.size();
	for(int i = 1; i < frame_group.size(); ++i) 
		num_particles[i] = frame_group[i].particles.size() + num_particles[i-1];

    for (int cam_a = 0; cam_a < frame_group.size() - 3; ++cam_a) {
		int a_start = (cam_a !=0 ? num_particles[cam_a - 1] : 0);
		for (int a = 0; a < frame_group[cam_a].particles.size(); ++a) {
            lpt::ParticleImage::Ptr Pa = frame_group[cam_a].particles[a];
			if( ! Pa->is_4way_matched ) 
            for (int cam_b = cam_a + 1; cam_b < frame_group.size() - 2; ++cam_b) {
				int b_start = (cam_b !=0 ? num_particles[cam_b-1] : 0);
                for(int match_ab = 0; match_ab < NUM_MATCHES; ++match_ab) { //loop through all A,B matches
					int b = matchmap[a + a_start][cam_b][match_ab]; 
					if (b < 0)
						break;
					lpt::ParticleImage::Ptr Pb = frame_group[cam_b].particles[b];
						
					if( ! Pb->is_4way_matched ) 
					for (int cam_c = cam_b + 1; cam_c < frame_group.size() - 1; ++cam_c) {
                        int c_start = (cam_c !=0 ? num_particles[cam_c-1] : 0);
						for (int match_bc = 0; match_bc < NUM_MATCHES; ++match_bc) {
                            int c = matchmap[b + b_start][cam_c][match_bc];
							if (c < 0) 
								break;
								
							lpt::ParticleImage::Ptr Pc = frame_group[cam_c].particles[c];

							if( ! Pc->is_4way_matched && std::count(matchmap[a + a_start][cam_c].begin(), matchmap[a + a_start][cam_c].end(), c) )  
                            for (int cam_d = cam_c + 1; cam_d < frame_group.size(); ++cam_d) {
								vector<lpt::Match::Ptr> matches4way;
                                int d_start = (cam_d !=0 ? num_particles[cam_d-1] : 0);
								for (int match_cd = 0; match_cd < NUM_MATCHES; ++match_cd) {
									int d = matchmap[c + c_start][cam_d][match_cd];
									if (d < 0)
										break;
									lpt::ParticleImage::Ptr Pd = frame_group[cam_d].particles[d];
									if( ! Pd->is_4way_matched && std::count(matchmap[a + a_start][cam_d].begin(), matchmap[a + a_start][cam_d].end(), d)  && std::count(matchmap[b + b_start][cam_d].begin(), matchmap[b+b_start][cam_d].end(), d)  ) {
										if(! Pa->is_4way_matched && ! Pb->is_4way_matched && ! Pc->is_4way_matched && ! Pd->is_4way_matched) { 
											lpt::Match::Ptr newmatch = lpt::Match::create();
											newmatch->addParticle(Pa,cam_a);
											newmatch->addParticle(Pb,cam_b);
											newmatch->addParticle(Pc,cam_c);
											newmatch->addParticle(Pd,cam_d);
											matches4way.push_back(std::move(newmatch));
											Pa->is_4way_matched = true;
											Pb->is_4way_matched = true;
											Pc->is_4way_matched = true;
											Pd->is_4way_matched = true;
											match_ab = NUM_MATCHES;
											match_bc = NUM_MATCHES;
											match_cd = NUM_MATCHES;
												
										}
                                    } 
                                }
								std::move(matches4way.begin(), matches4way.end(), std::back_inserter(matches) );
                            }
                        }
                    }
                }
            }
		}
    }
}

void PointMatcherCUDA::find3WayMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches) {
	vector<int> num_particles(frame_group.size());
	num_particles[0] = frame_group[0].particles.size();
	for(int i = 1; i < frame_group.size(); ++i) 
		num_particles[i] = frame_group[i].particles.size() + num_particles[i-1];
	
	int num_cameras = frame_group.size();
	matches.clear();

	for (int cam_a = 0; cam_a < frame_group.size() - 2; ++cam_a) {
		int a_start = (cam_a !=0 ? num_particles[cam_a - 1] : 0);
		for (int a = 0; a < frame_group[cam_a].particles.size(); ++a) {
            lpt::ParticleImage::Ptr Pa = frame_group[cam_a].particles[a];
			if( ! Pa->is_4way_matched ) {
                for (int cam_b = cam_a + 1; cam_b < frame_group.size() - 1; ++cam_b) {
					int b_start = (cam_b !=0 ? num_particles[cam_b-1] : 0);
                    for(int match_ab = 0; match_ab < NUM_MATCHES; ++match_ab) { //loop through all A,B matches
						int b = matchmap[a + a_start][cam_b][match_ab]; 
						if (b < 0)
							break;
						lpt::ParticleImage::Ptr Pb = frame_group[cam_b].particles[b];
						
						if( ! Pb->is_4way_matched ) {
							for (int cam_c = cam_b + 1; cam_c < frame_group.size(); ++cam_c) {
								int c_start = (cam_c !=0 ? num_particles[cam_c-1] : 0);
								for (int match_bc = 0; match_bc < NUM_MATCHES; ++match_bc) {
									int c = matchmap[b + b_start][cam_c][match_bc];
									if (c < 0) 
										break;
								
									lpt::ParticleImage::Ptr Pc = frame_group[cam_c].particles[c];

									if( ! Pc->is_4way_matched && std::count(matchmap[a + a_start][cam_c].begin(), matchmap[a + a_start][cam_c].end(), c) ) {
										lpt::Match::Ptr newmatch = lpt::Match::create();
										newmatch->addParticle(Pa,cam_a);
										newmatch->addParticle(Pb,cam_b);
										newmatch->addParticle(Pc,cam_c);

										matches.push_back(std::move(newmatch));
												
										Pa->is_4way_matched = true;
										Pb->is_4way_matched = true;
										Pc->is_4way_matched = true;

										match_ab = NUM_MATCHES;
										match_bc = NUM_MATCHES;
												
										cam_b = num_cameras;
										cam_c = num_cameras;
									}
                                }
                            }
                        }
                    }
                }
			}
		}
	}
	//cout << matches.size() << endl;
}

void PointMatcherCUDA::findEpipolarMatchesStreams(lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap) {
	thrust::fill(matches2way_d.begin(), matches2way_d.end(), match_initializer);
	num_particles_h[0] = frame_group[0].particles.size();
	for(int p = 0; p < frame_group[0].particles.size(); ++p) {
			particles_x_h[p] = static_cast<float>(frame_group[0].particles[p]->x);
			particles_y_h[p] = static_cast<float>(frame_group[0].particles[p]->y);
	}

	int max_particles = num_particles_h[0];
	for(int i = 1; i < frame_group.size(); ++i) {
		num_particles_h[i] = frame_group[i].particles.size() + num_particles_h[i-1];
		for(int p = 0; p < frame_group[i].particles.size(); ++p) {
			particles_x_h[ num_particles_h[i-1] + p] = static_cast<float>(frame_group[i].particles[p]->x);
			particles_y_h[ num_particles_h[i-1] + p] = static_cast<float>(frame_group[i].particles[p]->y);
		}
		if (frame_group[i].particles.size() > max_particles)
			max_particles = frame_group[i].particles.size();
	}
	
	cudaMemcpyAsync(thrust::raw_pointer_cast(&num_particles_d[0]), thrust::raw_pointer_cast(&num_particles_h[0]), num_particles_h.size() * sizeof(int), cudaMemcpyHostToDevice, streams[0]);
	cudaMemcpyAsync(thrust::raw_pointer_cast(&particles_x_d[0]), thrust::raw_pointer_cast(&particles_x_h[0]), *num_particles_h.rbegin() * sizeof(float), cudaMemcpyHostToDevice, streams[0]);
	cudaMemcpyAsync(thrust::raw_pointer_cast(&particles_y_d[0]), thrust::raw_pointer_cast(&particles_y_h[0]), *num_particles_h.rbegin() * sizeof(float), cudaMemcpyHostToDevice, streams[0]);

	int num_pairs = camera_pairs_h.size();
	dim3 dimblock(128,1,1);
	dim3 dimgrid( static_cast<unsigned int>(num_pairs / streams.size()), (static_cast<unsigned int>(max_particles) / dimblock.x ) + 1, 1 );
	 
	num_matches_h[0] = frame_group[camera_pairs_h[0].cam_b_id].particles.size(); 
	for(int i = 1; i < this->camera_pairs_h.size(); ++i ) {
		num_matches_h[i] = frame_group[camera_pairs_h[i].cam_b_id].particles.size() + num_matches_h[i-1]; 
	}

	cudaMemcpyAsync(thrust::raw_pointer_cast(&num_matches_d[0]), thrust::raw_pointer_cast(&num_matches_h[0]), num_matches_h.size() * sizeof(int), cudaMemcpyHostToDevice, streams[0]);
	
    for(int i = 0; i < streams.size(); i++) 
		calcEpipolarResidualAllInOneStreams_kernel <<< dimgrid, dimblock, 0, streams[i] >>> (i * dimgrid.x , params.match_threshold, particles_x_d, particles_y_d, num_particles_d, matches2way_d, num_matches_d);

	for(int i = 0; i < streams.size(); i++) {
		 int index = (i == 0 ? 0 : num_matches_h[ i * dimgrid.x - 1]);
		 int nbytes = ( i == 0 ?  num_matches_h[dimgrid.x - 1] :  num_matches_h[(i+1)*dimgrid.x - 1] -  num_matches_h[i*dimgrid.x - 1]) * sizeof(MatchIDs);
		 cudaMemcpyAsync(thrust::raw_pointer_cast(&matches2way_h[index]), thrust::raw_pointer_cast(&matches2way_d[index]), nbytes, cudaMemcpyDeviceToHost, streams[i]);
	}
	
	int match_overload = 0;
	for (unsigned int i = 0; i < streams.size(); ++i) {
		cudaStreamSynchronize(streams[i]);
		for (unsigned int p = i * dimgrid.x; p < (i + 1) * dimgrid.x; ++p) {
			int match_id = (p == 0 ? 0 : num_matches_h[p-1]);
			int cam_b = camera_pairs_h[p].cam_b_id; 
			int cam_a = camera_pairs_h[p].cam_a_id;
			int b_end = num_particles_h[cam_b];
			int b_start = (cam_b !=0 ? num_particles_h[cam_b-1] : 0);
			int a_start = (cam_a !=0 ? num_particles_h[cam_a-1] : 0);
			for (int b_id = b_start; b_id < b_end; ++b_id, ++match_id) {
				for (int m = 0; m < NUM_MATCHES; ++m) {
					int a_id = matches2way_h[match_id].ids[m];
					if (a_id >= 0) {
						//matchmap[b_id][cam_a][m] = a_id - a_start;
						//matchmap[a_id][cam_b][m] = b_id - b_start; 
						auto itb = std::find(matchmap[b_id][cam_a].begin(), matchmap[b_id][cam_a].end(), -1);
						auto ita = std::find(matchmap[a_id][cam_b].begin(), matchmap[a_id][cam_b].end(), -1);
						if (itb != matchmap[b_id][cam_a].end() && ita != matchmap[a_id][cam_b].end() ) {
							*itb = a_id - static_cast<int>(a_start);
							*ita = static_cast<int>(b_id - b_start);
						} else
							match_overload++;
					}
					else				
						break;
				}				
			}
		}
	}

	if (match_overload > 0)
		;//cout << "WARNING: MORE MATCHES THAN ARRAY SIZE NUM_MATCHES: total overload = " << match_overload << endl;
}

void PointMatcherCUDA::findEpipolarMatchesManyThreads(lpt::ImageFrameGroup& frame_group) {
	
	num_particles_h[0] = frame_group[0].particles.size();
	for(int p = 0; p < frame_group[0].particles.size(); ++p) {
			particles_x_h[p] = static_cast<float>(frame_group[0].particles[p]->x);
			particles_y_h[p] = static_cast<float>(frame_group[0].particles[p]->y);
	}

	int max_particles = num_particles_h[0];
	for(int i = 1; i < frame_group.size(); ++i) {
		num_particles_h[i] = frame_group[i].particles.size() + num_particles_h[i-1];
		for(int p = 0; p < frame_group[i].particles.size(); ++p) {
			particles_x_h[ num_particles_h[i-1] + p] = static_cast<float>(frame_group[i].particles[p]->x);
			particles_y_h[ num_particles_h[i-1] + p] = static_cast<float>(frame_group[i].particles[p]->y);
		}
		if (frame_group[i].particles.size() > max_particles)
			max_particles = frame_group[i].particles.size();
	}
		
	cudaMemcpyAsync(thrust::raw_pointer_cast(&num_particles_d[0]), thrust::raw_pointer_cast(&num_particles_h[0]), num_particles_h.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpyAsync(thrust::raw_pointer_cast(&particles_x_d[0]), thrust::raw_pointer_cast(&particles_x_h[0]), *num_particles_h.rbegin() * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpyAsync(thrust::raw_pointer_cast(&particles_y_d[0]), thrust::raw_pointer_cast(&particles_y_h[0]), *num_particles_h.rbegin() * sizeof(float), cudaMemcpyHostToDevice);

	int num_pairs = camera_pairs_h.size();
	dim3 dimblock(256);
	dim3 dimgrid( static_cast<unsigned int>(num_pairs), (static_cast<unsigned int>(max_particles) / dimblock.x ) + 1 );

	thrust::host_vector<int> num_lines_h(camera_pairs_h.size(), 0);
	num_lines_h[0] = frame_group[camera_pairs_h[0].cam_b_id].particles.size();
	for(int i = 1; i < this->camera_pairs_h.size(); ++i ) 
		num_lines_h[i] = frame_group[camera_pairs_h[i].cam_b_id].particles.size() + num_lines_h[i-1];
	
	thrust::device_vector<int> num_lines_d = num_lines_h;
	thrust::device_vector<float> lines_x( *num_lines_h.rbegin(), 0.f );
	thrust::device_vector<float> lines_y( *num_lines_h.rbegin(), 0.f );
	thrust::device_vector<float> lines_z( *num_lines_h.rbegin(), 0.f );
	
//	lines_x.resize( *num_lines_h.rbegin(), 0.f );
//	lines_y.resize( *num_lines_h.rbegin(), 0.f );
//	lines_z.resize( *num_lines_h.rbegin(), 0.f );	
	
	calcEpipolarLines_kernel <<< dimgrid, dimblock >>> (particles_x_d, particles_y_d, num_particles_d, lines_x, lines_y, lines_z, num_lines_d);
	
	num_matches_h[0] = frame_group[camera_pairs_h[0].cam_b_id].particles.size(); 
	for(int i = 1; i < this->camera_pairs_h.size(); ++i ) {
		num_matches_h[i] = frame_group[camera_pairs_h[i].cam_b_id].particles.size() + num_matches_h[i-1]; 
	}

	cudaMemcpyAsync(thrust::raw_pointer_cast(&num_matches_d[0]), thrust::raw_pointer_cast(&num_matches_h[0]), num_matches_h.size() * sizeof(int), cudaMemcpyHostToDevice);
	cudaStreamSynchronize(0);
	
	dim3 dimblock2(512,1,1);
	dim3 dimgrid2( static_cast<unsigned int>(*num_lines_h.rbegin()), ( static_cast<unsigned int>(max_particles) / dimblock2.x ) + 1, 1 );
	//cout <<"K2 Grid = " << dimgrid2.x << " x " << dimgrid2.y << " x " << dimgrid2.z << endl;
	//cout <<"K2 Block = " << dimblock2.x << " x " << dimblock2.y << " x " << dimblock2.z << endl;
	
	calcEpipolarResiduals_kernel <<< dimgrid2, dimblock2 >>> (params.match_threshold, particles_x_d, particles_y_d, num_particles_d, lines_x, lines_y, lines_z, num_lines_d, matches2way_d, num_matches_d);

	thrust::copy(matches2way_d.begin(), matches2way_d.begin() + *num_matches_h.rbegin(), matches2way_h.begin() );

}

}
