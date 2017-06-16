
#include "correspond.hpp"

namespace lpt {

using namespace std;


Correspondence::Correspondence()
    : initial_max_particles_per_image(2048), map_storage_size(500)
{
}

Correspondence::~Correspondence() {}

void Correspondence::findMatches(const ImageFrameGroup &frame_group, vector<Match::Ptr> &matches)
{
    resetMatchMap(current_matchmap);
    findEpipolarMatches(frame_group, current_matchmap);
    findUniqueMatches(frame_group, current_matchmap, matches);
}

void  Correspondence::run(
		lpt::concurrent_queue<lpt::ImageFrameGroup>* in_queue, 
		lpt::concurrent_queue<std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > >* out_queue,
		int number_epipolarmatching_threads,
        int number_uniquematching_threads )
{
	cout << "Running correspondence with " << number_epipolarmatching_threads<< " epi threads, and " << number_uniquematching_threads << " unique match threads " << endl; 
	//for (int i = 0; i < number_2way_threads; ++i)
	//	epipolar_match_threads.create_thread( boost::bind( &Correspondence::runEpipolarMatching, this, 1, in_queue ) );
	//for (int i = 0; i < number_4way_threads; ++i)
	//	unique_match_threads.create_thread( boost::bind( &Correspondence::runUniqueMatching, this, 1) );

	boost::thread epipolar_matcher_th = boost::thread( boost::bind(&Correspondence::runEpipolarMatching, this, 1, in_queue ) );

	runUniqueMatching(1, out_queue);

	//while (true) {
	//	// order output from 4way threads and send to out_queue
	//	lpt::ImageFrameGroup frame_group;
	//	vector<lpt::Match::Ptr> matches;
	//	auto pair = std::make_pair(frame_group, matches);
	//	final_match_queue.wait_and_pop(pair);
	//	// wait_pop from match_queue containing frame_group, matches pairs
	//	out_queue->push( pair );
	//	boost::this_thread::interruption_point();
	//}
}

void Correspondence::stop()
{
	epipolar_match_threads.interrupt_all();
	epipolar_match_threads.join_all();

	unique_match_threads.interrupt_all();
	unique_match_threads.join_all();
}

void Correspondence::runEpipolarMatching(int thread_id, lpt::concurrent_queue<lpt::ImageFrameGroup>* in_queue)
{
	this->initializeEpipolarMatchThread(thread_id);
	while (true) {
		lpt::ImageFrameGroup frame_group;
		
		in_queue->wait_and_pop(frame_group);
		
		vector<pair<lpt::ImageFrameGroup, lpt::MatchMap>>::iterator iter;
		
		this->empty_maps_queue.wait_and_pop( iter );
		
		lpt::MatchMap& current_map = iter->second;
		resetMatchMap(current_map);
		findEpipolarMatches(frame_group, current_map);
		//printMatchMap(frame_group, current_map, "matchmap.txt");
		iter->first = frame_group;
		
		this->full_maps_queue.push(iter);

		boost::this_thread::interruption_point();
	}
}

void Correspondence::runUniqueMatching(int thread_id, lpt::concurrent_queue<std::pair<lpt::ImageFrameGroup, vector<lpt::Match::Ptr> > >* out_queue)
{
	while (true) {
		vector<lpt::Match::Ptr> matches;
		vector<pair<lpt::ImageFrameGroup, lpt::MatchMap>>::iterator iter;
		this->full_maps_queue.wait_and_pop(iter);
		
		lpt::MatchMap& current_map = iter->second;
		lpt::ImageFrameGroup& group = iter->first;
		findUniqueMatches(group, current_map, matches);

		// push frame_group and matches to internal match_queue for reordering and pushing to out_queue
		
		//this->final_match_queue.push(std::make_pair(group, matches));
		out_queue->push(std::make_pair(group, matches));
		this->empty_maps_queue.push(iter);
		boost::this_thread::interruption_point();
	}
}

double Correspondence::calculateEpipolarResidual(double line[3], const lpt::ParticleImage& point) const
{
    return fabs( point.x * line[0] + point.y * line[1] + line[2] );
}

void Correspondence::calculateEpiline(const lpt::ParticleImage& point, const double F[3][3], double line[3]) const
{
    double x = point.x;
    double y = point.y;
    double z = 1.0;

    line[0] = F[0][0] * x + F[1][0] * y + F[2][0] * z;
    line[1] = F[0][1] * x + F[1][1] * y + F[2][1] * z;
    line[2] = F[0][2] * x + F[1][2] * y + F[2][2] * z;

    double factor = line[0] * line[0] + line[1] * line[1];
    factor = factor ? 1./sqrt(factor) : 1.;
    line[0] *= factor;
    line[1] *= factor;
    line[2] *= factor;
}

void Correspondence::testMatches(const ImageFrameGroup &cameragroup, const vector<Match::Ptr> &matches) const
{
    size_t N_correct = 0;
    size_t N_matched = 0;
    size_t N_total = 100000000;
    for(int c = 0; c < cameragroup.size(); c++) {
        size_t temp = cameragroup[c].particles.size();
        if (temp < N_total)
            N_total = temp;
    }

    for (int m = 0; m < matches.size(); m++){
        if ( matches[m] ) {
            N_matched++;
            auto& particles = matches[m]->particles;
            bool matched = true;
            for (int j = 0; j < particles.size() - 1; j++){
                matched &= particles[j].first->id == particles[j+1].first->id ? true : false;
            }
            if ( matched )
                N_correct++;
        }
    }
    cout << "\t Correct Ratio: " << N_correct << " / " << N_matched << " = " <<
        setprecision(2) << (N_matched > 0 ? (double) N_correct / N_matched : 0.0 ) << endl;
    cout << "\t Cover Ratio: " << N_correct << " / " << N_total << " = " <<
        setprecision(2) <<  (N_total > 0 ? (double) N_correct / N_total : 0.0 ) << endl;

}

void Correspondence::printMatchMap(const lpt::ImageFrameGroup& frame_group, string output_file_name) const
{
	cout << "printing matchmap..." << endl;
    ofstream fout(output_file_name.c_str());
    int p_id = 0;
    for(int i = 0; i < frame_group.size(); ++i) {
        fout << "------------------------------------------" << endl;
        for (int p = 0; p < frame_group[i].particles.size();  ++p, ++p_id) {
            fout << "Cam " << i << ", " << p << ": " << endl;
            for (int c = 0; c < frame_group.size(); ++c) {
                fout << "\tC" << c << " -- ";
                for (int d = 0; d < NUM_MATCHES; ++d) {
                    int index = current_matchmap[p_id][c][d];
                    if (index >=0 )
                        fout << "[" << index << ", " << frame_group[c].particles[index]->id << "]  ";
                    else
                        fout << index << "  ";
                }
                fout << endl;
            }
        }
    }
    fout.close();
}

void Correspondence::printMatchMap(const lpt::ImageFrameGroup& frame_group, const lpt::MatchMap& match_map, string output_file_name) const
{
	cout << "printing matchmap..." << endl;
    ofstream fout(output_file_name.c_str());
    int p_id = 0;
    for(int i = 0; i < frame_group.size(); ++i) {
        fout << "------------------------------------------" << endl;
        for (int p = 0; p < frame_group[i].particles.size();  ++p, ++p_id) {
            fout << "Cam " << i << ", " << p << ": " << endl;
            for (int c = 0; c < frame_group.size(); ++c) {
                fout << "\tC" << c << " -- ";
                for (int d = 0; d < NUM_MATCHES; ++d) {
                    int index = match_map[p_id][c][d];
                    if (index >=0 )
                        fout << "[" << index << ", " << frame_group[c].particles[index]->id << "]  ";
                    else
                        fout << index << "  ";
                }
                fout << endl;
            }
        }
    }
    fout.close();
}

void Correspondence::resetMatchMap(lpt::MatchMap& matchmap)
{
	for (int a = 0; a < matchmap.size(); ++a) { 
        for (int b = 0; b < matchmap[b].size(); ++b) {
            for (int c = 0; c < matchmap[a][b].size(); ++c) {
				matchmap[a][b][c] = -1;
            }
        }
	}
}

void Correspondence::initializeMatchMap()
{
	auto& cameras = shared_objects->cameras;
    init_matchmap.resize( cameras.size() * this->initial_max_particles_per_image );
	array<int, NUM_MATCHES> temp;
	for (int z = 0; z < temp.size(); ++z)
		temp[z] = -1;
	for (int a = 0; a < init_matchmap.size(); ++a)  
		init_matchmap[a].resize(cameras.size(), temp);
	current_matchmap = init_matchmap;
	lpt::ImageFrameGroup init_framegroup(cameras.size());
	auto pair = std::make_pair(init_framegroup, init_matchmap);
	matchmap_storage.resize(map_storage_size, pair);

	empty_maps_queue.setCapacity(map_storage_size);
	full_maps_queue.setCapacity(map_storage_size);

	for (auto iter = matchmap_storage.begin(); iter != matchmap_storage.end(); ++iter)
		empty_maps_queue.push(iter);
}

PointMatcher::PointMatcher()
{
    cout << "Epipolar Point matcher constructed" << endl;
}

PointMatcher::~PointMatcher() {}

void PointMatcher::initialize()
{
    this->initializeMatchMap();
}

void PointMatcher::addControls()
{
    void* matcher_void_ptr = static_cast<void*> (this);
    cv::createTrackbar("Match Threshold", string(), &params.match_thresh_level, 100, callbackMatchThresh, matcher_void_ptr);
}

void PointMatcher::findEpipolarMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap )
{
	auto& camera_pairs = shared_objects->camera_pairs;
	size_t match_overload = 0;
	vector<size_t> num_particles(frame_group.size());
	num_particles[0] = frame_group[0].particles.size();
	for(int i = 1; i < frame_group.size(); ++i) 
		num_particles[i] = frame_group[i].particles.size() + num_particles[i-1];

	for (vector<lpt::CameraPair>::iterator pair = camera_pairs.begin(); pair != camera_pairs.end(); ++pair){

		int cam_a = pair->cam_A.id;
		int cam_b = pair->cam_B.id;
		const vector<ParticleImage::Ptr>& cam_A_particles = frame_group[cam_a].particles;
		const vector<ParticleImage::Ptr>& cam_B_particles = frame_group[cam_b].particles;	

	
		size_t b_start = (cam_b !=0 ? num_particles[cam_b-1] : 0);
		size_t b_end = num_particles[cam_b];

		size_t a_start = (cam_a !=0 ? num_particles[cam_a-1] : 0);
		size_t a_end = num_particles[cam_a];

		for (size_t b_id = b_start; b_id < b_end; ++b_id) {
			double lineA[3];
			calculateEpiline(*cam_B_particles[b_id - b_start], pair->F, lineA);

			for (size_t a_id = a_start; a_id < a_end; ++a_id) {
				double residual = calculateEpipolarResidual(lineA, *cam_A_particles[a_id - a_start]);

				if (residual <= params.match_threshold){
					auto itb = std::find(matchmap[b_id][cam_a].begin(), matchmap[b_id][cam_a].end(), -1);
					auto ita = std::find(matchmap[a_id][cam_b].begin(), matchmap[a_id][cam_b].end(), -1);
					
					if (itb != matchmap[b_id][cam_a].end() && ita != matchmap[a_id][cam_b].end() ) {
						*itb = static_cast<int>(a_id - a_start);
						*ita = static_cast<int>(b_id - b_start);
					} else
						match_overload++;
				}
			}
		}
	}
    //cout << match_overload << endl;
}

void PointMatcher::findUniqueMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches)
{

   vector<size_t> num_particles(frame_group.size());
   num_particles[0] = frame_group[0].particles.size();
   for(int i = 1; i < frame_group.size(); ++i)
      num_particles[i] = frame_group[i].particles.size() + num_particles[i-1];

   size_t num_cameras = frame_group.size();
   matches.clear();

   for (size_t cam_a = 0; cam_a < num_cameras - 3; ++cam_a)
   {
      size_t a_start = (cam_a !=0 ? num_particles[cam_a - 1] : 0);
      for (int a = 0; a < frame_group[cam_a].particles.size(); ++a)
      {
         lpt::ParticleImage::Ptr Pa = frame_group[cam_a].particles[a];
         if( ! Pa->is_4way_matched )
         {
            for (size_t cam_b = cam_a + 1; cam_b < num_cameras - 2; ++cam_b)
            {
               size_t b_start = (cam_b !=0 ? num_particles[cam_b-1] : 0);
               for(size_t match_ab = 0; match_ab < NUM_MATCHES; ++match_ab) { //loop through all A,B matches
                  int b = matchmap[a + a_start][cam_b][match_ab];
                  if (b < 0)
                     break;
                  lpt::ParticleImage::Ptr Pb = frame_group[cam_b].particles[b];

                  if( ! Pb->is_4way_matched )
                  {
                     for (size_t cam_c = cam_b + 1; cam_c < num_cameras - 1; ++cam_c)
                     {
                        size_t c_start = (cam_c !=0 ? num_particles[cam_c-1] : 0);
                        for (int match_bc = 0; match_bc < NUM_MATCHES; ++match_bc)
                        {
                           int c = matchmap[b + b_start][cam_c][match_bc];
                           if (c < 0)
                              break;

                           lpt::ParticleImage::Ptr Pc = frame_group[cam_c].particles[c];

                           if( ! Pc->is_4way_matched && std::count(matchmap[a + a_start][cam_c].begin(), matchmap[a + a_start][cam_c].end(), c) )
                           {
                              for (size_t cam_d = cam_c + 1; cam_d < num_cameras; ++cam_d)
                              {
                                 size_t d_start = (cam_d !=0 ? num_particles[cam_d-1] : 0);
                                 for (int match_cd = 0; match_cd < NUM_MATCHES; ++match_cd)
                                 {
                                    int d = matchmap[c + c_start][cam_d][match_cd];
                                    if (d < 0)
                                       break;
                                    lpt::ParticleImage::Ptr Pd = frame_group[cam_d].particles[d];
                                    if( ! Pd->is_4way_matched && std::count(matchmap[a + a_start][cam_d].begin(), matchmap[a + a_start][cam_d].end(), d)  && std::count(matchmap[b + b_start][cam_d].begin(), matchmap[b+b_start][cam_d].end(), d)  )
                                    {
                                       lpt::Match::Ptr newmatch = lpt::Match::create();
                                       newmatch->addParticle(Pa,cam_a);
                                       newmatch->addParticle(Pb,cam_b);
                                       newmatch->addParticle(Pc,cam_c);
                                       newmatch->addParticle(Pd,cam_d);
                                       matches.push_back(std::move(newmatch));
                                       Pa->is_4way_matched = true;
                                       Pb->is_4way_matched = true;
                                       Pc->is_4way_matched = true;
                                       Pd->is_4way_matched = true;
                                       match_ab = NUM_MATCHES;
                                       match_bc = NUM_MATCHES;
                                       match_cd = NUM_MATCHES;
                                       cam_b = num_cameras;
                                       cam_c = num_cameras;
                                       cam_d = num_cameras;
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

void PointMatcher::find3WayMatches(const lpt::ImageFrameGroup& frame_group, lpt::MatchMap& matchmap, vector<lpt::Match::Ptr>& matches) {
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
}

void PointMatcher::removeNonUniqueMatches(vector<Match::Ptr>& matches)
{
	for (int i = 0; i < matches.size(); ++i) {
		if (matches[i]->isUnique() == false) {
			matches[i].reset();
		}
	}
}

//void PointMatcher::find3wayMatches(tuple<int,int,int>& triple, lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches) {
//	int camA = std::get<0>(triple);
//	int camB = std::get<1>(triple);
//	int camC = std::get<2>(triple);
//	vector<lpt::Match::Ptr> matches3way;
//	for (int i = 0; i < cameragroup[ camA ].particles.size(); ++i) {
//		lpt::ParticleImage* Pa = cameragroup[camA].particles[i].get();
//		set<lpt::ParticleImage*>& matchesAB = Pa->matches2way[camB];
//		set<lpt::ParticleImage*>::iterator matchAB;
//		if( ! Pa->is_3way_matched ) 
//		for (matchAB = matchesAB.begin(); matchAB != matchesAB.end(); ++matchAB){
//			ParticleImage* Pb = (*matchAB);
//			set<ParticleImage*>& matchesBC = Pb->matches2way[camC];
//			set<ParticleImage*>& matchesAC = Pa->matches2way[camC];
//			set<ParticleImage*>::iterator matchAC;
//			if( ! Pb->is_3way_matched ) 
//			for (matchAC = matchesAC.begin(); matchAC != matchesAC.end(); ++matchAC) {
//				ParticleImage* Pc = (*matchAC);
//				if( ! Pc->is_3way_matched ) 
//				if( matchesBC.count(Pc) ) { 
//					lpt::Match::Ptr newmatch = lpt::Match::create();
//					newmatch->addParticle(Pa,camA);
//					newmatch->addParticle(Pb,camB);
//					newmatch->addParticle(Pc,camC);
//					matches3way.push_back(newmatch);
//					Pa->is_3way_matched = true;
//					Pb->is_3way_matched = true;
//					Pc->is_3way_matched = true;
//					matchAB = matchesAB.end();
//					matchAC = matchesAC.end();
//				}
//			}
//		}
//	}
//	removeNonUniqueMatches( matches3way );
//	std::move(matches3way.begin(), matches3way.end(), std::back_inserter(matches) );
//}

//void PointMatcher::find4WayMatches(tuple<int,int,int,int>& quad, lpt::ImageFrameGroup& cameragroup, vector<lpt::Match::Ptr>& matches) {
//	int camA = std::get<0>(quad);
//	int camB = std::get<1>(quad);
//	int camC = std::get<2>(quad);
//	int camD = std::get<3>(quad);
//	vector<lpt::Match::Ptr> matches4way;
//	for (int i = 0; i < cameragroup[ camA ].particles.size(); ++i) {
//		lpt::ParticleImage* Pa = cameragroup[camA].particles[i].get();
//		set<lpt::ParticleImage*>& matchesAB = Pa->matches2way[camB];
//		set<lpt::ParticleImage*>::iterator matchAB;
//		if( ! Pa->is_4way_matched ) 
//		for(matchAB = matchesAB.begin(); matchAB != matchesAB.end(); ++matchAB) { //loop through all A,B matches
//			lpt::ParticleImage* Pb = (*matchAB);
//			set<lpt::ParticleImage*>& matchesBC = Pb->matches2way[camC];
//			set<lpt::ParticleImage*>::iterator matchBC;
//			if( ! Pb->is_4way_matched ) 
//			for (matchBC = matchesBC.begin(); matchBC != matchesBC.end(); ++matchBC) {
//				lpt::ParticleImage* Pc = (*matchBC);
//				set<lpt::ParticleImage*>& matchesCD = Pc->matches2way[camD];
//				set<lpt::ParticleImage*>::iterator matchCD;
//				if( ! Pc->is_4way_matched ) 
//				for (matchCD = matchesCD.begin(); matchCD != matchesCD.end(); ++matchCD) {
//					lpt::ParticleImage* Pd = (*matchCD);
//					if( ! Pd->is_4way_matched ) 
//					if(Pd->matches2way[camA].count(Pa) && Pb->matches2way[camD].count(Pd) && Pa->matches2way[camC].count(Pc)) {
//						if(! Pa->is_4way_matched && ! Pb->is_4way_matched && ! Pc->is_4way_matched && ! Pd->is_4way_matched) { 
//							lpt::Match::Ptr newmatch = lpt::Match::create();
//							newmatch->addParticle(Pa,camA);
//							newmatch->addParticle(Pb,camB);
//							newmatch->addParticle(Pc,camC);
//							newmatch->addParticle(Pd,camD);
//							matches4way.push_back(std::move(newmatch));
//							Pa->is_4way_matched = true;
//							Pb->is_4way_matched = true;
//							Pc->is_4way_matched = true;
//							Pd->is_4way_matched = true;
//							matchAB = matchesAB.end();
//							matchBC = matchesBC.end();
//							matchCD = matchesCD.end();
//						}
//					} 
//				}
//			}
//		}
//	}
//	removeNonUniqueMatches( matches4way );
//	std::move(matches4way.begin(), matches4way.end(), std::back_inserter(matches) );
//}

Reconstruct3D::Reconstruct3D() : frame_count(0) {}

Reconstruct3D::~Reconstruct3D() 
{
	for (size_t i=0; i<positions.size(); i++) {
		double x = positions[i][0] / frame_count;
		double y = positions[i][1] / frame_count;
		double z = positions[i][2] / frame_count;
		axis << i << "\t" << x << "\t" << y << "\t" << z << endl;
	}
	cout << "Reconstruct3D destructed" << endl;
	axis.close();
}

void Reconstruct3D::reconstruct3DFrame(vector<lpt::Match::Ptr>& matches, lpt::Frame3d& frame)
{
	auto& cameras = shared_objects->cameras;
	for (int m = 0; m < matches.size(); ++m) {
		if ( matches[m] ) {
			size_t number_of_cams = matches[m]->particles.size();
			if ( number_of_cams >= 2 ) {
				double empty[3] = {0};
				double empty1 = 0;
				double** A = new double* [ 2 * number_of_cams ];
				double* B = new double [ 2 * number_of_cams ];
				double* X = new double[3];

				int id = -1;
				for(int i = 0; i < number_of_cams; ++i) {
					int s = i * 2;
					int e = s + 1;
					A[s] = new double[3];
					A[e] = new double[3];
					lpt::ParticleImage::Ptr P = matches[m]->particles[i].first;
                    size_t camID = matches[m]->particles[i].second;
					double x = ( P->x - cameras[camID].c[0] ) / (1. * cameras[camID].f[0]);    // FIXME: May need -1 multiplier for pixel coordinate system (upper left corner origin)
					double y = ( P->y - cameras[camID].c[1] ) / (1. * cameras[camID].f[1]);    // Convert P.x and P.y to normalized coordinates through intrinsic parameters
					id = P->id;
					A[s][0] = x * cameras[camID].R[2][0] - cameras[camID].R[0][0];   
					A[s][1] = x * cameras[camID].R[2][1] - cameras[camID].R[0][1];
					A[s][2] = x * cameras[camID].R[2][2] - cameras[camID].R[0][2];
					A[e][0] = y * cameras[camID].R[2][0] - cameras[camID].R[1][0];
					A[e][1] = y * cameras[camID].R[2][1] - cameras[camID].R[1][1];
					A[e][2] = y * cameras[camID].R[2][2] - cameras[camID].R[1][2];

					B[s] = cameras[camID].T[0] - x * cameras[camID].T[2];
					B[e] = cameras[camID].T[1] - y * cameras[camID].T[2];
				}

				solver.Householder(A, B, number_of_cams * 2, 3);   //Transform A into upper triangular form
				solver.Backsub(A, B, 3, X);
				//array<lpt::Particle3d::float_type, lpt::Particle3d::dim> coords = {{X[0], X[1], X[2]}};
				lpt::Particle3d_Ptr newparticle = lpt::Particle3d::create();
				newparticle->id = id;
				if (this->shared_objects->isRotation_Correction) {
					auto& S = this->shared_objects->S;
					auto& P = this->shared_objects->P;
					newparticle->X[0] = S[0][0] * X[0] + S[0][1] * X[1] + S[0][2] * X[2] - P[0];
					newparticle->X[1] = S[1][0] * X[0] + S[1][1] * X[1] + S[1][2] * X[2] - P[1];
					newparticle->X[2] = S[2][0] * X[0] + S[2][1] * X[1] + S[2][2] * X[2] - P[2];
				}
				else {
					newparticle->X[0] = X[0];
					newparticle->X[1] = X[1];
					newparticle->X[2] = X[2];
				}
				newparticle->frame_index = frame.frame_index; 

				frame.objects.push_back(std::move(newparticle));

				for(int i = 0; i < number_of_cams * 2; ++i)
					delete [] A[i];
				delete [] A;
				delete [] B;
				delete [] X;
			}
		}
	}
}

struct myclass{
	bool operator() (lpt::Particle3d_Ptr p1, lpt::Particle3d_Ptr p2)
	{
		return (p1->X[1] < p2->X[1]);
	}
} mycomp;

void Reconstruct3D::draw(lpt::Frame3d& frame)
{
	auto& cameras = shared_objects->cameras;
	vector<cv::Point3d> object_points(frame.objects.size());
	vector<cv::Point2d> image_points(frame.objects.size());

	auto particles = frame.objects;

	if (positions.size() != frame.objects.size() )
		positions.resize(frame.objects.size());

	std::sort(particles.begin(), particles.end(), mycomp);

	for (int p = 0; p < frame.objects.size(); ++p) {
		object_points[p].x = frame.objects[p]->X[0];
		object_points[p].y = frame.objects[p]->X[1];
		object_points[p].z = frame.objects[p]->X[2];

		for (int i=0; i<frame.objects[p]->X.size(); i++) {
			positions[p][i] += particles[p]->X[i];
		}
	}

	frame_count++;
	
	for (int i = 0; i < frame.camera_frames.size(); ++i) {
		//if (frame.camera_frames[i].image.channels() == 1)
		//	cv::cvtColor(frame.camera_frames[i].image, frame.camera_frames[i].image, CV_GRAY2BGR);
		cv::Mat R = cv::Mat(3, 3, CV_64F, cameras[i].R);
		cv::Mat t_vec = cv::Mat(3, 1, CV_64F, cameras[i].T);
		cv::Mat r_vec = cv::Mat::zeros(3,1, CV_64F);
		cv::Rodrigues(R, r_vec);

		cv::projectPoints(cv::Mat(object_points), r_vec, t_vec, cameras[i].getCameraMatrix(), cameras[i].getDistCoeffs(), image_points);

		for (int p = 0; p < image_points.size(); ++p) 
			cv::circle( frame.camera_frames[i].image, image_points[p], 3, 255, -1 );	
		
	}
}

Reconstruct3DwithSVD::Reconstruct3DwithSVD() : Reconstruct3D() {}

Reconstruct3DwithSVD::~Reconstruct3DwithSVD() {}

void Reconstruct3DwithSVD::reconstruct3DFrame(vector<lpt::Match::Ptr>& matches, lpt::Frame3d& frame)
{
	auto& cameras = shared_objects->cameras;
	for (int m = 0; m < matches.size(); ++m) {
		if ( matches[m] ) {
            size_t number_of_cams = matches[m]->particles.size();
			if ( number_of_cams >= 2 ) {
				cv::Mat A = cv::Mat::zeros( 2 * static_cast<int>( cameras.size() ), 3, CV_64F );
				cv::Mat B = cv::Mat::zeros( 2 * static_cast<int>( cameras.size() ), 1, CV_64F );
				cv::Mat X = cv::Mat::zeros( 3, 1, CV_64F );
				lpt::Particle3d_Ptr newparticle = lpt::Particle3d::create();

				int id = -1;
				for(int i = 0; i < number_of_cams; ++i) {
					int s = i * 2;
					int e = s + 1;
				
					lpt::ParticleImage::Ptr P = matches[m]->particles[i].first;
					size_t cam_id = matches[m]->particles[i].second;              
					newparticle->camera_ids[i] = cam_id;								//WARNING this will cause a segfault if number_of_cameras is greater than four;
					id = P->id;

					double x = ( P->x - cameras[cam_id].c[0] ) / cameras[cam_id].f[0];
					double y = ( P->y - cameras[cam_id].c[1] ) / cameras[cam_id].f[1];

					A.at<double>(s, 0) = x * cameras[cam_id].R[2][0] - cameras[cam_id].R[0][0];   
					A.at<double>(s, 1) = x * cameras[cam_id].R[2][1] - cameras[cam_id].R[0][1];
					A.at<double>(s, 2) = x * cameras[cam_id].R[2][2] - cameras[cam_id].R[0][2];

					A.at<double>(e, 0) = y * cameras[cam_id].R[2][0] - cameras[cam_id].R[1][0];
					A.at<double>(e, 1) = y * cameras[cam_id].R[2][1] - cameras[cam_id].R[1][1];
					A.at<double>(e, 2) = y * cameras[cam_id].R[2][2] - cameras[cam_id].R[1][2];

					B.at<double>(s) = cameras[cam_id].T[0] - x * cameras[cam_id].T[2];
					B.at<double>(e) = cameras[cam_id].T[1] - y * cameras[cam_id].T[2];
				}

				svd(A);
				svd.backSubst(B,X);
						
				newparticle->id = id;
				newparticle->X[0] = X.at<double>(0);
				newparticle->X[1] = X.at<double>(1);
				newparticle->X[2] = X.at<double>(2);
				newparticle->frame_index = frame.frame_index; 
				
				double max_svalue = 0, min_svalue = 0; 
				cv::minMaxIdx(svd.w, &min_svalue, &max_svalue);
				double cond_A = max_svalue / min_svalue;
				double normL2_A = max_svalue;
				double normL2_B = cv::norm(B);
				double normL2_X = cv::norm(X);

				cv::Mat AX = A * X;
				cv::Mat R = B - AX;
				double normL2_R = cv::norm(R); 
				double theta =  atan(normL2_R / cv::norm(AX) );
				
				
				//for(int i = 0; i < cameras.size()+10000; ++i) {
				//	
				//	/*double x_u2 = 0; 
				//	double y_u2 = 0;

				//	lpt::ParticleImage* P = matches[m]->particles[i].first;
				//	int cam_id = matches[m]->particles[i].second;
				//	
				//	lpt::Camera& cam = cameras[cam_id];
				//	cv::Mat camera_matrix = cam.getCameraMatrix();
				//	cv::Mat rotation_matrix(3,3, CV_64F, cam.R); 		
				//	cv::Mat translation_vec(3,1, CV_64F, cam.T);
				//	cv::Mat dist_coeffs = cam.getDistCoeffs();

				//	cv::Mat rotation_vec;
				//	cv::Rodrigues(rotation_matrix, rotation_vec);

				//	cv::Point3d point;
				//	point.x = newparticle->X[0];
				//	point.x = newparticle->X[1];
				//	point.x = newparticle->X[2];

				//	vector<cv::Point3d> point3D; 
				//	point3D.push_back(point);
				//	vector<cv::Point2d> image_point_distorted;
				//	vector<cv::Point2d> image_point;*/
				//	/*cv::projectPoints( cv::Mat(point3D), rotation_vec, translation_vec, camera_matrix, dist_coeffs, image_point_distorted );
				//	cv::undistortPoints( image_point_distorted, image_point, camera_matrix, dist_coeffs, cv::Mat(), camera_matrix);

				//	for(int dist_id = 0; dist_id < 4; ++dist_id) {

				//		cv::Mat dist_coeffs_perturbed = dist_coeffs.clone();
				//		dist_coeffs_perturbed.at<double>(dist_id) += cam.dist_coeffs_u[dist_id];

				//		vector<cv::Point2d> image_point_perturbed;

				//		cv::undistortPoints( image_point_distorted, image_point_perturbed, camera_matrix, dist_coeffs_perturbed, cv::Mat(), camera_matrix);

				//		double x_dxprime = 0;
				//		double y_dyprime = 0;

				//		if ( cam.dist_coeffs_u[dist_id] != 0 ) { 
				//			x_dxprime = ( image_point[0].x - image_point_perturbed[0].x ); 
				//			y_dyprime = ( image_point[0].y - image_point_perturbed[0].y ); 
				//		}

				//		x_u2 += x_dxprime * x_dxprime;
				//		y_u2 += y_dyprime * y_dyprime;

				//	}*/

				//	//x_u2 += cameras[cam_id].centriod_loc_uncertainty * cameras[cam_id].centriod_loc_uncertainty; // 0.25 = 0.5 * 0.5 = squared uncertainty of centriod pixel coordinate
				//	//y_u2 += cameras[cam_id].centriod_loc_uncertainty * cameras[cam_id].centriod_loc_uncertainty; // 0.25 = 0.5 * 0.5 = squared uncertainty of centriod pixel coordinate

				//	/*double x = ( P->x - cameras[cam_id].c[0] ) / cameras[cam_id].f[0];
				//	double y = ( P->y - cameras[cam_id].c[1] ) / cameras[cam_id].f[1];

				//	double dximdx = 1.0 / cameras[cam_id].f[0];
				//	double dximdc = -1.0 / cameras[cam_id].f[0];
				//	double dximdf = ( cameras[cam_id].c[0] - P->x ) / ( cameras[cam_id].f[0] * cameras[cam_id].f[0] );

				//	cameras[cam_id].X_u[0] = 
				//		sqrt( 
				//		dximdx * dximdx * x_u2
				//		+ dximdc * dximdc * cameras[cam_id].c_u[0] * cameras[cam_id].c_u[0]
				//		+ dximdf * dximdf * cameras[cam_id].f_u[0] * cameras[cam_id].f_u[0]
				//	);

				//	double dyimdy =  1.0 / cameras[cam_id].f[1];
				//	double dyimdc = -1.0 / cameras[cam_id].f[1];
				//	double dyimdf = ( cameras[cam_id].c[1] - P->y ) / ( cameras[cam_id].f[1] * cameras[cam_id].f[1] );

				//	cameras[cam_id].X_u[1] = 
				//		sqrt( 
				//		dyimdy * dyimdy * y_u2 
				//		+ dyimdc * dyimdc * cameras[cam_id].c_u[1] * cameras[cam_id].c_u[1]
				//		+ dyimdf * dyimdf * cameras[cam_id].f_u[1] * cameras[cam_id].f_u[1]
				//		);*/

				//	//// 3D position uncertainty due to x' and y' image coordinates
				//	/*if ( cameras[cam_id].X_u[0] != 0 && cameras[cam_id].X_u[1] != 0) {
				//		for (int dim = 0; dim < 2; ++dim) {

				//			cv::Mat delta_B = cv::Mat::zeros( 2 * cameras.size(), 1, CV_64F );
				//			cv::Mat E = cv::Mat::zeros( 2 * cameras.size(), 3, CV_64F );

				//			int s = (dim == 0 ? cam_id * 2 : cam_id * 2 + 1);

				//			E.at<double>(s, 0) = cameras[cam_id].X_u[dim] * cameras[cam_id].R[2][0];   
				//			E.at<double>(s, 1) = cameras[cam_id].X_u[dim] * cameras[cam_id].R[2][1];
				//			E.at<double>(s, 2) = cameras[cam_id].X_u[dim] * cameras[cam_id].R[2][2];
				//			delta_B.at<double>(s) = -1.0 * cameras[cam_id].X_u[dim] * cameras[cam_id].T[2];

				//			svd(E);
				//			double max_s = 0, min_s = 0; 
				//			cv::minMaxIdx(svd.w, &min_s, &max_s);
				//			double normL2_E = max_s;
				//			double normL2_deltaB = cv::norm(delta_B);
				//			double E_sen = (cond_A * cond_A * tan(theta) + cond_A ) * normL2_E / normL2_A * normL2_X / cameras[cam_id].X_u[dim];
				//			double deltaB_sen = cond_A / cos(theta) * normL2_deltaB / normL2_B * normL2_X / cameras[cam_id].X_u[dim];
				//			double squared_uncertainty_homogenious_coordinate =	
				//				( E_sen * E_sen + deltaB_sen * deltaB_sen) * cameras[cam_id].X_u[dim] * cameras[cam_id].X_u[dim];

				//			newparticle->uncertainty += squared_uncertainty_homogenious_coordinate;
				//		}
				//	}*/
				//	//// 3D position uncertainty to translation vector T
				//	//for (int t = 0; t < 3; ++t) {
				//	//	double T_u[3] = {0};
				//	//	T_u[t] = cameras[cam_id].T_u[t];

				//	//	cv::Mat delta_B = cv::Mat::zeros( 2 * cameras.size(), 1, CV_64F );

				//	//	int s = cam_id * 2;
				//	//	int e = s + 1;

				//	//	delta_B.at<double>(s) = T_u[0] - x * T_u[2];
				//	//	delta_B.at<double>(e) = T_u[1] - y * T_u[2];

				//	//	double normL2_deltaB = cv::norm(delta_B);
				//	//	double deltaB_sen = 0;
				//	//	if ( T_u[t] != 0)
				//	//		deltaB_sen = cond_A / cos(theta) * normL2_deltaB / normL2_B * normL2_X / cameras[cam_id].T_u[t];
				//	//	double squared_uncertainty_tvec =	deltaB_sen * deltaB_sen * T_u[t] * T_u[t];

				//	//	newparticle->uncertainty += squared_uncertainty_tvec;		
				//	//}

				//	//// 3D position uncertainty due to Rotation vector r_vec


				//	//cv::Mat A_u = cv::Mat::zeros( 2 * cameras.size(), 3, CV_64F );
				//	//cv::Mat B_u = cv::Mat::zeros( 2 * cameras.size(), 1, CV_64F );
				//	//cv::Mat X_u = cv::Mat::zeros( 3, 1, CV_64F );
				//	////cout << "\n" << i << "\t";
				//	//for(int c = 0; c < cameras.size(); ++c) {
				//	//	lpt::ParticleImage* Pc = matches[m]->particles[c].first;
				//	//	int cam_id = matches[m]->particles[c].second;

				//	//	double xc = ( Pc->x - cameras[c].c[0] ) / cameras[c].f[0];
				//	//	double yc = ( Pc->y - cameras[c].c[1] ) / cameras[c].f[1];

				//	//	cv::Mat R_u = cv::Mat(3, 3, CV_64F, cameras[c].R).clone(); 		
				//	//	cv::Mat r_vec(3, 1, CV_64F);
				//	//	cv::Rodrigues(R_u, r_vec);

				//	//	r_vec.at<double>(0) += cameras[c].r_vec_u[0];
				//	//	r_vec.at<double>(1) += cameras[c].r_vec_u[1];
				//	//	r_vec.at<double>(2) += cameras[c].r_vec_u[2];

				//	//	cv::Rodrigues(r_vec, R_u);

				//	//	int s = c * 2;
				//	//	int e = s + 1;

				//	//	A_u.at<double>(s, 0) = xc * R_u.at<double>(2,0) - R_u.at<double>(0,0);   
				//	//	A_u.at<double>(s, 1) = xc * R_u.at<double>(2,1) - R_u.at<double>(0,1);
				//	//	A_u.at<double>(s, 2) = xc * R_u.at<double>(2,2) - R_u.at<double>(0,2);

				//	//	A_u.at<double>(e, 0) = yc * R_u.at<double>(2,0) - R_u.at<double>(1,0);
				//	//	A_u.at<double>(e, 1) = yc * R_u.at<double>(2,1) - R_u.at<double>(1,1);
				//	//	A_u.at<double>(e, 2) = yc * R_u.at<double>(2,2) - R_u.at<double>(1,2);

				//	//	B_u.at<double>(s) = cameras[c].T[0] - xc * cameras[c].T[2];
				//	//	B_u.at<double>(e) = cameras[c].T[1] - yc * cameras[c].T[2];
				//	//}

				//	//svd(A_u);
				//	//svd.backSubst(B_u, X_u);
				//	//double dX = 
				//	//	sqrt(
				//	//	(newparticle->X[0] - X_u.at<double>(0)) * (newparticle->X[0] - X_u.at<double>(0)) +
				//	//	(newparticle->X[1] - X_u.at<double>(1)) * (newparticle->X[1] - X_u.at<double>(1)) +
				//	//	(newparticle->X[2] - X_u.at<double>(2)) * (newparticle->X[2] - X_u.at<double>(2)) 
				//	//	);

				//	//double squared_uncertainty_rvec = dX*dX;

				//	//newparticle->uncertainty += squared_uncertainty_rvec;					

				//}
				//newparticle->uncertainty = ( newparticle->uncertainty > 0 ? sqrt(newparticle->uncertainty) : 0 );
				
				frame.objects.push_back(std::move(newparticle));
		}
	}
}
}


////------OLD CODE---------------------------------------------------------------------------
//void SpatialMatcher::findMatches() {
//	cout << "Starting correspondence search: Total Frames = " << (int)cameras[0].frames.size() <<  endl;
//    globalmatches.resize( (int)cameras[0].frames.size() );
//    for (int frame_index = 0; frame_index < (int)cameras[0].frames.size(); ++frame_index) {
////    	cout << "Frame " << frame_index << ":" << endl;
////    	cout << "\t Finding two way matches: ";
//    	for (vector<CameraPair>::iterator pair = stereo_pairs.begin(); pair != stereo_pairs.end(); ++pair){
//        	find2WayMatches(frame_index, *pair);
//        }
////        cout<<"\t--Complete" <<endl << "\t Finding multicam matches:";
//
//        for (int camA = 0; camA < cameras.size() - 1; ++camA) {
//            // printf("%d,", cameras_[camA]->ID_);
//        	bool matched4way = false;
//            for (int i = 0; i < cameras[camA].frames[frame_index].particles.size(); ++i) {
//                ParticleImage* Pa = cameras[camA].frames[frame_index].particles[i].get();
//                for (int camB = camA + 1; camB < cameras.size() - 1; ++camB) {
//                    // printf("%d,", cameras_[camB]->ID_);
//                    set<ParticleImage*>& matchesAB = Pa->matches2way[camB];
//                    set<ParticleImage*>::iterator matchAB;
//                    for(matchAB = matchesAB.begin(); matchAB != matchesAB.end(); ++matchAB) { //loop through all A,B matches
//                        ParticleImage* Pb = (*matchAB);
//
//                        for (int camC = camB + 1; camC < cameras.size(); ++camC) {
//                            // printf("%d\n", cameras_[camC]->ID_);
//                            set<ParticleImage*> &matchesBC = Pb->matches2way[camC];
//                            set<ParticleImage*>::iterator matchBC;
//
//                            for (matchBC = matchesBC.begin(); matchBC != matchesBC.end(); ++matchBC) {
//                                ParticleImage* Pc = (*matchBC);
//                                for (int camD = camC + 1; camD < cameras.size(); ++camD) {
//                                    // printf("%d\n", cameras_[camD]->ID_);
//                                    set<ParticleImage*> &matchesCD = Pc->matches2way[camD];
//                                    set<ParticleImage*>::iterator matchCD;
//                                    for (matchCD = matchesCD.begin(); matchCD != matchesCD.end(); ++matchCD) {
//                                        ParticleImage* Pd = (*matchCD);
//                                        if(Pd->matches2way[camA].count(Pa) && Pb->matches2way[camD].count(Pd) && Pa->matches2way[camC].count(Pc)) {
//                                           //printf("\t SuperMATCH!!!!");
//                                            Match* newmatch = new Match();
//                                            newmatch->addParticle(Pa,camA);
//                                            newmatch->addParticle(Pb,camB);
//                                            newmatch->addParticle(Pc,camC);
//                                            newmatch->addParticle(Pd,camD);
//                                            globalmatches[frame_index].push_back(newmatch);
////                                            printf("4Matching: %d-%d, %d-%d, %d-%d, %d-%d\n", camA, Pa->id, camB, Pb->id, camC, Pc->id, camD, Pd->id);
//                                            matched4way = true;
//                                        } //TODO: Add code to identify simple 4-way matches 1234, 1243 and 1324.
//                                    }
//                                }
//                                //TODO: add code to store the match residual in the Match object
//                                if (matched4way == false)
//
//                                    if( Pc->matches2way[camA].count(Pa) ) {
//                                    	;
////                                        Match* newmatch = new Match();
////                                        newmatch->addParticle(Pa,camA);
////                                        newmatch->addParticle(Pb,camB);
////                                        newmatch->addParticle(Pc,camC);
////                                        globalmatches[frame_index].push_back(newmatch);
////                                        printf("3Matching: %d-%d, %d-%d, %d-%d\n", camA, Pa->id, camB, Pb->id, camC, Pc->id);
//                                    }
//                            }
//                        }
//                    }
//                }
//            }
//        }
////        cout<< "\t--Complete"<< endl << "\t Removing non unique matches";
//        removeNonUniqueMatches( globalmatches[frame_index] );
////        cout<< "\t--Complete"<<endl<<endl;
//    }
//}
//
///*void SpatialMatcher::InitializeMatchMap(int camID, ParticleImage* P) {
//
//    for (int c = 0; c < number_of_cameras_; ++c) {
//        matchmaps_[c];
//    }
//}*/
//
//void SpatialMatcher::find2WayMatches(int frame_index, CameraPair &pair) {
//    vector<ParticleImage::Ptr>& cam_A_particles = pair.cam_A.frames[frame_index].particles;
//    vector<ParticleImage::Ptr>& cam_B_particles = pair.cam_B.frames[frame_index].particles;
//
//    for (int b = 0; b < (int)cam_B_particles.size(); ++b) {
//        double lineA[3];
//        calculateEpiline(cam_B_particles[b], pair.F, lineA);
//        for (int a = 0; a < (int)cam_A_particles.size(); ++a) {
//            double residual = calculateEpipolarResidual(lineA, cam_A_particles[a]);
////            if (i == j){
////                printf("Should match C%d_%d and C%d_%d, r = %f\n", pair.cam_A_id, cam_A_particles[a]->id, pair.cam_B_id, cam_B_particles[b]->id, residual);
////                //printf("Cam%d_%d: [%E, %E]\n", A_->ID_, Cam_A_particles[a]->id_,Cam_A_particles[a]->px_, Cam_A_particles[a]->py_);
////                //printf("Cam%d_%d: [%E, %E]\n", B_->ID_, Cam_B_particles[b]->id_,Cam_B_particles[b]->px_, Cam_B_particles[b]->py_);
////            }
//            if (residual <= match_residual_threshold){
//                cam_A_particles[a]->add2WayMatch(pair.cam_B.id, cam_B_particles[b].get());
//                cam_B_particles[b]->add2WayMatch(pair.cam_A.id, cam_A_particles[a].get());
//                //matches_[ Cam_A_particles[i] ].insert( Cam_B_particles[j] );
//                //FIXME: Find a way to store the residual value and later add 3 and 4 cam matches
////                printf("matched C%d_%d and C%d_%d, r = %f\n", pair.cam_A_id, cam_A_particles[i]->id, pair.cam_B_id, cam_B_particles[j]->id, residual);
//            }
//        }
//    }
//
//}
//
//double SpatialMatcher::calculateEpipolarResidual(double line[3], ParticleImage::Ptr point)  {
//    double residual = fabs( point->x * line[0] + point->y * line[1] + line[2] );
//    return residual;
//}
//
//void SpatialMatcher::calculateEpiline(ParticleImage::Ptr point, const double F[3][3], double line[3]) {
//    double x = point->x;
//    double y = point->y;
//    double z = 1.0;
//
//    line[0] = F[0][0] * x + F[1][0] * y + F[2][0] * z;
//    line[1] = F[0][1] * x + F[1][1] * y + F[2][1] * z;
//    line[2] = F[0][2] * x + F[1][2] * y + F[2][2] * z;
//    double factor = line[0] * line[0] + line[1] * line[1];
//    factor = factor ? 1./sqrt(factor) : 1.;
//    line[0] *= factor;
//    line[1] *= factor;
//    line[2] *= factor;
//}
//
//void SpatialMatcher::removeNonUniqueMatches(vector<Match*> &matches) {
//    for (int i = 0; i < matches.size(); ++i)
//        if (matches[i]->isUnique() == false) {
//            delete matches[i];//printf(" not unique\n");
//            matches[i] = 0;
//            // TODO: Erase null matches
//        }
//}
//
//void SpatialMatcher::get3wayMatch(int camB, int camC, ParticleImage::Ptr Pa) {
//    set<ParticleImage*>::iterator matchAB;
//    for (matchAB = Pa->matches2way[camB].begin(); matchAB != Pa->matches2way[camB].end(); ++matchAB){
//        ParticleImage* Pb = (*matchAB);
//        set<ParticleImage*> &matchesBC = Pb->matches2way[camC];
//        set<ParticleImage*>::iterator matchAC;
//        for (matchAC = Pa->matches2way[camC].begin(); matchAC != Pa->matches2way[camC].end(); ++matchAC) {
//            ParticleImage* particle_matched_AC = (*matchAC);
//            if(matchesBC.count(particle_matched_AC)) //FIXME: this compares two pointers and my error if parallelized
//                printf("matched: %d,%d,%d \n", Pa->id, Pb->id, particle_matched_AC->id);
//        }
//    }
//}
//
//void SpatialMatcher::testMatches(){
//	for (int f = 0; f < globalmatches.size(); f++) {
//	    	int N_correct = 0;
//	    	int N_matched = 0;
//	    	int N_total = 100000000;
//	    	for(int c = 0; c < cameras.size(); c++) {
//	    			int temp = cameras[c].frames[f].particles.size();
//	    			if (temp < N_total)
//	    				N_total = temp;
//	    	}
//
//	    	for (int m = 0; m < globalmatches[f].size(); m++){
//	    		if ( globalmatches[f][m] ) {
//	    			N_matched++;
//	    			vector<pair < ParticleImage*, int > >& particles = globalmatches[f][m]->particles;
//	    			bool matched = true;
//	    			for (int j = 0; j < particles.size() - 1; j++){
//	    				matched &= particles[j].first->id == particles[j+1].first->id ? true : false;
//	    			}
//	    			if ( matched )
//	    				N_correct++;
//	    		}
//	    	}
//	    	cout << "Frame "<< f << ":";
//	        cout << "\t Correct Ratio: " << N_correct << " / " << N_matched << " = " <<
//	        		setprecision(2) << (N_matched > 0 ? (double) N_correct / N_matched : 0.0 );
//	        cout << "\t Cover Ratio: " << N_correct << " / " << N_total << " = " <<
//	        		setprecision(2) <<  (N_total > 0 ? (double) N_correct / N_total : 0.0 ) << endl;
//	    }
//}
//
//SpatialMatcher::~SpatialMatcher() {
//  //TODO: Auto generated destructor
//}
//
//void reconstruct3D(vector<lpt::Camera>& cameras, vector<vector<lpt::Match*> >& matches, vector<lpt::Frame*>& frames) {
//	lpt::Regression solve;
//	for (int f = 0; f < matches.size(); ++f) {
//		cout << "Frame " << f <<":  number of matches = " << matches[f].size() << endl;
//		lpt::Frame* newframe = new lpt::Frame(f);     // FIXME: Assumes frame_index is sequential starting at zero!!!
//		for (int m = 0; m < matches[f].size(); ++m) {
//			if (matches[f][m]) { //if valid pointer
//				int number_of_cams = matches[f][m]->particles.size();
//				//cout << "\t match " << m << ":  number of cameras = " << number_of_cams <<  endl;
//				double empty[3] = {0};
//				double empty1 = 0;
//				double** A = new double* [ 2 * number_of_cams ];
//				double* B = new double [ 2 * number_of_cams ];
//
//				double* X = new double[3];
//				int id = -1;
//				for(int i = 0; i < number_of_cams; ++i) {
//					int s = i * 2;
//					int e = s + 1;
//					A[s] = new double[3];
//					A[e] = new double[3];
//					lpt::ParticleImage* P = matches[f][m]->particles[i].first;
//					int camID = matches[f][m]->particles[i].second;
//					double x = ( P->x - cameras[camID].c[0] ) / (1. * cameras[camID].f[0]);    // FIXME: May need -1 multiplier for pixel coordinate system (upper left corner origin)
//					double y = ( P->y - cameras[camID].c[1] ) / (1. * cameras[camID].f[1]);
//					id = P->id;
//					A[s][0] = x * cameras[camID].R[2][0] - cameras[camID].R[0][0];   // Convert P.x and P.y to normalized coordinates through intrinsic parameters
//					A[s][1] = x * cameras[camID].R[2][1] - cameras[camID].R[0][1];
//					A[s][2] = x * cameras[camID].R[2][2] - cameras[camID].R[0][2];
//					A[e][0] = y * cameras[camID].R[2][0] - cameras[camID].R[1][0];
//					A[e][1] = y * cameras[camID].R[2][1] - cameras[camID].R[1][1];
//					A[e][2] = y * cameras[camID].R[2][2] - cameras[camID].R[1][2];
//
//					B[s] = cameras[camID].T[0] - x * cameras[camID].T[2];
//					B[e] = cameras[camID].T[1] - y * cameras[camID].T[2];
//				}
//
//				solve.Householder(A, B, number_of_cams * 2, 3);   //Transform A into upper triangular form
////				cout << "ready to solve: " << m << endl;
////				for (int a = 0; a < number_of_cams * 2; a++){
////					for (int b = 0; b < 3; b++){
////						cout << A[a][b] << " ";
////					}
////					cout << endl;
////				}
////				cout << endl << "B = " << endl;
////				for (int a = 0; a < number_of_cams * 2; a++)
////					cout << B[a] << endl;
//				solve.Backsub(A, B, 3, X);
//				//cout << "X = " << X[0] << ", "<< X[1] << ", " << X[2]<< endl;
//
//				lpt::Particle::Ptr  newparticle = lpt::Particle::create(id, X[0], X[1], X[2]);
//				newframe->particles.push_back(newparticle);
//				for(int i = 0; i < number_of_cams * 2; ++i)
//					delete [] A[i];
//				delete [] A;
//				delete [] B;
//				delete [] X;
//			}
//		}
//		frames.push_back(newframe);
//	}
//}

} /* namespace lpt */
