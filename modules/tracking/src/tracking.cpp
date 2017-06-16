/**
 * @file tracking.cpp
 * The tracking module definition
 */
#include "tracking.hpp"

namespace lpt {
void callbackSetAlpha(int state, void* data)
{
    Tracker* tracker = static_cast<Tracker*>(data);
    tracker->params.alpha = tracker->params.alpha_level / 100.0;
}

void callbackSetRMax(int state, void* data)
{
    Tracker* tracker = static_cast<Tracker*>(data);
    tracker->params.max_radius = tracker->params.max_radius_level / 1.0;
}

void callbackSetRMin(int state, void* data)
{
    Tracker* tracker = static_cast<Tracker*>(data);
    tracker->params.min_radius = tracker->params.min_radius_level / 1.0;
}

template<typename Object_T>
void callbackSetBetaT(int state, void* data)
{
    Tracker_<Object_T>* tracker = static_cast<Tracker_<Object_T>*>(data);
    tracker->params.beta = tracker->params.beta_level / 100.0;
}

template<typename Object_T>
void callbackSetRMaxT(int state, void* data)
{
    Tracker_<Object_T>* tracker = static_cast<Tracker_<Object_T>*>(data);
    tracker->params.r_max = tracker->params.r_max_level / 1.0;
}
	
CostCalculator::~CostCalculator()
{
    cout << "CostCalculator destructed" << endl;
}

CostNearestNeighbor::CostNearestNeighbor()
{
    cout << "Cost Nearest Neighbor constructed" << endl;
}

CostNearestNeighbor::~CostNearestNeighbor()
{
    cout << "Cost Nearest Neighbor destructed" << endl;
}

double CostNearestNeighbor::calculate(const Trajectory3d &traj, const Particle3d &particle) const
{
    return lpt::calculateDistance(*(traj.objects.back().first), particle);
}

CostMinimumAcceleration::CostMinimumAcceleration()
{
    cout << "Cost Minimum Acceleration constructed" << endl;
}

CostMinimumAcceleration::~CostMinimumAcceleration()
{
    cout << "Cost Minimum Acceleration destructed" << endl;
}

double CostMinimumAcceleration::calculate(const Trajectory3d &traj, const Particle3d &particle) const
{
    return lpt::calculateDistance(traj.extrap_object, particle);
}

Tracker::Tracker() : clear_drawings(false), particle_count(0)
{
    cout << "Tracker constructed" << endl;
}

Tracker::~Tracker()
{
    cout << "Tracker destructed" << endl;
	cout << "particles detected: " << particle_count << endl;
	
    cout << "Save trajectories? (1 for yes, 0 for no)" << endl;
    cin >> isSave;
	//isSave = false;
    if (isSave) {
		stringstream file_name;
		file_name << shared_objects->output_path << "Trajectories.txt";
		ofstream output(file_name.str());
		output << trajectories.size() + active_trajs.size() << endl;
		output << endl;

        int id = 0;
		size_t link_count = 0;
        for (auto traj_it = trajectories.begin(); traj_it != trajectories.end(); traj_it++) {
			
			if ( (*traj_it)->objects.size() >= 10 ) {
				link_count += (*traj_it)->objects.size();
				saveTrajectory(output, *traj_it, id);
				id++;
			}
        }
		for (auto traj_it = active_trajs.begin(); traj_it != active_trajs.end(); traj_it++) {
			
			if ( (*traj_it)->objects.size() >= 10 ) {
				link_count == (*traj_it)->objects.size();
				saveTrajectory(output, *traj_it, id);
				id++;
			}
		}
		cout << id << endl;
		cout << "links tracked: " << link_count << endl;
    }
}

void Tracker::saveTrajectory(ofstream& output, Trajectory3d_Ptr traj, int id)
{
    //stringstream file_name;
    //file_name << shared_objects->output_path << "trajectories/Trajectory_" << id << ".txt";
    //ofstream output(file_name.str());

    for (size_t i=0; i<traj->objects.size(); i++) {
        output << traj->objects[i].first->frame_index << " ";
		for (size_t j=0; j<traj->objects[i].first->X.size(); j++) {
			output << traj->objects[i].first->X[j] << " ";
		}
        for (size_t j=traj->objects[i].first->X.size(); j<traj->objects[i].second.size(); j++) {
            output << traj->objects[i].second[j] << " ";
		}
		output << endl;
	}
	output << endl;
}

void Tracker::trackFrames(lpt::Frame3d& frame1, lpt::Frame3d& frame2 )
{
	particle_count += frame1.objects.size() + frame2.objects.size();
    list<lpt::Particle3d_Ptr> next_unmatched_particles;
    list<lpt::Trajectory3d_Ptr> new_trajs;

    auto& frame_rate = this->shared_objects->frame_rate;

    if ( active_trajs.empty() && unmatched_particles.empty() )
    {
        std::move(frame1.objects.begin(), frame1.objects.end(), std::back_inserter(unmatched_particles));
    }
	      
    // TODO: Maybe try parallelizing this loop for multi-core system
    // Extrapolate all active trajectories
    for (list<lpt::Trajectory3d_Ptr>::iterator traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter)
    {
        (*traj_iter)->candidates_costs.clear();
        lpt::extrapolatePosition( *(*traj_iter) );
    }
			
    // Loop over all particles in frame2
    for (int p2 = 0; p2 < frame2.objects.size(); ++p2) {
        lpt::Particle3d& particle2 = *frame2.objects[p2];
        double lowest_cost = std::numeric_limits<double>::max();
        lpt::Trajectory3d* traj_match = 0;
		
        // Search among active trajectories
        for (list<lpt::Trajectory3d_Ptr>::iterator traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter) {
            lpt::Trajectory3d& traj = *(*traj_iter);
            double distance = lpt::calculateDistance(particle2, traj.extrap_object);
            double radius = (traj.search_radius > params.min_radius) ? traj.search_radius : params.min_radius;

            if (distance <= params.alpha * radius) {
                double cost = cost_calculator->calculate(traj, particle2);
                if (cost < lowest_cost) {
                    lowest_cost = cost;
                    traj_match = &traj;
                }
            }
        }		
			
        // Candidate exists within search radius among active trajectories
        if ( traj_match ) {
            traj_match->candidates_costs.push_back( std::make_pair(p2, lowest_cost) );
            std::push_heap( traj_match->candidates_costs.begin(), traj_match->candidates_costs.end(), lpt::compareCost );
        }

        // Candidate does not exist within search radius. Loop over unmatched particles in frame 1 to find candidates within maximum radius
        else {
            lowest_cost = std::numeric_limits<double>::max();
            list<lpt::Particle3d_Ptr>::iterator matched_particle_iter = unmatched_particles.end();
            for (list<lpt::Particle3d_Ptr>::iterator p1_iter = unmatched_particles.begin(); p1_iter != unmatched_particles.end(); ++p1_iter) {
                double cost = lpt::calculateDistance(*(*p1_iter), *frame2.objects[p2]);
                if (cost <= params.max_radius && cost < lowest_cost ) {
                    lowest_cost = cost;
                    matched_particle_iter = p1_iter;
                }
            }
            // Candidate found within maximum radius. Initialize new trajectory
            if ( matched_particle_iter != unmatched_particles.end() ) {
                // Create a new trajectory
                lpt::Trajectory3d_Ptr newtraj = lpt::Trajectory3d::create();
                // Push particles to the trajectory object vector
				auto particle1 = *matched_particle_iter;
				auto particle2 = frame2.objects[p2];
                array<double, 9> temp1, temp2;

                if (this->shared_objects->KF_isOn) {
                    // Create Kalman filter
                    cv::Mat stateIni = cv::Mat::zeros(9, 1, CV_64F);

                    for (int i=0; i<3; i++) {
						stateIni.at<double>(i) = particle1.get()->X[i];
                    }

                    newtraj->initializeKF(stateIni, params.KF_sigma_a, params.KF_sigma_z);
                    auto statePre = newtraj->kalman_filter.predict();

                    // Take measurement
                    cv::Mat measurement(9, 1, CV_64F);
                    for (int i=0; i<9; i++) {
                        if (i < 3) {
                            measurement.at<double>(i) = particle2.get()->X[i];
                        }
                        else if (i < 6) {
                            measurement.at<double>(i) = (particle2.get()->X[i-3] - particle1.get()->X[i-3]) * frame_rate;
                        }
                        else {
                            measurement.at<double>(i) = 0;
                        }
                    }

                    auto statePost = newtraj->kalman_filter.correct(measurement);
                    for (int i=0; i<9; i++){
                        temp1[i] = stateIni.at<double>(i);
                        temp2[i] = statePost.at<double>(i);
                    }

                }

                newtraj->objects.push_back( std::make_pair(particle1, std::move(temp1)) );
                newtraj->objects.push_back( std::make_pair(particle2, std::move(temp2)) );

                //newtraj->obj_state.push_back(std::move(temp1));
                //newtraj->obj_state.push_back(std::move(temp2));

                new_trajs.push_back(std::move(newtraj));
            }

            // Candidate not found within maximum radius. Add p2 to next unmatched particles list
            else {
                next_unmatched_particles.push_back(frame2.objects[p2]);
            }
        }
    }

    unmatched_particles.clear();
    std::move( next_unmatched_particles.begin(), next_unmatched_particles.end(), std::back_inserter(unmatched_particles) );

    // Loop over all active trajectories
    for (list<lpt::Trajectory3d_Ptr>::iterator traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter) {
        // If has a candidate, add the candidate to trajectory
        if ( ! (*traj_iter)->candidates_costs.empty() ) {
            int p2 = (*traj_iter)->candidates_costs.front().first;
			auto curr_obj = frame2.objects[p2];
            array<double, 9> temp;

            if (this->shared_objects->KF_isOn) {
                auto obj_it = (*traj_iter)->objects.end();
                auto X0 = curr_obj->X;
                auto Xm1 = (obj_it-1)->first->X;
                auto Xm2 = (obj_it-2)->first->X;

                cv::Mat measurement(9, 1, CV_64F);
                for (int i=0; i<9; i++) {
                    if (i < 3) {
                        measurement.at<double>(i) = X0[i];
                    }
                    else if (i < 6) {
                        measurement.at<double>(i) = (X0[i-3] - Xm1[i-3]) * frame_rate;
                    }
                    else {
                        measurement.at<double>(i) = (X0[i-6] - 2*Xm1[i-6] + Xm2[i-6]) * frame_rate * frame_rate;
                    }
                }
                auto statePre = (*traj_iter)->kalman_filter.predict();
                auto statePost = (*traj_iter)->kalman_filter.correct(measurement);


                for (int i=0; i<9; i++)
                    temp[i] = statePost.at<double>(i);
            }
            (*traj_iter)->objects.push_back( std::make_pair(curr_obj, std::move(temp)) );
            //(*traj_iter)->obj_state.push_back(std::move(temp));
        }
        // If no candidate, remove trajectory from active list
        else {
            trajectories.push_back(*traj_iter);  //TODO: If need to save trajectories then uncomment this!!!!!!!!!!!!!!
            traj_iter = --(active_trajs.erase(traj_iter));
        }
    }
    // Add newly created trajectories to the active list
    std::move( new_trajs.begin(), new_trajs.end(), std::back_inserter(active_trajs) );
}

void Tracker::addControls()
{
    void* tracker_void_ptr = static_cast<void*> ( this );
    cv::createTrackbar("Alpha", string() , &params.alpha_level, 200, callbackSetAlpha, tracker_void_ptr);
    cv::createTrackbar("MaxRadius", string() , &params.max_radius_level, 100, callbackSetRMax, tracker_void_ptr);
    cv::createTrackbar("MinRadius", string() , &params.min_radius_level, 50, callbackSetRMin, tracker_void_ptr);
    cout << "Added Tracker Controls to Window" << endl;
}
	
void Tracker::setTrajectoryViews(vector<lpt::Camera>& cameras, cv::Mat image_type)
{
    trajectory_views.clear();
    trajectory_views.resize( cameras.size() );
    for( int i = 0; i < cameras.size(); ++i ) {
        trajectory_views[i] = cv::Mat::zeros(image_type.size(), CV_8UC3);
    }
}
		
void Tracker::drawTrajectories(vector<lpt::Camera>& cameras)
{
    vector<cv::Point3d> object_points2( active_trajs.size() );
    vector<cv::Point3d> object_points1( active_trajs.size() );
    vector<cv::Point2d> image_points1( active_trajs.size() );
    vector<cv::Point2d> image_points2( active_trajs.size() );
    vector<pair<lpt::Particle3d_Ptr, array<double, 9>>>::iterator object_iter;

    int id = 0;
    for (list<lpt::Trajectory3d_Ptr>::iterator traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter) {
        object_iter = (*traj_iter)->objects.end();
        --object_iter;
        object_points2[id].x = (*object_iter).first->X[0];
        object_points2[id].y = (*object_iter).first->X[1];
        object_points2[id].z = (*object_iter).first->X[2];

        --object_iter;
        object_points1[id].x = (*object_iter).first->X[0];
        object_points1[id].y = (*object_iter).first->X[1];
        object_points1[id].z = (*object_iter).first->X[2];

        ++id;
    }
		
    for (int i = 0; i < trajectory_views.size(); ++i) {
        cv::Mat R = cv::Mat(3, 3, CV_64F, cameras[i].R);
        cv::Mat t_vec = cv::Mat(3, 1, CV_64F, cameras[i].T);
        cv::Mat r_vec = cv::Mat::zeros(3,1, CV_64F);
        cv::Rodrigues(R, r_vec);
        if ( ! object_points1.empty() && ! object_points2.empty() ) {
            cv::projectPoints(cv::Mat(object_points1), r_vec, t_vec, cameras[i].getCameraMatrix(), cameras[i].getDistCoeffs(), image_points1);
            cv::projectPoints(cv::Mat(object_points2), r_vec, t_vec, cameras[i].getCameraMatrix(), cameras[i].getDistCoeffs(), image_points2);
        }

        id = 0;
        for (list<lpt::Trajectory3d_Ptr>::iterator traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter) {
            cv::circle(trajectory_views[i], image_points2[id], 2, cv::Scalar(255,0,0) , -1 );
            cv::line(trajectory_views[i], image_points1[id], image_points2[id],  cv::Scalar(255,0,0) );
            ++id;
        }
        cv::imshow(cameras[i].name, trajectory_views[i] );
    }
}

TestTracker::TestTracker(vector<Trajectory3d_Ptr> &trajs, vector<Trajectory3d_Ptr> &gold_trajs)
  : trajectories(trajs), gold_trajectories(gold_trajs),
    correct_ratio(0.0), cover_ratio(0.0),
    links_created(0), links_total(0), links_correct(0)
{
}

void TestTracker::testTrajectories(size_t maxframes)
{
    cout << endl <<"-------Running gold test------" << endl;
	vector<vector<bool>> matches( gold_trajectories.size() );
	for(int i = 0; i < gold_trajectories.size(); i++){
		vector<bool> temp(gold_trajectories[i]->objects.size(), 0);
		matches[i] = temp;
	}

    for (int f = 0; f < maxframes - 1; f++){                  ///<FIXME Assumes that no gaps exist in the trajectorys which start at Traj->startframe
		for (int i = 0; i < gold_trajectories.size(); i++){

			int sg = gold_trajectories[i]->objects[0].first->frame_index;
			if (sg <= f   &&  f + 1 - sg < (int)(gold_trajectories[i]->objects.size()) ){
				auto& P_g1 = gold_trajectories[i]->objects[f - sg].first;
				auto& P_g2 = gold_trajectories[i]->objects[f + 1 - sg].first;

				for (int j = 0; j < trajectories.size(); j++){
					int st = trajectories[j]->objects[0].first->frame_index;
					if (st <= f  &&  f + 1 - st < (int)(trajectories[j]->objects.size()) ){
						auto& P_a1 = trajectories[j]->objects[f - st].first;
						auto& P_a2 = trajectories[j]->objects[f + 1 - st].first;

						if ( P_a1->id == P_g1->id  &&  P_a2->id == P_g2->id ){
							matches[i][f - sg] = 1;
						}
					}
				}
			}
		}
	}

	for (int j = 0; j < trajectories.size(); j++)
		links_created += trajectories[j]->objects.size()-1;

	for (int i = 0; i < gold_trajectories.size(); i++){
		if (gold_trajectories[i]->objects.size() >= 6){
			links_total += gold_trajectories[i]->objects.size()-1;
			int temp = 0;
			for (int p = 0; p < gold_trajectories[i]->objects.size(); p++){
				temp += matches[i][p];
			}
			links_correct += temp;
#if(EVALRESULT)
			gold_trajectories[i]->Eff = (double)temp / (double)(gold_trajectories[i]->particles.size() - 1);
#endif
		}
	}
#if(EVALRESULT)
	AnalyzeInput eval;
	Input in;
	Output out;
	string trajstat = "trajstats";
	string framestat = "framestats";
	vector<Frame*> frames = in.trajs2frames(gold_trajectories);
	eval.evalframes(frames);
	eval.evaltrajresult(gold_trajectories,frames);

	out.fprintTrajsStats(trajstat, gold_trajectories);
	out.fprintFramesStats(framestat, frames);

	string eval_results = "eval_results";

	ofstream fout(eval_results);
	// Print out the gold trajectories with match indicators to show which trajs could be correctly reconstructed
	if( fout.isOpen() ){
		for (int i = 0; i < gold_trajectories.size(); i++) {
			if(gold_trajectories[i]->particles.size() > 1) {
				for (int p = 0; p < gold_trajectories[i]->particles.size()-1; p++){
					Particle* PD = gold_trajectories[i]->particles[p];
					fout << PD->id << "\t" << PD->frame_index << "/t" << (PD->ds)/(PD->dr) << "\t" << (int)gold_trajectories[i]->matches[p] << endl;
				}
			}
		}
	}
	fout.close();
#endif
}
void TestTracker::printTestResults()
{

	cout << "Results:" << endl << endl;
	cout << "Total number of Trajectories Tracked:\t= " << trajectories.size() << endl;
	cout << "Total number of Gold Trajectories:\t= " << gold_trajectories.size() << endl << endl;
	cout << "Average Trajectory Length tracked:\t= " << (int) ((double)links_created / trajectories.size() + 1) << " frames" << endl;
	cout << "Average Trajectory Length (gold_trajectories)\t= " << (int) ((double)links_total / gold_trajectories.size() + 1) << " frames" << endl;

	cout << "Total number of links created:\t= " << links_created << endl;
	cout << "Total number of total true links:\t= " << links_total << endl;
	cout << "Total number of correct links:\t= " << links_correct << endl;

	cover_ratio = 100. * (double)links_correct / links_total;
	correct_ratio = 100. * (double)links_correct / links_created;

	cout << endl;
	cout << "Coverage Ratio (links_correct / links_total)\t= " << cover_ratio << " percent" << endl;
	cout << "Quality Tracking Efficiency (links_correct / links_created)\t= " << correct_ratio << endl;
}

} /* NAMESPACE_PT */
