#include <core.hpp>
#include <visualization.hpp>
#include "vtkConeSource.h"
#include "vtkSphereSource.h"
#include "vtkArrowSource.h"
#include <vtkRibbonFilter.h>
#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkProperty.h"
#include "vtkCommand.h"
#include <vtkCellArray.h>
#include <vtkMath.h>

using namespace std;

class RenderParticlesCallback : public vtkCommand
{
public:
	RenderParticlesCallback(vtkRenderer& renderer, vector<lpt::TrajectoryPathVTK::Ptr>& trajs) : renderer(renderer), trajs(trajs), timer_count(0), count(0), add(0) {}
	
	virtual void Execute(vtkObject *caller, unsigned long eventId, void * vtkNotUsed(callData))
	{
		vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
		if (vtkCommand::TimerEvent == eventId)
			++this->timer_count;
		for (int i = 0; i < trajs.size(); ++i) {
			lpt::ParticleVectors current_particle;
			current_particle[0][0] = i*2 + 10.0 * std::cos(this->timer_count * 0.01);
			current_particle[0][1] = i - 5.0 * ( std::cos(this->timer_count* 0.01) - i * std::sin(this->timer_count* 0.01));
			current_particle[0][2] = -i * std::cos(this->timer_count * .1);
			
			current_particle[1][0] = i*2 + 10.0 * std::cos(this->timer_count * 0.01);
			current_particle[1][1] = i - 5.0 * ( std::cos(this->timer_count* 0.01) - i * std::sin(this->timer_count* 0.01));
			current_particle[1][2] = -i * std::cos(this->timer_count * .1);
		
			trajs[i]->addNextPoint(current_particle);
		}
		if (timer_count %100 == 0 && !add ) {
			trajs.pop_back();
			if (trajs.size() == 0 )
				add = 1;
		} else if ( timer_count % 100 == 0 && add ) {
			trajs.push_back(lpt::TrajectoryPathVTK::create(&renderer));
			if (trajs.size() == 10)
				add = 0;
		}
		iren->GetRenderWindow()->Render();
	}

private:
	bool add;
	int timer_count;
	int count;
	vtkRenderer& renderer;
	vector<lpt::TrajectoryPathVTK::Ptr>& trajs;
};


int main(int argc, char** argv)
{
	vtkRenderer* ren1= vtkRenderer::New();
	ren1->SetBackground( .2, .2, .2 );
	lpt::TrajectoryPathVTK::initializeLookupTable();
	vector<lpt::TrajectoryPathVTK::Ptr> trajs;
	for (int i = 0; i < 10; ++i) {
		trajs.push_back( std::move( lpt::TrajectoryPathVTK::create(ren1) ) );
	}
	lpt::TrajectoryPathVTK::max_points = 100;
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer( ren1 );
	renWin->SetSize( 640, 480 );

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(renWin);

	vtkInteractorStyleTrackballCamera *style = 
		vtkInteractorStyleTrackballCamera::New();
	iren->SetInteractorStyle(style);
	iren->Initialize();

	RenderParticlesCallback particles_callback(*ren1, trajs);
	iren->AddObserver(vtkCommand::TimerEvent, &particles_callback);
	int timerId = iren->CreateRepeatingTimer(1);

	iren->Start();
	
	trajs.clear();
	ren1->Delete();
	renWin->Delete();
	iren->Delete();
	style->Delete();

	cout << "Finished" << endl;
	return 0;
}


