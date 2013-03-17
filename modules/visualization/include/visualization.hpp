#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_

#include <core.hpp>  //LPT core module

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/numeric/functional.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/error_of.hpp>
#include <boost/accumulators/statistics/error_of_mean.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_sum.hpp>
#include <boost/accumulators/statistics/weighted_covariance.hpp>
#include <boost/accumulators/framework/features.hpp>

#include <vtkObjectFactory.h>
#include <vtkConeSource.h>
#include <vtkArrowSource.h>
#include <vtkSphereSource.h>
#include <vtkLineSource.h>
#include <vtkPlaneSource.h>
#include <vtkRibbonFilter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkActor.h>
#include <vtkGlyph3DMapper.h>
#include <vtkRenderer.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkCommand.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkCell.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkMath.h>
#include <vtkAxesActor.h>
#include <vtkSmartPointer.h>
#include <vtkPropPicker.h>
#include <vtkVectorText.h>
#include <vtkFollower.h>
#include <vtkImageData.h>
#include <vtkStructuredGrid.h>
#include <vtkCellCenters.h>
#include <vtkExtractEdges.h>
#include <vtkOutlineFilter.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkPlaneWidget.h>
#include <vtkImplicitPlaneWidget.h>
#include <vtkImagePlaneWidget.h>
#include <vtkScalarBarWidget.h>
#include <vtkImageReslice.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkAlgorithmOutput.h>
#include <vtkFrustumSource.h>
#include <vtkPlanes.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkDataSetAttributes.h>
#include <vtkInformation.h>
#include "vtkAxis.h"
#include "vtkChartXY.h"
#include "vtkContextView.h"
#include "vtkContextScene.h"
#include "vtkContextActor.h"
#include "vtkTable.h"
#include "vtkPen.h"
#include "vtkTextProperty.h"
#include "vtkPlotPoints.h"
#include "vtkBrush.h"
#include <vtkStreamLine.h>

namespace lpt {

using namespace std;
using namespace boost::accumulators;

typedef accumulator_set< double, features<tag::mean, tag::variance, tag::count, tag::max, tag::min > > boost_accumulator;
typedef accumulator_set< double, features<tag::weighted_mean, tag::weighted_variance, tag::count, tag::max, tag::min >, double > boost_weighted_accumulator;
typedef accumulator_set< double, features<tag::weighted_mean, tag::weighted_covariance<double, tag::covariate1> >, double > boost_weighted_covariance_accumulator;
typedef vector < array <double, 3> > VectorField;
typedef vector < array <double, 6> > ReynoldsStressField;
typedef array < array <double, 3> , 4 > ParticleVectors;


enum ViewMode {TRAJECTORIES, VECTORGRID, VECTORGRID_AND_TRAJECTORIES, NONE};
enum VectorMode {VELOCITY, CORRECTED_VELOCITY, ACCELERATION, TURBULENT_KINETIC_ENERGY, TURBULENCE_DISSIPATION_RATE, MASS_RESIDUAL, VELOCITY_SD, ACCELERATION_SD, COUNT, PRESSURE, PRESSURE_LAGRANGIAN, VORTICITY};

class PickDim;
class CoordinateArrows;
class CameraVTK;
class TrajectoryPathVTK;
class ParticlesVTK;
class FluidProperties;
class TrajectoryHandler;
class FiniteVolumeGrid;
class VisualizerInteractorStyle;
class Visualizer;

void getVectorField(vector<array<lpt::boost_accumulator, 3> >& accumulators, lpt::VectorField& vector_field);
void getVectorField(vector<array<lpt::boost_weighted_accumulator, 3> >& weighted_accumulators, lpt::VectorField& vector_field);
void getReynoldsStressField(vector< array <lpt::boost_weighted_accumulator, 3 > >& velocity_accumulators, vector< array <lpt::boost_weighted_covariance_accumulator, 3 > >& reynolds_stress_accumulators, lpt::ReynoldsStressField& stress_field);

array<double, 3>  calcPressureGradiantNavierStokes( array<double,3>& U0, array<double,3>& U1, array<double,3>& U2, array<double,3>& A1, lpt::FluidProperties& fluid_props, array<double,3>& dX);
array<double, 3>  calcPressureGradiantBernoulli( array<double,3>& U0, array<double,3>& U2, lpt::FluidProperties& fluid_props, array<double,3>& dX);

class FluidProperties {
public:
	typedef	shared_ptr<FluidProperties> Ptr;
	static inline FluidProperties::Ptr create() { return FluidProperties::Ptr( new FluidProperties() ); }

	FluidProperties() : rho(1.225), Td(15), Tw(15), mu(1.78E-5), P_ref(101325) {}
	// Default air at 15 C
	double rho;  //density (kg/m3)
	double Td;   //dry bulb temperature (C)
	double Tw;   //wet bulb temperature (C)
	double mu;   // dynamic viscocity ( kg/(m·s) )
	double P_ref; // Reference Pressure (Pa)
};

class CoordinateArrows {
public:
	typedef std::shared_ptr<CoordinateArrows> Ptr; 
	static inline CoordinateArrows::Ptr create(	vtkSmartPointer<vtkRenderWindowInteractor> iren, double scale = 1) { return CoordinateArrows::Ptr(new CoordinateArrows(iren, scale) ); }

	CoordinateArrows( vtkSmartPointer<vtkRenderWindowInteractor> iren, double scale);
	virtual ~CoordinateArrows() { }
private:
	double scale;
	vtkSmartPointer<vtkRenderWindowInteractor> interactor;
	vtkSmartPointer < vtkArrowSource > arrow;
	vector < vtkSmartPointer < vtkActor >> actors;
	vtkSmartPointer < vtkPolyDataMapper > mapper;

	vtkSmartPointer<vtkAxesActor> axes;
	vtkSmartPointer<vtkOrientationMarkerWidget> widget; 
};

class CameraVTK {
public:
	
	static array<double, 3> getRandomColor() { 
		array<double, 3> color;
		color[0] = vtkMath::Random(0.5, 1.0);
		color[1] = vtkMath::Random(0.5, 1.0);
		color[2] = vtkMath::Random(0.5, 1.0);
		return color;
	}

	CameraVTK(lpt::Camera& camera);

	void calcPlaneEq(array<double,3>& A, array<double,3>& B, array<double,3>& C, array<double,4>& eq );

	void camToWorld(lpt::Camera& cam, array<double, 3>& pc, array<double, 3>& pw);

	void addToRenderer(vtkSmartPointer<vtkRenderer> renderer); 

	void removeFromRenderer(vtkSmartPointer<vtkRenderer> renderer);

private:

	vtkSmartPointer<vtkActor> lineactor;
	vtkSmartPointer<vtkLineSource> line;
	vtkSmartPointer<vtkPolyDataMapper> linemapper;

	vtkSmartPointer<vtkFrustumSource> frustrum;
	vtkSmartPointer<vtkPlanes> planes;
	vtkSmartPointer<vtkActor> actor;
	vtkSmartPointer<vtkPolyDataMapper> mapper;

	vtkSmartPointer<vtkVectorText> text;
	vtkSmartPointer<vtkPolyDataMapper> textmapper;
	vtkSmartPointer<vtkFollower> follower;
};

class TrajectoryPathVTK : public lpt::TrajectoryVTKBase {
public:
	typedef std::shared_ptr<TrajectoryPathVTK> Ptr; 
	
	static vtkSmartPointer<vtkLookupTable> lookuptable;
	static int max_points;
	static bool show_paths;
	static int scalar_index;

	static inline void initializeLookupTable() {
		lookuptable = vtkSmartPointer<vtkLookupTable>::New();
		lookuptable->SetRange( 0, 2.0 );
		lookuptable->SetNumberOfColors(256);
		lookuptable->SetHueRange(0.667, 0.0);
		lookuptable->Build();
	}

	static inline void setLookupTable(vtkSmartPointer<vtkLookupTable> lut) { lookuptable = lut; }

	static inline void setVectorMode(lpt::VectorMode mode) {	
		switch (mode) {
		case lpt::VELOCITY:
			scalar_index = 1;
			break;
		case lpt::ACCELERATION:
			scalar_index = 2;
			break;
		default:
			break;
		}
	}

	static inline TrajectoryPathVTK::Ptr create(vtkRenderer* renderer) { return TrajectoryPathVTK::Ptr(new TrajectoryPathVTK(renderer)); }
	
	TrajectoryPathVTK(vtkRenderer* renderer);
	
	void addNextPoint(lpt::ParticleVectors& current_particle);

	void removeFromRenderer();

	virtual ~TrajectoryPathVTK() { removeFromRenderer(); }

private:
	vector<vtkIdType> point_ids;
	deque<double> scalars_queue;
	deque<array<double, 3> > position_queue;
	vtkSmartPointer<vtkRibbonFilter> ribbon_filter;
	vtkSmartPointer<vtkTubeFilter> tube_filter;
	vtkSmartPointer<vtkDoubleArray> scalars;
	vtkSmartPointer<vtkPolyDataMapper> traj_mapper;
	vtkSmartPointer<vtkActor> trajactor;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkCellArray> lines;
	vtkSmartPointer<vtkPolyData> polydata;
	vtkRenderer* renderer;
};

class PressureFieldSolver {
public:
	typedef	std::shared_ptr<PressureFieldSolver> Ptr;

	PressureFieldSolver() : dt(0), dx(0), dy(0), dz(0), initial_pressure(0) {
		grid_cell_counts[0] = 0; 
		grid_cell_counts[1] = 0;
		grid_cell_counts[2] = 0;
		fluid_props = lpt::FluidProperties::create();  // Default to air at 15 C;
	}
	
	virtual void solve(lpt::VectorField& velocity_field, lpt::VectorField& acceleration_field, lpt::ReynoldsStressField& stress_field)=0;
	
	virtual void addControls(){}
	
	virtual void setFluidProperties(FluidProperties::Ptr props) {
		fluid_props = props;
	}

	virtual void setGridProperties(array<int, 3>& grid_cell_counts, array<double, 3>& cell_dimensions) {
		this->grid_cell_counts[0] = grid_cell_counts[0];
		this->grid_cell_counts[1] = grid_cell_counts[1];
		this->grid_cell_counts[2] = grid_cell_counts[2];
		dx = cell_dimensions[0] / 1000.0; // mm
		dy = cell_dimensions[1] / 1000.0; // mm
		dz = cell_dimensions[2] / 1000.0; // mm
		this->resetPressureField();
	}
	
	inline void setSharedObjects( std::shared_ptr < lpt::SharedObjects > new_shared_objects ) { 
		shared_objects = new_shared_objects; 
	}
	 
	inline const vector<double>& getPressureField() { return pressure; }

	inline void resetPressureField() {
		pressure.clear();
		int total_cell_count = grid_cell_counts[0] * grid_cell_counts[1] * grid_cell_counts[2]; 
		pressure.resize(total_cell_count, initial_pressure);
	}

protected:
	//index of grid cell (i,j,k)
	inline int getGridIndex(int i, int j, int k) {
		return ( i + j * grid_cell_counts[0] + k * (grid_cell_counts[0] * grid_cell_counts[1]) );     
	}

	double initial_pressure;
	FluidProperties::Ptr fluid_props;
	std::shared_ptr < lpt::SharedObjects > shared_objects;
	vector<double> pressure;
	double dx, dy, dz, dt;
	array<int, 3> grid_cell_counts;
	lpt::VectorField* velocity_field;
	lpt::VectorField* acceleration_field;
	lpt::ReynoldsStressField* stress_field;
};

class PressurePoissonSolver2Step : public PressureFieldSolver {
public:
	typedef	std::shared_ptr<PressurePoissonSolver2Step> Ptr;
	static inline PressurePoissonSolver2Step::Ptr create() { return PressurePoissonSolver2Step::Ptr(new PressurePoissonSolver2Step() ); }
	
	class Parameters {
	public:
		Parameters() : sor_weight(1.2), max_iterations(20000), convergence_tolerance(1E-6) {}
		double sor_weight;
		double max_iterations;
		double convergence_tolerance;
	};

	Parameters params;

	PressurePoissonSolver2Step() {}

	virtual void solve(lpt::VectorField& velocities, lpt::VectorField& accelerations, lpt::ReynoldsStressField& stresses);

	virtual void addControls();

private:
	void calcCellPressure(int i, int j, int k);
	void calcCellPressureGradientsSteadyRANS(int i, int j, int k);

	lpt::VectorField pressure_gradients;
};

class LagrangianPressureFieldSolver : public PressureFieldSolver {
public:
	typedef	shared_ptr<LagrangianPressureFieldSolver> Ptr;
	static inline LagrangianPressureFieldSolver::Ptr create() { return LagrangianPressureFieldSolver::Ptr(new LagrangianPressureFieldSolver() ); }
			
	class Parameters {
	public:
		Parameters() : sor_weight(1.2), max_iterations(20000), convergence_tolerance(1E-6) {}
		double sor_weight;
		double max_iterations;
		double convergence_tolerance;
	} params;

	LagrangianPressureFieldSolver() {}
	virtual void solve(lpt::VectorField& velocities, lpt::VectorField& accelerations, lpt::ReynoldsStressField& stresses);
	virtual void solve(lpt::VectorField& pressure_grads);
	virtual void addControls();

private:
	void calcCellPressure(int i, int j, int k);
	
	lpt::VectorField pressure_gradients;
};

class FiniteVolumeGrid {
public:
	typedef	shared_ptr<FiniteVolumeGrid> Ptr;
	static inline FiniteVolumeGrid::Ptr create(vtkSmartPointer < vtkRenderer > renderer) { return FiniteVolumeGrid::Ptr(new FiniteVolumeGrid(renderer) ); }
	
	class Parameters {
	public:
		Parameters() { 
			scalar_range[0] = 0; 
			scalar_range[1] = 2.5;

			grid_cell_counts[0] = 27;//35;//65;// 
			grid_cell_counts[1] = 37;//25;//31;// 
			grid_cell_counts[2] = 37;//25;//31;// 
			
			double length_size = 300.0;//160; ////
			grid_dimensions[0] = grid_cell_counts[0] * length_size / grid_cell_counts[1]; // mm
			grid_dimensions[1] = length_size; // mm
			grid_dimensions[2] = length_size; // mm

			cell_dimensions[0] = grid_dimensions[0] / static_cast<double>(grid_cell_counts[0]); //mm
			cell_dimensions[1] = grid_dimensions[1] / static_cast<double>(grid_cell_counts[1]); //mm
			cell_dimensions[2] = grid_dimensions[2] / static_cast<double>(grid_cell_counts[2]); //mm
			
			grid_origin[0] = 137; //-1.0*length_size/2.0; //  //mm 
			grid_origin[1] = -1.0*length_size/2.0 + 83; //mm
			grid_origin[2] = -1.0*length_size/2.0; //mm

			fixed_scale_factor = 0.65 * (*std::max_element(cell_dimensions.begin(), cell_dimensions.end())); 
		}

		array<double, 3> grid_dimensions;
		array<int, 3> grid_cell_counts;
		array<double, 3> cell_dimensions;
		array<double, 3> grid_origin;
		
		double fixed_scale_factor;
		double scalar_range[2];
	} params;

	FiniteVolumeGrid (vtkSmartPointer < vtkRenderer > renderer);
	void initialize();
	void setPressureFieldSolver(lpt::PressureFieldSolver::Ptr solver);
	void setGridOrigin(double x, double y, double z) {
		this->params.grid_origin[0] = x;
		this->params.grid_origin[1] = y;
		this->params.grid_origin[2] = z;
	}

	void setGridCellCounts(int nx, int ny, int nz) {
		this->params.grid_cell_counts[0] = nx;
		this->params.grid_cell_counts[1] = ny;
		this->params.grid_cell_counts[2] = nz;

		this->params.cell_dimensions[0] = this->params.grid_dimensions[0] / static_cast<double>(this->params.grid_cell_counts[0]); //mm
		this->params.cell_dimensions[1] = this->params.grid_dimensions[1] / static_cast<double>(this->params.grid_cell_counts[1]); //mm
		this->params.cell_dimensions[2] = this->params.grid_dimensions[2] / static_cast<double>(this->params.grid_cell_counts[2]); //mm
	}

	void setGridDimensions(double length_x, double length_y, double length_z) {
		this->params.grid_dimensions[0] = length_x;
		this->params.grid_dimensions[1] = length_y;
		this->params.grid_dimensions[2] = length_z;

		this->params.cell_dimensions[0] = this->params.grid_dimensions[0] / static_cast<double>(this->params.grid_cell_counts[0]); //mm
		this->params.cell_dimensions[1] = this->params.grid_dimensions[1] / static_cast<double>(this->params.grid_cell_counts[1]); //mm
		this->params.cell_dimensions[2] = this->params.grid_dimensions[2] / static_cast<double>(this->params.grid_cell_counts[2]); //mm
	}
	
	void updateAccumulators(lpt::Trajectory3d* traj_ptr, lpt::ParticleVectors& current_particle );
	void resetGrid();
	void savePlaneData();

	void updateGrid();
	bool gridIsFull();

	void calcPressures();
	void calcPressuresLagrangian();

	void calcMassResiduals();
	void calcVorticity();

	vtkSmartPointer<vtkLookupTable> getLookupTable() {
		return lookuptable;
	}

	vtkSmartPointer<vtkScalarBarWidget> getScalarBarWidget() {
		return scalarbar;
	}

	inline void setCellNeighborhood(vector<array<int,3>>& neighborhood) {
		cell_neighborhood.resize(neighborhood.size());
		std::copy(neighborhood.begin(), neighborhood.end(), cell_neighborhood.begin());
	}

	inline void addImplicitPlane() {
		renderer->AddActor(planeactor);
		plane_widget->On();
		renderer->RemoveActor(glyphactor);
	}
	inline void removeImplicitPlane() {
		renderer->RemoveActor(planeactor);
		plane_widget->Off();
		renderer->AddActor(glyphactor);
	}
	inline void addToRenderer() {
		renderer->AddActor(outlineactor);
		if ( !getImplicitPlane()->GetEnabled() )
			renderer->AddActor(glyphactor);
		scalarbar->On();
	}
	inline void removeFromRenderer() {
		renderer->RemoveActor(outlineactor);
		if ( getImplicitPlane()->GetEnabled() )
			removeImplicitPlane();
		renderer->RemoveActor(glyphactor);
		scalarbar->Off();
	}

	void setVectorMode(lpt::VectorMode mode);

	inline lpt::VectorMode getVectorMode() {return vector_mode;}
	
	inline void setFluidProperties(FluidProperties::Ptr props) { 
		fluid_props = props; 
		pressure_solver->setFluidProperties(fluid_props);
		lagrangian_pressure_solver->setFluidProperties(fluid_props);
	}

	vtkSmartPointer<vtkImplicitPlaneWidget> getImplicitPlane() {return plane_widget;}
	
	vtkSmartPointer<vtkActor> getImplicitPlaneActor() { return planeactor;}
	
	inline double* getBounds() { return grid->GetBounds(); }
	
	inline void setSharedObjects( std::shared_ptr < lpt::SharedObjects > new_shared_objects ) { 
		shared_objects = new_shared_objects; 
		pressure_solver->setSharedObjects(new_shared_objects);
		lagrangian_pressure_solver->setSharedObjects(new_shared_objects);
	}	
	
	~FiniteVolumeGrid() { }

private:
	inline int getGridIndex(int i, int j, int k) {
		return ( i + j * params.grid_cell_counts[0] + k * (params.grid_cell_counts[0] * params.grid_cell_counts[1]) );     //index of grid cell (i,j,k)
	}
	
	inline void computeCellCoordsIJK(const int cell_id, array<int,3>& ijk)
	{
		int ni  = params.grid_cell_counts[0]; 
		int nj  = params.grid_cell_counts[1]; 
		int nij = ni * nj;

		int k = cell_id / nij + 1;
		int j = (cell_id - (k - 1) * nij ) / ni + 1;
		int i = cell_id - (k - 1) * nij - (j - 1) * ni + 1;
		ijk[0] = i - 1;
		ijk[1] = j - 1;
		ijk[2] = k - 1;
	}

	vtkSmartPointer < vtkRenderer > renderer;
	vtkSmartPointer<vtkImageData> grid;
	vtkSmartPointer<vtkCellCenters> cellcentersfilter;
	vtkSmartPointer<vtkOutlineFilter> outline;
	vtkSmartPointer<vtkGlyph3DMapper> glyph3Dmapper;
	vtkSmartPointer<vtkArrowSource> arrow_source;
	vtkSmartPointer<vtkSphereSource> sphere_source;

	vtkSmartPointer<vtkActor> glyphactor;
	vtkSmartPointer<vtkActor> outlineactor;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
	vtkSmartPointer<vtkScalarBarWidget> scalarbar;
	vtkSmartPointer<vtkImplicitPlaneWidget> plane_widget;
	vtkSmartPointer<vtkImagePlaneWidget> planeWidgetX; 
	vtkSmartPointer<vtkPolyDataMapper> planemapper;
	vtkSmartPointer<vtkActor> planeactor;

	vtkSmartPointer<vtkLookupTable> lookuptable; 
	vtkSmartPointer<vtkDoubleArray> vectors;
	vtkSmartPointer<vtkDoubleArray> magnitudes;
	vtkSmartPointer<vtkIntArray> arrow_source_ids;
	vtkSmartPointer<vtkIntArray> sphere_source_ids;
	vector< array <boost_weighted_accumulator, 3> > acceleration_accumulators;
	vector< array <boost_weighted_accumulator, 3> > velocity_accumulators;
	vector< array <boost_weighted_accumulator, 3> > pressure_grad_accumulators;

	// Reynolds Stress Auccumulators: Each contains 3 of the 6 symmetric off diagonal covariance terms from the reynolds stress tensor; 
	// index [0](u,v), [1](u,w), [2](v,w),
	vector< array <boost_weighted_covariance_accumulator, 3 > > reynolds_stress_accumulators;  
	vector<double> pressure;
	vector<double> mass_residuals;
	lpt::VectorField vorticity;
	vector< array <int ,3> > cell_neighborhood;
	lpt::VectorMode vector_mode;
	lpt::PressureFieldSolver::Ptr pressure_solver;
	lpt::LagrangianPressureFieldSolver::Ptr lagrangian_pressure_solver;
	lpt::FluidProperties::Ptr fluid_props;
	
	std::shared_ptr < lpt::SharedObjects > shared_objects;
};


class ParticlesVTK {
public:
	typedef std::shared_ptr<ParticlesVTK> Ptr; 
	static inline ParticlesVTK::Ptr create(vtkSmartPointer<vtkRenderer> renderer) { return ParticlesVTK::Ptr(new ParticlesVTK(renderer) ); }
	
	class Parameters {
	public:
		Parameters() : scale(1.0) { 
			velocity_range[0] = 0; velocity_range[1] = 2.5;
			acceleration_range[0] = 0; acceleration_range[1] = 30;
		}
		double scale;
		array<double,2> velocity_range;
		array<double,2> acceleration_range;
	};
	
	Parameters params;

	ParticlesVTK(vtkSmartPointer<vtkRenderer> renderer);

	void resizeGlyphArrays(int size) {
		points->SetNumberOfPoints( size );
		velocity_magnitudes->SetNumberOfTuples( size );
		acceleration_magnitudes->SetNumberOfTuples( size );
	}

	void setScalarsModified() {
		points->Modified();
		velocity_magnitudes->Modified();
		acceleration_magnitudes->Modified();
	}

	void updateParticle( int id, lpt::ParticleVectors& current_particle);

	inline void setLookupTable(vtkSmartPointer<vtkLookupTable> lut) {
		lookuptable = lut;
		glyph3Dmapper->SetLookupTable( lookuptable );
		glyph3Dmapper->Update();
	}

	inline void setScalarBarWidget(vtkSmartPointer<vtkScalarBarWidget> bar) {
		scalarbar = bar;
	}

	inline void addToRenderer() {
		renderer->AddActor(glyphactor);
		//scalarbar->On();
	}

	inline void removeFromRenderer() {
		renderer->RemoveActor(glyphactor);
		//scalarbar->Off();
	}

	inline void setVectorMode(lpt::VectorMode mode);

	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkSphereSource> sphere;
	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkPolyData> polydata;	
	vtkSmartPointer<vtkPolyDataAlgorithm> source;
	vtkSmartPointer<vtkLookupTable> lookuptable; 
	vtkSmartPointer<vtkDoubleArray> velocity_magnitudes;
	vtkSmartPointer<vtkDoubleArray> acceleration_magnitudes;
	vtkSmartPointer<vtkGlyph3DMapper> glyph3Dmapper;
	vtkSmartPointer<vtkScalarBarWidget> scalarbar;
	vtkSmartPointer<vtkActor> glyphactor;
	lpt::VectorMode vector_mode;
};

class TrajectoryHandler : public vtkCommand
{
public:
	typedef std::shared_ptr<TrajectoryHandler> Ptr; 
	static inline TrajectoryHandler::Ptr create(vtkSmartPointer<vtkRenderer> renderer) { return TrajectoryHandler::Ptr(new TrajectoryHandler(renderer) ); }
	
	class Parameters {
	public:
		Parameters() : traj_update_stride(0), grid_update_stride(100), queue_capacity(500) {}
		int traj_update_stride;
		int grid_update_stride;
		int queue_capacity;
	} params;

	TrajectoryHandler(vtkSmartPointer<vtkRenderer> renderer) : 
		renderer(renderer), tick_count1(0), 
		view_paths(false), clear_trajs(false), reset_volume_grid(false), update_volume_grid(false), 
		add_vector_view(false), add_traj_view(true), remove_vector_view(false), 
		remove_traj_view(false), add_cameras_view(true), remove_cameras_view(false), save_plane(false), view_mode(lpt::TRAJECTORIES) {}

	virtual void Execute(vtkObject* caller, unsigned long eventId, void * vtkNotUsed(callData));
	virtual inline void addControls();

	void clearTrajectoryPaths();
	void updateTrajectoryVisualization(  pair< vector <lpt::Trajectory3d*>, vector <lpt::ParticleVectors > >& traj_updates);
	
	inline void pushToRenderQueue( pair< vector < lpt::Trajectory3d* >, vector < ParticleVectors > >& traj_update ){
		if (render_queue.size() < params.queue_capacity) 
			render_queue.push(traj_update);
	}

	inline void setViewMode(ViewMode mode);
	inline void setFiniteVolumeGrid(lpt::FiniteVolumeGrid* grid) { volume_grid = grid; }
	inline void setTrajectoryGlyphs(lpt::ParticlesVTK* glyphs) { traj_glyphs = glyphs; }
	inline int getQueueSize() { return render_queue.size(); }
	inline void setViewPaths(bool state) {view_paths = state;}
	inline bool getViewPaths(){return view_paths;}
	inline void setClearView() { clear_trajs = true; }
	inline void setResetVolumeGrid() { reset_volume_grid = true; }
	inline void setUpdateVolumeGrid() { update_volume_grid = true; }
	inline void setClearQueue() { clear_queue = true; }
	inline void setCamerasVTK(vector<lpt::CameraVTK>& cameras) { camerasvtk = &cameras; }
	inline lpt::ViewMode getViewMode(){return view_mode;}
	
	inline void setVectorMode(lpt::VectorMode mode) {
		traj_glyphs->setVectorMode(mode);
		volume_grid->setVectorMode(mode);
		lpt::TrajectoryPathVTK::setVectorMode(mode);
	}
		
	inline void setSharedObjects( std::shared_ptr < lpt::SharedObjects > new_shared_objects ) { 
		shared_objects = new_shared_objects; 
	}
	
	inline void setQueueCapacity(int cap) { 
		params.queue_capacity = cap;
		render_queue.setCapacity(cap);
	}
	
	inline void setSavePlane(bool state) { 
		save_plane = state;
		cout << "Save plane data: " << ( state ? " ON " : " OFF ") << endl;
	}	
	
	inline void setDisplayCameras(bool state) { 
		state ? add_cameras_view = true : remove_cameras_view = true;	
		cout << "Camera Display: " << ( state ? " ON " : " OFF ") << endl;
	}
	
	void addCamerasToRenderer() {
		for (int c = 0; c < camerasvtk->size(); ++c)
				camerasvtk->operator[](c).addToRenderer(renderer);
		add_cameras_view = false;
	}
	
	void removeCamerasFromRenderer() {
		for (int c = 0; c < camerasvtk->size(); ++c)
				camerasvtk->operator[](c).removeFromRenderer(renderer);
		remove_cameras_view = false;
	}

	virtual ~TrajectoryHandler() {}
	
	friend void callbackResetVolumeGrid( int state, void* data ) {
		TrajectoryHandler* trajhandler = static_cast<TrajectoryHandler*>(data);
		trajhandler->setResetVolumeGrid();
	}
	friend void callbackUpdateVolumeGrid( int state, void* data ) {
		TrajectoryHandler* trajhandler = static_cast<TrajectoryHandler*>(data);
		trajhandler->setUpdateVolumeGrid();
	}
	friend void callbackClearTrajView( int state, void* data ) {
		TrajectoryHandler* trajhandler = static_cast<TrajectoryHandler*>(data);
		trajhandler->setClearView();
	}
	
	friend void callbackSetDisplayCameras( int state, void* data ) {
		TrajectoryHandler* trajhandler = static_cast<TrajectoryHandler*>(data);
		trajhandler->setDisplayCameras(state);
	}

	friend void callbackSavePlane( int state, void* data ) {
		TrajectoryHandler* trajhandler = static_cast<TrajectoryHandler*>(data);
		trajhandler->setSavePlane(true);
	}
	
private:
	vtkSmartPointer < vtkRenderer > renderer;
	lpt::ParticlesVTK* traj_glyphs;
	lpt::FiniteVolumeGrid* volume_grid;
	vector<lpt::CameraVTK>* camerasvtk;
	list<lpt::Trajectory3d*> current_traj_list;

	bool save_plane;
	bool clear_queue;
	bool add_traj_view, add_vector_view;
	bool remove_traj_view, remove_vector_view;
	bool view_paths, clear_trajs;
	bool update_volume_grid, reset_volume_grid;
	bool add_cameras_view, remove_cameras_view;
		
	std::shared_ptr < lpt::SharedObjects > shared_objects;

	int tick_count1;

	lpt::ViewMode view_mode;
	lpt::concurrent_queue< pair< vector< lpt::Trajectory3d*>, vector< ParticleVectors > > > render_queue;
};


class VisualizerInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
	static VisualizerInteractorStyle* New();
	vtkTypeMacro(VisualizerInteractorStyle, vtkInteractorStyleTrackballCamera);

	VisualizerInteractorStyle() { 
		vector_modes.push_back(lpt::VELOCITY);
		vector_modes.push_back(lpt::ACCELERATION);
		vector_modes.push_back(lpt::TURBULENT_KINETIC_ENERGY);
		vector_modes.push_back(lpt::VORTICITY);
		vector_modes.push_back(lpt::MASS_RESIDUAL);
		vector_modes.push_back(lpt::VELOCITY_SD);
		vector_modes.push_back(lpt::ACCELERATION_SD);
		vector_modes.push_back(lpt::COUNT);
		vector_modes.push_back(lpt::PRESSURE);
		vector_modes.push_back(lpt::PRESSURE_LAGRANGIAN);
		vector_mode_iter = vector_modes.begin();
	}

	virtual void OnKeyPress();

	lpt::FiniteVolumeGrid* grid;
	lpt::TrajectoryHandler* traj_handler;

private:
	vector<lpt::VectorMode> vector_modes;
	vector<lpt::VectorMode>::iterator vector_mode_iter;
};

class Visualizer { 
public: 
	typedef std::shared_ptr<Visualizer> Ptr;
	static inline Visualizer::Ptr create() { return std::make_shared<lpt::Visualizer>(); }
	
	class Parameters {
	public:
		Parameters() : mode(0), stride(2), scale(30), timer_duration(1), queue_capacity(500), min_traj_size(6) {
			array<int, 2> size = {{ 640, 480 }};
			window_size = size;
			
			array<double,3> color = {{.1,.1,.1}};
			background_color = color;
		}

		int mode;
		array<int,2> window_size;
		array<double,3> background_color;
		int scale;
		int stride;
		int timer_duration;
		int queue_capacity;
		int min_traj_size;
	} params;
	
	Visualizer();

	virtual inline void initialize(); 

	virtual inline void start();

	virtual inline void stop();

	virtual inline void addControls();

	void manageQueue();

	void addTrajectoriesToQueue( list<lpt::Trajectory3d_Ptr>& active_trajs );
	
	void takeMeasurement(const vector < lpt::ParticleVectors >& particle_vectors);

	void accumulateCentroidDetectionUncertainty( vector<lpt::Match::Ptr>& matches );
		
	void setCameras(vector<lpt::Camera>& cameras) { 
		for (int c = 0; c < cameras.size(); ++c)
			camerasvtk.emplace_back(cameras[c]);
		handler->setCamerasVTK(camerasvtk);
		centroid_uncertainty_accumulators.resize(cameras.size());
	}

	inline bool getTakeMeasurement() const { return take_measurement; }
	inline bool getTakeImageMeasurement() const { return image_measurement; }
	inline bool getVisualizationStatus() const { return visualization_status; }
	inline void setAccumulate(bool status) { accumulate_status = status;}
	inline int getQueueSize() { return traj_queue.size(); }
	inline int getRenderQueueSize() { return handler->getQueueSize(); }
	inline lpt::FiniteVolumeGrid::Ptr getVolumeGrid() { return volume_grid; }

	inline void setSharedObjects( std::shared_ptr < lpt::SharedObjects > new_shared_objects ) { 
		cout << "Visualizer Setting Shared Objects " << endl;
		shared_objects = new_shared_objects; 
		handler->setSharedObjects(new_shared_objects);
		volume_grid->setSharedObjects(new_shared_objects);
	}
	
	inline void setVisualizationStatus(bool state) { 
		visualization_status = state;
		cout << "Visualization Status: " << ( state ? " ON " : " OFF ") << endl;
	}

	inline void setTakeMeasurement(bool state) { 
		take_measurement = state;
		image_measurement = state;
		cout << "Take measurement Status: " << ( state ? " ON " : " OFF ") << endl;
	}	
		
	~Visualizer() {	cout << "Visualizer closed" << endl; }

	friend void callbackSetVisualizationStatus( int state, void* data ) {
		Visualizer* visualizer = static_cast<Visualizer*>(data);
		visualizer->setVisualizationStatus(state);
	}

	friend void callbackSetStride(int state, void* data) {
		Visualizer* visualizer = static_cast<Visualizer*>(data);
	}

	friend void callbackSetAccumulate( int state, void* data ) {
		Visualizer* visualizer = static_cast<Visualizer*>(data);
		visualizer->setAccumulate(state);
	}
	
	friend void callbackTakeMeasurement( int state, void* data ) {
		Visualizer* visualizer = static_cast<Visualizer*>(data);
		visualizer->setTakeMeasurement(true);
	}

	
private:
	lpt::concurrent_queue< vector< pair<lpt::Trajectory3d*, vector<lpt::Particle3d_Ptr>::iterator > > > traj_queue;
	boost::thread queue_manager;
	
	lpt::TrajectoryHandler::Ptr handler;
	lpt::FiniteVolumeGrid::Ptr volume_grid;
	lpt::ParticlesVTK::Ptr traj_glyphs;
	lpt::CoordinateArrows::Ptr coordinates;
	lpt::FluidProperties::Ptr fluid_props;
	vector<lpt::CameraVTK> camerasvtk;

	vtkSmartPointer < VisualizerInteractorStyle > style;
	vtkSmartPointer < vtkRenderWindowInteractor > window_interactor;	
	vtkSmartPointer < vtkRenderWindow > render_window;
	vtkSmartPointer < vtkRenderer > renderer;

	std::shared_ptr < lpt::SharedObjects > shared_objects; 
	
	int timer_id;

	bool visualization_status;	
	bool accumulate_status;
	bool take_measurement;
	bool image_measurement;

	lpt::boost_accumulator distance_accumulator;
	array<lpt::boost_accumulator, 2> scalar_accumulator;
	vector<array<lpt::boost_accumulator, 3>> position_accumulators;
	vector<vector<array<lpt::boost_accumulator, 2>>> centroid_uncertainty_accumulators; // [cam_id][particle_id][pixel dimension]
	int measurement_count;
	int image_measurement_count;
};

// Potential new functionality 

class HistogramVTK : vtkCommand {
public:
	typedef	shared_ptr<HistogramVTK> Ptr;
	static HistogramVTK::Ptr create() { return HistogramVTK::Ptr(new HistogramVTK() ); }

	HistogramVTK() : number_of_bins(100) {
		view = vtkSmartPointer<vtkContextView>::New();
		histogram = vtkSmartPointer<vtkChartXY>::New();
		table = vtkSmartPointer<vtkTable>::New();
		counts = vtkSmartPointer<vtkIntArray>::New();
		bin_values = vtkSmartPointer<vtkDoubleArray>::New();
		
		counts->SetName("Counts");
		bin_values->SetName("Bin Range");
		table->AddColumn(bin_values);
		table->AddColumn(counts);
		bin_size = 1; //default value
		view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
		view->GetRenderWindow()->SetSize(400, 300);
		view->GetScene()->AddItem(histogram);
		histogram->SetTitle("Histogram");
		histogram->GetTitleProperties()->BoldOn();

		histogram->GetAxis(1)->SetTitle("Residual level");
		histogram->GetAxis(0)->SetTitle("Count");
		vtkPlot* plot = histogram->AddPlot(vtkChart::BAR);
		plot->SetInput(table, 0, 1);
		plot->SetColor(0, 0, 255, 255);
		//plot->GetYAxis()->SetLogScale(true);
		histogram->SetAutoSize(false);
		histogram->SetSize(vtkRectf(350, 0, 350, 250));

		chart_scene = vtkSmartPointer<vtkContextScene>::New();
		chart_actor = vtkSmartPointer<vtkContextActor>::New();
		chart_scene->AddItem( histogram.GetPointer() );
		chart_actor->SetScene( chart_scene.GetPointer() );
		range[0] = 0;
		range[1] = 100;
		setBins(number_of_bins, range);
		view->AddObserver(vtkCommand::UserEvent, this);
	} 
	
	virtual void Execute(vtkObject* caller, unsigned long eventId, void * callData)
	{
		vtkContextView* view_caller = static_cast<vtkContextView*>(caller);
		if (eventId == vtkCommand::UserEvent) {
			view_caller->Update();
			view_caller->Render();
		}
	}

	void setBins(int num_bins, double bin_range[2]) {
		number_of_bins = num_bins;
		range[0] = bin_range[0];
		range[1] = bin_range[1];
		bins.resize(number_of_bins, 0);
		count_data.resize(number_of_bins, 0);
		bin_size = (range[1] - range[0]) / static_cast<double>(number_of_bins);
		for (int i = 0; i < bins.size(); ++i) 
			bins[i] = bin_size * (i+1);
	}

	void setAxisColor(int color[3]) {
		histogram->GetTitleProperties()->SetColor(color[0]/255.0, color[1]/255.0, color[2]/255.0);
		histogram->GetAxis(0)->GetPen()->SetColor(color[0], color[1], color[2]);
		histogram->GetAxis(1)->GetPen()->SetColor(color[0], color[1], color[2]);
		histogram->GetAxis(0)->GetTitleProperties()->SetColor(color[0]/255.0, color[1]/255.0, color[2]/255.0);
		histogram->GetAxis(1)->GetTitleProperties()->SetColor(color[0]/255.0, color[1]/255.0, color[2]/255.0);
		histogram->GetAxis(1)->GetGridPen()->SetColor(color[0], color[1], color[2],255);
		histogram->GetAxis(0)->GetGridPen()->SetColor(color[0], color[1], color[2],255);
		histogram->GetAxis(1)->GetGridPen()->SetColor(color[0], color[1], color[2],255);
		histogram->GetAxis(0)->GetLabelProperties()->SetColor(color[0]/255.0, color[1]/255.0, color[2]/255.0);
		histogram->GetAxis(1)->GetLabelProperties()->SetColor(color[0]/255.0, color[1]/255.0, color[2]/255.0);
	}

	void setTextStrings(string& title, string& x_title, string& y_title) {
		histogram->SetTitle( title.c_str() );
		histogram->GetAxis(1)->SetTitle( x_title.c_str() );
		histogram->GetAxis(0)->SetTitle( y_title.c_str() );
	}

	void updateData(vector<double>& data) {
		for (int i = 0; i < data.size(); i++) 
			count_data[ getBinID( data[i] ) ] = count_data[ getBinID( data[i] ) ] + 1 ;
		
		table->SetNumberOfRows( bins.size() );
		for (int i = 0; i < bins.size(); i++) {
			table->SetValue(i, 0, bins[i] - bin_size / 2);
			table->SetValue(i, 1, count_data[i] );
		}
		counts->Modified();
		bin_values->Modified();
		table->Modified();		
		histogram->RecalculateBounds();
		histogram->Modified();
		//view->InvokeEvent(vtkCommand::UserEvent);
	}

	int getBinID(double value) {
		int id;
		for ( id = 0; id < bins.size(); id++) {
			if ( bins[id] > value ) 
				break;
		}
		return id;
	}
	
	void addToRenderer(vtkSmartPointer<vtkRenderer> renderer) {
		renderer->AddActor( chart_actor.GetPointer() );
		chart_scene->SetRenderer(renderer);
	}

	void startChartWindow() {
		if (! view->GetInteractor()->GetEnabled() ) {
			view->GetInteractor()->Initialize();
			view->GetInteractor()->Start();
		}
	}
	void stopChartWindow() {
		if ( view->GetInteractor()->GetEnabled() ) {
			cout <<" find a way to close the window without terminating the app" << endl;
		}
	}
	void removeFromRenderer(vtkSmartPointer<vtkRenderer> renderer) {
		renderer->RemoveActor( chart_actor.GetPointer() );
	}

private:
	double bin_size;
	double range[2];
	int number_of_bins;
	vector<double> bins;
	vector<int> count_data;
	vtkSmartPointer<vtkContextView> view;
	vtkSmartPointer<vtkChartXY> histogram;
	vtkSmartPointer<vtkTable> table;
	vtkSmartPointer<vtkDoubleArray> bin_values;
	vtkSmartPointer<vtkIntArray> counts;
	vtkSmartPointer<vtkContextScene> chart_scene;
	vtkSmartPointer<vtkContextActor> chart_actor;
};

class StreamLines {
public:
	typedef std::shared_ptr<StreamLines> Ptr; 
	static inline StreamLines::Ptr create(vtkAlgorithmOutput* output_port) { return StreamLines::Ptr(new StreamLines(output_port) ); }
	StreamLines(vtkAlgorithmOutput* output_port) {
		// Source of the streamlines
		seeds = vtkSmartPointer<vtkPlaneSource>::New();
		seeds->SetXResolution(8);
		seeds->SetYResolution(8);
		
		// Streamline itself
		streamline = vtkSmartPointer<vtkStreamLine>::New();
		//streamline->SetSource(seeds->GetOutput());
		streamline->SetInputConnection(output_port);
		
		//streamLine->SetStartPosition(2,-2,30);
		// as alternative to the SetSource(), which can handle multiple
		// streamlines, you can set a SINGLE streamline from
		// SetStartPosition()
		streamline->SetStartPosition(0,0,0);
		streamline->SetMaximumPropagationTime(1000);
		streamline->SetIntegrationStepLength(0.5);
		streamline->SetStepLength(0.5);
		streamline->SetIntegrationDirectionToIntegrateBothDirections(); 
		streamline->SetNumberOfThreads(1);
		streamline->VorticityOn();
		streamline->Update();

		streamline_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		streamline_mapper->SetInputConnection( streamline->GetOutputPort() );

		streamline_actor = vtkSmartPointer<vtkActor>::New();
		streamline_actor->SetMapper(streamline_mapper);
		streamline_actor->VisibilityOn();
	}
	inline void setSeedPlane(double origin[3], double point1[3], double point2[3]) {
		seeds->SetOrigin(origin);
		seeds->SetPoint1(point1);
		seeds->SetPoint2(point2);
		seeds->Modified();
		seeds->Update();
	}
	inline void addToRenderer(vtkSmartPointer<vtkRenderer> renderer) {
		streamline->Modified();
		streamline->Update();
		streamline_actor->Modified();
		renderer->AddActor(streamline_actor);
		
	}
	inline void removeFromRenderer(vtkSmartPointer<vtkRenderer> renderer) {
		renderer->RemoveActor(streamline_actor);
	}
private:
	vtkSmartPointer<vtkPlaneSource> seeds;
	vtkSmartPointer<vtkStreamLine> streamline;
	vtkSmartPointer<vtkPolyDataMapper> streamline_mapper;
	vtkSmartPointer<vtkActor> streamline_actor; 
};

class PickDim : public vtkCommand {
public:
	typedef std::shared_ptr<PickDim> Ptr; 
	static inline PickDim::Ptr create() { return PickDim::Ptr(new PickDim() ); }
	
	PickDim() {
		index = 0;
		initial = true;
		text = vtkSmartPointer<vtkVectorText>::New();
		mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		mapper->SetInputConnection( text->GetOutputPort() );
		follower = vtkSmartPointer<vtkFollower>::New();
		follower->SetMapper( mapper );
		follower->GetProperty()->SetColor( 1, 0, 0 ); // red
		points.resize(2);
	}
	virtual void Execute(vtkObject* caller, unsigned long eventId, void * vtkNotUsed(callData))
	{
		switch (eventId) {
		case vtkCommand::EndPickEvent: 
			{
				vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
				double* clickPos = interactor->GetPicker()->GetSelectionPoint();
				vtkSmartPointer<vtkPropPicker>  picker =
					vtkSmartPointer<vtkPropPicker>::New();
				picker->Pick(clickPos, interactor->GetInteractorStyle()->GetDefaultRenderer() );
				vtkActor* actor = picker->GetActor();
				if (actor) { 
					points[index] = actor;
					//double* pos = points[index]->GetPosition();
					//textstream.str("");
					//textstream << "["<< pos[0] << " " << pos[1] << " " << pos[2] << "]";
					//cout << textstream.str() << endl;
					//text->SetText( textstream.str().c_str() );					
					if (index == 0) {
						if (initial == true ) {
							follower->SetPosition( 100, 100 , 0 );//points[index]->GetPosition()[0] + 10, points[index]->GetPosition()[1] + 10, points[index]->GetPosition()[2]); 
							follower->SetCamera( interactor->GetInteractorStyle()->GetDefaultRenderer()->GetActiveCamera() );
						}
						follower->SetScale( 5 * (*actor->GetScale()) );
						follower->Modified();
						++index;
					}
					else {
						interactor->GetInteractorStyle()->GetDefaultRenderer()->AddActor( follower );
						index = 0;
						if (initial == true ) {
							interactor->AddObserver(vtkCommand::TimerEvent, this);
							initial = false;
						}
					}
					//interactor->Render();
				}
			break;
			} 
		case vtkCommand::TimerEvent: {
			vtkRenderWindowInteractor* interactor = static_cast<vtkRenderWindowInteractor*>(caller);
			if (points[0]->GetReferenceCount() > 1 && points[1]->GetReferenceCount() > 1) {
				double* pos1 = points[0]->GetPosition();
				double* pos2 = points[1]->GetPosition();
				double dist = sqrt( (pos1[0] - pos2[0]) * (pos1[0] - pos2[0]) + (pos1[1] - pos2[1]) * (pos1[1] - pos2[1]) + (pos1[2] - pos2[2]) * (pos1[2] - pos2[2]) ); 
				textstream.str("");
				textstream << "Distance = " << dist << " mm";
				//cout << textstream.str() << endl;
				text->SetText( textstream.str().c_str() );
				//follower->SetPosition( pos1[0] + 10, pos1[1] + 10, pos1[2]); 
			} else {
				textstream.str("");
				//interactor->RemoveObserver(this);
				interactor->GetInteractorStyle()->GetDefaultRenderer()->RemoveActor( follower );
			}
									  }
		break;
		default:
			break;
		}
	}
private:
	int index;
	bool initial;
	std::stringstream textstream;
	vector<vtkActor*> points;
	vtkSmartPointer<vtkVectorText> text;
	vtkSmartPointer<vtkFollower> follower;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
};

}/* NAMESPACE_PT */

#endif /*VISUALIZATION_H_*/

