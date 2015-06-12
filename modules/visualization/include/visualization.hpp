/**
 * @file visualization.hpp
 * The visualization module declaration
 * Contains classes for fluid properties, flow variables, and particle/trajectory display
 */
#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_

#include <core.hpp>

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
enum VectorMode {VELOCITY, ACCELERATION, TURBULENT_KINETIC_ENERGY, VORTICITY, MASS_RESIDUAL, VELOCITY_SD, ACCELERATION_SD, COUNT, PRESSURE, PRESSURE_LAGRANGIAN, CORRECTED_VELOCITY, TURBULENCE_DISSIPATION_RATE};

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

/**
 * @brief The FluidProperties class
 */
class FluidProperties {
public:
	typedef	shared_ptr<FluidProperties> Ptr;
    static inline FluidProperties::Ptr create() { return std::make_shared<FluidProperties>(); }

    ///< Default air at 15 C
    FluidProperties();
    FluidProperties(double n_rho, double n_Td, double n_Tw, double n_mu, double n_P_ref);

	double rho;  //density (kg/m3)
	double Td;   //dry bulb temperature (C)
	double Tw;   //wet bulb temperature (C)
	double mu;   // dynamic viscocity ( kg/(m·s) )
	double P_ref; // Reference Pressure (Pa)
};

/**
 * @brief The CoordinateArrows class
 */
class CoordinateArrows {
public:
	typedef std::shared_ptr<CoordinateArrows> Ptr; 
    static inline CoordinateArrows::Ptr create(	vtkSmartPointer<vtkRenderWindowInteractor> iren, double scale = 1) { return std::make_shared<CoordinateArrows>(iren, scale); }

	CoordinateArrows( vtkSmartPointer<vtkRenderWindowInteractor> iren, double scale);
    virtual ~CoordinateArrows();
private:
	double scale;
	vtkSmartPointer<vtkRenderWindowInteractor> interactor;
	vtkSmartPointer < vtkArrowSource > arrow;
	vector < vtkSmartPointer < vtkActor >> actors;
	vtkSmartPointer < vtkPolyDataMapper > mapper;

	vtkSmartPointer<vtkAxesActor> axes;
	vtkSmartPointer<vtkOrientationMarkerWidget> widget; 
};

/**
 * @brief The CameraVTK class
 */
class CameraVTK {
public:
    /**
     * @brief Get Random Color
     * @return A 3 by 1 array containing the random color
     */
    static array<double, 3> getRandomColor();

    /**
     * @brief CameraVTK constructor
     * @param camera
     */
	CameraVTK(lpt::Camera& camera);

    /**
     * @brief Calculate Plane Equation
     * @param A
     * @param B
     * @param C
     * @param eq The resulting Plane equation
     */
	void calcPlaneEq(array<double,3>& A, array<double,3>& B, array<double,3>& C, array<double,4>& eq );

    /**
     * @brief Project from camera coordinate to world coordinate
     * @param cam Camera
     * @param pc Camera coordinate
     * @param pw World coordinate
     */
	void camToWorld(lpt::Camera& cam, array<double, 3>& pc, array<double, 3>& pw);

    /**
     * @brief Add actors to renderer
     * @param renderer
     */
	void addToRenderer(vtkSmartPointer<vtkRenderer> renderer); 

    /**
     * @brief Remove actors from renderer
     * @param renderer
     */
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

/**
 * @brief The TrajectoryPathVTK class
 * Derived from the TrajectoryVTKBase class
 */
class TrajectoryPathVTK : public lpt::TrajectoryVTKBase {
public:
	typedef std::shared_ptr<TrajectoryPathVTK> Ptr; 
    static inline TrajectoryPathVTK::Ptr create(vtkRenderer* renderer) { return std::make_shared<lpt::TrajectoryPathVTK>(renderer); }
	
	static vtkSmartPointer<vtkLookupTable> lookuptable;
	static int max_points;
	static bool show_paths;
	static int scalar_index;

    /**
     * @brief TrajectoryPathVTK constructor
     * @param renderer
     */
    TrajectoryPathVTK(vtkRenderer* renderer);
    /**
     * @brief TrajectoryPathVTK destructor
     */
    virtual ~TrajectoryPathVTK();

    /**
     * @brief Initialize Lookup Table
     */
	static inline void initializeLookupTable() {
		lookuptable = vtkSmartPointer<vtkLookupTable>::New();
		lookuptable->SetRange( 0, 2.0 );
		lookuptable->SetNumberOfColors(256);
		lookuptable->SetHueRange(0.667, 0.0);
		lookuptable->Build();
	}

    /**
     * @brief Set Lookup Table
     * @param lut Lookup Table
     */
	static inline void setLookupTable(vtkSmartPointer<vtkLookupTable> lut) { lookuptable = lut; }

    /**
     * @brief Set Vector Mode
     * @param mode
     */
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

    /**
     * @brief Add next point
     * @param current_particle
     */
	void addNextPoint(lpt::ParticleVectors& current_particle);

    /**
     * @brief Remove actors from renderer
     */
	void removeFromRenderer();

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

/**
 * @brief The PressureFieldSolver class
 * Base class for solving pressure field
 */
class PressureFieldSolver {
public:
	typedef	std::shared_ptr<PressureFieldSolver> Ptr;

    /**
     * @brief PressureFieldSolver constructor
     */
    PressureFieldSolver();

    /**
     * @brief PressureFieldSolver destructor
     */
    virtual ~PressureFieldSolver();

    /**
     * @brief Solve pressure field
     * @param velocity_field Velocity
     * @param acceleration_field Acceleration
     * @param stress_field Reynolds stress
     */
	virtual void solve(lpt::VectorField& velocity_field, lpt::VectorField& acceleration_field, lpt::ReynoldsStressField& stress_field)=0;
	
    /**
     * @brief addControls
     */
    virtual void addControls();
	
    /**
     * @brief Set Fluid Properties
     * @param props Fluid properties
     */
    virtual void setFluidProperties(FluidProperties::Ptr props);

    /**
     * @brief Set Grid Properties
     * @param grid_cell_counts Number of grid cells in each direction
     * @param cell_dimensions Size of each grid cell in meter
     */
    virtual void setGridProperties(array<int, 3>& grid_cell_counts, array<double, 3>& cell_dimensions);
	
    /**
     * @brief Set SharedObjects
     * @param new_shared_objects
     */
    inline void setSharedObjects( lpt::SharedObjects::Ptr new_shared_objects ) {
		shared_objects = new_shared_objects; 
	}

    /**
     * @brief Get Pressure Field
     * @return Pressure field
     */
	inline const vector<double>& getPressureField() { return pressure; }

    /**
     * @brief Reset Pressure Field
     */
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

/**
 * @brief The PressurePoissonSolver2Step class
 * Derived from the PressureFieldSolver class
 * Solves pressure poisson equation using a 2-step method
 * First calculates cell pressure gradient from the RANS equation
 * Then calculates cell pressure by solving the pressure Poisson equation with pressure gradient as source terms
 */
class PressurePoissonSolver2Step : public PressureFieldSolver {
public:
	typedef	std::shared_ptr<PressurePoissonSolver2Step> Ptr;
    static inline PressurePoissonSolver2Step::Ptr create() { return std::make_shared<PressurePoissonSolver2Step>(); }
	
    /**
     * @brief The Parameters class
     * Parameters for successvie over relaxation solver
     */
	class Parameters {
	public:
		Parameters() : sor_weight(1.2), max_iterations(20000), convergence_tolerance(1E-6) {}
		double sor_weight;
		double max_iterations;
		double convergence_tolerance;
    } params;

    /**
     * @brief PressurePoissonSolver2Step constructor
     */
    PressurePoissonSolver2Step();

    /**
     * @brief PressurePoissonSolver2Step destructor
     */
    virtual ~PressurePoissonSolver2Step();

    /**
     * @brief Solve for pressure field
     * Iteratively calls the calcCellPressure function over the entire domain
     * @param velocities Velocity
     * @param accelerations Acceleraton
     * @param stresses Reynolds stress
     */
	virtual void solve(lpt::VectorField& velocities, lpt::VectorField& accelerations, lpt::ReynoldsStressField& stresses);

    /**
     * @brief addControls
     */
	virtual void addControls();

private:
    /**
     * @brief Helper function. Calculate cell pressure by solving pressure Poisson equation
     * @param i Cell index in X-direction
     * @param j Cell index in Y-direction
     * @param k Cell index in Z-direction
     */
	void calcCellPressure(int i, int j, int k);

    /**
     * @brief Helper function. Calculate cell pressure gradients using Steady RANS method
     * @param i Cell index in X-direction
     * @param j Cell index in Y-direction
     * @param k Cell index in Z-direction
     */
	void calcCellPressureGradientsSteadyRANS(int i, int j, int k);

	lpt::VectorField pressure_gradients;
};

/**
 * @brief The LagrangianPressureFieldSolver class
 * Derived from the PressureFieldSolver class
 * Acquires cell pressure gradient from the accumulator grid
 * Then calculates cell pressure following the same procedure as the PressurePoissonSolver2Step
 */
class LagrangianPressureFieldSolver : public PressureFieldSolver {
public:
	typedef	shared_ptr<LagrangianPressureFieldSolver> Ptr;
    static inline LagrangianPressureFieldSolver::Ptr create() { return std::make_shared<LagrangianPressureFieldSolver>(); }
			
    /**
     * @brief The Parameters class
     * Parameters for successive over relaxation solver
     */
	class Parameters {
	public:
		Parameters() : sor_weight(1.2), max_iterations(20000), convergence_tolerance(1E-6) {}
		double sor_weight;
		double max_iterations;
		double convergence_tolerance;
	} params;

    /**
     * @brief LagrangianPressureFieldSolver constructor
     */
    LagrangianPressureFieldSolver();

    /**
     * @brief ~LagrangianPressureFieldSolver destructor
     */
    virtual ~LagrangianPressureFieldSolver();

    /**
     * @brief Solve pressure field (not actually used)
     * @param velocities
     * @param accelerations
     * @param stresses
     */
	virtual void solve(lpt::VectorField& velocities, lpt::VectorField& accelerations, lpt::ReynoldsStressField& stresses);

    /**
     * @brief Solve pressure field
     * Iteratively calls the calcCellPressure function over the entire domain
     * @param pressure_grads Pressure gradients obtained from Lagrangian pressure gradients
     */
	virtual void solve(lpt::VectorField& pressure_grads);

    /**
     * @brief addControls
     */
	virtual void addControls();

private:
    /**
     * @brief Helper function. Calculate cell pressure by solving pressure Poisson equation
     * @param i Cell index in X-direction
     * @param j Cell index in Y-direction
     * @param k Cell index in Z-direction
     */
	void calcCellPressure(int i, int j, int k);
	
	lpt::VectorField pressure_gradients;
};

/**
 * @brief The FiniteVolumeGrid class
 * Manages the finite volume grid
 */
class FiniteVolumeGrid {
public:
	typedef	shared_ptr<FiniteVolumeGrid> Ptr;
    static inline FiniteVolumeGrid::Ptr create(vtkSmartPointer < vtkRenderer > renderer) { return std::make_shared<FiniteVolumeGrid>(renderer); }
	
	class Parameters {
	public:
		Parameters() { 
			scalar_range[0] = 0; 
			scalar_range[1] = 2.5;

			grid_cell_counts[0] = 50;
			grid_cell_counts[1] = 50;
			grid_cell_counts[2] = 50; 
			
			double length_size = 500;			// mm
			grid_dimensions[0] = length_size;	// mm
			grid_dimensions[1] = length_size;	// mm
			grid_dimensions[2] = length_size;	// mm

			cell_dimensions[0] = grid_dimensions[0] / static_cast<double>(grid_cell_counts[0]);		//mm
			cell_dimensions[1] = grid_dimensions[1] / static_cast<double>(grid_cell_counts[1]);		//mm
			cell_dimensions[2] = grid_dimensions[2] / static_cast<double>(grid_cell_counts[2]);		//mm
			
			grid_origin[0] = 0;		//mm 
			grid_origin[1] = 0;		//mm
			grid_origin[2] = 0;		//mm

			fixed_scale_factor = 0.65 * (*std::max_element(cell_dimensions.begin(), cell_dimensions.end())); 
		}

		array<double, 3> grid_dimensions;
		array<int, 3> grid_cell_counts;
		array<double, 3> cell_dimensions;
		array<double, 3> grid_origin;
		
		double fixed_scale_factor;
		double scalar_range[2];
	} params;

    /**
     * @brief FiniteVolumeGrid constructor
     * @param renderer
     */
	FiniteVolumeGrid (vtkSmartPointer < vtkRenderer > renderer);

    /**
     * @brief FiniteVolumeGrid destructor
     */
    ~FiniteVolumeGrid();

    /**
     * @brief initialize
     */
	void initialize();

    /**
     * @brief Set pressure field solver
     * @param solver Desired pressure solver
     */
	void setPressureFieldSolver(lpt::PressureFieldSolver::Ptr solver);

    /**
     * @brief Set Grid Origin
     * @param x X origin
     * @param y Y origin
     * @param z Z origin
     */
    void setGridOrigin(double x, double y, double z);

    /**
     * @brief Set Grid Cell Counts
     * @param nx X cell counts
     * @param ny Y cell counts
     * @param nz Z cell counts
     */
    void setGridCellCounts(int nx, int ny, int nz);

    /**
     * @brief Set Grid Dimensions
     * @param length_x X dimension
     * @param length_y Y dimension
     * @param length_z Z dimension
     */
    void setGridDimensions(double length_x, double length_y, double length_z);

    /**
     * @brief Update Grid
     */
    void updateGrid();

    /**
     * @brief Reset Grid
     */
    void resetGrid();

    /**
     * @brief Grid Is Full, i.e. every cell has at least one measurement
     * @return True if grid is full; false if not full
     */
    bool gridIsFull();

    /**
     * @brief updateAccumulators
     * @param traj_ptr
     * @param current_particle
     */
	void updateAccumulators(lpt::Trajectory3d* traj_ptr, lpt::ParticleVectors& current_particle );

    /**
     * @brief Save Plane Data
     * Save flow variables on 2-D planes (XY and XZ mid planes) to text files
     * TODO: take parameters to determine which plane to save
     */
	void savePlaneData();

    /**
     * @brief Calculates Pressure using 2-step pressure poisson
     */
	void calcPressures();

    /**
     * @brief Calculates Pressure using Lagrangian pressure gradient
     */
	void calcPressuresLagrangian();

    /**
     * @brief Calculates Mass Residuals
     */
	void calcMassResiduals();

    /**
     * @brief Calculates Vorticity
     */
	void calcVorticity();

    /**
     * @brief Get Lookup Table
     * @return Lookup Table
     */
    vtkSmartPointer<vtkLookupTable> getLookupTable();
    /**
     * @brief Get Scalar Bar Widget
     * @return Scalar bar
     */
    vtkSmartPointer<vtkScalarBarWidget> getScalarBarWidget();
    /**
     * @brief Get Implicit Plane
     * @return Implicit plane
     */
    vtkSmartPointer<vtkImplicitPlaneWidget> getImplicitPlane();
    /**
     * @brief Get Implicit Plane Actor
     * @return Implicit plane actor
     */
    vtkSmartPointer<vtkActor> getImplicitPlaneActor();
    /**
     * @brief Set Cell Neighborhood
     * @param neighborhood
     */
	inline void setCellNeighborhood(vector<array<int,3>>& neighborhood) {
		cell_neighborhood.resize(neighborhood.size());
		std::copy(neighborhood.begin(), neighborhood.end(), cell_neighborhood.begin());
	}

    /**
     * @brief Add Implicit Plane
     */
	inline void addImplicitPlane() {
		renderer->AddActor(planeactor);
		plane_widget->On();
		renderer->RemoveActor(glyphactor);
	}
    /**
     * @brief Remove Implicit Plane
     */
	inline void removeImplicitPlane() {
		renderer->RemoveActor(planeactor);
		plane_widget->Off();
		renderer->AddActor(glyphactor);
	}

    /**
     * @brief Add To Renderer
     */
	inline void addToRenderer() {
		renderer->AddActor(outlineactor);
		if ( !getImplicitPlane()->GetEnabled() )
			renderer->AddActor(glyphactor);
		scalarbar->On();
	}
    /**
     * @brief Remove From Renderer
     */
	inline void removeFromRenderer() {
		renderer->RemoveActor(outlineactor);
		if ( getImplicitPlane()->GetEnabled() )
			removeImplicitPlane();
		renderer->RemoveActor(glyphactor);
		scalarbar->Off();
	}
    /**
     * @brief Set Vector Mode
     * Invoked when switched between vector modes (left and right button)
     * @param mode
     */
	void setVectorMode(lpt::VectorMode mode);
    /**
     * @brief Get Vector Mode
     * @return Vecotr mode
     */
	inline lpt::VectorMode getVectorMode() {return vector_mode;}
	
    /**
     * @brief Set Fluid Properties
     * @param props Fluid properties
     */
	inline void setFluidProperties(FluidProperties::Ptr props) { 
		fluid_props = props; 
		pressure_solver->setFluidProperties(fluid_props);
		lagrangian_pressure_solver->setFluidProperties(fluid_props);
	}
    /**
     * @brief Get Bounds
     * @return
     */
	inline double* getBounds() { return grid->GetBounds(); }
	
    /**
     * @brief Set Shared Objects
     * @param new_shared_objects
     */
	inline void setSharedObjects( std::shared_ptr < lpt::SharedObjects > new_shared_objects ) { 
		shared_objects = new_shared_objects; 
		pressure_solver->setSharedObjects(new_shared_objects);
		lagrangian_pressure_solver->setSharedObjects(new_shared_objects);
	}	
	
private:
    /**
     * @brief Get 1-D cell index from 3-D index
     * @param i X cell index
     * @param j Y cell index
     * @param k Z cell index
     * @return 1-D cell id
     */
	inline int getGridIndex(int i, int j, int k) {
		return ( i + j * params.grid_cell_counts[0] + k * (params.grid_cell_counts[0] * params.grid_cell_counts[1]) );     //index of grid cell (i,j,k)
	}
	
    /**
     * @brief Compute 3-D Cell Coordinate from 1-D index
     * @param cell_id 1-D cell index
     * @param ijk 3-D array of cell coordinates
     */
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

/**
 * @brief The ParticlesVTK class
 * Manages particles display
 */
class ParticlesVTK {
public:
	typedef std::shared_ptr<ParticlesVTK> Ptr; 
    static inline ParticlesVTK::Ptr create(vtkSmartPointer<vtkRenderer> renderer) { return std::make_shared<ParticlesVTK>(renderer); }
	
	class Parameters {
	public:
		Parameters() : scale(1.0) { 
			velocity_range[0] = 0; velocity_range[1] = 2.5;
			acceleration_range[0] = 0; acceleration_range[1] = 30;
		}
		double scale;
		array<double,2> velocity_range;
		array<double,2> acceleration_range;
    } params;

    /**
     * @brief ParticlesVTK constructor
     * @param renderer
     */
	ParticlesVTK(vtkSmartPointer<vtkRenderer> renderer);
    /**
     * @brief ParticlesVTK destructor
     */
    ~ParticlesVTK();

    /**
     * @brief Resize Glyph Arrays
     * @param size New size
     */
    void resizeGlyphArrays(size_t size);

    /**
     * @brief Set Scalars when modified
     */
    void setScalarsModified();

    /**
     * @brief Updates Particle velocity and acceleration
     * @param id The particle id
     * @param current_particle The particle to be updated
     */
    void updateParticle(int id, lpt::ParticleVectors& current_particle);

    /**
     * @brief Set Lookup Table
     * @param lut Lookup Table
     */
	inline void setLookupTable(vtkSmartPointer<vtkLookupTable> lut) {
		lookuptable = lut;
		glyph3Dmapper->SetLookupTable( lookuptable );
		glyph3Dmapper->Update();
	}

    /**
     * @brief Set Scalar Bar Widget
     * @param bar Scalar bar
     */
	inline void setScalarBarWidget(vtkSmartPointer<vtkScalarBarWidget> bar) {
		scalarbar = bar;
	}

    /**
     * @brief Add To Renderer
     */
	inline void addToRenderer() {
		renderer->AddActor(glyphactor);
        //scalarbar->On();
	}

    /**
     * @brief Remove From Renderer
     */
	inline void removeFromRenderer() {
		renderer->RemoveActor(glyphactor);
        //scalarbar->Off();
	}

    /**
     * @brief Set Vector Mode by which particles are colored
     * @param mode Velocity or acceleration
     */
    void setVectorMode(lpt::VectorMode mode);

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

/**
 * @brief The TrajectoryHandler class
 * Manages Trajectories display
 */
class TrajectoryHandler : public vtkCommand
{
public:
	typedef std::shared_ptr<TrajectoryHandler> Ptr; 
    static inline TrajectoryHandler::Ptr create(vtkSmartPointer<vtkRenderer> renderer) { return TrajectoryHandler::Ptr(new TrajectoryHandler(renderer)); }
	
	class Parameters {
	public:
		Parameters() : traj_update_stride(0), grid_update_stride(100), queue_capacity(500) {}
		int traj_update_stride;
		int grid_update_stride;
		int queue_capacity;
	} params;

    /**
     * @brief TrajectoryHandler constructor
     * @param renderer
     */
    TrajectoryHandler(vtkSmartPointer<vtkRenderer> renderer);
    /**
     * @brief TrajectoryHandler destructor
     */
    virtual ~TrajectoryHandler();

    /**
     * @brief Execute
     * @param caller
     * @param eventId
     * @param vtkNotUsed
     */
	virtual void Execute(vtkObject* caller, unsigned long eventId, void * vtkNotUsed(callData));

    /**
     * @brief addControls
     */
    virtual void addControls();

    /**
     * @brief Clear Trajectory Paths
     * Rests each trajectory and Clear the current trajectory list
     */
	void clearTrajectoryPaths();

    /**
     * @brief Update Trajectory Visualization
     * @param traj_updates
     */
    void updateTrajectoryVisualization(pair< vector <lpt::Trajectory3d*>, vector <lpt::ParticleVectors > >& traj_updates);
	
    /**
     * @brief Push To Render Queue
     * @param traj_update
     */
	inline void pushToRenderQueue( pair< vector < lpt::Trajectory3d* >, vector < ParticleVectors > >& traj_update ){
		if (render_queue.size() < params.queue_capacity) 
			render_queue.push(traj_update);
	}

    /**
     * @brief Set View Mode
     * @param mode New view mode
     */
    void setViewMode(ViewMode mode);
	inline void setFiniteVolumeGrid(lpt::FiniteVolumeGrid* grid) { volume_grid = grid; }
	inline void setTrajectoryGlyphs(lpt::ParticlesVTK* glyphs) { traj_glyphs = glyphs; }
	inline size_t getQueueSize() { return render_queue.size(); }
	inline void setViewPaths(bool state) {view_paths = state;}
	inline bool getViewPaths(){return view_paths;}
	inline void setClearView() { clear_trajs = true; }
	inline void setResetVolumeGrid() { reset_volume_grid = true; }
	inline void setUpdateVolumeGrid() { update_volume_grid = true; }
	inline void setClearQueue() { clear_queue = true; }
	inline void setCamerasVTK(vector<lpt::CameraVTK>& cameras) { camerasvtk = &cameras; }
    inline lpt::ViewMode getViewMode() { return view_mode; }
	

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
	
    void addCamerasToRenderer();
	
    void removeCamerasFromRenderer();
	
    friend void callbackResetVolumeGrid( int state, void* data );
    friend void callbackUpdateVolumeGrid( int state, void* data );
    friend void callbackClearTrajView( int state, void* data );
    friend void callbackSetDisplayCameras( int state, void* data );
    friend void callbackSavePlane( int state, void* data );
	
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

/**
 * @brief The VisualizerInteractorStyle class
 * Defines keyboard interaction
 */
class VisualizerInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
    static VisualizerInteractorStyle* New();
	vtkTypeMacro(VisualizerInteractorStyle, vtkInteractorStyleTrackballCamera);

    /**
     * @brief VisualizerInteractorStyle constructor
     */
    VisualizerInteractorStyle();
    /**
     * @brief VisualizerInteractorStyle destructor
     */
    virtual ~VisualizerInteractorStyle();

    /**
     * @brief On Key Press
     * Defines the respond of each keyboard pressing
     */
	virtual void OnKeyPress();

	lpt::FiniteVolumeGrid* grid;
	lpt::TrajectoryHandler* traj_handler;

private:
	vector<lpt::VectorMode> vector_modes;
	vector<lpt::VectorMode>::iterator vector_mode_iter;
};

/**
 * @brief The Visualizer class
 */
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
	
    /**
     * @brief Visualizer constructor
     */
	Visualizer();
    /**
     * @brief Visualizer destructor
     */
    virtual ~Visualizer();

    /**
     * @brief initialize
     */
	virtual inline void initialize(); 
    /**
     * @brief start
     */
	virtual inline void start();
    /**
     * @brief stop
     */
	virtual inline void stop();
    /**
     * @brief addControls
     */
	virtual inline void addControls();

    /**
     * @brief manage trajectory queue
     * Calculates Lagrangian flow variables: velocity, accleration and pressure gradient
     */
	void manageQueue();

    /**
     * @brief Add Trajectories To Queue
     * @param active_trajs Trajectory to be added
     */
	void addTrajectoriesToQueue( list<lpt::Trajectory3d_Ptr>& active_trajs );

    /**
     * @brief Take Measurement
     * Accumulate position, displacement, velocity and acceleration magnitude
     * @param particle_vectors Vector of particles in the form of ParticleVectors
     */
	void takeMeasurement(const vector < lpt::ParticleVectors >& particle_vectors);

    /**
     * @brief Accumulate Centroid Detection Uncertainty
     * @param matches
     */
	void accumulateCentroidDetectionUncertainty( vector<lpt::Match::Ptr>& matches );
		
    /**
     * @brief Set Cameras
     * @param cameras
     */
    void setCameras(vector<lpt::Camera>& cameras);

	inline bool getTakeMeasurement() const { return take_measurement; }
	inline bool getTakeImageMeasurement() const { return image_measurement; }
	inline bool getVisualizationStatus() const { return visualization_status; }
	inline void setAccumulate(bool status) { accumulate_status = status;}
	inline size_t getQueueSize() { return traj_queue.size(); }
	inline size_t getRenderQueueSize() { return handler->getQueueSize(); }
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
		
    friend void callbackSetVisualizationStatus( int state, void* data );
    friend void callbackSetStride(int state, void* data);
    friend void callbackSetAccumulate( int state, void* data );
    friend void callbackTakeMeasurement( int state, void* data );

private:
	//lpt::concurrent_queue< vector< pair<lpt::Trajectory3d*, vector<array<double, 9>>::iterator > > > traj_queue;
    lpt::concurrent_queue< vector< pair< lpt::Trajectory3d*, vector< pair< lpt::Particle3d_Ptr, array<double, 9> > >::iterator > > > traj_queue;
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

class HistogramVTK : public vtkCommand {
public:
	typedef	shared_ptr<HistogramVTK> Ptr;
    static inline HistogramVTK::Ptr create() { return HistogramVTK::Ptr( new HistogramVTK() ); }

    HistogramVTK();
    virtual ~HistogramVTK();

    virtual void Execute(vtkObject* caller, unsigned long eventId, void * callData);

    void setBins(int num_bins, double bin_range[2]);

    void setAxisColor(int color[3]);

    void setTextStrings(string& title, string& x_title, string& y_title);

    void updateData(vector<double>& data);

    int getBinID(double value);
	
    void addToRenderer(vtkSmartPointer<vtkRenderer> renderer);

    void startChartWindow();

    void stopChartWindow();

    void removeFromRenderer(vtkSmartPointer<vtkRenderer> renderer);

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
    static inline StreamLines::Ptr create(vtkAlgorithmOutput* output_port) { return std::make_shared<StreamLines>(output_port); }

    StreamLines(vtkAlgorithmOutput* output_port);

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
    static inline PickDim::Ptr create() { return PickDim::Ptr( new PickDim() ); }
	
    PickDim();
    virtual ~PickDim();

    virtual void Execute(vtkObject* caller, unsigned long eventId, void * vtkNotUsed(callData));

private:
	int index;
	bool initial;
	std::stringstream textstream;
	vector<vtkActor*> points;
	vtkSmartPointer<vtkVectorText> text;
	vtkSmartPointer<vtkFollower> follower;
	vtkSmartPointer<vtkPolyDataMapper> mapper;
};

/*
class KalmanFilter {
public:
    typedef std::shared_ptr<KalmanFilter> Ptr;
    static inline KalmanFilter::Ptr create(std::shared_ptr<SharedObjects> new_shared_objects)
    { return std::make_shared<KalmanFilter>(new_shared_objects); }

    KalmanFilter(std::shared_ptr<SharedObjects> new_shared_objects);
    ~KalmanFilter();

    void filter();

    Eigen::Matrix<double, 6, 1> getState() const;

    void setState(Eigen::Matrix<double, 6, 1> new_state);

    void setObservation(Eigen::Matrix<double, 6, 1> new_observation);

    void setSharedObjects(std::shared_ptr<SharedObjects> new_shared_objects);

private:
    Eigen::Matrix<double, 6, 1> s;      //state vector
    Eigen::Matrix<double, 6, 1> z;      //observation
    Eigen::Matrix<double, 6, 6> F;      //transition model
    Eigen::Matrix<double, 6, 6> H;      //observation model
    Eigen::Matrix<double, 6, 6> P;      //error matrix
    Eigen::Matrix<double, 6, 6> Q;      //prediction uncertanty
    Eigen::Matrix<double, 6, 6> R;      //observation uncertanty

    std::shared_ptr<SharedObjects> shared_objects;
}; */

}/* NAMESPACE_LPT */

#endif /*VISUALIZATION_H_*/

