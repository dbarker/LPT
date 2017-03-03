/**
 * @file visualization.cpp
 * The visualization module definition
 */
#include "visualization.hpp"

namespace lpt {

using namespace std;

vtkSmartPointer<vtkLookupTable> TrajectoryPathVTK::lookuptable = 0;
int TrajectoryPathVTK::max_points = 10;
int TrajectoryPathVTK::scalar_index = 1;
bool TrajectoryPathVTK::show_paths = false;

vtkStandardNewMacro(VisualizerInteractorStyle);

/***** Definitions of global functions *****/

void getVectorField(vector<array<lpt::boost_accumulator, 3> >& accumulators, lpt::VectorField& vector_field)
{
    vector_field.resize( accumulators.size() );

    for (int p = 0; p < accumulators.size() ; ++p) {
        for (int d = 0; d < 3; ++d) {
            vector_field[p][d] = extract_result<tag::mean>( accumulators[p][d] );
        }
    }
}

void getVectorField(vector<array<lpt::boost_weighted_accumulator, 3> >& weighted_accumulators, lpt::VectorField& vector_field)
{
    vector_field.resize( weighted_accumulators.size() );

    for (int p = 0; p < weighted_accumulators.size() ; ++p) {
        for (int d = 0; d < 3; ++d) {
            vector_field[p][d] = extract_result<tag::weighted_mean>( weighted_accumulators[p][d] );
        }
    }
}

void getReynoldsStressField(vector< array <lpt::boost_weighted_accumulator, 3 > >& velocity_accumulators,
                            vector< array <lpt::boost_weighted_covariance_accumulator, 3 > >& reynolds_stress_accumulators,
                            lpt::ReynoldsStressField& stress_field)
{
    stress_field.resize( velocity_accumulators.size() );

    for (int p = 0; p < velocity_accumulators.size(); ++p) {
        stress_field[p][0] = extract_result<tag::weighted_variance>( velocity_accumulators[p][0] );  // uu
        stress_field[p][1] = extract_result<tag::weighted_variance>( velocity_accumulators[p][1] );  // vv
        stress_field[p][2] = extract_result<tag::weighted_variance>( velocity_accumulators[p][2] );  // ww
        stress_field[p][3] = extract_result<tag::weighted_covariance<double, tag::covariate1>>( reynolds_stress_accumulators[p][0] );  // uv
        stress_field[p][4] = extract_result<tag::weighted_covariance<double, tag::covariate1>>( reynolds_stress_accumulators[p][1] );  // uw
        stress_field[p][5] = extract_result<tag::weighted_covariance<double, tag::covariate1>>( reynolds_stress_accumulators[p][2] );  // vw
    }
}

array<double, 3>  calcPressureGradiantNavierStokes(array<double,3>& U0, array<double,3>& U1, array<double,3>& U2,
                                                   array<double,3>& A1, lpt::FluidProperties& fluid_props, array<double,3>& dX)
{
    array<double, 3> dP_dX;

    for(int d = 0; d < 3; ++d) {
        double dU_dXX = 0;

        dU_dXX = ( U2[d] -2.0 * U1[d] + U0[d] ) / ( dX[0] * dX[0] ) +
                 ( U2[d] -2.0 * U1[d] + U0[d] ) / ( dX[1] * dX[1] ) +
                 ( U2[d] -2.0 * U1[d] + U0[d] ) / ( dX[2] * dX[2] ) ;

        dP_dX[d] = -1.0 * fluid_props.rho * ( A1[d] ) + fluid_props.mu * dU_dXX;
    }
		return dP_dX;
}

array<double, 3>  calcPressureGradiantBernoulli( array<double,3>& U0, array<double,3>& U2, lpt::FluidProperties& fluid_props, array<double,3>& dX)
{
    array<double, 3> dP_dX;

    dP_dX[0] = ( fluid_props.rho / 2.0 * ( U0[0] * U0[0] - U2[0] * U2[0] ) ) / ( 2.0 * dX[0] );

    dP_dX[1] = ( fluid_props.rho / 2.0 * ( U0[1] * U0[1] - U2[1] * U2[1] ) ) / ( 2.0 * dX[1] );

    dP_dX[2] = ( fluid_props.rho / 2.0 * ( U0[2] * U0[2] - U2[2] * U2[2] ) ) / ( 2.0 * dX[2] );

    return dP_dX;
}

/***** Definition of FluidProperties class *****/

FluidProperties::FluidProperties()
  : rho(1.225), Td(15), Tw(15), mu(1.78E-5), P_ref(101325) {}

FluidProperties::FluidProperties(double n_rho, double n_Td, double n_Tw, double n_mu, double n_P_ref)
  : rho(n_rho), Td(n_Td), Tw(n_Tw), mu(n_mu), P_ref(n_P_ref) {}

/***** Definition of CoordinateArrows class *****/

CoordinateArrows::CoordinateArrows( vtkSmartPointer<vtkRenderWindowInteractor> iren, double scale) : interactor(iren), scale(scale)
{
    arrow = vtkSmartPointer<vtkArrowSource>::New();
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection( arrow->GetOutputPort() );
    actors.resize(3);
    for (int i = 0; i < actors.size(); ++i) {
        actors[i] = vtkSmartPointer<vtkActor>::New();
        actors[i]->SetMapper(mapper);
        actors[i]->SetScale(scale);
        //actors[i]->GetProperty()->SetOpacity(0.995);
        interactor->GetInteractorStyle()->GetDefaultRenderer()->AddActor(actors[i]);
    }

    actors[1]->RotateWXYZ(90, 0, 0, 1);
    actors[2]->RotateWXYZ(-90, 0, 1, 0);
    actors[0]->GetProperty()->SetColor(1, 0, 0);
    actors[1]->GetProperty()->SetColor(0, 1, 0);
    actors[2]->GetProperty()->SetColor(0, 0, 1);

    axes = vtkSmartPointer<vtkAxesActor>::New();

    widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( interactor );
    widget->SetViewport( 0.0, 0.0, 0.25, 0.25 );
    widget->SetEnabled( 1 );
    widget->InteractiveOff();
    widget->KeyPressActivationOff();
}

CoordinateArrows::~CoordinateArrows()
{
}

/***** Definition of CameraVTK class *****/

array<double, 3> CameraVTK::getRandomColor()
{
    array<double, 3> color;
    color[0] = vtkMath::Random(0.5, 1.0);
    color[1] = vtkMath::Random(0.5, 1.0);
    color[2] = vtkMath::Random(0.5, 1.0);
    return color;
}

CameraVTK::CameraVTK(lpt::Camera& camera)
{
    frustrum = vtkSmartPointer<vtkFrustumSource>::New();
    planes = vtkSmartPointer<vtkPlanes>::New();
    double pixel_size = camera.pixel_size[0]; //mm
    double sensor_width = camera.sensor_size[0]; //mm
    double sensor_height = camera.sensor_size[1]; //mm
    double projection_depth = 150; //mm
    double f = camera.f[0] * pixel_size; //mm
    double alpha = 2.0 * atan( sensor_width / 2 / f );
    double beta = 2.0 * atan( sensor_height / 2 / f );
		
    vector<array<double, 3>> Ac(4);
    vector<array<double, 3>> Bc(4);
    vector<array<double, 3>> Aw(4);
    vector<array<double, 3>> Bw(4);

    Ac[0][0] = 0;
    Ac[0][1] = 0;
    Ac[0][2] = 0;

    Ac[1][0] = 0;
    Ac[1][1] = sensor_height;
    Ac[1][2] = 0;

    Ac[2][0] = sensor_width;
    Ac[2][1] = sensor_height;
    Ac[2][2] = 0;

    Ac[3][0] = sensor_width;
    Ac[3][1] = 0;
    Ac[3][2] = 0;
		
    double lw = 2.0 * (f + projection_depth) * tan(alpha / 2.0);
    double lh = 2.0 * (f + projection_depth) * tan(beta / 2.0);

    Bc[0][0] = -1.0 * lw / 2.0 - sensor_width / 2.0;
    Bc[0][1] = -1.0 * lh / 2.0 - sensor_height / 2.0;
    Bc[0][2] = f + projection_depth;

    Bc[1][0] = -1.0 * lw / 2.0 - sensor_width / 2.0;
    Bc[1][1] = lh / 2.0 - sensor_height / 2.0;
    Bc[1][2] = f + projection_depth;

    Bc[2][0] = lw / 2.0 - sensor_width / 2.0;
    Bc[2][1] = lh / 2.0 - sensor_height / 2.0;
    Bc[2][2] = f + projection_depth;

    Bc[3][0] = lw / 2.0 - sensor_width / 2.0;
    Bc[3][1] = -1.0 * lh / 2.0 - sensor_height / 2.0;
    Bc[3][2] = f + projection_depth;

    for( int i = 0; i < 4; ++i) {
        lpt::convertCameraCoordinatesToWorld(camera, Ac[i], Aw[i]);//camToWorld(camera, Ac[i], Aw[i]);
        lpt::convertCameraCoordinatesToWorld(camera, Bc[i], Bw[i]);//camToWorld(camera, Bc[i], Bw[i]);
    }

    vector<array<double, 4>> plainarray(6);

    calcPlaneEq(Aw[0], Aw[1], Bw[0], plainarray[0]);
    calcPlaneEq(Aw[2], Aw[3], Bw[2], plainarray[1]);
    calcPlaneEq(Aw[0], Aw[3], Bw[0], plainarray[2]);
    calcPlaneEq(Aw[1], Aw[2], Bw[1], plainarray[3]);
    calcPlaneEq(Aw[0], Aw[1], Aw[2], plainarray[4]);
    calcPlaneEq(Bw[0], Bw[1], Bw[2], plainarray[5]);

    double plain_eqs[24];
		
    for (int i = 0; i < plainarray.size(); ++i) {
        for (int j = 0; j < plainarray[i].size(); ++j) {
            plain_eqs[i*plainarray[i].size() + j] = plainarray[i][j];
        }
    }

    planes->SetFrustumPlanes(plain_eqs);
    frustrum->SetPlanes(planes);
    frustrum->ShowLinesOff();
    frustrum->Update();
		
    actor = vtkSmartPointer<vtkActor>::New();
    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->AddInputConnection(frustrum->GetOutputPort());
    actor->SetMapper(mapper);
    actor->SetScale(1);
    actor->PickableOff();
    mapper->Update();
		
    double position[3];
    position[0] = -1 * camera.R[0][0]*camera.T[0] - camera.R[1][0]*camera.T[1] - camera.R[2][0]*camera.T[2];
    position[1] = -1 *camera.R[0][1]*camera.T[0] - camera.R[1][1]*camera.T[1] - camera.R[2][1]*camera.T[2];
    position[2] = -1 *camera.R[0][2]*camera.T[0] - camera.R[1][2]*camera.T[1] - camera.R[2][2]*camera.T[2];

    double direction[3];
    direction[0] = camera.R[0][0] * 0.0 + camera.R[1][0] * 0.0 + camera.R[2][0] * 1.0;
    direction[1] = camera.R[0][1] * 0.0 + camera.R[1][1] * 0.0 + camera.R[2][1] * 1.0;
    direction[2] = camera.R[0][2] * 0.0 + camera.R[1][2] * 0.0 + camera.R[2][2] * 1.0;
			
    double mag = sqrt( direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2] );
		
    line = vtkSmartPointer<vtkLineSource>::New();
    line->SetPoint1(position[0], position[1], position[2]);
    line->SetPoint2(position[0] + 1000*direction[0], position[1] + 1000*direction[1], position[2] + 1000* direction[2]);
		
    lineactor = vtkSmartPointer<vtkActor>::New();
    linemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    linemapper->AddInputConnection(line->GetOutputPort());
    lineactor->SetMapper(linemapper);

    array<double, 3> color = CameraVTK::getRandomColor();
    //actor->GetProperty()->SetColor(color.data());
    actor->GetProperty()->BackfaceCullingOn();
    actor->GetProperty()->SetOpacity(0.50);
    actor->Modified();

    stringstream textstream;
    textstream << camera.id; // << ": " << camera.name;
    text = vtkSmartPointer<vtkVectorText>::New();
    text->SetText( textstream.str().c_str() );
    text->Update();
    textmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    textmapper->SetInputConnection( text->GetOutputPort() );
    follower = vtkSmartPointer<vtkFollower>::New();
    follower->SetMapper( textmapper );
    follower->GetProperty()->SetColor( 1, 0, 0 ); // red
    follower->SetPosition(position[0] + 50, position[1] + 50, position[2] - 20 );
    follower->SetScale( 20 * (*actor->GetScale()) );
    follower->PickableOff();
    follower->Modified();
}

void CameraVTK::calcPlaneEq(array<double,3>& A, array<double,3>& B, array<double,3>& C, array<double,4>& eq )
{
    eq[0] = (B[1] - A[1]) * (C[2] - A[2]) - (C[1] - A[1]) * (B[2] - A[2]);
    eq[1] = (B[2] - A[2]) * (C[0] - A[0]) - (C[2] - A[2]) * (B[0] - A[0]);
    eq[2] = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1]);
    eq[3] = -1.0 * (eq[0] * A[0] + eq[1] * A[1] + eq[2] * A[2]);
}

void CameraVTK::camToWorld(lpt::Camera& cam, array<double, 3>& pc, array<double, 3>& pw)
{
    pw[0] = cam.R[0][0] * (pc[0] - cam.T[0]) + cam.R[1][0] * (pc[1] - cam.T[1]) + cam.R[2][0] * (pc[2] - cam.T[2]);
    pw[1] = cam.R[0][1] * (pc[0] - cam.T[0]) + cam.R[1][1] * (pc[1] - cam.T[1]) + cam.R[2][1] * (pc[2] - cam.T[2]);
    pw[2] = cam.R[0][2] * (pc[0] - cam.T[0]) + cam.R[1][2] * (pc[1] - cam.T[1]) + cam.R[2][2] * (pc[2] - cam.T[2]);
}

void CameraVTK::addToRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
    renderer->AddActor(actor);
    //renderer->AddActor(lineactor);
    follower->SetCamera(renderer->GetActiveCamera() );
    follower->Modified();
    renderer->AddActor(follower);
}

void CameraVTK::removeFromRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
    renderer->RemoveActor(actor);
    //renderer->RemoveActor(lineactor);
    renderer->RemoveActor(follower);
}

/***** Definition of TrajectoryPathVTK class *****/

TrajectoryPathVTK::TrajectoryPathVTK(vtkRenderer* renderer)
    : renderer(renderer)
{
    points = vtkSmartPointer<vtkPoints>::New();
    lines = vtkSmartPointer<vtkCellArray>::New();
    polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);

/*
    ribbon_filter = vtkSmartPointer<vtkRibbonFilter>::New();
    ribbon_filter->SetInput(polydata);
    ribbon_filter->SetWidth(1);
    ribbon_filter->SetVaryWidth(1);
*/
		
    scalars = vtkSmartPointer<vtkDoubleArray>::New();
    scalars->SetName("scalars_data");
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(this->max_points);

    polydata->GetPointData()->AddArray(scalars);
    polydata->GetPointData()->SetActiveScalars("scalars_data");

    tube_filter = vtkSmartPointer<vtkTubeFilter>::New();
    tube_filter->SetInput(polydata);
    tube_filter->SetRadius(1.0);
    tube_filter->SetNumberOfSides(4);
    tube_filter->Update();
		
    traj_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    traj_mapper->SetInputConnection( tube_filter->GetOutputPort() );

    traj_mapper->ScalarVisibilityOn();
    traj_mapper->SetScalarModeToUsePointFieldData();
    traj_mapper->SelectColorArray("scalars_data");
    traj_mapper->UseLookupTableScalarRangeOn();
    traj_mapper->SetColorModeToMapScalars();
					
    traj_mapper->SetLookupTable(lookuptable);

    trajactor = vtkSmartPointer<vtkActor>::New();
    trajactor->SetMapper(traj_mapper);
    renderer->AddActor(trajactor);
}
	
TrajectoryPathVTK::~TrajectoryPathVTK()
{
    removeFromRenderer();
}

void TrajectoryPathVTK::addNextPoint(lpt::ParticleVectors& current_particle)
{
    position_queue.push_back( current_particle[0] );
    scalars_queue.push_back(
        sqrt(
            current_particle[this->scalar_index][0] * current_particle[this->scalar_index][0]
          + current_particle[this->scalar_index][1] * current_particle[this->scalar_index][1]
          + current_particle[this->scalar_index][2] * current_particle[this->scalar_index][2])
    );
		
    if (position_queue.size() < this->max_points) {
        point_ids.push_back( position_queue.size() - 1 );
        points->SetNumberOfPoints( position_queue.size() );
        scalars->SetNumberOfTuples( position_queue.size() );

        if (point_ids.size() > 1) {
            lines->Reset();
            lines->InsertNextCell(point_ids.size(), point_ids.data() );
        }
    } else {
        position_queue.pop_front();
        scalars_queue.pop_front();
    }
		
    auto position_iter = position_queue.begin();
    auto scalars_iter = scalars_queue.begin();

    for ( vtkIdType point_id = point_ids.front(); point_id <= point_ids.back(); ++position_iter, ++scalars_iter, ++point_id) {
        points->SetPoint(point_id, position_iter->data() );
        scalars->SetTuple1(point_id, *scalars_iter);
    }
			
    lines->Modified();
    points->Modified();
    scalars->Modified();
    traj_mapper->Update();
}

void TrajectoryPathVTK::removeFromRenderer()
{
    if (renderer)
        renderer->RemoveActor(this->trajactor);
}

//Pressure Calculation
/***** Definition of PressureFieldSolver *****/

PressureFieldSolver::PressureFieldSolver()
  : dt(0), dx(0), dy(0), dz(0), initial_pressure(0)
{
    for (size_t i=0; i<grid_cell_counts.size(); i++)
        grid_cell_counts[i] = 0;

    fluid_props = lpt::FluidProperties::create(); // default fluid: air at 15 degree C
}

PressureFieldSolver::~PressureFieldSolver()
{
}

void PressureFieldSolver::addControls()
{
}

void PressureFieldSolver::setFluidProperties(FluidProperties::Ptr props)
{
    fluid_props = props;
}

void PressureFieldSolver::setGridProperties(array<int, 3> &grid_cell_counts, array<double, 3> &cell_dimensions)
{
    for (size_t i=0; i<grid_cell_counts.size(); i++) {
        this->grid_cell_counts[i] = grid_cell_counts[i];
    }

    dx = cell_dimensions[0] / 1000.0; // mm
    dy = cell_dimensions[1] / 1000.0;// mm
    dz = cell_dimensions[2] / 1000.0; // mm

    this->resetPressureField();
}

/***** Definition of PressurePoissonSolver2Step class *****/

PressurePoissonSolver2Step::PressurePoissonSolver2Step()
{
}

PressurePoissonSolver2Step::~PressurePoissonSolver2Step()
{
}

void PressurePoissonSolver2Step::solve(lpt::VectorField& velocities, lpt::VectorField& accelerations, lpt::ReynoldsStressField& stresses)
{
    velocity_field = &velocities;
    acceleration_field = &accelerations;
    stress_field = &stresses;
    this->dt = 1.0 / static_cast<double>(shared_objects->frame_rate);  // seconds
    cout << "Pressure solver\n"<< "dx, dy, dz, dt = " << dx << " " << dy << " " << dz << " " << dt << " rho = " << fluid_props->rho << " sor_weight = " << params.sor_weight << endl;
    vector<double> pressure_old( pressure.begin(), pressure.end() );

    cout << "Solving Pressure Field with SOR:" << endl;
    pressure_gradients.resize( pressure.size() );
		
    for ( int i = 0; i < grid_cell_counts[0]; ++i ) {
        for ( int j = 0; j < grid_cell_counts[1]; ++j ) {
            for ( int k = 0; k < grid_cell_counts[2]; ++k ) {
                calcCellPressureGradientsSteadyRANS(i,j,k);
            }
        }
    }

    for (int r = 0; r < params.max_iterations; ++r)  {
        for ( int i = 0; i < grid_cell_counts[0]; ++i ) {
            for ( int j = 0; j < grid_cell_counts[1]; ++j ) {
                for ( int k = 0; k < grid_cell_counts[2]; ++k ) {
                    calcCellPressure(i,j,k);
                }
            }
        }

        // set reference pressure
        int ref_id = getGridIndex(2,2,2);
        for (int p = 0; p < pressure.size(); ++p)
            pressure[p] -= pressure[ref_id];

        if(r % 200 == 0) {
            double norm2 = 0;
            for (int p = 0; p < pressure.size(); ++p) {
                norm2 += (pressure[p] - pressure_old[p]) * (pressure[p] - pressure_old[p]);
            }
            norm2 = sqrt(norm2);

            cout << "\t" << r <<": L2 residual norm = " <<  norm2 << endl;
            if (norm2 <= params.convergence_tolerance)
                break;
        }
        std::copy( pressure.begin(), pressure.end(), pressure_old.begin() );
    }
}

void PressurePoissonSolver2Step::addControls()
{
}

void PressurePoissonSolver2Step::calcCellPressure(int i, int j, int k)
{
    int p;
    array<int,2> x, y, z;  // define the points surrounding point ijk in the finite differnece stencile

    double delta_x = 2.0 * dx ;
    double delta_y = 2.0 * dy ;
    double delta_z = 2.0 * dz ;

    array<double,2> ax, ay, az;

    //Source terms for Nuemman boundary conditions
    double S_bc_x = 0;
    double S_bc_y = 0;
    double S_bc_z = 0;

    p = getGridIndex( i  , j  , k   );     //index p = currnet point, (i,j,k)
				
    ax[0] = dy*dy * dz*dz;
    ay[0] = dx*dx * dz*dz;
    az[0] = dx*dx * dy*dy;
		
    ax[1] = ax[0];
    ay[1] = ay[0];
    az[1] = az[0];

    // X-Boundary Check
    if (i > 0 && i < grid_cell_counts[0] - 1) {
        x[0] = getGridIndex( i + 1 , j  , k    );
        x[1] = getGridIndex( i - 1 , j  , k    );
			
        S_bc_x = 0;
    }
    else {
        delta_x = dx;
        if (i == 0) {
            x[0] = getGridIndex( i + 1 , j  , k    );
            x[1] = getGridIndex( i     , j  , k    );
				
            ax[1] = 0;
            S_bc_x = -1.0 * ax[0] * dx * pressure_gradients[ p ][0];
        }
        else if( i == grid_cell_counts[0] - 1) {
            x[0] = getGridIndex( i     , j  , k    );
            x[1] = getGridIndex( i - 1 , j  , k    );
					
            ax[0] = 0;
            S_bc_x = ax[1] * dx * pressure_gradients[ p ][0];
        }
    }
    // Y-Boundary Check
    if (j > 0 && j < grid_cell_counts[1] - 1) {
        y[0] = getGridIndex( i     , j + 1 , k    );
        y[1] = getGridIndex( i     , j - 1 , k    );
		
        S_bc_y = 0;
    }
    else  {
        delta_y = dy;
        if (j == 0) {
            y[0] = getGridIndex( i    , j + 1, k    );
            y[1] = getGridIndex( i    , j    , k    );
				
            ay[1] = 0;
            S_bc_y = -1.0 * ay[0] * dy * pressure_gradients[ p ][1];
        }
        else if( j == grid_cell_counts[1] - 1) {
            y[0] = getGridIndex( i    , j    , k    );
            y[1] = getGridIndex( i    , j - 1, k    );

            ay[0] = 0;
            S_bc_y = ay[1] * dy * pressure_gradients[ p ][1];
        }
    }
    // Z-Boundary Check
    if (k > 0 && k < grid_cell_counts[2] - 1) {
        z[0] = getGridIndex( i     , j     , k + 1 );
        z[1] = getGridIndex( i     , j     , k - 1 );

        S_bc_z = 0;
    }
    else {
        delta_z = dz;
        if (k == 0) {
            z[0] = getGridIndex( i    , j    , k + 1 );
            z[1] = getGridIndex( i    , j    , k     );

            az[1] = 0;
            S_bc_z = -1.0 * az[0] * dz * pressure_gradients[ p ][2];
        }
        else if( k == grid_cell_counts[2] - 1) {
            z[0] = getGridIndex( i    , j    , k     );
            z[1] = getGridIndex( i    , j    , k - 1 );

            az[0] = 0;
            S_bc_z = az[1] * dz * pressure_gradients[ p ][2];
        }
    }
						
    double S =
            ( pressure_gradients[x[0]][0] - pressure_gradients[x[1]][0] ) / delta_x
        +	( pressure_gradients[y[0]][1] - pressure_gradients[y[1]][1] ) / delta_y
        +	( pressure_gradients[z[0]][2] - pressure_gradients[z[1]][2] ) / delta_z
        ;
				
    double as = dx*dx * dy*dy * dz*dz;
    double ap = ax[0] + ax[1] + ay[0] + ay[1] + az[0] + az[1];

    double pressure_temp;
		
    pressure_temp =
            (
                ax[0] * pressure[ x[0] ] + ax[1] * pressure[ x[1] ] +

				ay[0] * pressure[ y[0] ] + ay[1] * pressure[ y[1] ] +

				az[0] * pressure[ z[0] ] + az[1] * pressure[ z[1] ] -							

				as * S + 
				
				S_bc_x + 
				
				S_bc_y + 
				
				S_bc_z

            ) / ap;
		
    pressure[p] = params.sor_weight * pressure_temp + (1 - params.sor_weight) * pressure[p];
}

void  PressurePoissonSolver2Step::calcCellPressureGradientsSteadyRANS(int i, int j, int k)
{
    int p;
    array<int,2> x, y, z;  // define the points surrounding point ijk in the finite differnece stencile
    bool x_boundary = false, y_boundary = false, z_boundary = false;
		
    double delta_x = 2.0 * dx ;
    double delta_y = 2.0 * dy ;
    double delta_z = 2.0 * dz ;

    p = getGridIndex( i  , j  , k   );     //index p = current point, (i,j,k)

    // X-Boundary Check
    if (i > 0 && i < grid_cell_counts[0] - 1) {
        x[0] = getGridIndex( i + 1 , j  , k    );
        x[1] = getGridIndex( i - 1 , j  , k    );
    }
    else {
        x_boundary = true;
        delta_x = dx;
        if (i == 0) {
            x[0] = getGridIndex( i + 1 , j  , k    );
            x[1] = getGridIndex( i     , j  , k    );
        }
        else if( i == grid_cell_counts[0] - 1) {
            x[0] = getGridIndex( i     , j  , k    );
            x[1] = getGridIndex( i - 1 , j  , k    );
        }
    }
    // Y-Boundary Check
    if (j > 0 && j < grid_cell_counts[1] - 1) {
        y[0] = getGridIndex( i     , j + 1 , k    );
        y[1] = getGridIndex( i     , j - 1 , k    );
    }
    else  {
        y_boundary = true;
        delta_y = dy;
        if (j == 0) {
            y[0] = getGridIndex( i    , j + 1, k    );
            y[1] = getGridIndex( i    , j    , k    );
        }
        else if( j == grid_cell_counts[1] - 1) {
            y[0] = getGridIndex( i    , j    , k    );
            y[1] = getGridIndex( i    , j - 1, k    );
        }
    }
    // Z-Boundary Check
    if (k > 0 && k < grid_cell_counts[2] - 1) {
        z[0] = getGridIndex( i     , j     , k + 1 );
        z[1] = getGridIndex( i     , j     , k - 1 );
    }
    else {
        z_boundary = true;
        delta_z = dz;
        if (k == 0) {
            z[0] = getGridIndex( i    , j    , k + 1 );
            z[1] = getGridIndex( i    , j    , k     );
        }
        else if( k == grid_cell_counts[2] - 1) {
            z[0] = getGridIndex( i    , j    , k     );
            z[1] = getGridIndex( i    , j    , k - 1 );
        }
    }
	
    array<array<double,3>,2> Ui, Uj, Uk;
    array<double, 3> Up;
    array<array<double, 6>,2> ti, tj, tk;

    auto& velocities = *velocity_field;
    auto& stresses = *stress_field;

    for (int d = 0; d < 3; ++d) {
        Ui[0][d] = velocities[ x[0] ][d];
        Ui[1][d] = velocities[ x[1] ][d];

        Uj[0][d] = velocities[ y[0] ][d];
        Uj[1][d] = velocities[ y[1] ][d];

        Uk[0][d] = velocities[ z[0] ][d];
        Uk[1][d] = velocities[ z[1] ][d];
			
        Up[d] = velocities[p][d];
    }
 	
    for (int s = 0; s < 6; ++s) {
        ti[0][s] = stresses[x[0] ][s];
        tj[0][s] = stresses[y[0] ][s];
        tk[0][s] = stresses[z[0] ][s];

        ti[1][s] = stresses[x[1] ][s];
        tj[1][s] = stresses[y[1] ][s];
        tk[1][s] = stresses[z[1] ][s];
    }
		
    // velocity array index key: [0,1] = [X1,X2],  [2,3] = [Y1,Y2],  [4,5] = [Z1, Z2], [6] = center
			
    double U_grad_u;
    double U_grad_v;
    double U_grad_w;

    double du_dXX = 0;
    double dv_dXX = 0;
    double dw_dXX = 0;

    double dt1_dX = 0;
    double dt2_dX = 0;
    double dt3_dX = 0;
		
    U_grad_u =
        Up[0] * ( Ui[0][0] - Ui[1][0] ) / delta_x +
			
        Up[1] * ( Uj[0][0] - Uj[1][0] ) / delta_y +

        Up[2] * ( Uk[0][0] - Uk[1][0] ) / delta_z;

    U_grad_v =
        Up[0] * ( Ui[0][1] - Ui[1][1] ) / delta_x +
			
        Up[1] * ( Uj[0][1] - Uj[1][1] ) / delta_y +
			
        Up[2] * ( Uk[0][1] - Uk[1][1] ) / delta_z;

    U_grad_w =
        Up[0] * ( Ui[0][2] - Ui[1][2] ) / delta_x +
			
        Up[1] * ( Uj[0][2] - Uj[1][2] ) / delta_y +
			
        Up[2] * ( Uk[0][2] - Uk[1][2] ) / delta_z;


    if( x_boundary ) {
        du_dXX += ( Ui[0][0] -2.0 * Ui[1][0] + Up[0] ) / ( dx * dx );
        dv_dXX += ( Ui[0][1] -2.0 * Ui[1][1] + Up[1] ) / ( dx * dx );
        dw_dXX += ( Ui[0][2] -2.0 * Ui[1][2] + Up[2] ) / ( dx * dx );
    } else {
        du_dXX += ( Ui[0][0] -2.0 * Up[0] + Ui[1][0] ) / ( dx * dx );
        dv_dXX += ( Ui[0][1] -2.0 * Up[1] + Ui[1][1] ) / ( dx * dx );
        dw_dXX += ( Ui[0][2] -2.0 * Up[2] + Ui[1][2] ) / ( dx * dx );
    }

    if( y_boundary ) {
        du_dXX += ( Uj[0][0] -2.0 * Uj[1][0] + Up[0] ) / ( dy * dy );
        dv_dXX += ( Uj[0][1] -2.0 * Uj[1][1] + Up[1] ) / ( dy * dy );
        dw_dXX += ( Uj[0][2] -2.0 * Uj[1][2] + Up[2] ) / ( dy * dy );

    } else {
        du_dXX += ( Uj[0][0] -2.0 * Up[0] + Uj[1][0] ) / ( dy * dy );
        dv_dXX += ( Uj[0][1] -2.0 * Up[1] + Uj[1][1] ) / ( dy * dy );
        dw_dXX += ( Uj[0][2] -2.0 * Up[2] + Uj[1][2] ) / ( dy * dy );

    }

    if( z_boundary ) {
        du_dXX += ( Uk[0][0] -2.0 * Uk[1][0] + Up[0] ) / ( dz * dz );
        dv_dXX += ( Uk[0][1] -2.0 * Uk[1][1] + Up[1] ) / ( dz * dz );
        dw_dXX += ( Uk[0][2] -2.0 * Uk[1][2] + Up[2] ) / ( dz * dz );
    } else {
        du_dXX += ( Uk[0][0] -2.0 * Up[0] + Uk[1][0] ) / ( dz * dz );
        dv_dXX += ( Uk[0][1] -2.0 * Up[1] + Uk[1][1] ) / ( dz * dz );
        dw_dXX += ( Uk[0][2] -2.0 * Up[2] + Uk[1][2] ) / ( dz * dz );
    }


    //reynolds stress indicies: 0 = u2, 1 = v2, 2 = w2, 3 = uv, 4 = uw, 5 = vw
    dt1_dX =
            (
            ( ti[0][0] + ti[0][3] + ti[0][4] ) -
            ( ti[1][0] + ti[1][3] + ti[1][4] )
            ) / (delta_x)

            +
            (
            ( tj[0][0] + tj[0][3] + tj[0][4] ) -
            ( tj[1][0] + tj[1][3] + tj[1][4] )
            ) / (delta_y)

            +
            (
            ( tk[0][0] + tk[0][3] + tk[0][4] ) -
            ( tk[1][0] + tk[1][3] + tk[1][4] )
            ) / (delta_z)
            ;

    dt2_dX =
            (
            ( ti[0][1] + ti[0][3] + ti[0][5] ) -
            ( ti[1][1] + ti[1][3] + ti[1][5] )
            ) / (delta_x)

            +
            (
            ( tj[0][1] + tj[0][3] + tj[0][5] ) -
            ( tj[1][1] + tj[1][3] + tj[1][5] )
            ) / (delta_y)

            +
            (
            ( tk[0][1] + tk[0][3] + tk[0][5] ) -
            ( tk[1][1] + tk[1][3] + tk[1][5] )
            ) / (delta_z)
            ;
		
    dt3_dX =
            (
            ( ti[0][2] + ti[0][4] + ti[0][5] ) -
            ( ti[1][2] + ti[1][4] + ti[1][5] )
            ) / (delta_x)

            +
            (
            ( tj[0][2] + tj[0][4] + tj[0][5] ) -
            ( tj[1][2] + tj[1][4] + tj[1][5] )
            ) / (delta_y)

            +
            (
            ( tk[0][2] + tk[0][4] + tk[0][5] ) -
            ( tk[1][2] + tk[1][4] + tk[1][5] )
            ) / (delta_z)
            ;

    pressure_gradients[p][0] = -1.0 * fluid_props->rho *  U_grad_u + fluid_props->mu * du_dXX -  fluid_props->rho * dt1_dX;
    pressure_gradients[p][1] = -1.0 * fluid_props->rho *  U_grad_v + fluid_props->mu * dv_dXX -  fluid_props->rho * dt2_dX;
    pressure_gradients[p][2] = -1.0 * fluid_props->rho *  U_grad_w + fluid_props->mu * dw_dXX -  fluid_props->rho * dt3_dX;
}

/***** Definition of LagrangianPressureFieldSolver class *****/

LagrangianPressureFieldSolver::LagrangianPressureFieldSolver()
{
}

LagrangianPressureFieldSolver::~LagrangianPressureFieldSolver()
{
}

void LagrangianPressureFieldSolver::solve(lpt::VectorField& velocities, lpt::VectorField& accelerations, lpt::ReynoldsStressField& stresses)
{
}

void LagrangianPressureFieldSolver::solve(lpt::VectorField& pressure_grads)
{
    this->dt = 1.0 / static_cast<double>(shared_objects->frame_rate);  // seconds
    cout << "Pressure solver\n"<< "dx, dy, dz, dt = " << dx << " " << dy << " " << dz << " " << dt << " rho = " << fluid_props->rho << " sor_weight = " << params.sor_weight << endl;
    vector<double> pressure_old( pressure.begin(), pressure.end() );

    cout << "Solving Pressure Field with SOR:" << endl;
    pressure_gradients.resize( pressure_grads.size() );
    std::copy( pressure_grads.begin(), pressure_grads.end(), pressure_gradients.begin() );
				
    for (int r = 0; r < params.max_iterations; ++r)  {
        for ( int i = 0; i < grid_cell_counts[0]; ++i ) {
            for ( int j = 0; j < grid_cell_counts[1]; ++j ) {
                for ( int k = 0; k < grid_cell_counts[2]; ++k ) {
                    calcCellPressure(i,j,k);
                }
            }
        }

        // set reference pressure
        int ref_id = getGridIndex(2,2,2);
        for (int p = 0; p < pressure.size(); ++p)
            pressure[p] -= pressure[ref_id];

        if(r % 200 == 0) {
            double norm2 = 0;
            for (int p = 0; p < pressure.size(); ++p)
                norm2 += (pressure[p] - pressure_old[p]) * (pressure[p] - pressure_old[p]);
            norm2 = sqrt(norm2);

            cout << "\t" << r <<": L2 residual norm = " <<  norm2 << endl;
            if (norm2 <= params.convergence_tolerance)
                break;
        }

        std::copy( pressure.begin(), pressure.end(), pressure_old.begin() );
    }
}

void LagrangianPressureFieldSolver::calcCellPressure(int i, int j, int k)
{
    int p;
    array<int,2> x, y, z;  // define the points surrounding point ijk in the finite differnece stencile

    double delta_x = 2.0 * dx ;
    double delta_y = 2.0 * dy ;
    double delta_z = 2.0 * dz ;

    array<double,2> ax, ay, az;

    //Source terms for Nuemman boundary conditions
    double S_bc_x = 0;
    double S_bc_y = 0;
    double S_bc_z = 0;

    p = getGridIndex( i  , j  , k   );     //index p = currnet point, (i,j,k)
				
    ax[0] = dy*dy * dz*dz;
    ay[0] = dx*dx * dz*dz;
    az[0] = dx*dx * dy*dy;
		
    ax[1] = ax[0];
    ay[1] = ay[0];
    az[1] = az[0];

    // X-Boundary Check
    if (i > 0 && i < grid_cell_counts[0] - 1) {
        x[0] = getGridIndex( i + 1 , j  , k    );
        x[1] = getGridIndex( i - 1 , j  , k    );
			
        S_bc_x = 0;
    }
    else {
        delta_x = dx;
        if (i == 0) {
            x[0] = getGridIndex( i + 1 , j  , k    );
            x[1] = getGridIndex( i     , j  , k    );
				
            ax[1] = 0;
            S_bc_x = -1.0 * ax[0] * dx * pressure_gradients[ p ][0];
        }
        else if( i == grid_cell_counts[0] - 1) {
            x[0] = getGridIndex( i     , j  , k    );
            x[1] = getGridIndex( i - 1 , j  , k    );
					
            ax[0] = 0;
            S_bc_x = ax[1] * dx * pressure_gradients[ p ][0];
        }
    }
    // Y-Boundary Check
    if (j > 0 && j < grid_cell_counts[1] - 1) {
        y[0] = getGridIndex( i     , j + 1 , k    );
        y[1] = getGridIndex( i     , j - 1 , k    );
		
        S_bc_y = 0;
    }
    else  {
        delta_y = dy;
        if (j == 0) {
            y[0] = getGridIndex( i    , j + 1, k    );
            y[1] = getGridIndex( i    , j    , k    );
				
            ay[1] = 0;
            S_bc_y = -1.0 * ay[0] * dy * pressure_gradients[ p ][1];
        }
        else if( j == grid_cell_counts[1] - 1) {
            y[0] = getGridIndex( i    , j    , k    );
            y[1] = getGridIndex( i    , j - 1, k    );
				
            ay[0] = 0;
            S_bc_y = ay[1] * dy * pressure_gradients[ p ][1];
        }
    }
    // Z-Boundary Check
    if (k > 0 && k < grid_cell_counts[2] - 1) {
        z[0] = getGridIndex( i     , j     , k + 1 );
        z[1] = getGridIndex( i     , j     , k - 1 );

        S_bc_z = 0;
    }
    else {
        delta_z = dz;
        if (k == 0) {
            z[0] = getGridIndex( i    , j    , k + 1 );
            z[1] = getGridIndex( i    , j    , k     );

            az[1] = 0;
            S_bc_z = -1.0 * az[0] * dz * pressure_gradients[ p ][2];
        }
        else if( k == grid_cell_counts[2] - 1) {
            z[0] = getGridIndex( i    , j    , k     );
            z[1] = getGridIndex( i    , j    , k - 1 );

            az[0] = 0;
            S_bc_z = az[1] * dz * pressure_gradients[ p ][2];
        }
    }
						
    double S =
            ( pressure_gradients[x[0]][0] - pressure_gradients[x[1]][0] ) / delta_x
        +	( pressure_gradients[y[0]][1] - pressure_gradients[y[1]][1] ) / delta_y
        +	( pressure_gradients[z[0]][2] - pressure_gradients[z[1]][2] ) / delta_z
        ;
				
    double as = dx*dx * dy*dy * dz*dz;
    double ap = ax[0] + ax[1] + ay[0] + ay[1] + az[0] + az[1];

    double pressure_temp;
		
    pressure_temp =
            (
                ax[0] * pressure[ x[0] ] + ax[1] * pressure[ x[1] ] +

                ay[0] * pressure[ y[0] ] + ay[1] * pressure[ y[1] ] +

                az[0] * pressure[ z[0] ] + az[1] * pressure[ z[1] ] -

                as * S +
				
                S_bc_x +

                S_bc_y +
				
                S_bc_z

            ) / ap;
		
    pressure[p] = params.sor_weight * pressure_temp + (1 - params.sor_weight) * pressure[p];
}

void LagrangianPressureFieldSolver::addControls()
{
}
	
/***** Definition of FiniteVolumeGrid class *****/

FiniteVolumeGrid::FiniteVolumeGrid (vtkSmartPointer < vtkRenderer > renderer) : renderer(renderer), vector_mode(lpt::VELOCITY)
{
    fluid_props = lpt::FluidProperties::create();
    pressure_solver = lpt::PressurePoissonSolver2Step::create();
    lagrangian_pressure_solver = lpt::LagrangianPressureFieldSolver::create();
}

FiniteVolumeGrid::~FiniteVolumeGrid()
{
}

void FiniteVolumeGrid::initialize()
{
    grid = vtkSmartPointer<vtkImageData>::New();
    grid->SetExtent(0, params.grid_cell_counts[0],
                    0, params.grid_cell_counts[1],
                    0, params.grid_cell_counts[2]);

    grid->SetSpacing( params.cell_dimensions.data() );
    grid->SetOrigin( params.grid_origin.data() );
    grid->SetNumberOfScalarComponents(1);
    grid->SetScalarTypeToInt();
		
    pressure_solver->setGridProperties(params.grid_cell_counts, params.cell_dimensions);
    pressure_solver->setSharedObjects(this->shared_objects);
    pressure_solver->setFluidProperties(fluid_props);
		
    lagrangian_pressure_solver->setGridProperties(params.grid_cell_counts, params.cell_dimensions);
    lagrangian_pressure_solver->setSharedObjects(this->shared_objects);
    lagrangian_pressure_solver->setFluidProperties(fluid_props);

    int id = 0;
    //loop over all points in the grid;  there is one more point in each dimension than number of cells
    for (int i = 0; i < params.grid_cell_counts[0] + 1; i++)
        for (int j = 0; j < params.grid_cell_counts[1] + 1; j++)
            for (int k = 0; k < params.grid_cell_counts[2] + 1; k++, ++id)
                grid->SetScalarComponentFromDouble(i,j,k,0, id);

    cellcentersfilter =	vtkSmartPointer<vtkCellCenters>::New();
    cellcentersfilter->SetInputConnection(grid->GetProducerPort());
    cellcentersfilter->VertexCellsOn();
    cellcentersfilter->Update();

    outline = vtkSmartPointer<vtkOutlineFilter>::New();
    outline->SetInputConnection(grid->GetProducerPort());
    magnitudes = vtkSmartPointer<vtkDoubleArray>::New();

    arrow_source_ids = vtkSmartPointer<vtkIntArray>::New();
    arrow_source_ids->SetName("arrow");
    arrow_source_ids->SetNumberOfComponents(1);
    arrow_source_ids->SetNumberOfTuples( grid->GetNumberOfCells() );

    sphere_source_ids = vtkSmartPointer<vtkIntArray>::New();
    sphere_source_ids->SetName("sphere");
    sphere_source_ids->SetNumberOfComponents(1);
    sphere_source_ids->SetNumberOfTuples( grid->GetNumberOfCells() );

    magnitudes->SetName("magnitudes");
    magnitudes->SetNumberOfComponents(1);
    magnitudes->SetNumberOfTuples( grid->GetNumberOfCells() );

    vectors = vtkSmartPointer<vtkDoubleArray>::New();
    vectors->SetName("vectors");
    vectors->SetNumberOfComponents(3);
    vectors->SetNumberOfTuples( grid->GetNumberOfCells() );

    velocity_accumulators.resize( grid->GetNumberOfCells() );
    acceleration_accumulators.resize( grid->GetNumberOfCells() );
    reynolds_stress_accumulators.resize( grid->GetNumberOfCells() );
    pressure_grad_accumulators.resize( grid->GetNumberOfCells() );

    mass_residuals.resize( grid->GetNumberOfCells(), 0 );

    for(vtkIdType i = 0; i < cellcentersfilter->GetOutput()->GetNumberOfCells(); i++) {
        magnitudes->SetTuple1(i, 0 );
        vectors->SetTuple3(i, 0, 0, 0);
        arrow_source_ids->SetTuple1(i,1);
        sphere_source_ids->SetTuple1(i,1);
    }

    grid->GetCellData()->AddArray(magnitudes);
    grid->GetCellData()->AddArray(vectors);
    grid->GetCellData()->AddArray(arrow_source_ids);
    grid->GetCellData()->AddArray(sphere_source_ids);

    grid->GetCellData()->SetActiveScalars("magnitudes");
    grid->GetCellData()->SetActiveVectors("vectors");
    grid->SetActiveAttribute(grid->GetInformation(), vtkDataObject::FIELD_ASSOCIATION_CELLS, "magnitudes", vtkDataSetAttributes::SCALARS);
    grid->SetActiveAttribute(grid->GetInformation(), vtkDataObject::FIELD_ASSOCIATION_CELLS, "vectors", vtkDataSetAttributes::VECTORS);
    grid->Modified();

    std::cout << "There are "
        << grid->GetNumberOfPoints() << " points."
        << std::endl;
    std::cout << "There are "
        << grid->GetNumberOfCells() << " cells."
        << std::endl;

    cout << "Number of cells = " << params.grid_cell_counts[0] * params.grid_cell_counts[1] * params.grid_cell_counts[2] << endl;
    arrow_source = vtkSmartPointer<vtkArrowSource>::New();
    sphere_source = vtkSmartPointer<vtkSphereSource>::New();
    sphere_source->SetRadius(0.5);

    glyph3Dmapper = vtkSmartPointer<vtkGlyph3DMapper>::New();
    glyph3Dmapper->SetSourceConnection( 0, arrow_source->GetOutputPort() );
    glyph3Dmapper->SetSourceConnection( 1, sphere_source->GetOutputPort() );
    glyph3Dmapper->SetSourceIndexing(true);
    glyph3Dmapper->SetSourceIndexArray("arrow");
    glyph3Dmapper->SetInputConnection( cellcentersfilter->GetOutputPort() );
    //glyph3Dmapper->SetScalarModeToUseCellData();
    glyph3Dmapper->SetOrientationModeToDirection();
    glyph3Dmapper->SetScaleFactor(params.fixed_scale_factor);
    glyph3Dmapper->UseLookupTableScalarRangeOn();
    glyph3Dmapper->SetColorModeToMapScalars();
    glyph3Dmapper->SetScalarVisibility(true);

    glyph3Dmapper->SetScaleArray("magnitudes");
    //glyph3Dmapper->SetScaleModeToScaleByMagnitude();
    glyph3Dmapper->SetScaleModeToNoDataScaling();
	
    lookuptable = vtkSmartPointer<vtkLookupTable>::New();
    lookuptable->SetNumberOfColors(256);
    lookuptable->SetHueRange(0.667, 0.0);
    lookuptable->SetRange(params.scalar_range);
    lookuptable->Build();

    glyph3Dmapper->SetLookupTable( lookuptable );
    glyph3Dmapper->SelectColorArray("magnitudes");
    glyph3Dmapper->Update();
    outlineactor = vtkSmartPointer<vtkActor>::New();

    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(outline->GetOutput());
    outlineactor->SetMapper(mapper);

    glyphactor = vtkSmartPointer<vtkActor>::New();
    glyphactor->SetMapper(glyph3Dmapper);

    plane_widget = vtkSmartPointer<vtkImplicitPlaneWidget>::New();
    plane_widget->SetInteractor( renderer->GetRenderWindow()->GetInteractor() );
    plane_widget->SetPlaceFactor(1);
    plane_widget->SetOutlineTranslation(false);
    //plane_widget->SetOriginTranslation(false);
    plane_widget->PlaceWidget( grid->GetBounds() );
    plane_widget->UpdatePlacement();
    plane_widget->SetDrawPlane(false);
    plane_widget->GetPolyDataAlgorithm()->SetInput(grid);
    plane_widget->SetKeyPressActivation(false);
    plane_widget->SetScaleEnabled(false);
    plane_widget->On();

    renderer->GetActors()->GetLastActor()->GetMapper()->SetUseLookupTableScalarRange(1);
    renderer->GetActors()->GetLastActor()->GetMapper()->SetScalarModeToUseCellData();
    planemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    planemapper->SetInputConnection( plane_widget->GetPolyDataAlgorithm()->GetOutputPort() );
    planemapper->SetUseLookupTableScalarRange(true);
    planemapper->SetScalarModeToUseCellData();
    planemapper->SetColorModeToMapScalars();
    //planemapper->InterpolateScalarsBeforeMappingOff();
    planemapper->SetLookupTable(lookuptable);
    planemapper->Update();
    planeactor = vtkSmartPointer<vtkActor>::New();
    planeactor->SetMapper( planemapper);
    plane_widget->Off();

    //vtkSmartPointer<vtkCellPicker> picker =	vtkSmartPointer<vtkCellPicker>::New();
    //picker->SetTolerance(0.005);

    vtkSmartPointer<vtkProperty> ipwProp = vtkSmartPointer<vtkProperty>::New();
    planeWidgetX =	vtkSmartPointer<vtkImagePlaneWidget>::New();
    planeWidgetX->SetInteractor( renderer->GetRenderWindow()->GetInteractor() );
    planeWidgetX->SetKeyPressActivationValue('x');
    //planeWidgetX->SetPicker(picker);
    planeWidgetX->RestrictPlaneToVolumeOn();
    planeWidgetX->GetPlaneProperty()->SetColor(1,0,0);
    planeWidgetX->SetTexturePlaneProperty(ipwProp);
    planeWidgetX->TextureInterpolateOff();
    planeWidgetX->SetResliceInterpolateToNearestNeighbour();
    planeWidgetX->SetInput(grid);
    planeWidgetX->SetPlaneOrientationToXAxes();
    planeWidgetX->SetSliceIndex(1);
    planeWidgetX->DisplayTextOn();
    planeWidgetX->On();
    //planeWidgetX->InteractionOff();
    planeWidgetX->InteractionOn();
    planeWidgetX->Off();
    //planeWidgetX->GetReslice()->Update();
		
    scalarbar = vtkSmartPointer<vtkScalarBarWidget>::New();
    scalarbar->SetInteractor( renderer->GetRenderWindow()->GetInteractor() );
    scalarbar->GetScalarBarActor()->SetTitle("Velocity Magnitude (m/s)");
    scalarbar->GetScalarBarActor()->SetLookupTable( glyph3Dmapper->GetLookupTable() );
    scalarbar->SetRepositionable(false);
    scalarbar->SelectableOff();
    scalarbar->SetManagesCursor(false);
    scalarbar->SetProcessEvents(false);

    for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
            for (int k = -1; k <=1; ++k) {
                array<int,3> ijk = {{i,j,k}};
                cell_neighborhood.push_back(ijk);
            }

    cout << "cell neighborhood size = " << cell_neighborhood.size() << endl;
	
    //streamlines = lpt::StreamLines::create( cellcentersfilter->GetOutputPort() );
    /*double point1[3] = {
            params.grid_origin[0],
            params.grid_origin[1] + params.cell_dimensions[1] / 2,
            params.grid_origin[2] + params.cell_dimensions[2] / 2,
            };
    double point2[3] = {
            params.grid_origin[0] + params.grid_dimensions[0] / 2,
            params.grid_origin[1] + params.grid_dimensions[1] / 2,
            params.grid_origin[2] + params.grid_dimensions[2] / 2,
            };

    streamlines->setSeedPlane(params.grid_origin.data(), point1, point2);*/
    //histogram = lpt::HistogramVTK::create();
    //int color[3] = {255,255,255};
    //histogram->setAxisColor(color);
    //histogram->startChartWindow();
    //histogram->addToRenderer(&renderer);
}

void FiniteVolumeGrid::setPressureFieldSolver(lpt::PressureFieldSolver::Ptr solver) {
    pressure_solver = solver;
    pressure_solver->setGridProperties(params.grid_cell_counts, params.cell_dimensions);
    pressure_solver->setSharedObjects(this->shared_objects);
    pressure_solver->setFluidProperties(fluid_props);
}

void FiniteVolumeGrid::setGridOrigin(double x, double y, double z)
{
    this->params.grid_origin[0] = x;
    this->params.grid_origin[1] = y;
    this->params.grid_origin[2] = z;
}

void FiniteVolumeGrid::setGridCellCounts(int nx, int ny, int nz)
{
    this->params.grid_cell_counts[0] = nx;
    this->params.grid_cell_counts[1] = ny;
    this->params.grid_cell_counts[2] = nz;

    for(int i=0; i<3; i++)
        this->params.cell_dimensions[i] = this->params.grid_dimensions[i] / static_cast<double>(this->params.grid_cell_counts[i]); //mm
}

void FiniteVolumeGrid::setGridDimensions(double length_x, double length_y, double length_z)
{
    this->params.grid_dimensions[0] = length_x;
    this->params.grid_dimensions[1] = length_y;
    this->params.grid_dimensions[2] = length_z;

    for (int i=0; i<3; i++)
        this->params.cell_dimensions[i] = this->params.grid_dimensions[i] / static_cast<double>(this->params.grid_cell_counts[i]); //mm
}

void FiniteVolumeGrid::updateGrid()
{
    cout << "Updating grid..." << endl;

    renderer->RemoveActor(glyphactor);
    double vec[3];
    switch(vector_mode) {
    case lpt::VELOCITY:
        {
            boost_accumulator velocity_stats;

            for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
            {
                size_t total_count = extract_result<tag::count>( velocity_accumulators[(int)i][0] );
                if ( total_count > 0)
                        arrow_source_ids->SetTuple1(i, 0);
                vec[0] = extract_result<tag::weighted_mean>( velocity_accumulators[(int)i][0] );
                vec[1] = extract_result<tag::weighted_mean>( velocity_accumulators[(int)i][1] );
                vec[2] = extract_result<tag::weighted_mean>( velocity_accumulators[(int)i][2] );
                if (vec[0] == vec[0] && vec[1] == vec[1] && vec[2] == vec[2] ) { // check for NaN values
                    double mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                    magnitudes->SetTuple1(i, mag );
                    vectors->SetTuple3(i, vec[0], vec[1], vec[2] );
                    velocity_stats(mag);
                } else {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            }
            double min = extract_result<tag::min>( velocity_stats );
            double max = extract_result<tag::max>( velocity_stats );
            double mean = extract_result<tag::mean>( velocity_stats );
            double stdev = extract_result<tag::variance>( velocity_stats );
            stdev =  (stdev > 0 ? sqrt(stdev) : stdev);

            cout << "Velocity Magnitude (m/s):  Mean = " << mean <<" +- "<< stdev << ", Max = "<< max << ", Min = " << min << endl;
            double avg_cell_size = (params.cell_dimensions[0] + params.cell_dimensions[0] + params.cell_dimensions[0] ) / 3.0;
            //glyph3Dmapper->SetScaleFactor( avg_cell_size / max );
            lookuptable->SetRange(min, max);
            break;
        }
    case lpt::MASS_RESIDUAL:
        {
            if ( gridIsFull() ) {
                boost_accumulator residual_stats;
                calcMassResiduals();

                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                        residual_stats(mass_residuals[i]);
                        magnitudes->SetTuple1(i, mass_residuals[i] );
                        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                double min = extract_result<tag::min>( residual_stats );
                double max = extract_result<tag::max>( residual_stats );
                double mean = extract_result<tag::mean>( residual_stats );
                double stdev = extract_result<tag::variance>( residual_stats );
                stdev =  (stdev > 0 ? sqrt(stdev) : stdev);

                double residual_range[2] =  {min, max };
                //histogram->setBins(100, residual_range );
                //histogram->updateData(mass_residuals);
                lookuptable->SetRange(residual_range);

                cout << "Relative Mass Residual:  Mean = " << mean <<" +- "<< stdev << ", Max = "<< max << ", Min = " << min << endl;
            } else{
                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                cout << "Grid Not ready" << endl;
            }

            break;
        }
    case lpt::VORTICITY:
        {
            if ( gridIsFull() ) {
                boost_accumulator vorticity_stats;
                calcVorticity();

                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    vec[0] = vorticity[i][0];
                    vec[1] = vorticity[i][1];
                    vec[2] = vorticity[i][2];
                    if (vec[0] == vec[0] && vec[1] == vec[1] && vec[2] == vec[2] ) { // check for NaN values
                        double mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                        vorticity_stats(mag);
                        magnitudes->SetTuple1(i, mag );
                        vectors->SetTuple3(i, vec[0], vec[1], vec[2] );
                    } else {
                        magnitudes->SetTuple1(i, 0.0 );
                        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                    }

                }
                double min = extract_result<tag::min>( vorticity_stats );
                double max = extract_result<tag::max>( vorticity_stats );
                double mean = extract_result<tag::mean>( vorticity_stats );
                double stdev = extract_result<tag::variance>( vorticity_stats );
                stdev =  (stdev > 0 ? sqrt(stdev) : stdev);

                double residual_range[2] =  {min, max };
                //histogram->setBins(100, residual_range );
                //histogram->updateData(mass_residuals);
                lookuptable->SetRange(residual_range);

                cout << "Vorticity :  Mean = " << mean <<" +- "<< stdev << ", Max = "<< max << ", Min = " << min << endl;
            } else{
                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                cout << "Grid Not ready" << endl;
            }

            break;
        }
    case lpt::TURBULENT_KINETIC_ENERGY:
        {
            if ( gridIsFull() ) {
                boost_accumulator turbulent_ke_stats;

                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    double u_variance = extract_result<tag::weighted_variance>( velocity_accumulators[(int)i][0] );
                    double v_variance = extract_result<tag::weighted_variance>( velocity_accumulators[(int)i][1] );
                    double w_variance = extract_result<tag::weighted_variance>( velocity_accumulators[(int)i][2] );

                    double turbulent_kinetic_energy = 1.0 / 2.0 * (u_variance + v_variance + w_variance);
                    turbulent_ke_stats(turbulent_kinetic_energy);
                    magnitudes->SetTuple1(i, turbulent_kinetic_energy );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                double min = extract_result<tag::min>( turbulent_ke_stats );
                double max = extract_result<tag::max>( turbulent_ke_stats );
                double mean = extract_result<tag::mean>( turbulent_ke_stats );
                double stdev = extract_result<tag::variance>( turbulent_ke_stats );
                stdev =  (stdev > 0 ? sqrt(stdev) : stdev);

                double residual_range[2] =  {min, max };

                lookuptable->SetRange(residual_range);

                cout << "Turbulent KE m2/s2:  Mean = " << mean << " +- " << stdev << ", Max = "<< max << ", Min = " << min << endl;
            } else{
                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                cout << "Grid Not ready" << endl;
            }
            break;
        }
    case lpt::ACCELERATION:
        {
            boost_accumulator acceleration_stats;

            for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
            {
                size_t total_count = extract_result<tag::count>( acceleration_accumulators[(int)i][0] );
                if ( total_count > 0)
                        arrow_source_ids->SetTuple1(i, 0);
                vec[0] = extract_result<tag::weighted_mean>( acceleration_accumulators[(int)i][0] );
                vec[1] = extract_result<tag::weighted_mean>( acceleration_accumulators[(int)i][1] );
                vec[2] = extract_result<tag::weighted_mean>( acceleration_accumulators[(int)i][2] );
                if (vec[0] == vec[0] && vec[1] == vec[1] && vec[2] == vec[2] ) { // check for NaN values
                    double mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                    acceleration_stats(mag);
                    magnitudes->SetTuple1(i, mag );
                    vectors->SetTuple3(i, vec[0], vec[1], vec[2] );
                } else {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            }
            double min = extract_result< tag::min >( acceleration_stats );
            double max = extract_result< tag::max >( acceleration_stats );
            double mean = extract_result< tag::mean >( acceleration_stats );
            double stdev = extract_result<tag::variance>( acceleration_stats );
            stdev =  (stdev > 0 ? sqrt(stdev) : stdev);
            lookuptable->SetRange( min, max );
            double avg_cell_size = (params.cell_dimensions[0] + params.cell_dimensions[0] + params.cell_dimensions[0] ) / 3.0;
            //glyph3Dmapper->SetScaleFactor( avg_cell_size / max );
            cout << "Acceleration magnitude (m/s2):  Mean = " << mean <<" +- "<< stdev << ", Max = "<< max <<  ", Min = "<< min << endl;
            break;
        }
    case lpt::VELOCITY_SD:
        {
            boost_accumulator velocity_stdev_stats;

            for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
            {
                size_t total_count = extract_result<tag::count>( acceleration_accumulators[(int)i][0] );
                if (total_count > 10) {
                    vec[0] = extract_result<tag::weighted_variance>( velocity_accumulators[(int)i][0] );
                    vec[1] = extract_result<tag::weighted_variance>( velocity_accumulators[(int)i][1] );
                    vec[2] = extract_result<tag::weighted_variance>( velocity_accumulators[(int)i][2] );
                    if (vec[0] == vec[0] && vec[1] == vec[1] && vec[2] == vec[2] ) { // check for NaN values
                        double variance =  sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                        double standard_deviation = (variance > 0 ? sqrt(variance) : variance);
                        velocity_stdev_stats( standard_deviation );
                        magnitudes->SetTuple1(i, standard_deviation);
                        vectors->SetTuple3(i, vec[0], vec[1], vec[2] );
                    } else {
                        magnitudes->SetTuple1(i, 0.0 );
                        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                    }
                }else {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            }

            double min = extract_result< tag::min >( velocity_stdev_stats );
            double max = extract_result< tag::max >( velocity_stdev_stats );
            double mean = extract_result< tag::mean >( velocity_stdev_stats );
            double stdev = extract_result<tag::variance>( velocity_stdev_stats );
            stdev =  (stdev > 0 ? sqrt(stdev) : stdev);
            lookuptable->SetRange(min, max);
            cout << "Velocity Standard Deviation m/s: Mean = " << mean <<" +- "<< stdev << ", Max = "<< max <<  ", Min = "<< min << endl;
            break;
        }
    case lpt::ACCELERATION_SD:
        {
            boost_accumulator acceleration_stdev_stats;

            for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
            {
                size_t total_count = extract_result<tag::count>( acceleration_accumulators[(int)i][0] );
                if (total_count > 10) {
                    vec[0] = extract_result<tag::weighted_variance>( acceleration_accumulators[(int)i][0] );
                    vec[1] = extract_result<tag::weighted_variance>( acceleration_accumulators[(int)i][1] );
                    vec[2] = extract_result<tag::weighted_variance>( acceleration_accumulators[(int)i][2] );
                    if (vec[0] == vec[0] && vec[1] == vec[1] && vec[2] == vec[2] ) { // check for NaN values
                        double variance =  sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
                        double standard_deviation = (variance > 0 ? sqrt(variance) : variance);
                        acceleration_stdev_stats( standard_deviation );
                        magnitudes->SetTuple1(i, standard_deviation);
                        vectors->SetTuple3(i, vec[0], vec[1], vec[2] );
                    } else {
                        magnitudes->SetTuple1(i, 0.0 );
                        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                    }
                }else {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            }
            double min = extract_result< tag::min >( acceleration_stdev_stats );
            double max = extract_result< tag::max >( acceleration_stdev_stats );
            double mean = extract_result< tag::mean >( acceleration_stdev_stats );
            double stdev = extract_result<tag::variance>( acceleration_stdev_stats );
            stdev =  (stdev > 0 ? sqrt(stdev) : stdev);
            lookuptable->SetRange(min, max);
            cout << "Acceleration Standard Deviation m/s2:  Mean = " << mean <<" +- "<< stdev << ", Max = "<< max <<  ", Min = "<< min << endl;
            break;
        }
    case lpt::COUNT:
        {
            boost_accumulator count_stats;

            for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
            {
                size_t total_count = extract_result<tag::count>( acceleration_accumulators[(int)i][0] );
                if (total_count == total_count ) { // check for NaN values
                    count_stats(total_count);
                    magnitudes->SetTuple1(i, total_count );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                } else {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            }

            double min = extract_result<tag::min>( count_stats );
            double max = extract_result<tag::max>( count_stats );
            double mean = extract_result<tag::mean>( count_stats );
            double stdev = extract_result<tag::variance>( count_stats );
            stdev =  (stdev > 0 ? sqrt(stdev) : stdev);
            double range[2] =  {min, max};//{ mean - 4.0 * error_of_mean, mean + 4.0 * error_of_mean };
            lookuptable->SetRange(range);

            cout << "Count: Mean = " << mean <<" +- "<< stdev << ", Max = "<< max <<  ", Min = "<< min << endl;
            break;
        }
    case lpt::PRESSURE:
        {
            if ( gridIsFull() ) {

                boost::thread pressure_solver_thread(&FiniteVolumeGrid::calcPressures, this );
                const vector<double>& pressure = pressure_solver->getPressureField();

                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, pressure[i] );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            } else{
                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                cout << "Grid Not ready" << endl;
            }

            break;
        }
    case lpt::PRESSURE_LAGRANGIAN:
        {
            bool ready = true;
            int counter = 0;
            for (int p = 0; p < pressure_grad_accumulators.size() ; ++p) {
                ready &= ( extract_result<tag::count>( pressure_grad_accumulators[p][0] ) > 0 ? true : false);
                if (!ready)
                    ++counter;
            }
            if ( ready ) {
                boost::thread pressure_solver_thread(&FiniteVolumeGrid::calcPressuresLagrangian, this );
                const vector<double>& pressure = lagrangian_pressure_solver->getPressureField();

                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, pressure[i] );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
            } else{
                for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
                {
                    magnitudes->SetTuple1(i, 0.0 );
                    vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
                }
                cout << " Number empty = " << counter << endl;
                cout << "Grid Not ready" << endl;
            }

            break;
        }
    default:
        break;
    }

    magnitudes->Modified();
    vectors->Modified();
    grid->GetCellData()->Modified();
    lookuptable->Modified();
    grid->GetProducerPort()->Modified();
    grid->Modified();
    glyph3Dmapper->Update();
    if (! plane_widget->GetEnabled() )
        renderer->AddActor(glyphactor);
    cout << "Grid update complete" << endl;
}

void FiniteVolumeGrid::resetGrid()
{
    renderer->RemoveActor(glyphactor);
    for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
    {
        //double vec[3];
        //boost_accumulator newx, newy, newz;
        velocity_accumulators[(int)i][0] = boost_weighted_accumulator();
        velocity_accumulators[(int)i][1] = boost_weighted_accumulator();
        velocity_accumulators[(int)i][2] = boost_weighted_accumulator();

        acceleration_accumulators[(int)i][0] = boost_weighted_accumulator();
        acceleration_accumulators[(int)i][1] = boost_weighted_accumulator();
        acceleration_accumulators[(int)i][2] = boost_weighted_accumulator();

        pressure_grad_accumulators[(int)i][0] = boost_weighted_accumulator();
        pressure_grad_accumulators[(int)i][1] = boost_weighted_accumulator();
        pressure_grad_accumulators[(int)i][2] = boost_weighted_accumulator();

        for (int n = 0; n < 6; ++n)
            reynolds_stress_accumulators[(int)i][n] = boost_weighted_covariance_accumulator();

        magnitudes->SetTuple1(i, 0.0 );
        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
        arrow_source_ids->SetTuple1(i, 1);
    }
    pressure_solver->resetPressureField();
    magnitudes->Modified();
    vectors->Modified();
    grid->Modified();
    if ( ! plane_widget->GetEnabled() )
        renderer->AddActor(glyphactor);
}

bool FiniteVolumeGrid::gridIsFull()
{
    bool ready = true;
    for (int p = 0; p < velocity_accumulators.size() ; ++p)
        ready &= ( extract_result<tag::count>( velocity_accumulators[p][0] ) > 0 ? true : false);
    return ready;
}

void FiniteVolumeGrid::updateAccumulators(lpt::Trajectory3d* traj_ptr, lpt::ParticleVectors& current_particle )
{
    int subIdNotUsed;
    double pcoordsNotUsed[3];
    double weightsNotUsed[8]; //needs to be the same length as the number of points in the cell

    array<double, 3>& X = current_particle[0]; //position
    array<double, 3>& U = current_particle[1]; //velocity
    array<double, 3>& A = current_particle[2]; //acceleration
    array<double, 3>& dPdX = current_particle[3]; //pressure gradient

    //vtkIdType last_cell_id = ( traj_ptr->last_cell_id >= 0 && traj_ptr->last_cell_id < grid->GetNumberOfCells() ? traj_ptr->last_cell_id : 0);
    //vtkSmartPointer<vtkGenericCell> gen_cell = vtkSmartPointer<vtkGenericCell>::New();

    //vtkCell* last_cell = grid->GetCell( last_cell_id );

    vtkIdType cell_id = grid->FindCell(
        current_particle[0].data(),
        NULL,
        0,
        1e-6, subIdNotUsed, pcoordsNotUsed, weightsNotUsed
        );

    if (cell_id >= 0 && cell_id < grid->GetNumberOfCells() ) {

        array<int,3> ijk_base;
        this->computeCellCoordsIJK(cell_id, ijk_base);
        array<double,3> base_cell_center;

        base_cell_center[0] = static_cast<double>( ijk_base[0] ) * params.cell_dimensions[0] + 0.5 * params.cell_dimensions[0] + params.grid_origin[0];
        base_cell_center[1] = static_cast<double>( ijk_base[1] ) * params.cell_dimensions[1] + 0.5 * params.cell_dimensions[1] + params.grid_origin[1];
        base_cell_center[2] = static_cast<double>( ijk_base[2] ) * params.cell_dimensions[2] + 0.5 * params.cell_dimensions[2] + params.grid_origin[2];

        for (int n = 0; n < cell_neighborhood.size(); ++n) {
            array<int, 3> ijk = cell_neighborhood[n];
            ijk[0] += ijk_base[0];
            ijk[1] += ijk_base[1];
            ijk[2] += ijk_base[2];

            // X Boundary Check
            if (ijk[0] < 0 || ijk[0] > params.grid_cell_counts[0] - 1 )
                continue;

            // Y Boundary Check
            if (ijk[1] < 0 || ijk[1] > params.grid_cell_counts[1] - 1 )
                continue;

            // Z Boundary Check
            if (ijk[2] < 0 || ijk[2] > params.grid_cell_counts[2] - 1 )
                continue;

            array<double, 3> X2;

            X2[0] = base_cell_center[0] + cell_neighborhood[n][0] * params.cell_dimensions[0];
            X2[1] = base_cell_center[1] + cell_neighborhood[n][1] * params.cell_dimensions[1];
            X2[2] = base_cell_center[2] + cell_neighborhood[n][2] * params.cell_dimensions[2];

            double distance_weight =
                1.0 / (
                sqrt(
                ( X[0] - X2[0] ) * ( X[0] - X2[0] )
                +   ( X[1] - X2[1] ) * ( X[1] - X2[1] )
                +	( X[2] - X2[2] ) * ( X[2] - X2[2] )
                )
                );

            cell_id = getGridIndex( ijk[0], ijk[1], ijk[2] );

            velocity_accumulators[cell_id][0]( U[0], weight = distance_weight ); // m/s
            velocity_accumulators[cell_id][1]( U[1], weight = distance_weight );
            velocity_accumulators[cell_id][2]( U[2], weight = distance_weight );

            acceleration_accumulators[cell_id][0]( A[0],  weight = distance_weight ); // m/s2
            acceleration_accumulators[cell_id][1]( A[1],  weight = distance_weight );
            acceleration_accumulators[cell_id][2]( A[2],  weight = distance_weight );

            reynolds_stress_accumulators[cell_id][0]( U[0], weight = distance_weight, covariate1 = U[1] );  // Covariance of (u, v) = Cov(v, u)
            reynolds_stress_accumulators[cell_id][1]( U[0], weight = distance_weight, covariate1 = U[2] );  // Covariance of (u, w) = Cov(w, u)
            reynolds_stress_accumulators[cell_id][2]( U[1], weight = distance_weight, covariate1 = U[2] );	 // Covariance of (v, w) = Cov(w, v)

            if (dPdX[0] != 0 || dPdX[1] != 0 || dPdX[2] != 0) {
                pressure_grad_accumulators[cell_id][0]( dPdX[0],  weight = distance_weight );
                pressure_grad_accumulators[cell_id][1]( dPdX[1],  weight = distance_weight );
                pressure_grad_accumulators[cell_id][2]( dPdX[2],  weight = distance_weight );
            }
        }

        // update cell counts here!
        traj_ptr->last_cell_id = static_cast<int>(cell_id);
    }
}

void FiniteVolumeGrid::savePlaneData()
{
    lpt::VectorField velocities;
    lpt::VectorField accelerations;

    lpt::getVectorField(velocity_accumulators, velocities);
    lpt::getVectorField(acceleration_accumulators, accelerations);

    lpt::ReynoldsStressField stress_field;
    lpt::getReynoldsStressField(velocity_accumulators, reynolds_stress_accumulators, stress_field);


    const vector<double>& pressures = pressure_solver->getPressureField();
    const vector<double>& lagrangian_pressures = lagrangian_pressure_solver->getPressureField();

    int k = this->params.grid_cell_counts[2] / 2 + 1;
    cout << "Saving XY plane data: k = " << k << endl;

    stringstream header;

    header << "\tGrid count (X,Y,Z)\t"
        << this->params.grid_cell_counts[0] << "\t" << this->params.grid_cell_counts[1] << "\t" << this->params.grid_cell_counts[2]
        << "\tcell size\t"  << this->params.cell_dimensions[0] << "\t"  << this->params.cell_dimensions[1] << "\t"  << this->params.cell_dimensions[2]
        << "\torigin\t" << this->params.grid_origin[0] << "\t"<< this->params.grid_origin[1] << "\t"<< this->params.grid_origin[2] << "\t" <<  endl;


    ofstream U_out(this->shared_objects->output_path + "U_XY.txt");
    ofstream V_out(this->shared_objects->output_path + "V_XY.txt");
    ofstream W_out(this->shared_objects->output_path + "W_XY.txt");

    ofstream u_var_out(this->shared_objects->output_path + "u_var_XY.txt");
    ofstream v_var_out(this->shared_objects->output_path + "v_var_XY.txt");
    ofstream w_var_out(this->shared_objects->output_path + "w_var_XY.txt");

    ofstream tke_out(this->shared_objects->output_path + "TKE_XY.txt");
    ofstream uv_stress_out(this->shared_objects->output_path + "uv_stress_planeXY.txt");
    ofstream uw_stress_out(this->shared_objects->output_path + "uw_stress_planeXY.txt");
    ofstream vw_stress_out(this->shared_objects->output_path + "vw_stress_planeXY.txt");

    ofstream P_e_out(this->shared_objects->output_path + "P_e_XY.txt");
    ofstream P_l_out(this->shared_objects->output_path + "P_l_XY.txt");
    ofstream all_data_out(this->shared_objects->output_path + "all_data.txt");
	ofstream axis(this->shared_objects->output_path + "axis.txt");

    U_out << header.str();
    U_out << "XY Plane Z = \t" << k << endl;

    V_out << header.str();
    V_out << "XY Plane Z = \t" << k << endl;

    W_out << header.str();
    W_out << "XY Plane Z = \t" << k << endl;

    u_var_out << header.str();
    u_var_out << "XY Plane Z = \t" << k << endl;

    v_var_out << header.str();
    v_var_out << "XY Plane Z = \t" << k << endl;

    w_var_out << header.str();
    w_var_out << "XY Plane Z = \t" << k << endl;

    tke_out << header.str();
    tke_out << "XY Plane Z = \t" << k << endl;

    uv_stress_out << header.str();
    uv_stress_out << "XY Plane Z = \t" << k << endl;

    uw_stress_out << header.str();
    uw_stress_out << "XY Plane Z = \t" << k << endl;

    vw_stress_out << header.str();
    vw_stress_out << "XY Plane Z = \t" << k << endl;

    P_e_out << header.str();
    P_e_out << "XY Plane Z = \t" << k << endl;

    P_l_out << header.str();
    P_l_out << "XY Plane Z = \t" << k << endl;

    for (int i = 0; i < this->params.grid_cell_counts[0]; ++i) {
        for (int j = 0; j < this->params.grid_cell_counts[1]; ++j) {
            int cell_id = this->getGridIndex(i,j,k);

            double u_variance = stress_field[cell_id][0];
            double v_variance = stress_field[cell_id][1];
            double w_variance = stress_field[cell_id][2];

            double uv_covariance = stress_field[cell_id][3];
            double uw_covariance = stress_field[cell_id][4];
            double vw_covariance = stress_field[cell_id][5];

            double turbulent_kinetic_energy = 1.0 / 2.0 * (u_variance + v_variance + w_variance);

            U_out << velocities[cell_id][0] << "\t";
            V_out << velocities[cell_id][1] << "\t";
            W_out << velocities[cell_id][2] << "\t";

            tke_out << turbulent_kinetic_energy << "\t";

            u_var_out << u_variance << "\t";
            v_var_out << v_variance << "\t";
            w_var_out << w_variance << "\t";

            uv_stress_out << uv_covariance << "\t";
            uw_stress_out << uw_covariance << "\t";
            vw_stress_out << vw_covariance << "\t";

            P_e_out << pressures[cell_id] << "\t";
            P_l_out << lagrangian_pressures[cell_id] << "\t";
        }

        U_out << endl;
        V_out << endl;
        W_out << endl;

        u_var_out << endl;
        v_var_out << endl;
        w_var_out << endl;

        tke_out << endl;

        uv_stress_out << endl;
        uw_stress_out << endl;
        vw_stress_out << endl;

        P_e_out << endl;
        P_l_out << endl;

    }
    U_out.close();
    V_out.close();
    W_out.close();

    u_var_out.close();
    v_var_out.close();
    w_var_out.close();

    tke_out.close();

    uv_stress_out.close();
    uw_stress_out.close();
    vw_stress_out.close();

    P_e_out.close();
    P_l_out.close();

    int j = this->params.grid_cell_counts[1] / 2 + 1;
    cout << "Saving XZ plane data: j = " << j << endl;

    U_out.open(this->shared_objects->output_path + "U_XZ.txt");
    V_out.open(this->shared_objects->output_path + "V_XZ.txt");
    W_out.open(this->shared_objects->output_path + "W_XZ.txt");

    u_var_out.open(this->shared_objects->output_path + "u_var_XZ.txt");
    v_var_out.open(this->shared_objects->output_path + "v_var_XZ.txt");
    w_var_out.open(this->shared_objects->output_path + "w_var_XZ.txt");

    tke_out.open(this->shared_objects->output_path + "TKE_XZ.txt");
    uv_stress_out.open(this->shared_objects->output_path + "uv_stress_planeXZ.txt");
    uw_stress_out.open(this->shared_objects->output_path + "uw_stress_planeXZ.txt");
    vw_stress_out.open(this->shared_objects->output_path + "vw_stress_planeXZ.txt");


    P_e_out.open(this->shared_objects->output_path + "P_e_XZ.txt");
    P_l_out.open(this->shared_objects->output_path + "P_l_XZ.txt");

    U_out << header.str();
    U_out << "XZ Plane Y = \t" << j << endl;

    V_out << header.str();
    V_out << "XZ Plane Y = \t" << j << endl;

    W_out << header.str();
    W_out << "XZ Plane Y = \t" << j << endl;

    u_var_out << header.str();
    u_var_out << "XZ Plane Y = \t" << j << endl;

    v_var_out << header.str();
    v_var_out << "XZ Plane Y = \t" << j << endl;

    w_var_out << header.str();
    w_var_out << "XZ Plane Y = \t" << j << endl;

    tke_out << header.str();
    tke_out << "XZ Plane Y = \t" << j << endl;

    uv_stress_out << header.str();
    uv_stress_out << "XZ Plane Y = \t" << j << endl;

    uw_stress_out << header.str();
    uw_stress_out <<  "XZ Plane Y = \t" << j << endl;

    vw_stress_out << header.str();
    vw_stress_out <<  "XZ Plane Y = \t" << j << endl;

    P_e_out << header.str();
    P_e_out << "XZ Plane Y = \t" << j << endl;

    P_l_out << header.str();
    P_l_out << "XZ Plane Y = \t" << j << endl;

    for (int i = 0; i < this->params.grid_cell_counts[0]; ++i) {
        for (int k = 0; k < this->params.grid_cell_counts[2]; ++k) {
            int cell_id = this->getGridIndex(i,j,k);

            double u_variance = stress_field[cell_id][0];
            double v_variance = stress_field[cell_id][1];
            double w_variance = stress_field[cell_id][2];

            double uv_covariance = stress_field[cell_id][3];
            double uw_covariance = stress_field[cell_id][4];
            double vw_covariance = stress_field[cell_id][5];

            double turbulent_kinetic_energy = 1.0 / 2.0 * (u_variance + v_variance + w_variance);

            U_out << velocities[cell_id][0] << "\t";
            V_out << velocities[cell_id][1] << "\t";
            W_out << velocities[cell_id][2] << "\t";

            tke_out << turbulent_kinetic_energy << "\t";

            u_var_out << u_variance << "\t";
            v_var_out << v_variance << "\t";
            w_var_out << w_variance << "\t";

            uv_stress_out << uv_covariance << "\t";
            uw_stress_out << uw_covariance << "\t";
            vw_stress_out << vw_covariance << "\t";

            P_e_out << pressures[cell_id] << "\t";
            P_l_out << lagrangian_pressures[cell_id] << "\t";
        }
        U_out << endl;
        V_out << endl;
        W_out << endl;

        u_var_out << endl;
        v_var_out << endl;
        w_var_out << endl;

        tke_out << endl;

        uv_stress_out << endl;
        uw_stress_out << endl;
        vw_stress_out << endl;

        P_e_out << endl;
        P_l_out << endl;
    }
    U_out.close();
    V_out.close();
    W_out.close();

    u_var_out.close();
    v_var_out.close();
    w_var_out.close();

    tke_out.close();

    uv_stress_out.close();
    uw_stress_out.close();
    vw_stress_out.close();

    all_data_out << header.str();
    for (int cell_id = 0; cell_id < velocities.size(); cell_id++) {
        double u_variance = stress_field[cell_id][0];
        double v_variance = stress_field[cell_id][1];
        double w_variance = stress_field[cell_id][2];

        double uv_covariance = stress_field[cell_id][3];
        double uw_covariance = stress_field[cell_id][4];
        double vw_covariance = stress_field[cell_id][5];

        double turbulent_kinetic_energy = 1.0 / 2.0 * (u_variance + v_variance + w_variance);
        array<int,3> ijk;
        this->computeCellCoordsIJK(cell_id, ijk);

        all_data_out << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << "\t";
        all_data_out << velocities[cell_id][0] << "\t" << velocities[cell_id][1] << "\t" << velocities[cell_id][2] << "\t";

        all_data_out << u_variance << "\t" << v_variance << "\t" << w_variance << "\t";

        all_data_out << uv_covariance << "\t" << uw_covariance << "\t" << vw_covariance << "\t";

        all_data_out << turbulent_kinetic_energy <<"\t";

        all_data_out << pressures[cell_id] << "\t" << lagrangian_pressures[cell_id] << "\t";

        all_data_out << endl;
    }
    all_data_out.close();

	array<int,3> ijk;
	for (int i=0; i<this->params.grid_cell_counts[0]; i++) {
		double max = std::numeric_limits<double>::lowest();
		for (int j=0; j<this->params.grid_cell_counts[1]; j++) {
			for (int k=0; k<this->params.grid_cell_counts[2]; k++) {
				int cell_id = this->getGridIndex(i,j,k);
				if (velocities[cell_id][0] > max) {
					ijk[0] = i;
					ijk[1] = j;
					ijk[2] = k;
					max = velocities[cell_id][0];
				}
			}
		}
		axis << ijk[0] << "\t" << ijk[1] << "\t" << ijk[2] << endl;
	}
	axis.close();
	
    cout << "Done Saving data" << endl;
}

void FiniteVolumeGrid::calcPressures()
{
    cout << "Thread " << boost::this_thread::get_id() << "  solving pressure field" << endl;
    lpt::VectorField velocities;
    lpt::VectorField accelerations;
    lpt::ReynoldsStressField stress_field;

    lpt::getVectorField(velocity_accumulators, velocities);
    lpt::getVectorField(acceleration_accumulators, accelerations);
    lpt::getReynoldsStressField(velocity_accumulators, reynolds_stress_accumulators, stress_field);

    pressure_solver->solve(velocities, accelerations, stress_field);
    cout << "Pressure solution complete" << endl;

    boost_accumulator pressure_stats;
    const vector<double>& pressure = pressure_solver->getPressureField();

    for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
    {
        pressure_stats(pressure[i]);
        magnitudes->SetTuple1(i, pressure[i] );
        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
    }

    double min = extract_result<tag::min>( pressure_stats );
    double max = extract_result<tag::max>( pressure_stats );
    double mean = extract_result<tag::mean>( pressure_stats );
    double stdev = extract_result<tag::variance>( pressure_stats );
    stdev =  (stdev > 0 ? sqrt(stdev) : stdev);
    double pressure_range[2] =  {min, max};//{ mean - 4.0 * error_of_mean, mean + 4.0 * error_of_mean };
    lookuptable->SetRange(pressure_range);

    magnitudes->Modified();
    vectors->Modified();
    grid->GetCellData()->Modified();
    lookuptable->Modified();
    grid->GetProducerPort()->Modified();
    grid->Modified();


    cout << "Pressure (Pa): Mean = " << mean <<" +- "<< stdev << ", Max = "<< max <<  ", Min = "<< min << endl;
}

void FiniteVolumeGrid::calcPressuresLagrangian()
{
    cout << "Thread " << boost::this_thread::get_id() << "  solving pressure field" << endl;
    lpt::VectorField pressure_grads;
    getVectorField(pressure_grad_accumulators, pressure_grads);

    lagrangian_pressure_solver->solve(pressure_grads);
    cout << "Pressure solution complete" << endl;

    boost_accumulator pressure_stats;
    const vector<double>& pressure = lagrangian_pressure_solver->getPressureField();

    for(vtkIdType i = 0; i < grid->GetNumberOfCells(); i++)
    {
        pressure_stats(pressure[i]);
        magnitudes->SetTuple1(i, pressure[i] );
        vectors->SetTuple3(i, 0.0, 0.0, 0.0 );
    }

    double min = extract_result<tag::min>( pressure_stats );
    double max = extract_result<tag::max>( pressure_stats );
    double mean = extract_result<tag::mean>( pressure_stats );
    double stdev = extract_result<tag::variance>( pressure_stats );
    stdev =  (stdev > 0 ? sqrt(stdev) : stdev);
    double pressure_range[2] =  {min, max};//{ mean - 4.0 * error_of_mean, mean + 4.0 * error_of_mean };
    lookuptable->SetRange(pressure_range);

    magnitudes->Modified();
    vectors->Modified();
    grid->GetCellData()->Modified();
    lookuptable->Modified();
    grid->GetProducerPort()->Modified();
    grid->Modified();


    cout << "Pressure (Pa): Mean = " << mean <<" +- "<< stdev << ", Max = "<< max <<  ", Min = "<< min << endl;

}

void FiniteVolumeGrid::calcMassResiduals()
{
    mass_residuals.resize( grid->GetNumberOfCells(), 0 );
    cout << "Calculating mass residual " << endl;

    int p, u1, u2, v1, v2, w1, w2;

    double dx = params.cell_dimensions[0] / 1000.0;
    double dy = params.cell_dimensions[1] / 1000.0;
    double dz = params.cell_dimensions[2] / 1000.0;

    double norm2 = 0;
    for ( int i = 1; i < params.grid_cell_counts[0] - 1; ++i ) {
        for ( int j = 1; j < params.grid_cell_counts[1] - 1; ++j ) {
            for ( int k = 1; k < params.grid_cell_counts[2] - 1; ++k ) {

                p  = getGridIndex( i     , j     , k     );
                u1 = getGridIndex( i + 1 , j     , k     );
                u2 = getGridIndex( i - 1 , j     , k     );
                v1 = getGridIndex( i     , j + 1 , k     );
                v2 = getGridIndex( i     , j - 1 , k     );
                w1 = getGridIndex( i     , j     , k + 1 );
                w2 = getGridIndex( i     , j     , k - 1 );

                double U1 = extract_result<tag::weighted_mean>( velocity_accumulators[u1][0] );  // X velocity
                double U2 = extract_result<tag::weighted_mean>( velocity_accumulators[u2][0] );

                double V1 = extract_result<tag::weighted_mean>( velocity_accumulators[v1][1] );  // Y velocity
                double V2 = extract_result<tag::weighted_mean>( velocity_accumulators[v2][1] );

                double W1 = extract_result<tag::weighted_mean>( velocity_accumulators[w1][2] );  // Z velocity
                double W2 = extract_result<tag::weighted_mean>( velocity_accumulators[w2][2] );


                double U = extract_result<tag::weighted_mean>( velocity_accumulators[p][0] );  // X velocity
                double V = extract_result<tag::weighted_mean>( velocity_accumulators[p][1] );  // Y velocity
                double W = extract_result<tag::weighted_mean>( velocity_accumulators[p][2] );  // Z velocity

                mass_residuals[p] =  (U1 - U2) + (V1 - V2) + (W1 - W2); //( (U1 - U2) / (2.0 * dx) + (V1 - V2) / (2.0 * dy) + (W1 - W2) / (2.0 * dz) );
                mass_residuals[p] = sqrt( mass_residuals[p] * mass_residuals[p] ) / sqrt(U*U + V*V + W*W);
                norm2 += mass_residuals[p] * mass_residuals[p];
            }
        }
    }

    cout << "Mass residual calculation complete: mass residual L2 norm = " << (norm2 > 0 ? sqrt(norm2) : norm2) << endl;
}

void FiniteVolumeGrid::calcVorticity()
{
    array<double,3> temp = {{ 0, 0, 0 }};
    vorticity.resize( grid->GetNumberOfCells(), temp );

    cout << "Calculating vorticity " << endl;

    int p, u1, u2, v1, v2, w1, w2;

    double dx = params.cell_dimensions[0] / 1000.0;
    double dy = params.cell_dimensions[1] / 1000.0;
    double dz = params.cell_dimensions[2] / 1000.0;

    double norm2 = 0;
    for ( int i = 1; i < params.grid_cell_counts[0] - 1; ++i ) {
        for ( int j = 1; j < params.grid_cell_counts[1] - 1; ++j ) {
            for ( int k = 1; k < params.grid_cell_counts[2] - 1; ++k ) {

                p  = getGridIndex( i     , j     , k     );
                u1 = getGridIndex( i + 1 , j     , k     );
                u2 = getGridIndex( i - 1 , j     , k     );
                v1 = getGridIndex( i     , j + 1 , k     );
                v2 = getGridIndex( i     , j - 1 , k     );
                w1 = getGridIndex( i     , j     , k + 1 );
                w2 = getGridIndex( i     , j     , k - 1 );

                double U1 = extract_result<tag::weighted_mean>( velocity_accumulators[u1][0] );  // X velocity
                double U2 = extract_result<tag::weighted_mean>( velocity_accumulators[u2][0] );

                double V1 = extract_result<tag::weighted_mean>( velocity_accumulators[v1][1] );  // Y velocity
                double V2 = extract_result<tag::weighted_mean>( velocity_accumulators[v2][1] );

                double W1 = extract_result<tag::weighted_mean>( velocity_accumulators[w1][2] );  // Z velocity
                double W2 = extract_result<tag::weighted_mean>( velocity_accumulators[w2][2] );


                double U = extract_result<tag::weighted_mean>( velocity_accumulators[p][0] );  // X velocity
                double V = extract_result<tag::weighted_mean>( velocity_accumulators[p][1] );  // Y velocity
                double W = extract_result<tag::weighted_mean>( velocity_accumulators[p][2] );  // Z velocity

                vorticity[p][0] =  (W1 - W2) / (2.0 * dy) - (V1 - V2) / ( 2.0 * dz );  // dW/dy - dV/dz

                vorticity[p][1] =  (U1 - U2) / (2.0 * dz) - (W1 - W2) / ( 2.0 * dx );  // dU/dz - dW/dx

                vorticity[p][2] =  (V1 - V2) / (2.0 * dx) - (U1 - U2) / ( 2.0 * dy ); //( (U1 - U2) / (2.0 * dx) + (V1 - V2) / (2.0 * dy) + (W1 - W2) / (2.0 * dz) );

                norm2 += vorticity[p][0] * vorticity[p][0] + vorticity[p][1] * vorticity[p][1] + vorticity[p][2] * vorticity[p][2];
            }
        }
    }

    cout << "Vorticity calculation complete: L2 norm = " << (norm2 > 0 ? sqrt(norm2) : norm2) << endl;
}

void FiniteVolumeGrid::setVectorMode(lpt::VectorMode mode)
{
    vector_mode = mode;
    switch (mode) {
    case lpt::VELOCITY:
            scalarbar->GetScalarBarActor()->SetTitle("Velocity Magnitude (m/s)");
            //glyph3Dmapper->SetScaleModeToScaleByMagnitude();
            glyph3Dmapper->SetSourceIndexArray("arrow");
            break;
    case lpt::VORTICITY:
            scalarbar->GetScalarBarActor()->SetTitle("Vorticity Magnitude (1/s)");
            //glyph3Dmapper->SetScaleModeToScaleByMagnitude();
            glyph3Dmapper->SetSourceIndexArray("arrow");
            break;
    case lpt::TURBULENT_KINETIC_ENERGY:
            scalarbar->GetScalarBarActor()->SetTitle("Turbulent Kinetic Energy (m2/s2)");
            //glyph3Dmapper->SetScaleModeToScaleByMagnitude();
            glyph3Dmapper->SetScaleFactor(params.fixed_scale_factor);
            glyph3Dmapper->SetSourceIndexArray("sphere");
            break;
    case lpt::TURBULENCE_DISSIPATION_RATE:
            scalarbar->GetScalarBarActor()->SetTitle("Turbulence Dissipation Rate (m2/s3)");
            //glyph3Dmapper->SetScaleModeToScaleByMagnitude();
            glyph3Dmapper->SetScaleFactor(params.fixed_scale_factor);
            glyph3Dmapper->SetSourceIndexArray("sphere");
            break;
    case lpt::MASS_RESIDUAL:
            scalarbar->GetScalarBarActor()->SetTitle("Mass Residual");
            glyph3Dmapper->SetSourceIndexArray("sphere");
            glyph3Dmapper->SetScaleModeToNoDataScaling();
            glyph3Dmapper->SetScaleFactor(params.fixed_scale_factor * 0.7);
            break;
    case lpt::ACCELERATION:
            scalarbar->GetScalarBarActor()->SetTitle("Acceleration Magnitude (m/s2)");
            //glyph3Dmapper->SetScaleModeToScaleByMagnitude();
            glyph3Dmapper->SetSourceIndexArray("arrow");
            break;
    case lpt::VELOCITY_SD:
            scalarbar->GetScalarBarActor()->SetTitle("Velocity Std Dev (m/s)");
            glyph3Dmapper->SetSourceIndexArray("sphere");
            glyph3Dmapper->SetScaleModeToNoDataScaling();
            glyph3Dmapper->SetScaleFactor(params.fixed_scale_factor * 0.7);
            break;
    case lpt::ACCELERATION_SD:
            scalarbar->GetScalarBarActor()->SetTitle("Acceleration Std Dev (m/s2)");
            glyph3Dmapper->SetSourceIndexArray("sphere");
            glyph3Dmapper->SetScaleModeToNoDataScaling();
            glyph3Dmapper->SetScaleFactor(params.fixed_scale_factor * 0.7);
            break;
    case lpt::COUNT:
            scalarbar->GetScalarBarActor()->SetTitle("Count");
            glyph3Dmapper->SetSourceIndexArray("sphere");
            glyph3Dmapper->SetScaleModeToNoDataScaling();
            glyph3Dmapper->SetScaleFactor( params.fixed_scale_factor * 0.7 );
            break;
    case lpt::PRESSURE:
            scalarbar->GetScalarBarActor()->SetTitle("Pressure (Pa)");
            glyph3Dmapper->SetSourceIndexArray("sphere");
            glyph3Dmapper->SetScaleModeToNoDataScaling();
            glyph3Dmapper->SetScaleFactor( params.fixed_scale_factor * 0.7 );
            break;
    case lpt::PRESSURE_LAGRANGIAN:
            scalarbar->GetScalarBarActor()->SetTitle("Pressure Lagrangian (Pa)");
            glyph3Dmapper->SetSourceIndexArray("sphere");
            glyph3Dmapper->SetScaleModeToNoDataScaling();
            glyph3Dmapper->SetScaleFactor( params.fixed_scale_factor * 0.7 );
            break;
    default:
        break;
    }
    updateGrid();
}

vtkSmartPointer<vtkLookupTable> FiniteVolumeGrid::getLookupTable()
{
    return lookuptable;
}

vtkSmartPointer<vtkScalarBarWidget> FiniteVolumeGrid::getScalarBarWidget()
{
    return scalarbar;
}

vtkSmartPointer<vtkImplicitPlaneWidget> FiniteVolumeGrid::getImplicitPlane()
{
    return plane_widget;
}

vtkSmartPointer<vtkActor> FiniteVolumeGrid::getImplicitPlaneActor()
{
    return planeactor;
}

/***** Definition of ParticlesVTK class *****/

ParticlesVTK::ParticlesVTK(vtkSmartPointer<vtkRenderer> renderer) : renderer(renderer), vector_mode(lpt::VELOCITY)
{
    sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetRadius( 3 );
    sphere->SetPhiResolution(10);
    sphere->SetThetaResolution(10);
    this->source = sphere;
    points = vtkSmartPointer<vtkPoints>::New();
    polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    velocity_magnitudes = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_magnitudes->SetName("velocity_mag");
    velocity_magnitudes->SetNumberOfComponents(1);

    acceleration_magnitudes = vtkSmartPointer<vtkDoubleArray>::New();
    acceleration_magnitudes->SetName("acceleration_mag");
    acceleration_magnitudes->SetNumberOfComponents(1);

    polydata->GetPointData()->AddArray(velocity_magnitudes);
    polydata->GetPointData()->AddArray(acceleration_magnitudes);
    polydata->GetPointData()->SetActiveScalars("velocity_mag");

    glyph3Dmapper = vtkSmartPointer<vtkGlyph3DMapper>::New();
    glyph3Dmapper->SetSourceConnection(source->GetOutputPort());
    glyph3Dmapper->SetInputConnection(polydata->GetProducerPort());
    glyph3Dmapper->SetScaleFactor(params.scale);

    glyph3Dmapper->UseLookupTableScalarRangeOn();
    glyph3Dmapper->SelectColorArray("velocity_mag");
    glyph3Dmapper->SetColorModeToMapScalars();
    glyph3Dmapper->SetScaleModeToNoDataScaling();

    lookuptable = vtkSmartPointer<vtkLookupTable>::New();
    lookuptable->SetRange( params.velocity_range.data() );
    lookuptable->SetNumberOfColors(256);
    lookuptable->SetHueRange(0.667, 0.0);
    lookuptable->Build();

    glyph3Dmapper->SetLookupTable( lookuptable );
    glyph3Dmapper->Update();

    scalarbar = vtkSmartPointer<vtkScalarBarWidget>::New();
    scalarbar->SetInteractor( renderer->GetRenderWindow()->GetInteractor() );
    scalarbar->GetScalarBarActor()->SetTitle("Velocity Magnitude (m/s)");
    scalarbar->GetScalarBarActor()->SetLookupTable( glyph3Dmapper->GetLookupTable() );
    scalarbar->SetRepositionable(false);
    scalarbar->SelectableOff();
    scalarbar->SetManagesCursor(false);
    scalarbar->SetProcessEvents(false);

    glyphactor = vtkSmartPointer<vtkActor>::New();
    glyphactor->SetMapper(glyph3Dmapper);
}

ParticlesVTK::~ParticlesVTK()
{
}

void ParticlesVTK::resizeGlyphArrays(size_t size)
{
    points->SetNumberOfPoints(size);
    velocity_magnitudes->SetNumberOfTuples(size);
    acceleration_magnitudes->SetNumberOfTuples(size);
}

void ParticlesVTK::setScalarsModified()
{
    points->Modified();
    velocity_magnitudes->Modified();
    acceleration_magnitudes->Modified();
}

void ParticlesVTK::updateParticle(int id, lpt::ParticleVectors& current_particle)
{
    points->SetPoint( id, current_particle[0].data() );
    double vec[3];
    vec[0] = current_particle[1][0]; // m/s
    vec[1] = current_particle[1][1];
    vec[2] = current_particle[1][2];
    velocity_magnitudes->SetTuple1(id, sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) );

    vec[0] = current_particle[2][0]; // m/s^2
    vec[1] = current_particle[2][1];
    vec[2] = current_particle[2][2];
    acceleration_magnitudes->SetTuple1(id, sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) );
}

void ParticlesVTK::setVectorMode(lpt::VectorMode mode)
{
    vector_mode = mode;
    switch (mode) {
    case lpt::VELOCITY:
            scalarbar->GetScalarBarActor()->SetTitle("Velocity Magnitude (m/s)");
            polydata->GetPointData()->SetActiveScalars("velocity_mag");
            lookuptable->SetRange( params.velocity_range.data() );
            break;
    case lpt::ACCELERATION:
            scalarbar->GetScalarBarActor()->SetTitle("Acceleration Magnitude (m/s^2)");
            polydata->GetPointData()->SetActiveScalars("acceleration_mag");
            lookuptable->SetRange( params.acceleration_range.data() );
            break;
    default:
        break;
    }
    polydata->Modified();
}

/***** Definition of TrajectoryHandler class *****/

TrajectoryHandler::TrajectoryHandler(vtkSmartPointer<vtkRenderer> renderer)
  : renderer(renderer), tick_count1(0),
    view_paths(false), clear_trajs(false), reset_volume_grid(false), update_volume_grid(false),
    add_vector_view(false), add_traj_view(true), remove_vector_view(false),
    remove_traj_view(false), add_cameras_view(true), remove_cameras_view(false), save_plane(false), view_mode(lpt::TRAJECTORIES)
{
}

TrajectoryHandler::~TrajectoryHandler()
{
}

void TrajectoryHandler::Execute(vtkObject* caller, unsigned long eventId, void * vtkNotUsed(callData))
{
    vtkRenderWindowInteractor* window_interactor = static_cast<vtkRenderWindowInteractor*>(caller);

    if (add_cameras_view)
        this->addCamerasToRenderer();
    else if (remove_cameras_view)
        this->removeCamerasFromRenderer();

    pair< vector < lpt::Trajectory3d*>, vector < ParticleVectors > > traj_update;
    if( render_queue.try_pop( traj_update ) ) {

        if (view_mode == lpt::TRAJECTORIES || view_mode == lpt::VECTORGRID_AND_TRAJECTORIES) {
            if (add_traj_view) {
                traj_glyphs->addToRenderer();
                add_traj_view = false;
            }
            if (remove_vector_view) {
                volume_grid->removeFromRenderer();
                remove_vector_view = false;
            }
            if (clear_trajs) {
                clearTrajectoryPaths();
                clear_trajs = false;
            }

            updateTrajectoryVisualization( traj_update );
            if ( tick_count1 >= params.traj_update_stride ) {
                window_interactor->GetRenderWindow()->Render();
                tick_count1 = 0;
            }
            ++tick_count1;
        }
    }

    if (view_mode == lpt::VECTORGRID || view_mode == lpt::VECTORGRID_AND_TRAJECTORIES ) {
            if (add_vector_view) {
                volume_grid->addToRenderer();
                add_vector_view = false;
            }
            if (remove_traj_view) {
                traj_glyphs->removeFromRenderer();
                clearTrajectoryPaths();
                remove_traj_view = false;
            }
            if (reset_volume_grid) {
                volume_grid->resetGrid();
                reset_volume_grid = false;
            }
            if (update_volume_grid) {
                update_volume_grid = false;
                volume_grid->updateGrid();
                window_interactor->GetRenderWindow()->Render();
            }
            if (save_plane) {
                save_plane = false;
                volume_grid->savePlaneData();
            }
    }

    if (clear_queue) {
        render_queue.clear();
        clear_queue = false;
    }
}

void TrajectoryHandler::clearTrajectoryPaths()
{
    for ( auto traj_iter = current_traj_list.begin(); traj_iter != current_traj_list.end(); ++traj_iter ) {
        if ( (*traj_iter)->trajvtk_ptr.get() ) {
            (*traj_iter)->trajvtk_ptr.reset();
        }
    }
    current_traj_list.clear();
}

void TrajectoryHandler::updateTrajectoryVisualization(pair< vector <lpt::Trajectory3d*>, vector <lpt::ParticleVectors > >& traj_updates)
{
    // update active trajectory renderers and create new renderers for new trajectories

    traj_glyphs->resizeGlyphArrays( traj_updates.first.size() );

    for (int id = 0; id < traj_updates.first.size(); ++id) {
        auto& current_traj = traj_updates.first[id];
        auto& current_vectors = traj_updates.second[id];

        traj_glyphs->updateParticle(id, current_vectors);

        if (view_paths) {
            if ( ! current_traj->trajvtk_ptr.get() ) {
                current_traj->trajvtk_ptr = lpt::TrajectoryPathVTK::create( renderer.GetPointer() );
            }
            lpt::TrajectoryPathVTK* trajvtk = static_cast<lpt::TrajectoryPathVTK*>( current_traj->trajvtk_ptr.get() );
            trajvtk->addNextPoint(current_vectors);
        }
    }

    traj_glyphs->setScalarsModified();

    list< lpt::Trajectory3d*> removed_traj_list( traj_updates.first.size() + current_traj_list.size() );
    list< lpt::Trajectory3d*> new_traj_list( traj_updates.first.begin(), traj_updates.first.end() );

    auto iter_end = std::set_difference(current_traj_list.begin(), current_traj_list.end(), new_traj_list.begin(), new_traj_list.end(), removed_traj_list.begin());

    for ( auto list_iter = removed_traj_list.begin(); list_iter != iter_end; ++list_iter) {
        if ( (*list_iter)->trajvtk_ptr.get() ) {
            (*list_iter)->trajvtk_ptr.reset();
        }
    }

    current_traj_list.swap( std::move( new_traj_list ) );
}

void TrajectoryHandler::addControls()
{
    void* trajhandler_void_ptr = static_cast<void*> ( this );
    cv::createButton("Clear TrajView", callbackClearTrajView, trajhandler_void_ptr, CV_PUSH_BUTTON, 0 );
    cv::createTrackbar("Update Stride", string() , &this->params.traj_update_stride, 60, 0, 0);
    cv::createButton("Reset Grid", callbackResetVolumeGrid, trajhandler_void_ptr, CV_PUSH_BUTTON, 0 );
    cv::createButton("Update Grid", callbackUpdateVolumeGrid, trajhandler_void_ptr, CV_PUSH_BUTTON, 0 );
    //cv::createTrackbar("GridUpStride", string() , &this->params.grid_update_stride, 1000, 0, 0);
    cv::createTrackbar("TrajLength", string() , &lpt::TrajectoryPathVTK::max_points, 250, 0, 0);
    cv::createButton("Show Cameras", callbackSetDisplayCameras, trajhandler_void_ptr , CV_CHECKBOX, 1 );
    cv::createButton("Save Data", callbackSavePlane, trajhandler_void_ptr , CV_PUSH_BUTTON, 0 );
    cout << "Added Visualizer Handler Controls to Window" << endl;
}

void TrajectoryHandler::setViewMode(ViewMode mode)
{
    view_mode = mode;
    switch(mode) {
    case lpt::TRAJECTORIES:
        remove_vector_view = true;
        add_traj_view = true;
        break;
    case lpt::VECTORGRID:
        remove_traj_view = true;
        add_vector_view = true;
        break;
    case lpt::VECTORGRID_AND_TRAJECTORIES:
        add_traj_view = true;
        add_vector_view = true;
    default:
        break;
    }
}

void TrajectoryHandler::addCamerasToRenderer()
{
    for (size_t c = 0; c < camerasvtk->size(); ++c)
        camerasvtk->operator[](c).addToRenderer(renderer);
    add_cameras_view = false;
}

void TrajectoryHandler::removeCamerasFromRenderer()
{
    for (size_t c = 0; c < camerasvtk->size(); ++c)
        camerasvtk->operator[](c).removeFromRenderer(renderer);
    remove_cameras_view = false;
}

/***** call back functions of TrajectoryHandler class *****/

void callbackResetVolumeGrid(int state, void* data)
{
    TrajectoryHandler* traj_handler = static_cast<TrajectoryHandler*>(data);
    traj_handler->setResetVolumeGrid();
}

void callbackUpdateVolumeGrid(int state, void* data)
{
    TrajectoryHandler* traj_handler = static_cast<TrajectoryHandler*>(data);
    traj_handler->setUpdateVolumeGrid();
}

void callbackClearTrajView(int state, void* data)
{
    TrajectoryHandler* traj_handler = static_cast<TrajectoryHandler*>(data);
    traj_handler->setClearView();
}

void callbackSetDisplayCameras(int state, void* data)
{
    TrajectoryHandler* traj_handler = static_cast<TrajectoryHandler*>(data);
    traj_handler->setDisplayCameras( (state != 0) );
}

void callbackSavePlane(int state, void* data)
{
    TrajectoryHandler* traj_handler = static_cast<TrajectoryHandler*>(data);
    traj_handler->setSavePlane(true);
}

/***** Definition of VisualizerInteractorStyle class *****/

VisualizerInteractorStyle::VisualizerInteractorStyle()
{
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

VisualizerInteractorStyle::~VisualizerInteractorStyle()
{
}

void VisualizerInteractorStyle::OnKeyPress()
{
    // Get the keypress
    vtkRenderWindowInteractor* iren = this->Interactor;
    string key = iren->GetKeySym();

    // Output the key that was pressed
    cout << "Pressed " << key << std::endl;

    //  Turn Implicit Plane on and off
    if (key == "z") {
        if ( grid->getImplicitPlane()->GetEnabled() ) {
            grid->removeImplicitPlane();
        }
        else {
            grid->addImplicitPlane();
        }
    }

    if(key == "0") {
        std::cout << "Selected Trajectory And Vector View" << std::endl;
        traj_handler->setViewMode(lpt::VECTORGRID_AND_TRAJECTORIES);
    }

    if(key == "1") {
        std::cout << "Selected Trajectory View" << std::endl;
        traj_handler->setViewMode(lpt::TRAJECTORIES);
    }

    if(key == "2") {
        std::cout << "Selected Vector Field View" << std::endl;
        traj_handler->setViewMode(lpt::VECTORGRID);
    }

    if(key == "Up") {
        if (traj_handler->getViewMode() == lpt::TRAJECTORIES || traj_handler->getViewMode() == lpt::VECTORGRID_AND_TRAJECTORIES)
            traj_handler->getViewPaths() ? traj_handler->setViewPaths(false) : traj_handler->setViewPaths(true);
        if ( traj_handler->getViewPaths() == false )
            traj_handler->setClearView();
        std::cout << "Showing trajectories" << std::endl;
    }

    if(key == "Left") {
        if ( vector_mode_iter == vector_modes.begin() )
            vector_mode_iter = --vector_modes.end();
        else
            --vector_mode_iter;
        traj_handler->setVectorMode(*vector_mode_iter);
        cout << "Setting vector mode # " << *vector_mode_iter << endl;
    }

    if(key == "Right") {
        ++vector_mode_iter;
        if ( vector_mode_iter == vector_modes.end() )
            vector_mode_iter = vector_modes.begin();
        traj_handler->setVectorMode(*vector_mode_iter);
        cout << "Setting vector mode # " << *vector_mode_iter << endl;
    }

    if(key == "c") {
        traj_handler->setClearQueue();
        std::cout << "Cleared Visualization Queue." << std::endl;
    }

    // Forward events
    vtkInteractorStyleTrackballCamera::OnKeyPress();
}

/***** Definition of Visualizer class *****/

Visualizer::Visualizer()
  : visualization_status(false), accumulate_status(false), take_measurement(false), image_measurement(false), measurement_count(0)
{
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground( params.background_color.data() );
    render_window = vtkSmartPointer<vtkRenderWindow>::New();
    render_window->AddRenderer( renderer );
    render_window->SetSize( params.window_size.data() );
    window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    window_interactor->SetRenderWindow(render_window);

    fluid_props = lpt::FluidProperties::create();
    volume_grid =  lpt::FiniteVolumeGrid::create(renderer);
    traj_glyphs = lpt::ParticlesVTK::create(renderer);
    handler = lpt::TrajectoryHandler::create(renderer);
    handler->setFiniteVolumeGrid( volume_grid.get() );
    handler->setTrajectoryGlyphs( traj_glyphs.get() );
}

Visualizer::~Visualizer()
{
    cout << "Visualizer closed" << endl;
}

void Visualizer::initialize()
{
    volume_grid->setSharedObjects(this->shared_objects);
    volume_grid->setFluidProperties(this->fluid_props);
    volume_grid->initialize();

    lpt::TrajectoryPathVTK::setLookupTable(volume_grid->getLookupTable() );

    traj_glyphs->setLookupTable( volume_grid->getLookupTable() );
    traj_glyphs->setScalarBarWidget( volume_grid->getScalarBarWidget() );

    style = vtkSmartPointer < VisualizerInteractorStyle >::New();
    style->SetDefaultRenderer(renderer);
    style->traj_handler = handler.get();
    style->grid = volume_grid.get();
    style->Modified();

    window_interactor->SetInteractorStyle(style);
    coordinates = lpt::CoordinateArrows::create(window_interactor, params.scale);
}

void Visualizer::start()
{
    window_interactor->Initialize();
    window_interactor->AddObserver( vtkCommand::TimerEvent, handler.get() );
    timer_id = window_interactor->CreateRepeatingTimer(params.timer_duration);
    renderer->ResetCamera();
    queue_manager = boost::thread(&Visualizer::manageQueue, this);
    window_interactor->Start();
}

void Visualizer::stop()
{
    queue_manager.interrupt();
    queue_manager.join();
    window_interactor->DestroyTimer(timer_id);
    render_window->Finalize();
    window_interactor->TerminateApp();
    cout << "Closing VTK window..." << endl;
}

void Visualizer::addControls()
{
    void* visualizer_void_ptr = static_cast<void*> ( this );
    cv::createButton("Run Visualization", callbackSetVisualizationStatus, visualizer_void_ptr , CV_CHECKBOX, 0 );
    cv::createButton("Accumulate", callbackSetAccumulate, visualizer_void_ptr , CV_CHECKBOX, 0 );
    cv::createTrackbar("FrameStride", string() , &this->params.stride, 60, 0, 0);
    cv::createButton("Take Measurment", callbackTakeMeasurement,visualizer_void_ptr , CV_PUSH_BUTTON, 0 );
    handler->addControls();
    cout << "Added Visualizer Controls to Window" << endl;
}

void Visualizer::manageQueue()
{
    int count = 0;
    while(true) {
        vector < pair < lpt::Trajectory3d*, vector<pair< lpt::Particle3d_Ptr, array<double, 9> > >::iterator > > traj_updates;

        traj_queue.wait_and_pop( traj_updates );

        vector < lpt::Trajectory3d*> traj_pointers( traj_updates.size() );
        vector < lpt::ParticleVectors > particle_vectors ( traj_updates.size() );

        int id = 0;
        for ( auto traj_iter = traj_updates.begin(); traj_iter != traj_updates.end(); ++traj_iter, ++id) {
            auto current_traj = traj_iter->first;
            auto object_iter = traj_iter->second;

            auto iter_next = object_iter + 1;
            auto iter_previous = object_iter - 1;

            lpt::ParticleVectors new_vectors;

            array<double, 3> Xm1, X0, X1, X2, X3, U0, U1, U2, A, dX, dPdX;

            dX[0] = volume_grid->params.cell_dimensions[0] / 1000.0;
            dX[1] = volume_grid->params.cell_dimensions[1] / 1000.0;
            dX[2] = volume_grid->params.cell_dimensions[2] / 1000.0;

            auto& frame_rate = this->shared_objects->frame_rate;

            if (this->shared_objects->KF_isOn) {
                for (int i=0; i<X1.size(); i++) {
                    X1[i] = object_iter->second.at(i);
                    U0[i] = iter_previous->second.at(i+3) / 1000.0;
                    U1[i] = object_iter->second.at(i+3) / 1000.0;
                    U2[i] = iter_next->second.at(i+3) / 1000.0;
                    A[i] = object_iter->second.at(i+6) / 1000.0;
                }
            }

            else {
                X0 = iter_previous->first->X; // previous particle, index 0
                X1 = object_iter->first->X;   // current particle, index 1
                X2 = iter_next->first->X;     // next particle, index 2

                for (int d = 0; d < X0.size(); ++d) {
                    //object_iter->second.at(d) = X1[d];
                    U1[d] = ( ( X2[d] - X0[d] ) * frame_rate ) / 2.0 / 1000.0;                      // Velocity
                    object_iter->second.at(d+3) = U1[d] * 1000.0;
					A[d] = ( X2[d] - 2*X1[d] + X0[d] ) * frame_rate * frame_rate / 1000.0;			// Acceleration
                    object_iter->second.at(d+6) = A[d] * 1000.0;
                }

                if (iter_next != ( current_traj->objects.end() - 1 ) && iter_previous != current_traj->objects.begin() ) {
                    Xm1 = (iter_previous - 1)->first->X;  // previous^2 particle, index -1
                    X3 = (iter_next + 1)->first->X;       // next^2 particle, index 3

                    for (int d = 0; d < X0.size(); ++d) {
                        U0[d] =  ( X1[d] - Xm1[d] ) * frame_rate / 2.0 / 1000.0;
                        //iter_previous->second.at(d+3) = U0[d] * 1000.0;
                        U2[d] =  ( X3[d] - X1[d] ) * frame_rate / 2.0 / 1000.0;
                        //iter_next->second.at(d+3) = U2[d] * 1000.0;
                    }
                }

            }

            dPdX = lpt::calcPressureGradiantNavierStokes(U0, U1, U2, A, *fluid_props, dX);

            for (int d = 0; d < X1.size(); ++d) {
                new_vectors[0][d] = X1[d];          // Coordinate
                new_vectors[1][d] = U1[d];          // Velocity
                new_vectors[2][d] = A[d];           // Acceleration
                new_vectors[3][d] = dPdX[d];		// Pressure gradient
            }

            if (accumulate_status) {
                volume_grid->updateAccumulators(current_traj, new_vectors);
            }

            traj_pointers[id] = current_traj;
            particle_vectors[id] = std::move(new_vectors);
        }

        if (take_measurement)
            this->takeMeasurement(particle_vectors);

        if (count >= params.stride ) {
            handler->pushToRenderQueue(std::move( std::make_pair( traj_pointers, particle_vectors ) ) );
            count = 0;
        }
        count++;
        boost::this_thread::interruption_point();
    }
}

//void Visualizer::manageQueue()
//{
//    int count = 0;
//    while(true) {
//        vector < pair < lpt::Trajectory3d*, vector<array<double, 9>>::iterator > > traj_updates;
//
//        traj_queue.wait_and_pop( traj_updates );
//
//        vector < lpt::Trajectory3d*> traj_pointers( traj_updates.size() );
//        vector < lpt::ParticleVectors > particle_vectors ( traj_updates.size() );
//
//        int id = 0;
//        for ( auto traj_iter = traj_updates.begin(); traj_iter != traj_updates.end(); ++traj_iter, ++id) {
//            auto current_traj = traj_iter->first;
//            auto state_iter = traj_iter->second;
//
//            auto next_iter = state_iter + 1;
//            auto previous_iter = state_iter - 1;
//
//			array<double, 3> X, U0, U1, U2, A, dX;
//
//			for (int i=0; i<X.size(); i++) {
//				X[i] = state_iter->at(i);
//				U0[i] = previous_iter->at(i+3);
//				U1[i] = state_iter->at(i+3);
//				U2[i] = next_iter->at(i+3);
//				A[i] = state_iter->at(i+6);
//			}
//
//			dX[0] = volume_grid->params.cell_dimensions[0] / 1000.0;
//            dX[1] = volume_grid->params.cell_dimensions[1] / 1000.0;
//            dX[2] = volume_grid->params.cell_dimensions[2] / 1000.0;
//
//			array<double,3> dPdX = lpt::calcPressureGradiantNavierStokes(U0, U1, U2, A, *fluid_props, dX);
//
//            lpt::ParticleVectors new_vectors;
//
//            auto& frame_rate = this->shared_objects->frame_rate;
//
//            for (int d = 0; d < X.size(); ++d) {
//                new_vectors[0][d] = X[d];               // Coordinate
//                new_vectors[1][d] = U1[d] / 1000.0;     // Velocity
//                new_vectors[2][d] = A[d] / 1000.0;		// Acceleration
//                new_vectors[3][d] = dPdX[d];			// Pressure gradient
//            }
//          
//            if (accumulate_status) {
//                volume_grid->updateAccumulators(current_traj, new_vectors);
//            }
//
//            traj_pointers[id] = current_traj;
//            particle_vectors[id] = std::move(new_vectors);
//        }
//
//        if (take_measurement)
//            this->takeMeasurement(particle_vectors);
//
//        if (count >= params.stride ) {
//            handler->pushToRenderQueue(std::move( std::make_pair( traj_pointers, particle_vectors ) ) );
//            count = 0;
//        }
//        count++;
//        boost::this_thread::interruption_point();
//    }
//}

void Visualizer::addTrajectoriesToQueue( list<lpt::Trajectory3d_Ptr>& active_trajs )
{
    vector< pair<lpt::Trajectory3d*, vector<pair< lpt::Particle3d_Ptr, array<double, 9> > >::iterator > > traj_updates;

    if ( ! traj_queue.full() ) {
        list<lpt::Trajectory3d_Ptr>::iterator traj_iter;
        for (traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter) {
            auto current_traj = traj_iter->get();
            if (current_traj->objects.size() > 3) {
                if (current_traj->objects.size() >= params.min_traj_size ) {
                    auto current_object_iter = current_traj->objects.end() - 3;
                    traj_updates.push_back( std::move( make_pair( current_traj, current_object_iter ) ) );
                }
            }
        }
    }
    this->traj_queue.push( std::move( traj_updates ) );
}

//void Visualizer::addTrajectoriesToQueue( list<lpt::Trajectory3d_Ptr>& active_trajs )
//{
//    vector< pair<lpt::Trajectory3d*, vector<array<double, 9>>::iterator > > traj_updates;
//
//    if ( ! traj_queue.full() ) {
//        list<lpt::Trajectory3d_Ptr>::iterator traj_iter;
//        for (traj_iter = active_trajs.begin(); traj_iter != active_trajs.end(); ++traj_iter) {
//            auto current_traj = traj_iter->get();
//            if (current_traj->obj_state.size() > 3) {
//				auto current_state_iter = current_traj->obj_state.end() - 2;
//				traj_updates.push_back( std::move( make_pair( current_traj, current_state_iter ) ) );
//			}
//        }
//    }
//    this->traj_queue.push( std::move( traj_updates ) );
//}
	
void Visualizer::takeMeasurement(const vector < lpt::ParticleVectors >& particle_vectors)
{
    if ( particle_vectors.size() > 1 ) {
        int num_measured = 0;

        for (int i = 0; i < particle_vectors.size() - 1; ++i) {
            for (int j = i + 1; j < particle_vectors.size(); ++j) {
                auto& X1 = particle_vectors[i][0];
                auto& X2 = particle_vectors[j][0];
                double dist =
                    sqrt(
                    (X1[0] - X2[0]) * (X1[0] - X2[0])
                    + (X1[1] - X2[1]) * (X1[1] - X2[1])
                    + (X1[2] - X2[2]) * (X1[2] - X2[2])
                    );
                if ( dist > 15 && dist < 25) {
                    distance_accumulator(dist);
                    num_measured++;
                }
            }
        }

        if ( measurement_count == 0 ) {
            position_accumulators.resize( particle_vectors.size() );
            for (int i = 0; i < particle_vectors.size(); ++i) {
                auto& X1 = particle_vectors[i][0];
                position_accumulators[i][0]( X1[0] );
                position_accumulators[i][1]( X1[1] );
                position_accumulators[i][2]( X1[2] );
            }
        } else if ( particle_vectors.size() == position_accumulators.size() ) {
            for (int i = 0; i < particle_vectors.size(); ++i) {
                double x = extract_result<tag::mean>(position_accumulators[i][0]);
                double y = extract_result<tag::mean>(position_accumulators[i][1]);
                double z = extract_result<tag::mean>(position_accumulators[i][2]);
                for (int j = 0; j < particle_vectors.size(); ++j) {
                    auto& X1 = particle_vectors[j][0];
                    if ( abs(X1[0] - x) < 5.0 && abs(X1[1] - y) < 5.0 && abs(X1[2] - z) < 5.0 ) {
                        position_accumulators[i][0](X1[0]);
                        position_accumulators[i][1](X1[1]);
                        position_accumulators[i][2](X1[2]);
                        break;
                    }
                }
            }
        }

        cout << "Measurement count = " << measurement_count << " num measured = " << num_measured << endl;
        if ( measurement_count >= 100 ) {
            ofstream fout;
            fout.open("measurement.txt", ios_base::app);
            take_measurement = false;
            double mean = extract_result<tag::mean>(distance_accumulator);
            double stdev = extract_result<tag::variance>(distance_accumulator);
            stdev = stdev != 0 ? sqrt(stdev) : 0;
            double max = extract_result<tag::max>(distance_accumulator);
            double min = extract_result<tag::min>(distance_accumulator);
            size_t num = extract_result<tag::count>(distance_accumulator);
            cout << "\n\nMeasurement Complete: Mean distance = " << mean << " +- " << stdev << " max, min = " << max << ", " << min << endl;
            fout << this->shared_objects->frame_rate << "\t" << mean << "\t" << stdev << "\t" << max << "\t" << min << "\t" << num << "\t\n";
            fout.close();
            fout.open("positions.txt", ios_base::app);
            lpt::boost_accumulator uncertainty_accumulator;

            for (int i = 0; i < position_accumulators.size(); ++i) {
                fout << this->shared_objects->frame_rate << "\t" << i << "\t";
                fout << extract_result<tag::count>(position_accumulators[i][0]) << "\t";
                double uncertainty = 0;
                for (int d = 0; d < 3; ++d) {
                    double x_mean = extract_result<tag::mean>(position_accumulators[i][d]);
                    double stdev_x = extract_result<tag::variance>(position_accumulators[i][d]);
                    stdev_x = stdev_x != 0 ? sqrt(stdev_x) : 0;
                    //double max_x = extract_result<tag::max>(position_accumulators[i][d]);
                    //double min_x = extract_result<tag::min>(position_accumulators[i][d]);
                    fout << x_mean << "\t" << stdev_x << "\t";
                    uncertainty += stdev_x * stdev_x;
                }
                fout << endl;
                uncertainty = uncertainty != 0 ? sqrt(uncertainty) : 0;
                uncertainty_accumulator(uncertainty);
            }
            double mean_uncertainty = extract_result<tag::mean>(uncertainty_accumulator);
            size_t num_sample = extract_result<tag::count>(uncertainty_accumulator);
            cout << "3D Position uncertainty = " << mean_uncertainty << endl;

            fout << mean_uncertainty << "\t" << num_sample << endl;
            fout.close();
            measurement_count = -1;
            distance_accumulator = lpt::boost_accumulator();
            position_accumulators.clear();
        }

        measurement_count++;

    } else if (particle_vectors.size() == 1) {
        auto& V = particle_vectors[0][1];
        auto& A = particle_vectors[0][2];
        double velocity_magnitude =
                sqrt(
                  (V[0] * V[0])
                + (V[1] * V[1])
                + (V[2] * V[2])
                );
        double acceleration_magnitude =
                sqrt(
                  (A[0] * A[0])
                + (A[1] * A[1])
                + (A[2] * A[2])
                );
        scalar_accumulator[0](velocity_magnitude);
        scalar_accumulator[1](acceleration_magnitude);

        if (measurement_count >= 10 * this->shared_objects->frame_rate) {
            ofstream fout;
            fout.open("velocity_acceleration.txt", ios_base::app);
            take_measurement = false;
            double mean = extract_result<tag::mean>(scalar_accumulator[0]);
            double stdev = extract_result<tag::variance>(scalar_accumulator[0]);
            stdev = stdev != 0 ? sqrt(stdev) : 0;
            double max = extract_result<tag::max>(scalar_accumulator[0]);
            double min = extract_result<tag::min>(scalar_accumulator[0]);
            size_t num = extract_result<tag::count>(scalar_accumulator[0]);
            cout << "\n\nMeasurement Complete: Mean velocity = " << mean << " +- " << stdev << " max, min = " << max << ", " << min << endl;
            fout << mean << "\t" << stdev << "\t" << max << "\t" << min << "\t" << num << "\t\t";

            mean = extract_result<tag::mean>(scalar_accumulator[1]);
            stdev = extract_result<tag::variance>(scalar_accumulator[1]);
            stdev = stdev != 0 ? sqrt(stdev) : 0;
            max = extract_result<tag::max>(scalar_accumulator[1]);
            min = extract_result<tag::min>(scalar_accumulator[1]);
            num = extract_result<tag::count>(scalar_accumulator[1]);
            cout << "Measurement Complete: Mean acceleration = " << mean << " +- " << stdev << " max, min = " << max << ", " << min << endl;
            fout << mean << "\t" << stdev << "\t" << max << "\t" << min << "\t" << num << endl;

            measurement_count = -1;
            scalar_accumulator[0] = lpt::boost_accumulator();
            scalar_accumulator[1] = lpt::boost_accumulator();

            fout.close();
        }
        measurement_count++;
    }
}

void Visualizer::accumulateCentroidDetectionUncertainty( vector<lpt::Match::Ptr>& matches )
{
    if (matches.size() > 0) {
        if ( image_measurement_count == 0 ) {

            centroid_uncertainty_accumulators.resize( matches.size() );

            for (int match_id = 0; match_id < matches.size(); ++match_id) {
                centroid_uncertainty_accumulators[match_id].resize( matches[match_id]->particles.size() );

                for (int j = 0; j < matches[match_id]->particles.size(); ++j ) {
                    size_t cam_id = matches[match_id]->particles[j].second;
                    lpt::ParticleImage::Ptr particle = matches[match_id]->particles[j].first;

                    centroid_uncertainty_accumulators[match_id][j][0]( particle->x );
                    centroid_uncertainty_accumulators[match_id][j][1]( particle->y );
                }
            }
        } else if ( matches.size() == centroid_uncertainty_accumulators.size() ) {
            for (int match_id = 0; match_id < matches.size(); ++match_id) {
                for (int j = 0; j < matches[match_id]->particles.size(); ++j ) {
                    size_t cam_id = matches[match_id]->particles[j].second;
                    lpt::ParticleImage::Ptr particle = matches[match_id]->particles[j].first;

                    for (int p = 0; p < matches.size(); ++p) {
                        double x = extract_result<tag::mean>(centroid_uncertainty_accumulators[p][j][0]);
                        double y = extract_result<tag::mean>(centroid_uncertainty_accumulators[p][j][1]);
                        if ( abs(particle->x - x) < 4.0 && abs(particle->y - y) < 4.0 ) {
                            centroid_uncertainty_accumulators[p][j][0]( particle->x );
                            centroid_uncertainty_accumulators[p][j][1]( particle->y );
                            break;
                        }
                    }
                }
            }
        }

        cout << "Measurement count = " << image_measurement_count << " num matches = " << matches.size() << endl;
        if ( image_measurement_count >= 100 ) {
            ofstream fout;
            fout.open("imagemeasurement.txt", ios_base::app);
            image_measurement = false;

            lpt::boost_accumulator uncertainty_accumulator;

            for (int i = 0; i < centroid_uncertainty_accumulators.size(); ++i) {
                fout << this->shared_objects->frame_rate << "\t" << i << "\t";
                for (int c = 0; c < centroid_uncertainty_accumulators[i].size(); c++) {
                    size_t cam_id = matches[i]->particles[c].second;

                    fout << cam_id << "\t" << extract_result<tag::count>(centroid_uncertainty_accumulators[i][c][0]) << "\t";
                    double pixel_uncertainty = 0;
                    for (int d = 0; d < 2; ++d) {
                        double x_mean = extract_result<tag::mean>(centroid_uncertainty_accumulators[i][c][d]);
                        double stdev_x = extract_result<tag::variance>(centroid_uncertainty_accumulators[i][c][d]);
                        stdev_x = stdev_x != 0 ? sqrt(stdev_x) : 0;
                        //double max_x = extract_result<tag::max>(position_accumulators[i][d]);
                        //double min_x = extract_result<tag::min>(position_accumulators[i][d]);
                        fout << x_mean << "\t" << stdev_x << "\t";
                        pixel_uncertainty += stdev_x * stdev_x;
                    }
                    pixel_uncertainty = pixel_uncertainty != 0 ? sqrt(pixel_uncertainty) : 0;
                    uncertainty_accumulator(pixel_uncertainty);
                }
                fout << endl;
            }
            double mean_uncertainty = extract_result<tag::mean>(uncertainty_accumulator);
            size_t num_sample = extract_result<tag::count>(uncertainty_accumulator);
            cout << "2D centroid uncertainty = " << mean_uncertainty << endl;

            fout << mean_uncertainty << "\t" << num_sample << endl;
            fout.close();
            image_measurement_count = -1;
            centroid_uncertainty_accumulators.clear();
        }
        image_measurement_count++;
    }
}

void Visualizer::setCameras(vector<Camera> &cameras)
{
    for (size_t c = 0; c < cameras.size(); ++c)
        camerasvtk.emplace_back(cameras[c]);
    handler->setCamerasVTK(camerasvtk);
    centroid_uncertainty_accumulators.resize(cameras.size());
}

/***** call back functions of Visualizer class *****/

void callbackSetVisualizationStatus(int state, void* data)
{
    Visualizer* visualizer = static_cast<Visualizer*>(data);
    visualizer->setVisualizationStatus( (state != 0) );
}

void callbackSetStride(int state, void* data)
{
    Visualizer* visualizer = static_cast<Visualizer*>(data);
}

void callbackSetAccumulate(int state, void* data)
{
    Visualizer* visualizer = static_cast<Visualizer*>(data);
    visualizer->setAccumulate( (state != 0) );
}

void callbackTakeMeasurement(int state, void*data)
{
    Visualizer* visualizer = static_cast<Visualizer*>(data);
    visualizer->setTakeMeasurement(true);
}


/***** Definition of HistogramVTK class *****/

HistogramVTK::HistogramVTK() : number_of_bins(100)
{
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

HistogramVTK::~HistogramVTK()
{
}

void HistogramVTK::Execute(vtkObject *caller, unsigned long eventId, void *callData)
{
    vtkContextView* view_caller = static_cast<vtkContextView*>(caller);
    if (eventId == vtkCommand::UserEvent) {
        view_caller->Update();
        view_caller->Render();
    }
}

void HistogramVTK::setBins(int num_bins, double bin_range[])
{
    number_of_bins = num_bins;
    range[0] = bin_range[0];
    range[1] = bin_range[1];
    bins.resize(number_of_bins, 0);
    count_data.resize(number_of_bins, 0);
    bin_size = (range[1] - range[0]) / static_cast<double>(number_of_bins);
    for (int i = 0; i < bins.size(); ++i)
        bins[i] = bin_size * (i+1);
}

void HistogramVTK::setAxisColor(int color[])
{
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

void HistogramVTK::setTextStrings(string &title, string &x_title, string &y_title)
{
    histogram->SetTitle( title.c_str() );
    histogram->GetAxis(1)->SetTitle( x_title.c_str() );
    histogram->GetAxis(0)->SetTitle( y_title.c_str() );
}

void HistogramVTK::updateData(vector<double> &data)
{
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

int HistogramVTK::getBinID(double value)
{
    int id;
    for ( id = 0; id < bins.size(); id++) {
        if ( bins[id] > value )
            break;
    }
    return id;
}

void HistogramVTK::addToRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
    renderer->AddActor(chart_actor.GetPointer());
    chart_scene->SetRenderer(renderer);
}

void HistogramVTK::startChartWindow()
{
    if (! view->GetInteractor()->GetEnabled() ) {
        view->GetInteractor()->Initialize();
        view->GetInteractor()->Start();
    }
}

void HistogramVTK::stopChartWindow()
{
    if (view->GetInteractor()->GetEnabled()) {
        cout << " fina a way to close the window without terminating the app " << endl;
    }
}

void HistogramVTK::removeFromRenderer(vtkSmartPointer<vtkRenderer> renderer)
{
    renderer->RemoveActor(chart_actor.GetPointer());
}

/***** Definition of Streamlines class *****/

StreamLines::StreamLines(vtkAlgorithmOutput *output_port)
{
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

/***** Definition of PickDim class *****/

PickDim::PickDim()
{
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

PickDim::~PickDim()
{
}

void PickDim::Execute(vtkObject *caller, unsigned long eventId, void *)
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

/*
KalmanFilter::KalmanFilter(std::shared_ptr<SharedObjects> new_shared_objects)
  : shared_objects(new_shared_objects)
{
    F = Eigen::Matrix<double, 6, 6>::Identity();
    H = Eigen::Matrix<double, 6, 6>::Identity();
    P = 100 * Eigen::Matrix<double, 6, 6>::Identity();
    for (i=0; i<3; i++)
        F(i,i+3) = 1.0 / shared_objects->frame_rate;
    //Q,R
}

KalmanFilter::~KalmanFilter()
{
}

void KalmanFilter::filter()
{
    Eigen::Matrix<double, 6, 1> temp_s, residual;
    Eigen::Matrix<double, 6, 6> S, K;
    auto I = Eigen::Matrix<double, 6, 6>::Identity();

    //predict
    temp_s = F * s;
    P = F * P * F.transpose() + Q;

    //update
    residual = z - H * temp_s;
    S = H * P * H.transpose() + R;
    K = P * H.transpose() * S.fullPivLu().solve(I);
    s = temp_s + K * residual;
    P = (I - K * H) * P;
}

Eigen::Matrix<double, 6, 1> KalmanFilter::GetState() const
{
    return s;
}

void KalmanFilter::setState(Eigen::Matrix<double, 6, 1> new_state)
{
    s = new_state;
}

void KalmanFilter::setObservation(Eigen::Matrix<double, 6, 1> new_observation)
{
    z = new_observation;
}

void KalmanFilter::setSharedObjects(std::shared_ptr<SharedObjects> new_shared_objects)
{
    shared_objects = new_shared_objects;
} */

} /*NAMESPACE_LPT*/
