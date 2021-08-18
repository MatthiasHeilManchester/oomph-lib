//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1307 $
//LIC//
//LIC// $LastChangedDate: 2018-01-18 11:30:14 +0000 (Thu, 18 Jan 2018) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

#include<fenv.h>

//Generic routines
#include "generic.h"

// Poisson
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"
 
// Get the mesh
#include "meshes/tetgen_mesh.h" 
#include "meshes/refineable_tetgen_mesh.h"
#include "meshes/gmsh_tet_mesh.h"

// Get the faceted surfaces
#include "tetmesh_faceted_surfaces.h"

// Tetgen or Gmsh
//#define DO_TETGEN


using namespace std;

using namespace oomph;



///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////



//=============================================================
/// Namespace for problem parameters
//=============================================================
namespace Global_Parameters
{

 /// Pointer to Gmsh parameters
 GmshParameters* Gmsh_parameters_pt=0;

 /// Counter for generated meshes
 unsigned Count=0;

 /// Initial element volume
 double Initial_element_volume=1.0;

 /// Typical variation in element volume in initial mesh
 double Typical_var_in_element_volume_in_percent=0.0;

 /// Number of refinements for uniform refinement test
 unsigned Nrefinement=3;

 /// Number of refinements for uniform unrefinement test
 unsigned Nunrefinement=3;

 /// \short Factor by which element volume is (attempted to be) increased
 /// during uniform unrefiment
 double Element_volume_increase_factor_during_unrefinement=8.0;
 
 /// \short Factor by which element volume is (attempted to be) reduced
 /// during uniform refinement
 double Element_volume_reduction_factor_during_refinement=1.0/8.0;
 
 /// (Half-)width of the box
 double Box_half_width = 1.5;

 /// (Half)height of the box
 double Box_half_length = 1.0;


}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////




//====================================================================
/// Demo class solves Poisson problem using Gmsh mesh
//====================================================================
template<class ELEMENT> 
class TetmeshPoissonProblem : public Problem
{

public:

 /// Constructor
 TetmeshPoissonProblem();
  
 /// Destructor (empty)
 ~TetmeshPoissonProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }

      
 /// Actions before adapt (empty)
 void actions_before_adapt()
  {
   Global_Parameters::Gmsh_parameters_pt->target_size_file_name()=
    "latest_initial_target_sizes_refinement"+
    to_string(Global_Parameters::Count)+".dat";
   Global_Parameters::Gmsh_parameters_pt->geo_and_msh_file_stem()=
    ".geo_file_refinement"+to_string(Global_Parameters::Count);
   
   oomph_info << "Before adapt: target size filename: "
              << Global_Parameters::Gmsh_parameters_pt->target_size_file_name()
              << std::endl;
  }

 /// Totally new mesh; build elements and apply boundary conditions
 void actions_after_adapt()
  {
   // Complete problem setup
   complete_problem_setup();

   char filename[100];
   sprintf(filename,"mesh%i.dat",Global_Parameters::Count+1);
   ofstream mesh_file;
   mesh_file.open(filename);
   Bulk_mesh_pt->output(mesh_file,2);
   mesh_file.close();
   
   oomph_info << "After adapt: written mesh to: "
              << filename
              << std::endl;
   
   Global_Parameters::Count++;
   
  }
 
 /// Update the problem specs before solve: (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}
 
 /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);

private:
 
 /// Apply BCs and make elements functional
 void complete_problem_setup();

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();
 
#ifdef DO_TETGEN

 /// Bulk mesh
 RefineableTetgenMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Bulk mesh
 RefineableGmshTetMesh<ELEMENT>* Bulk_mesh_pt;

#endif

 /// Storage for the outer boundary object
 TetMeshFacetedClosedSurface* Outer_boundary_pt;

 /// Inner boundary
 Vector<TetMeshFacetedSurface*> Inner_boundary_pt;

 /// First boundary ID for outer boundary
 unsigned First_boundary_id_for_outer_boundary;
 
 /// First boundary ID for hollow cube
 unsigned First_hollow_cube_boundary_id;

 /// First boundary ID for cube region
 unsigned First_cube_region_boundary_id;

 /// Region ID for cube region
 unsigned Cube_region_id;

 /// Boundary ID for rectangular facet in main volume
 unsigned Internal_rectangle_boundary_id;

 /// Boundary ID for rectangular facet embedded in region
 unsigned Internal_embedded_rectangle_boundary_id;

 /// Sanity check: Exact bounded volume
 double Exact_bounded_volume;

 /// Sanity check: Exact volume of hole
 double Exact_volume_of_hole;

 /// Sanity check: Exact volume of internal region
 double Exact_volume_of_internal_region;

};



//========================================================================
/// Constructor
//========================================================================
template<class ELEMENT>
TetmeshPoissonProblem<ELEMENT>::TetmeshPoissonProblem()
{ 

 // OUTER BOUNDARY
 //===============

 // Start boundary IDs for outer boundary from some crazy offset
 // (just for testing). By default the one-based boundary IDs go from
 // 1 to 6; let's start from 1001.
 unsigned outer_boundary_id_offset=1000;

 //Make the outer boundary object
 Outer_boundary_pt = new CubicTetMeshFacetedSurface(
  Global_Parameters::Box_half_width,
  Global_Parameters::Box_half_length,
  outer_boundary_id_offset);

 // Look, we can visualise the faceted surface!
 Outer_boundary_pt->output("outer_faceted_surface.dat");

 // First oomph-lib (zero-based!) boundary ID for outer boundary
 First_boundary_id_for_outer_boundary=outer_boundary_id_offset;
 
 // For sanity check:
 Exact_bounded_volume=
  2.0*Global_Parameters::Box_half_width*
  2.0*Global_Parameters::Box_half_width*
  2.0*Global_Parameters::Box_half_length;

 
 // INTERNAL BOUNDARIES
 //====================

 // A hollow cube
 //--------------

 // Size
 double cube_half_width=0.2*Global_Parameters::Box_half_width;

 // Offset of boundary ID. By default the one-based boundary IDs go from
 // 1 to 6; let's start from 1454.
 unsigned inner_cube_boundary_id_offset=1453;

 // Create the faceted surface
 CubicTetMeshFacetedSurface* hollow_cube_pt = 
  new CubicTetMeshFacetedSurface(
   cube_half_width,
   cube_half_width,
   inner_cube_boundary_id_offset);

 // Add it
 Inner_boundary_pt.push_back(hollow_cube_pt);
 
 // Look, we can visualise the faceted surface!
 hollow_cube_pt->output("hollow_cube_faceted_surface.dat");
 
 // Zero-based oomph-lib enumeration of boundaries:
 First_hollow_cube_boundary_id=inner_cube_boundary_id_offset;

 // For sanity check:
 Exact_volume_of_hole=
  2.0*cube_half_width*
  2.0*cube_half_width*
  2.0*cube_half_width;


 // A cube inside the main domain, used to indicate a region
 //---------------------------------------------------------

 // Size
 double cube_region_half_width=0.2*Global_Parameters::Box_half_width;

 // Cube is offset by this vector
 Vector<double> region_offset(3);
 region_offset[0]=0.5*Global_Parameters::Box_half_width;
 region_offset[1]=0.5*Global_Parameters::Box_half_width;
 region_offset[2]=0.3*Global_Parameters::Box_half_width;

 // Offset of boundary ID. By default the one-based boundary IDs go from
 // 1 to 6; let's start from 2454.
 unsigned cube_region_boundary_id_offset=2453;

 // This closed faceted surface defines the following (one-based) region
 unsigned one_based_region_id=10;
 Cube_region_id=one_based_region_id-1; 

 // Create the faceted surface
 CubicTetMeshFacetedSurface* region_cube_pt = 
  new CubicTetMeshFacetedSurface(
   one_based_region_id,
   cube_region_half_width,
   cube_region_half_width,
   region_offset,
   cube_region_boundary_id_offset);

 // Add it
 Inner_boundary_pt.push_back(region_cube_pt);
 
 // Look, we can visualise the faceted surface!
 region_cube_pt->output("cube_region_faceted_surface.dat");
 
 // Zero-based oomph-lib enumeration of boundaries:
 First_cube_region_boundary_id=cube_region_boundary_id_offset;
 
 // For sanity check:
 Exact_volume_of_internal_region=
  2.0*cube_region_half_width*
  2.0*cube_region_half_width*
  2.0*cube_region_half_width;


 // Planar rectangular facet in main volume defines internal boundary
 //------------------------------------------------------------------

 // Dimensions:
 double half_x_width =0.2*Global_Parameters::Box_half_width;
 double half_y_length=0.3*Global_Parameters::Box_half_width;
 
 
 // Offset 
 Vector<double> offset(3);
 offset[0]= 0.5*Global_Parameters::Box_half_width;
 offset[1]= 0.5*Global_Parameters::Box_half_width;
 offset[2]=-0.3*Global_Parameters::Box_half_width;
 
 
 // One-based boundary ID
 unsigned one_based_boundary_id=143;
 
 // ...and its oomph-lib counterpart
 Internal_rectangle_boundary_id=one_based_boundary_id-1;

 // Build it
 RectangularTetMeshFacetedSurface* rectanglar_facet_pt=
  new RectangularTetMeshFacetedSurface(half_x_width, 
                                       half_y_length,
                                       offset,
                                       one_based_boundary_id);

 // Add it
 Inner_boundary_pt.push_back(rectanglar_facet_pt);
 
 // Look, we can visualise the faceted surface!
 rectanglar_facet_pt->output("rectangular_faceted_surface.dat");



 // Planar rectangular facet in region defines internal boundary
 //-------------------------------------------------------------

 // Dimensions:
 double embedded_half_x_width =0.1*Global_Parameters::Box_half_width;
 double embedded_half_y_length=0.1*Global_Parameters::Box_half_width;
 
  
 // One-based boundary ID
 one_based_boundary_id=243;
 
 // ...and its oomph-lib counterpart
 Internal_embedded_rectangle_boundary_id=one_based_boundary_id-1;

 // Build it
 RectangularTetMeshFacetedSurface* embedded_rectanglar_facet_pt=
  new RectangularTetMeshFacetedSurface(embedded_half_x_width, 
                                       embedded_half_y_length,
                                       region_offset,
                                       one_based_boundary_id);

 // Add it
 Inner_boundary_pt.push_back(embedded_rectanglar_facet_pt);
 
 // Look, we can visualise the faceted surface!
 embedded_rectanglar_facet_pt->output("embedded_rectangular_faceted_surface.dat");


 // hierher
 // // Suppress generation of internal faceted surfaces for target volume check
 // if ((CommandLineArgs::command_line_flag_has_been_set(
 //      "--desired_target_volume"))||
 //     (CommandLineArgs::command_line_flag_has_been_set(
 //      "--check_uniform_refinement"))||
 //     (CommandLineArgs::command_line_flag_has_been_set(
 //      "--check_uniform_unrefinement")))
 //  {
 //oomph_info << "No inner regions\n";
 //  Inner_boundary_pt.resize(0);
   //}
 

 
 // Suppress generation of internal faceted surfaces for target volume check
 if (CommandLineArgs::command_line_flag_has_been_set(
      "--suppress_inner_boundaries"))
  { 
   oomph_info << "No inner regions\n";
   Inner_boundary_pt.resize(0);
  }

 // Build the mesh
 //--------------- 

 // Specify how to call gmsh from the command line
 std::string gmsh_command_line_invocation="/usr/local/bin/gmsh"; ///home/mheil/gmesh/bin/bin/gmsh";

 // Setup parameters for gmsh
 Global_Parameters::Gmsh_parameters_pt=
  new GmshParameters(Outer_boundary_pt,
                     gmsh_command_line_invocation);

 // Element volume
 Global_Parameters::Gmsh_parameters_pt->element_volume()=
  Global_Parameters::Initial_element_volume;

// hierher
 // // Volume ratio during unrefinement
 // if (CommandLineArgs::command_line_flag_has_been_set(
 //      "--unrefinement_volume_increase_factor"))
 //  {
 //   Global_Parameters::Gmsh_parameters_pt->
 //    element_volume_increase_factor_during_unrefinement()=
 //    Global_Parameters::Element_volume_increase_factor_during_unrefinement;
 //  }

 // // Volume ratio during refinement
 // if (CommandLineArgs::command_line_flag_has_been_set(
 //      "--refinement_volume_decrease_factor"))
 //  {
 //   oomph_info << "hierher command line changes refinement factor from " 
 //              <<  Global_Parameters::Gmsh_parameters_pt->
 //    element_volume_reduction_factor_during_refinement();
   
 //   Global_Parameters::Gmsh_parameters_pt->
 //    element_volume_reduction_factor_during_refinement()=
 //    Global_Parameters::Element_volume_reduction_factor_during_refinement;

 //   oomph_info << "to " 
 //              <<  Global_Parameters::Gmsh_parameters_pt->
 //    element_volume_reduction_factor_during_refinement();
   

 //  }

 
 // Specify inner boundaries
 Global_Parameters::Gmsh_parameters_pt->
  internal_surface_pt()=Inner_boundary_pt;

 // Filename for file in which target element size is stored
 // (for disk-based operation of gmsh)
 Global_Parameters::Gmsh_parameters_pt->
  stem_for_filename_gmsh_size_transfer()=
  "target_size_on_grid";
 Global_Parameters::Gmsh_parameters_pt->
  counter_for_filename_gmsh_size_transfer()=0;

 // Problem is linear so we don't need to transfer the solution to the
 // new mesh; we keep it on for self-test purposes...
 // hierher
 Global_Parameters::Gmsh_parameters_pt->disable_projection();

 // Redirect gmsh on-screen output
 Global_Parameters::Gmsh_parameters_pt->gmsh_onscreen_output_file_name()=
  "RESLT/gmsh_on_screen_output.dat";


 // Not needed, of course, but here to test out the handling
 // of timesteppers
 add_time_stepper_pt(new Steady<1>);

 double target_size=0.0;

#ifdef DO_TETGEN

 // And now build it...
 Bulk_mesh_pt =
  new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
                                    Inner_boundary_pt,
                                    Global_Parameters::Initial_element_volume,
                                    this->time_stepper_pt());

 // Problem is linear so we don't need to transfer the solution to the
 // new mesh; we keep it on for self-test purposes...
 //Bulk_mesh_pt->disable_projection();

#else

 // And now build it...
 bool use_mesh_grading_from_file=false; // hierher true; 
 Global_Parameters::Gmsh_parameters_pt->target_size_file_name()=
  "my_initial_target_sizes.dat";
 Global_Parameters::Gmsh_parameters_pt->geo_and_msh_file_stem()=
  ".geo_file_initial";



 // if (CommandLineArgs::command_line_flag_has_been_set(
 //      "--specify_gmsh_target_volume_via_file"))
 //  {   
 //   target_size=pow(
 //    Global_Parameters::Gmsh_parameters_pt->target_volume_multiplier()*
 //    Global_Parameters::Initial_element_volume,1.0/3.0);
 //     ofstream target_sizes_file;
 //     oomph_info << "target size: " << target_size << std::endl;
 //     target_sizes_file.open("my_initial_target_sizes.dat");
 //     unsigned n=20;
 //     double dx=20.0/double(n);
 //     double typical_var_in_percent=
 //      Global_Parameters::Typical_var_in_element_volume_in_percent;
 //     oomph_info << "typical_var_in_percent: " << typical_var_in_percent 
 //                << std::endl;
 //     target_sizes_file 
 //      << "-10.0 -10.0 -10.0\n"
 //      << dx << " " << dx << " "<< dx << "\n" 
 //      << n+1 << " " << n+1 << " "<< n+1 << "\n"; 
 //     for (unsigned i=0;i<=n;i++)
 //      {
 //       for (unsigned j=0;j<=n;j++)
 //        {     
 //         for (unsigned k=0;k<=n;k++)
 //          {
 //           // this is between -1 and +1
 //           double var = (-1.0+2.0*((double)rand() / RAND_MAX));
 //           target_sizes_file << 
 //            target_size*(1.0+typical_var_in_percent/100.0*var) 
 //                             << " ";
 //          }
 //        }
 //      }
 //     target_sizes_file.close();
 //  }


 Bulk_mesh_pt = new RefineableGmshTetMesh<ELEMENT>(
  Global_Parameters::Gmsh_parameters_pt,
  use_mesh_grading_from_file,
  this->time_stepper_pt());


#endif
 
 // hierher
 Bulk_mesh_pt->min_element_size()=0.0;

 // Add sub-mesh
 Problem::mesh_pt()=Bulk_mesh_pt;


 // Doc average element size; useful information to 
 // relate actual element sizes to the input. Gmsh 
 // tends to make finer meshes.
 unsigned nel=Bulk_mesh_pt->nelement();
 double av_el_size=Exact_bounded_volume/double(nel);
 double rms_dev=0.0;
 for (unsigned e=0;e<nel;e++)
  {
   rms_dev+=pow(Bulk_mesh_pt->finite_element_pt(e)->size()-av_el_size,2);
  }
 oomph_info << "Av. el. size: "
            << av_el_size << " " << " nel: " << nel
            << " and standard deviation : " << sqrt(rms_dev/double(nel))
            << " for target " << Global_Parameters::Initial_element_volume
            << " target size parameter for gmsh: " << target_size 
            << std::endl;

 Bulk_mesh_pt->output("mesh0.dat",2);


 // Do uniform mesh refinement tests?
 if (CommandLineArgs::command_line_flag_has_been_set(
      "--check_uniform_refinement"))
  {
   for (unsigned i=0;i<Global_Parameters::Nrefinement;i++)
    {

     Global_Parameters::Gmsh_parameters_pt->target_size_file_name()=
      "latest_initial_target_sizes_refinement"+
      to_string(Global_Parameters::Count)+".dat";
     Global_Parameters::Gmsh_parameters_pt->geo_and_msh_file_stem()=
      ".geo_file_refinement"+to_string(Global_Parameters::Count);


     oomph_info << "shite: doing refinement. count  = " 
                << Global_Parameters::Count << "\n";
     Bulk_mesh_pt->refine_uniformly();

     // Doc average element size; useful information to 
     // relate actual element sizes to the input. Gmsh 
     // tends to make finer meshes.
     unsigned nel=Bulk_mesh_pt->nelement();
     double av_el_size=Exact_bounded_volume/double(nel);
     double rms_dev=0.0;
     for (unsigned e=0;e<nel;e++)
      {
       rms_dev+=pow(Bulk_mesh_pt->finite_element_pt(e)->size()-av_el_size,2);
      }
     oomph_info << "Av. el. size: "
                << av_el_size << " " << " nel: " << nel
                << " and standard deviation : " << sqrt(rms_dev/double(nel))
                << " for target " << Global_Parameters::Initial_element_volume
                << std::endl;
     
     char filename[100];
     sprintf(filename,"mesh%i.dat",Global_Parameters::Count+1);
     ofstream mesh_file;
     mesh_file.open(filename);
     Bulk_mesh_pt->output(mesh_file,2);
     mesh_file.close();
     
     oomph_info << "After uniform refiment: written mesh to: "
                << filename
                << std::endl;

     Global_Parameters::Count++;

    }
  }
 if (CommandLineArgs::command_line_flag_has_been_set(
      "--check_uniform_unrefinement"))
  {
   for (unsigned i=0;i<Global_Parameters::Nunrefinement;i++)
    {
     Global_Parameters::Gmsh_parameters_pt->target_size_file_name()=
      "latest_initial_target_sizes_refinement"+
      to_string(Global_Parameters::Count)+".dat";
     Global_Parameters::Gmsh_parameters_pt->geo_and_msh_file_stem()=
      ".geo_file_refinement"+to_string(Global_Parameters::Count);

     oomph_info << "shite: doing unrefinement. count  = " 
                << Global_Parameters::Count << "\n";
     Bulk_mesh_pt->unrefine_uniformly();

     // Doc average element size; useful information to 
     // relate actual element sizes to the input. Gmsh 
     // tends to make finer meshes.
     unsigned nel=Bulk_mesh_pt->nelement();
     double av_el_size=Exact_bounded_volume/double(nel);
     double rms_dev=0.0;
     for (unsigned e=0;e<nel;e++)
      {
       rms_dev+=pow(Bulk_mesh_pt->finite_element_pt(e)->size()-av_el_size,2);
      }
     oomph_info << "Av. el. size: "
                << av_el_size << " " << " nel: " << nel
                << " and standard deviation : " << sqrt(rms_dev/double(nel))
                << " for target " << Global_Parameters::Initial_element_volume
                << std::endl;
     
     char filename[100];
     sprintf(filename,"mesh%i.dat",Global_Parameters::Count+1);
     ofstream mesh_file;
     mesh_file.open(filename);
     Bulk_mesh_pt->output(mesh_file,2);
     mesh_file.close();

     oomph_info << "After uniform unrefiment: written mesh to: "
                << filename
                << std::endl;

     Global_Parameters::Count++;
     
    }
  }
   

 // hierher
 //exit(0);


 // Bounding box for fake error estimator -- this is where elements are
 // refined
 Vector<double> lower_left(3,0.0);
 Vector<double> upper_right(3,1.2);

 // Set error estimator for bulk mesh
 DummyErrorEstimator* error_estimator_pt= 
  new DummyErrorEstimator(lower_left, 
                          upper_right);

 
  //new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Bulk_mesh_pt->max_permitted_error()=0.0005; 
 Bulk_mesh_pt->min_permitted_error()=0.00001;
 
 // Complete problem setup
 complete_problem_setup();
 
   
#ifdef OOMPH_HAS_HYPRE

 // Create a new Hypre linear solver
 HypreSolver* hypre_linear_solver_pt = new HypreSolver;
 
 // Set the linear solver for problem
 linear_solver_pt() = hypre_linear_solver_pt;
 
 // Set some solver parameters
 hypre_linear_solver_pt->max_iter() = 100;
 hypre_linear_solver_pt->tolerance() = 1e-10;
 hypre_linear_solver_pt->amg_simple_smoother() = 1;
 hypre_linear_solver_pt->disable_doc_time();
 hypre_linear_solver_pt->enable_hypre_error_messages();
 hypre_linear_solver_pt->amg_print_level() = 0;
 hypre_linear_solver_pt->krylov_print_level() = 0;
 hypre_linear_solver_pt->hypre_method() = HypreSolver::BoomerAMG;
   
#endif

 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}


//========================================================================
/// Complete problem setup
//========================================================================
template<class ELEMENT>
void TetmeshPoissonProblem<ELEMENT>::complete_problem_setup()
{
 // Apply bcs
 apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void TetmeshPoissonProblem<ELEMENT>::apply_boundary_conditions()
{
 

 ofstream pin_file;
 pin_file.open("pinned_nodes.dat");

 // Identify boundary ids of pinned nodes
 Vector<unsigned> pinned_boundary_id;
 if (!CommandLineArgs::command_line_flag_has_been_set(
      "--suppress_inner_boundaries"))
  {    
   for (unsigned ibound=First_hollow_cube_boundary_id;
        ibound<First_hollow_cube_boundary_id+6;ibound++)
    {
     pinned_boundary_id.push_back(ibound);
    }
   pinned_boundary_id.push_back(Internal_rectangle_boundary_id);
   pinned_boundary_id.push_back(Internal_embedded_rectangle_boundary_id);
  }
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }

 // Loop over pinned boundaries
 unsigned num_pin_bnd=pinned_boundary_id.size();
 for (unsigned bnd=0;bnd<num_pin_bnd;bnd++)
  {
   unsigned ibound=pinned_boundary_id[bnd];
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   if (num_nod==0)
    {
     std::ostringstream error_message;
     error_message << "No boundary nodes on boundary " 
                   << ibound << "! Something's gone wrong!\n";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
     Vector<double> x(3);
     x[0]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(0);
     x[1]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
     x[2]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(2);

     // Zero BC on outer boundary; unit BC elsewhere
     double value=0.0;
     if (ibound>=First_hollow_cube_boundary_id) 
      {
       value=1.0;
      }
     // ...apart from rectangular internal boundary where set 
     // solution to 2.0.
     if (ibound==Internal_rectangle_boundary_id)
      {
       value=2.0;
      }
     // ...apart from embedded rectangular internal boundary where set 
     // solution to 3.0.
     if (ibound==Internal_embedded_rectangle_boundary_id)
      {    
       value=3.0;
      }
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,value);
     pin_file << x[0] << " " 
              << x[1] << " " 
              << x[2] << " " 
              << std::endl;
    }
  }
 
 pin_file.close();



} // end set bc



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TetmeshPoissonProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                                DocInfo& doc_info)
{ 


 oomph_info << "Soln, # elements, # nodes: " 
            << doc_info.number() << " "
            << Bulk_mesh_pt->nelement() << " "
            << Bulk_mesh_pt->nnode() << " "
            << std::endl;

 // Output can take a long time; don't do it if we're just running 
 // the self-test (which analyses different data)
 bool do_bulk_output=true;
 if (CommandLineArgs::command_line_flag_has_been_set("--suppress_bulk_output"))
  {
   do_bulk_output=false;
  }


 ofstream some_file;
 ofstream some_file2;
 ofstream face_some_file;
 ofstream coarse_some_file;
 char filename[100];

 // Doc mesh quality (Ratio of max. edge length to min. height,
 /// so if it's very large it's BAAAAAD)
 sprintf(filename,"%s/mesh_quality%i.dat",
         doc_info.directory().c_str(),
         doc_info.number()); 
 ofstream quality_file;
 quality_file.open(filename);
 if (do_bulk_output) Bulk_mesh_pt->assess_mesh_quality(quality_file);
 quality_file.close();

 
 // Output elements adjacent to outer boundary
 //-------------------------------------------
 sprintf(filename,"%s/elements_next_to_outer_boundary%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6;ibound++)
  {
   unsigned n_el=Bulk_mesh_pt->nboundary_element(ibound);
   for (unsigned e=0;e<n_el;e++)
    {
     if (do_bulk_output) 
      {
       Bulk_mesh_pt->boundary_element_pt(ibound,e)->
        output(some_file,nplot);
      }
    }
  }
 some_file.close();


 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_boundaries(some_file);
 some_file.close();


 // Output volumes and things
 //--------------------------
 std::ofstream volume_file;
 sprintf(filename,"%s/volumes%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 volume_file.open(filename);


 if (!CommandLineArgs::command_line_flag_has_been_set(
      "--suppress_inner_boundaries"))
  { 
   // Output bulk elements in cube region
   //------------------------------------
   double volume_in_internal_region=0.0;
   sprintf(filename,"%s/soln_in_cube_region%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   unsigned region_id=Cube_region_id;
   unsigned n_el=Bulk_mesh_pt->nregion_element(region_id);
   for (unsigned e=0;e<n_el;e++)
    {
     if (do_bulk_output) 
      {
       Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
      }
     volume_in_internal_region+=Bulk_mesh_pt->
      region_element_pt(region_id,e)->size();
    }
   some_file.close();
   
   // Check volume:
   oomph_info 
    << "Volume of region 1 [in mesh, exact, diff (%)]: "
    << volume_in_internal_region << " " 
    << Exact_volume_of_internal_region << " " 
    << abs(volume_in_internal_region-Exact_volume_of_internal_region)/
    Exact_volume_of_internal_region*100.0 << std::endl;
   
   
   // Output bulk elements in region 0
   //--------------------------------- 
   double volume_in_region0=0.0;
   sprintf(filename,"%s/soln_in_zero_region%i.dat",doc_info.directory().c_str(),
           doc_info.number());
   some_file.open(filename);
   region_id=0;
   n_el=Bulk_mesh_pt->nregion_element(region_id);
   for (unsigned e=0;e<n_el;e++)
    {
     if (do_bulk_output) 
      {
       Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
      }
     volume_in_region0+= Bulk_mesh_pt->region_element_pt(region_id,e)->size();
    }
   some_file.close();
   
   // Check volume:
   double exact_volume=Exact_bounded_volume-Exact_volume_of_hole-
    Exact_volume_of_internal_region;
   oomph_info 
    << "Volume of region 0 [in mesh, exact, diff (%)]: "
    << volume_in_region0 << " " 
    << exact_volume << " " 
    << abs(volume_in_region0-exact_volume)/
    exact_volume*100.0 << std::endl;
   
   
   
   // Get total mesh volume
   double total_mesh_volume=0.0;
   n_el=Bulk_mesh_pt->nelement();
   for (unsigned e=0;e<n_el;e++)
    {
     total_mesh_volume+=Bulk_mesh_pt->finite_element_pt(e)->size();
    }
   
   
   // Check volume:
   exact_volume=Exact_bounded_volume-Exact_volume_of_hole;
   oomph_info 
    << "Total volume      [in mesh, exact, diff (%)]: "
    << total_mesh_volume << " " 
    << exact_volume << " " 
    << abs(total_mesh_volume-exact_volume)/
    exact_volume*100.0 << std::endl;
   
   // Doc volumes
   volume_file << volume_in_internal_region << " " 
               << total_mesh_volume << " " 
               << volume_in_region0 << std::endl;
   volume_file.close();

  }
 
 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 if (do_bulk_output) 
  {
   Bulk_mesh_pt->output(some_file,nplot);
  }
 some_file.close();

 // Output solution showing element outlines
 //-----------------------------------------
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 if (do_bulk_output) 
  {
   Bulk_mesh_pt->output(some_file,2);
  }
 some_file.close();

 // Output solution for paraview
 //-----------------------------
 sprintf(filename,"%s/soln%i.vtu",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 if (do_bulk_output) 
  {
   Bulk_mesh_pt->output_paraview(some_file,nplot);
  }
 some_file.close();

 // Output solution showing element outlines for paraview
 //------------------------------------------------------
 sprintf(filename,"%s/coarse_soln%i.vtu",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 if (do_bulk_output) 
  {
   Bulk_mesh_pt->output_paraview(some_file,2);
  }
 some_file.close();


 // Get norm of solution
 //---------------------
 sprintf(filename,"%s/norm%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double norm_soln=0.0;
 Bulk_mesh_pt->compute_norm(norm_soln);  
 some_file << sqrt(norm_soln) << std::endl;
 oomph_info << "Norm of computed solution: "   << sqrt(norm_soln)  << endl;
 some_file.close();

} // end of doc




//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{

 MPI_Helpers::init(argc,argv);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Suppress bulk output
 CommandLineArgs::specify_command_line_flag("--suppress_bulk_output");

 // Assess actual element sizes for specified target volume
 CommandLineArgs::specify_command_line_flag(
  "--desired_target_volume",
  &Global_Parameters::Initial_element_volume);

 // Do it via file or via overall specification?
 CommandLineArgs::specify_command_line_flag(
  "--specify_gmsh_target_volume_via_file");

 /// Do uniform refinement tests? Specify number of uniform refinements
 CommandLineArgs::specify_command_line_flag(
  "--check_uniform_refinement",
  &Global_Parameters::Nrefinement);

 /// Do uniform unrefinement tests? Specify number of uniform unrefinements
 CommandLineArgs::specify_command_line_flag(
  "--check_uniform_unrefinement",
  &Global_Parameters::Nunrefinement);

 // Factor by which element volume is (attempted to be) increased
 // during uniform unrefiment
 CommandLineArgs::specify_command_line_flag(
  "--unrefinement_volume_increase_factor",
  &Global_Parameters::Element_volume_increase_factor_during_unrefinement);

 // Factor by which element volume is (attempted to be) reduced
 // during uniform refinement
 CommandLineArgs::specify_command_line_flag(
  "--refinement_volume_decrease_factor",
  &Global_Parameters::Element_volume_reduction_factor_during_refinement);

 // Typical variation in element volume in initial mesh
 CommandLineArgs::specify_command_line_flag(
  "--typical_var_in_element_volume_in_percent",
  &Global_Parameters::Typical_var_in_element_volume_in_percent);

 // // Base adaptation on average size?
 // CommandLineArgs::specify_command_line_flag(
 //  "--base_target_volume_during_adapt_on_av_size_with_var_in_percent",
 //  &hierher_shite::var_in_target_based_on_average_element_size);

 // No innner boundaries
 CommandLineArgs::specify_command_line_flag(
  "--suppress_inner_boundaries");

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // if (CommandLineArgs::command_line_flag_has_been_set(
 //      "--base_target_volume_during_adapt_on_av_size_with_var_in_percent"))
 //  {
 //   hierher_shite::hierher_base_target_on_average_element_size=true;
 //  }


 // // hierher
 // tetgen_hierher_shite::element_volume_reduction_factor_during_refinement*=2.00; // 25;
 // oomph_info 
 //  << "reduced during refinement factor to "
 //  << tetgen_hierher_shite::element_volume_reduction_factor_during_refinement
 //  << std::endl;


 // Note that this can make tetgen die!
 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // // hierher
 // SamplePointContainer::Percentage_offset=0.0;
 // oomph_info << "hierher: SamplePointContainer::Percentage_offset = "
 //            << SamplePointContainer::Percentage_offset << std::endl;

 // Shut up prefix
 oomph_info.output_modifier_pt()=&default_output_modifier;

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=5;

 // Build problem
 TetmeshPoissonProblem<ProjectablePoissonElement<
  TPoissonElement<3,3> > > problem;

 //Output initial guess
 problem.doc_solution(nplot,doc_info);
 doc_info.number()++;



 // // hierher
 // oomph_info << "hierher changing refinement factor from " 
 //            <<  Global_Parameters::Gmsh_parameters_pt->
 //  element_volume_reduction_factor_during_refinement();

 // Global_Parameters::Gmsh_parameters_pt->
 //  element_volume_reduction_factor_during_refinement()=0.6;

 // oomph_info << "to " 
 //            <<  Global_Parameters::Gmsh_parameters_pt->
 //  element_volume_reduction_factor_during_refinement() << std::endl;

 for (unsigned up_down=0;up_down<2;up_down++)
  {
   
   unsigned max_adapt=10; // 5; 
   if (up_down==1) max_adapt=20; 
   for (unsigned i=0;i<=max_adapt;i++)
    {
     // Solve the bastard!
     //problem.newton_solve();
     
     //Output solution
     //problem.doc_solution(nplot,doc_info);
     //doc_info.number()++;
     
     if (i!=max_adapt)
      {
       if (up_down==0)
        {
         oomph_info << "Genuine adaptation top right box\n";
        }
       else        
        {
         oomph_info << "Adaptive unrefinement...\n";
        }
       problem.adapt();
       problem.doc_solution(nplot,doc_info); 
       doc_info.number()++;
      }
    }

   // // Move refinement region
   // if (up_down==0)
   //  {   
   //   Vector<double> lower_left(3,-1.2);
   //   Vector<double> upper_right(3,0.0);
   //   DummyErrorEstimator* err_est_pt=
   //   dynamic_cast<DummyErrorEstimator*>(
   //    dynamic_cast<RefineableMeshBase*>(problem.mesh_pt())->
   //    spatial_error_estimator_pt());
   //   err_est_pt->set_bounding_box(lower_left, 
   //                                upper_right);

   //   unsigned max_adapt=10; 
   //   for (unsigned i=0;i<=max_adapt;i++)
   //    {
   //     if (i!=max_adapt)
   //      {           
   //       oomph_info << "Genuine adaptation lower left box\n";
   //       problem.adapt();
   //       problem.doc_solution(nplot,doc_info); 
   //       doc_info.number()++;
   //      }
   //    }

    // }

   // Switch off forced refinement -- should now slowly revert back
   dynamic_cast<DummyErrorEstimator*>(
    dynamic_cast<RefineableMeshBase*>(problem.mesh_pt())
    ->spatial_error_estimator_pt())->
    set_error_in_target_region(0.0);
  }

}



