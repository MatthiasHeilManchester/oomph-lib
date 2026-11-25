//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2026 Matthias Heil and Andrew Hazel
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

// Oomph-lib includes
#include "generic.h"

// Include Poisson elements/equations
#include "poisson.h"

// Include the mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

namespace
{
  //============================================================================
  /// Throw an oomph-lib error if a self-test condition is not satisfied.
  //============================================================================
  void require_mpi_test(const bool passed, const std::string& message)
  {
    if (!passed)
    {
      throw OomphLibError(
        message, "generic_mpi_test", OOMPH_EXCEPTION_LOCATION);
    }
  }

  //============================================================================
  /// Check a condition and fail consistently on every MPI rank.
  //============================================================================
  void require_mpi_test_on_all_ranks(const OomphCommunicator* const comm_pt,
                                     const bool passed,
                                     const std::string& message)
  {
    unsigned local_passed = passed ? 1 : 0;
    unsigned all_passed = 0;
    MPI_Allreduce(&local_passed,
                  &all_passed,
                  1,
                  MPI_UNSIGNED,
                  MPI_MIN,
                  comm_pt->mpi_comm());
    require_mpi_test(all_passed != 0, message);
  }

  //============================================================================
  /// Construct a distributed distribution that should fail validation.
  //============================================================================
  void expect_invalid_distribution(const OomphCommunicator* const comm_pt,
                                  const unsigned first_row,
                                  const unsigned nrow_local,
                                  const unsigned nrow)
  {
    bool caught_error = false;
    try
    {
      LinearAlgebraDistribution invalid_distribution(
        comm_pt, first_row, nrow_local, nrow);
    }
    catch (OomphLibError& error)
    {
      caught_error = true;
      error.disable_error_message();
    }
    require_mpi_test(caught_error,
                     "An invalid distributed row partition was accepted.");
  }
} // namespace

//==start_of_namespace================================================
/// Namespace for fish-shaped solution of 1D Poisson equation
//====================================================================
namespace FishSolnOneDPoisson
{

 /// Sign of the source function 
 /// (- gives the upper half of the fish, + the lower half)
 int Sign=-1;


 /// Exact, fish-shaped solution as a 1D vector
 void get_exact_u(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = double(Sign)*((sin(sqrt(30.0))-1.0)*x[0]-sin(sqrt(30.0)*x[0]));
 }


 /// Source function required to make the fish shape an exact solution 
 void source_function(const Vector<double>& x, double& source)
 {
  source = double(Sign)*30.0*sin(sqrt(30.0)*x[0]);
 }

} // end of namespace







//==start_of_problem_class============================================
/// 1D Poisson problem in unit interval.
//====================================================================
template<class ELEMENT> 
class OneDPoissonProblem : public Problem
{

public:

 /// Constructor: Pass number of elements and pointer to source function
 OneDPoissonProblem(const unsigned& n_element, 
                    PoissonEquations<1>::PoissonSourceFctPt source_fct_pt);

 /// Destructor (empty)
 ~OneDPoissonProblem()
  {
   delete mesh_pt();
  }

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(const unsigned& label);

private:

 /// Pointer to source function
 PoissonEquations<1>::PoissonSourceFctPt Source_fct_pt;

}; // end of problem class





//=====start_of_constructor===============================================
/// Constructor for 1D Poisson problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT>
OneDPoissonProblem<ELEMENT>::OneDPoissonProblem(const unsigned& n_element,
 PoissonEquations<1>::PoissonSourceFctPt source_fct_pt) : 
 Source_fct_pt(source_fct_pt)
{ 
 Problem::Sparse_assembly_method = Perform_assembly_using_two_arrays;


 // shut up
 disable_info_in_newton_solve();

// Problem::Problem_is_nonlinear = false;
 // Set domain length 
 double L=1.0;

 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L);

 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 mesh_pt()->boundary_node_pt(1,0)->pin(0);

 // Complete the setup of the 1D Poisson problem:

 // Loop over elements and set pointers to source function
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *elem_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
   
   //Set the source function pointer
   elem_pt->source_fct_pt() = Source_fct_pt;
  }

 // Setup equation numbering scheme
 assign_eqn_numbers();

} // end of constructor




//===start_of_actions_before_newton_solve========================================
/// Update the problem specs before solve: (Re)set boundary values
/// from the exact solution. 
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::actions_before_newton_solve()
{
 
 // Assign boundary values for this problem by reading them out
 // from the exact solution.

 // Left boundary is node 0 in the mesh:
 Node* left_node_pt=mesh_pt()->node_pt(0);

 // Determine the position of the boundary node (the exact solution
 // requires the coordinate in a 1D vector!)
 Vector<double> x(1);
 x[0]=left_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 Vector<double> u(1);
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 left_node_pt->set_value(0,u[0]);


 // Right boundary is last node in the mesh:
 unsigned last_node=mesh_pt()->nnode()-1;
 Node* right_node_pt=mesh_pt()->node_pt(last_node);

 // Determine the position of the boundary node
 x[0]=right_node_pt->x(0);
 
 // Boundary value (read in from exact solution which returns
 // the solution in a 1D vector)
 FishSolnOneDPoisson::get_exact_u(x,u);
 
 // Assign the boundary condition to one (and only) nodal value
 right_node_pt->set_value(0,u[0]);

 
} // end of actions before solve



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT>
void OneDPoissonProblem<ELEMENT>::doc_solution(const unsigned& label)
{ 
 using namespace StringConversion;

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution with specified number of plot points per element
 ofstream solution_file(("soln" + to_string(label) + ".dat").c_str());
 mesh_pt()->output(solution_file,npts);
 solution_file.close();

 // Output exact solution at much higher resolution (so we can
 // see how well the solutions agree between nodal points)
 ofstream exact_file(("exact_soln" + to_string(label) + ".dat").c_str());
 mesh_pt()->output_fct(exact_file,20*npts,FishSolnOneDPoisson::get_exact_u); 
 exact_file.close();

 // Doc pointwise error and compute norm of error and of the solution
 double error,norm;
 ofstream error_file(("error" + to_string(label) + ".dat").c_str());
 mesh_pt()->compute_error(error_file,FishSolnOneDPoisson::get_exact_u,
                          error,norm); 
 error_file.close();

 // Doc error norm:
 cout << "\nNorm of error    : " << sqrt(error) << std::endl; 
 cout << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
 cout << std::endl;

} // end of doc

 

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//===start_of_main=============================================================
/// Driver code for generic mpi stuff
//=============================================================================
int main(int argc, char* argv[])
{

  // Initialise MPI
  MPI_Helpers::init(argc,argv);

  // Switch off output modifier
  oomph_info.output_modifier_pt() = &default_output_modifier;
  
  // Define processor-labeled output file for all on-screen stuff
  std::ofstream output_stream;
  char filename[100];
  snprintf(filename, sizeof(filename), "OUTPUT.%i",MPI_Helpers::communicator_pt()->my_rank());
  output_stream.open(filename);
  oomph_info.stream_pt() = &output_stream;
  OomphLibWarning::set_stream_pt(&output_stream);
  OomphLibError::set_stream_pt(&output_stream);  
  
  // Get the global oomph-lib communicator 
  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // Get rank and total number of processors. 
  unsigned my_rank = comm_pt->my_rank();
  unsigned nproc = comm_pt->nproc();

  if (nproc > 1)
  {
    OomphCommunicator local_comm(comm_pt);
    Vector<int> empty_int_vector;
    Vector<double> empty_double_vector;
    local_comm.broadcast(0, empty_int_vector);
    local_comm.broadcast(0, empty_double_vector);
    require_mpi_test_on_all_ranks(
      comm_pt,
      empty_int_vector.empty() && empty_double_vector.empty(),
      "Broadcasting an empty vector changed its size.");
  }

  // Tell us who you are...
  oomph_info << "I'm rank " << my_rank << " on a total of "
             << nproc << " processors" << std::endl;



  
  // Create a uniformly distributed LinearAlgebraDistribution
  // 100 global rows are uniformly distributed accross the processes of 
  // comm_pt
  unsigned nrow_global = 100;
  LinearAlgebraDistribution distributed_distribution(comm_pt,
                                                     nrow_global); 

  if (nproc > 1)
  {
    const unsigned uniform_first_row = my_rank * nrow_global / nproc;
    const unsigned uniform_end_row = (my_rank + 1) * nrow_global / nproc;
    require_mpi_test(
      distributed_distribution.first_row() == uniform_first_row &&
        distributed_distribution.nrow_local() ==
          uniform_end_row - uniform_first_row,
      "A valid uniform distributed row partition was not preserved.");

    const unsigned unequal_nrow_local = my_rank + 1;
    const unsigned unequal_first_row = my_rank * (my_rank + 1) / 2;
    const unsigned unequal_nrow = nproc * (nproc + 1) / 2;
    LinearAlgebraDistribution unequal_distribution(
      comm_pt, unequal_first_row, unequal_nrow_local, unequal_nrow);
    require_mpi_test(
      unequal_distribution.nrow_local() == unequal_nrow_local,
      "A valid unequal distributed row partition was not preserved.");

    if (my_rank == 0)
    {
      expect_invalid_distribution(comm_pt, 0, 1, nproc + 1); // gap at rank 1
      expect_invalid_distribution(comm_pt,
                                  0,
                                  2,
                                  nproc + 1); // overlap at rank 1
      expect_invalid_distribution(comm_pt,
                                  my_rank,
                                  1,
                                  nproc + 1); // incorrect final endpoint
    }
    else
    {
      const unsigned rank = my_rank;
      const unsigned gap_first_row = rank + 1;
      const unsigned overlap_first_row = rank;
      const unsigned endpoint_first_row = rank;
      expect_invalid_distribution(
        comm_pt, gap_first_row, 1, nproc + 1);
      expect_invalid_distribution(
        comm_pt, overlap_first_row, 1, nproc + 1);
      expect_invalid_distribution(
        comm_pt, endpoint_first_row, 1, nproc + 1);
    }

    unsigned zero_row_first_row = 0;
    unsigned zero_row_nrow_local = 0;
    unsigned zero_row_nrow = 1;
    if (nproc >= 3)
    {
      zero_row_nrow_local = (my_rank < 2) ? 0 : 1;
      zero_row_first_row = (my_rank < 2) ? 0 : my_rank - 2;
      zero_row_nrow = nproc - 2;
    }
    else if (my_rank == 1)
    {
      zero_row_nrow_local = 1;
    }

    std::ostringstream zero_row_warning_stream;
    OomphLibWarning::set_stream_pt(&zero_row_warning_stream);
    LinearAlgebraDistribution::Suppress_warning_about_zero_row_distribution =
      true;
    LinearAlgebraDistribution zero_row_distribution(
      comm_pt, zero_row_first_row, zero_row_nrow_local, zero_row_nrow);
    require_mpi_test(
      zero_row_distribution.nrow_local() == zero_row_nrow_local,
      "A valid distributed partition with zero-row processors was rejected.");
    require_mpi_test(zero_row_distribution.nrow() == zero_row_nrow,
                     "A valid zero-row partition has an incorrect global size.");
    require_mpi_test_on_all_ranks(
      comm_pt,
      zero_row_warning_stream.str().empty(),
      "The zero-row distribution warning was not suppressed.");

    LinearAlgebraDistribution::Suppress_warning_about_zero_row_distribution =
      false;
    LinearAlgebraDistribution zero_row_warning_distribution(
      comm_pt, zero_row_first_row, zero_row_nrow_local, zero_row_nrow);
    require_mpi_test(
      zero_row_warning_distribution.nrow_local() == zero_row_nrow_local,
      "A valid zero-row distribution changed its local row count.");
    const bool warning_was_emitted =
      (my_rank == 0) ?
        (zero_row_warning_stream.str().find(
           "This distributed LinearAlgebraDistribution is legal") !=
         std::string::npos) :
        zero_row_warning_stream.str().empty();
    OomphLibWarning::set_stream_pt(&output_stream);
    require_mpi_test_on_all_ranks(
      comm_pt,
      warning_was_emitted,
      "The zero-row distribution warning was not emitted once on rank 0.");

    zero_row_warning_stream.str("");
    zero_row_warning_stream.clear();
    LinearAlgebraDistribution repeated_zero_row_distribution(
      comm_pt, zero_row_first_row, zero_row_nrow_local, zero_row_nrow);
    require_mpi_test_on_all_ranks(
      comm_pt,
      zero_row_warning_stream.str().empty(),
      "The zero-row distribution warning was repeated.");

    LinearAlgebraDistribution::Suppress_warning_about_zero_row_distribution =
      true;
    const unsigned trailing_zero_first_row = (my_rank == 0) ? 0 : 1;
    const unsigned trailing_zero_nrow_local = (my_rank == 0) ? 1 : 0;
    LinearAlgebraDistribution trailing_zero_distribution(
      comm_pt, trailing_zero_first_row, trailing_zero_nrow_local, 1);
    require_mpi_test_on_all_ranks(
      comm_pt,
      trailing_zero_distribution.nrow_local() == trailing_zero_nrow_local,
      "A valid partition with trailing zero-row processors was rejected.");

    if (nproc >= 3)
    {
      const unsigned interior_zero_first_row = (my_rank == 0) ? 0 : 1;
      const unsigned interior_zero_nrow_local =
        (my_rank == 0 || my_rank == nproc - 1) ? 1 : 0;
      LinearAlgebraDistribution interior_zero_distribution(
        comm_pt, interior_zero_first_row, interior_zero_nrow_local, 2);
      require_mpi_test_on_all_ranks(
        comm_pt,
        interior_zero_distribution.nrow_local() ==
          interior_zero_nrow_local,
        "A valid partition with interior zero-row processors was rejected.");
    }
    LinearAlgebraDistribution::Suppress_warning_about_zero_row_distribution =
      false;
  }

  // Show us how many rows this processor holds
  oomph_info << "distributed distribution: first_row and nrow_local: "
             << distributed_distribution.first_row() << " " 
             << distributed_distribution.nrow_local() 
             << std::endl;
  

  // Construct an empty distribution object (does not specify a distribution)
  LinearAlgebraDistribution locally_replicated_distribution;
  
  // Build a locally replicated distribution such that every row is available
  // on every process
  bool distributed=false;
  locally_replicated_distribution.build(comm_pt,nrow_global,distributed);
  
  // Show us how many rows this processor holds
  oomph_info << "locally replicated distribution: first_row and nrow_local: "
             << locally_replicated_distribution.first_row() << " " 
             << locally_replicated_distribution.nrow_local() 
             << std::endl;
  
  // Get the number of local rows on this process
  unsigned nrow_local = distributed_distribution.nrow_local();
  
  // Get the first row on this process
  unsigned first_row = distributed_distribution.first_row();
  
  // Get the number of global rows
  nrow_global = distributed_distribution.nrow();
  
  // Is this distributed (true) or locally replicated (false)
  distributed = distributed_distribution.distributed();
  
  // Does this object specify a distribution
  bool built = distributed_distribution.built();
  


  // Construct a uniformly distributed DoubleVector with unit elements
  DoubleVector my_vector(distributed_distribution,1.0);
  
  // Increment every element of my_vector on this process by 1
  nrow_local = my_vector.distribution_pt()->nrow_local();
  for (unsigned i = 0; i < nrow_local; i++)
   {
    my_vector[i]+=1.0;
   }
  
  // Document elements on this process in my_vector
  nrow_local = my_vector.distribution_pt()->nrow_local();
  first_row = my_vector.distribution_pt()->first_row();
  for (unsigned i = 0; i < nrow_local; i++)
   {
    oomph_info << "local row " << i 
               << " is global row " << first_row+i 
               << " and has value " << my_vector[i] << std::endl; 
   }

  // Redistribute my_vector such that it is locally replicated on all processes
  my_vector.redistribute(&locally_replicated_distribution);
  
  // (Re)build my_vector such that it is uniformly distributed over all
  // processes
  my_vector.build(distributed_distribution,1.0);

  // Construct an empty DoubleVector
  DoubleVector another_vector;
  
  // Clear all data from an existing DoubleVector
  my_vector.clear();

  
  // Construct an empty CRDoubleMatrix
  CRDoubleMatrix my_matrix;

  // Specify that the rows be uniformly distributed
  my_matrix.build(&distributed_distribution);

  // Vector of coefficient of value 1.0
  Vector<double> values(nrow_local,1.0);
  
  // Column indices corresponding to values
  Vector<int> column_indices(nrow_local);
  
  // Index of vectors values and column_indices where the i-th row starts
  // (each row contains one coefficient)
  Vector<int> row_start(nrow_local+1);
  
  // populate column_indices and row_start
  for (unsigned i = 0; i < nrow_local; ++i)
   {
    column_indices[i]=first_row+nrow_local;
    row_start[i]=i;
   }
  row_start[nrow_local]=nrow_local;
  
  // Build the (square) matrix
  unsigned ncol = nrow_global;
  my_matrix.build(ncol,values,column_indices,row_start);
  
  CRDoubleMatrix my_matrix2(&distributed_distribution,ncol,values,
                           column_indices,row_start);
  
  if (nproc > 1)
  {
    const unsigned target_first_row = my_rank;
    const unsigned target_nrow_local =
      (my_rank == nproc - 1) ? nrow_global - nproc + 1 : 1;
    LinearAlgebraDistribution target_distribution(
      comm_pt, target_first_row, target_nrow_local, nrow_global);
    bool caught_unsupported_redistribution = false;
    try
    {
      my_matrix2.redistribute(&target_distribution);
    }
    catch (OomphLibError& error)
    {
      caught_unsupported_redistribution = true;
      error.disable_error_message();
    }
    require_mpi_test(caught_unsupported_redistribution,
                     "Distributed-to-distributed CRDoubleMatrix "
                     "redistribution did not raise OomphLibError.");
  }

  // Get the first (global) row of my_matrix on this process
  first_row =  my_matrix2.distribution_pt()->first_row();
  
  // Get the first (global) row of my_matrix on this process
  first_row =  my_matrix2.first_row();





  // Set up a problem: 
  // Solve a 1D Poisson problem using a source function that generates
  // a fish shaped exact solution
  unsigned n_element=40; 
  OneDPoissonProblem<QPoissonElement<1,4> > 
   problem(n_element,FishSolnOneDPoisson::source_function);
    
  // Get the residual and Jacobian, by default both are uniformly distributed
  // over all available processes
  my_vector.clear();
  my_matrix.clear();
  
  // Get the Jacobian
  problem.get_jacobian(my_vector,my_matrix);
  
  oomph_info 
   << "Uniformly distributed residual vector: first_row and nrow_local: " 
   << my_vector.first_row() << " " 
   << my_vector.nrow_local() << std::endl
   << "Uniformly distributed jacobian matrix: first_row, nrow_local and nnz: " 
   << my_matrix.first_row() << " " 
   << my_matrix.nrow_local() << " " 
   << my_matrix.nnz() << " " 
   << std::endl;


  // Request locally replicated residual and Jacobian from the problem
  distributed=false;
  LinearAlgebraDistribution locally_replicated_distribution_for_jac(
   comm_pt,problem.ndof(),distributed);
  
  // Show us how many rows this processor holds
  oomph_info 
   << "locally replicated distribution_for_jac: first_row and nrow_local: "
   << locally_replicated_distribution_for_jac.first_row() << " " 
   << locally_replicated_distribution_for_jac.nrow_local() 
   << std::endl;
  
  my_vector.build(locally_replicated_distribution_for_jac);
  my_matrix.build(&locally_replicated_distribution_for_jac);

  // Get the Jacobian
  problem.get_jacobian(my_vector,my_matrix);

  oomph_info 
   << "Replicated residual vector: first_row and nrow_local: " 
   << my_vector.first_row() << " " 
   << my_vector.nrow_local() << std::endl
   << "Replicated jacobian matrix: first_row, nrow_local and nnz: " 
   << my_matrix.first_row() << " " 
   << my_matrix.nrow_local() << " " 
   << my_matrix.nnz() << " " 
   << std::endl;

  Data global_data(1);
  problem.add_global_data(&global_data);
  bool caught_unsupported_global_data = false;
  try
  {
    problem.assign_eqn_numbers();
  }
  catch (OomphLibError& error)
  {
    caught_unsupported_global_data = true;
    error.disable_error_message();
  }
  require_mpi_test(caught_unsupported_global_data == (nproc > 1),
                   "Free problem-level global data had an unexpected MPI "
                   "support result.");

  global_data.pin(0);
  problem.assign_eqn_numbers();
  problem.flush_global_data();


  // Finalize MPI
  MPI_Helpers::finalize();


  // dummy line to suppress warning about unused built variable
  if (built)
   {
    built=false;
   }
  
} // end_of_main
