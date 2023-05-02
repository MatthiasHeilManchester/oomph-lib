//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
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
// Header file for the Biharmonic Bell elements
#ifndef OOMPH_MYBIHARMONIC_ELEMENTS_HEADER
#define OOMPH_MYBIHARMONIC_ELEMENTS_HEADER


#include<sstream>

//OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"
//#include "../C1_basis/SubparametricTriangleElement.h"


namespace oomph
{
 //=============================================================
 /// A class for all subparametric elements that solve the 2D-
 /// Biharmonic equations.
 /// \f[
 /// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
 /// \f]
 /// This contains the generic maths. Shape functions, geometric
 /// mapping etc. must get implemented in derived class.
 //=============================================================
 class FoepplVonKarmanEquations : public virtual FiniteElement
 {

 public:
  /// A pointer to a scalar function of the position. Can be used for
  /// out-of-plane forcing, swelling, isotropic-prestrain, etc.
  typedef void (*ScalarFctPt)(const Vector<double>& x, double& f);
    
  /// A pointer to a vector function of the position. Can be used for
  /// in-of-plane forcing, anisotropic-prestrain, etc.
  typedef void (*VectorFctPt)(const Vector<double>& x,
				      Vector<double>& forcing);

  /// Function pointer to the Error Metric we are using
  /// e.g could be that we are just interested in error on w etc.
  typedef void (*ErrorMetricFctPt)(const Vector<double>& x,
				   const Vector<double>& u,
				   const Vector<double>& u_exact,
				   double& error,
				   double& norm);

  /// Function pointer to the Error Metric we are using if we want multiple
  /// errors.  e.g could be we want errors seperately on each displacment
  typedef void (*MultipleErrorMetricFctPt)(const Vector<double>& x,
					   const Vector<double>& u,
					   const Vector<double>& u_exact,
					   Vector<double>& error,
					   Vector<double>& norm);

  /// (pure virtual) interface to fixing the out of plane displacement
  virtual void fix_out_of_plane_displacement_dof(const unsigned& dof_number, 
						 const unsigned& boundary_number,
						 const ScalarFctPt& w) = 0;

  /// (pure virtual) interface to fixing the in plane displacement
  virtual void fix_in_plane_displacement_dof(const unsigned& dof_number,
					     const unsigned& boundary_number,
					     const ScalarFctPt& u)=0;

  /// (pure virtual) interface to return the number of nodes used by u
  virtual unsigned nu_node() const = 0;
  /// (pure virtual) interface to return the number of nodes used by w
  virtual unsigned nw_node() const = 0;
  
  /// (pure virtual) interface to get the local indices of the nodes used by u
  virtual Vector<unsigned> get_u_nodes() const = 0;
  /// (pure virtual) interface to get the local indices of the nodes used by w
  virtual Vector<unsigned> get_w_nodes() const = 0;
  
  /// (pure virtual) interface to get the number of basis types for u at node j
  virtual unsigned nu_type_at_node(const unsigned& j_node) const = 0;
  /// (pure virtual) interface to get the number of basis types for w at node j
  virtual unsigned nw_type_at_node(const unsigned& j_node) const = 0;
  
  /// (pure virtual) interface to retrieve the value of u_alpha at node j of
  /// type k
  virtual double get_u_alpha_value_at_node_of_type(const unsigned& alpha,
						   const unsigned& node_j,
						   const unsigned& type_k) const = 0;
  /// (pure virtual) interface to retrieve the value of w at node j of type k
  virtual double get_w_value_at_node_of_type(const unsigned& node_j,
					     const unsigned& type_k) const = 0;
  
  /// (pure virtual) interface to get the pointer to the internal data used to
  /// interpolate u (NOTE: assumes each u field has exactly one internal data)
  virtual Vector<Data*> u_internal_data_pts() const = 0;
  /// (pure virtual) interface to get the pointer to the internal data used to
  /// interpolate u (NOTE: assumes w field has exactly one internal data)
  virtual Data* w_internal_data_pt() const = 0;

  /// (pure virtual) interface to get the number of internal types for the u
  /// fields
  virtual unsigned nu_internal_type() const = 0;
  /// (pure virtual) interface to get the number of internal types for the w
  /// fields
  virtual unsigned nw_internal_type() const = 0;
  
  /// (pure virtual) interface to retrieve the value of u_alpha of internal
  /// type k
  virtual double get_u_alpha_internal_value_of_type(const unsigned& alpha,
						    const unsigned& type_k) const = 0;
  /// (pure virtual) interface to retrieve the value of w of internal type k
  virtual double get_w_internal_value_of_type(const unsigned& type_k) const = 0;
  
  // // ===================== New functions [zdec] ================================
  // /// (pure virtual) interface to  get the number of fields (unknowns)
  // /// interpolated by the element
  // virtual unsigned nfield() const = 0;
 
  // // ------------------------ Nodal data ---------------------------------------
  // // Nodal basis functions might vary between field and node, here is the
  // // interface for accessing the appropriate dof/basis in generality.
  
  // /// (pure virtual) interface to get the number of nodes that field i is
  // /// interpolated over
  // virtual unsigned nnode_for_field(const unsigned& i_field) const = 0;

  // /// (pure virtual) interface to get the nodes associated with interpolating
  // /// field i
  // virtual Vector<unsigned> nodes_for_field(const unsigned& i_field) const = 0;
  
  // /// (pure virtual) interface to get the number of basis type for field i at
  // /// node j
  // virtual unsigned ntype_for_field_at_node(const unsigned& i_field,
  // 					   const unsigned& j_node) const = 0;

  // /// Get the dof of the field i at node j of type k
  // virtual double nodal_value_for_field_at_node_of_type(const unsigned& i_field,
  // 						       const unsigned& j_node,
  // 						       const unsigned& k_type) const = 0;
  
  // /// Get the dof of the field i at node j of type k at time t
  // virtual double nodal_value_for_field_at_node_of_type(const unsigned& t,
  // 						       const unsigned& i_field,
  // 						       const unsigned& j_node,
  // 						       const unsigned& k_type) const = 0;
  
  // // ----------------------- Internal data -------------------------------------
  // // Each field has its own internal data which is resized according to the
  // // number of dof/basis types required
  
  // /// Get the number of internal data for field i
  // virtual unsigned ninternal_types_for_field(const unsigned& i_field) const = 0;

  // /// Get the index of the internal data for field i
  // virtual unsigned index_of_internal_data_for_field(const unsigned& i_field) const = 0;
  
  // /// Return the pointer to the internal data for field i
  // virtual Data* internal_data_for_field_pt(const unsigned& i_field) const = 0;
  
  // /// Return the value at the internal data for field i type k
  // virtual double internal_value_for_field_of_type(const unsigned& i_field,
  // 						  const unsigned& k_type) const = 0;
   
  // /// Return the value at the internal data for field i type k at time t
  // virtual double internal_value_for_field_of_type(const unsigned& t,
  // 						  const unsigned& i_field,
  // 						  const unsigned& k_type) const = 0;
  
  // /// Get the jth bubble dof at the lth internal point. Deliberately broken 
  // /// for case when there is no curved edge.
  // virtual int local_internal_equation(const unsigned& i_field,
  // 				      const unsigned& k_type) const = 0;
   
  // // ========================= End new [zdec] ==================================

  /// Is this element curved or not
  virtual bool element_is_curved() const = 0;

  /// Get the value of the damping flag for u
  virtual bool u_is_damped() const
  {
    return U_is_damped;
  }
  /// Get the value of the damping flag for w
  virtual bool w_is_damped() const
  {
    return W_is_damped;
  }
 
  // // Get the number of nodes of the out-of-plane functions, pure virtual
  // virtual unsigned nnode_outofplane() const =  0;
 
  // // Get the number of nodes of the out-of-plane functions, pure virtual
  // virtual unsigned nnode_inplane() const =  0;
 
  // // Get the number of basis functions, pure virtual
  // virtual unsigned nnodal_basis_type() const =  0;

  // // Get the number of internal basis functions, pure virtual
  // virtual unsigned nbubble_basis() const = 0;
 
  // // Get the number of internal basis functions, pure virtual
  // virtual unsigned nbubble_basis_type() const = 0;

  /// Flag to control damping of the in-plane variables
  bool U_is_damped = false;

  /// Flag to control damping of the out-of-plane variable
  bool W_is_damped = true;
  
 public:
 
  /// Get pointer to association matrix 
  DenseMatrix<double> *get_association_matrix_pt()const {return Association_matrix_pt;};

  /// Eta
  const double &eta() const {return *Eta_pt;}

  /// Pointer to eta
  const double* &eta_pt() {return Eta_pt;}

  /// Pure virtual function to pin all deflection dofs
  virtual void pin_all_deflection_dofs() const=0;

  /// Constructor (must initialise the Pressure_fct_pt to null)
  FoepplVonKarmanEquations() : Pressure_fct_pt(0),
			       In_plane_forcing_fct_pt(0),
			       Swelling_fct_pt(0),
			       Error_metric_fct_pt(0),
			       Multiple_error_metric_fct_pt(0),
			       Association_matrix_pt(0) 
  {
   Eta_pt = &Default_Eta_Value;
   Nu_pt = &Default_Nu_Value;
   Mu_pt = &Default_Mu_Value;
  }

  /// Broken copy constructor
  FoepplVonKarmanEquations(const FoepplVonKarmanEquations& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanEquations");
  }

  /// Broken assignment operator
  void operator=(const FoepplVonKarmanEquations&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanEquations");
  }

  /// Return the index at which the first u unknown value is stored.
  /// In derived multi-physics elements, this function should be overloaded
  /// to reflect the chosen storage scheme. Note that these equations require
  /// that the unknown is always stored at the same index at each node.
  /// Note these are stored before w_dofs otherwise at vertex nodes the 
  /// stored value would be at a different index
  virtual inline unsigned u_nodal_index_foeppl_von_karman() const {return 0;}

  /// Return the index at which the first w unknown value is stored.
  /// In derived multi-physics elements, this function should be overloaded
  /// to reflect the chosen storage scheme. Note that these equations require
  /// that the unknown is always stored at the same index at each node.
  /// Note that at midside nodes there are no w dofs 
  virtual inline unsigned w_nodal_index_foeppl_von_karman() const {return u_nodal_index_foeppl_von_karman()+2;}

  /// Output with default number of plot points
  void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   FoepplVonKarmanEquations::output(outfile,n_plot);
  }
  
  /// \short Output FE representation of soln: x,y,u or x,y,z,u at
  /// n_plot^DIM plot points
  void output(std::ostream &outfile, const unsigned &n_plot);
  
  /// C_style output with default number of plot points
  void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   FoepplVonKarmanEquations::output(file_pt,n_plot);
  }

  /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
  /// n_plot^DIM plot points
  void output(FILE* file_pt, const unsigned &n_plot);

  /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot points
  void output_fct(std::ostream &outfile, const unsigned &n_plot,
		  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

  /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
  /// n_plot^DIM plot points (dummy time-dependent version to
  /// keep intel compiler happy)
  virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
			  const double& time,
			  FiniteElement::UnsteadyExactSolutionFctPt
			  exact_soln_pt)
  {
   throw OomphLibError(
		       "There is no time-dependent output_fct() for these elements ",
		       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


  /// Get error against and norm of exact solution
  void compute_error(std::ostream &outfile,
		     FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
		     double& error, double& norm);

  /// Get error against and norm of exact solution
  void compute_error_in_deflection(std::ostream &outfile,
				   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
				   double& error, double& norm);

  /// Dummy, time dependent error checker
  void compute_error(std::ostream &outfile,
		     FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
		     const double& time, double& error, double& norm)
  {
   throw OomphLibError(
		       "There is no time-dependent compute_error() for these elements",
		       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

  /// Fill in the strain tensor from displacement gradients
  void get_epsilon(DenseMatrix<double>& epsilon,
		   const DenseMatrix<double>& grad_u,
		   const DenseMatrix<double>& grad_w,
		   const double& c_swell)const
  {
   // Truncated Green Lagrange strain tensor
   DenseMatrix<double> dummy_epsilon(this->dim(),this->dim(),0.0);
   for(unsigned alpha=0;alpha<this->dim();++alpha)
    {
     for(unsigned beta=0;beta<this->dim();++beta)
      {
       // Truncated Green Lagrange strain tensor
       dummy_epsilon(alpha,beta) += 0.5* grad_u(alpha,beta)
	+ 0.5*grad_u(beta,alpha)
	+ 0.5*grad_w(0,alpha)*grad_w(0,beta);
      }
     // Swelling slack
     dummy_epsilon(alpha,alpha) -= c_swell;
    
    }
   epsilon=dummy_epsilon;
  }

  /// Fill in the stress tensor from displacement gradients
  void get_sigma(DenseMatrix<double>& sigma,
		 const DenseMatrix<double>& grad_u,
		 const DenseMatrix<double>& grad_w,
		 const double c_swell)const
  {
   // Get the Poisson ratio
   double nu(get_nu());

   // Truncated Green Lagrange strain tensor
   DenseMatrix<double> epsilon(this->dim(),this->dim(),0.0);
   get_epsilon(epsilon, grad_u, grad_w, c_swell);
   
   // Now construct the Stress
   for(unsigned alpha=0;alpha<this->dim();++alpha)
    {
     for(unsigned beta=0;beta<this->dim();++beta)
      {
       // The Laplacian term: Trace[ \epsilon ] I
       // \nu * \epsilon_{\alpha \beta} delta_{\gamma \gamma}
       sigma(alpha,alpha) += nu*epsilon(beta,beta)/(1-nu*nu);
      
       // The scalar transform term: \epsilon
       // (1-\nu) * \epsilon_{\alpha \beta}
       sigma(alpha,beta) += (1-nu)* epsilon(alpha,beta)/(1-nu*nu);
      }
    }
  }

  // TODO: Fix this function to use index loops
  /// Fill in the stress tensor using a precalculated strain tensor
  void get_sigma_from_epsilon(DenseMatrix<double>& sigma,
			      const DenseMatrix<double>& epsilon)const
  {
   // Get the Poisson ratio
   double nu(get_nu());
   
   // Now construct the Stress
   sigma(0,0) = (epsilon(0,0) + nu*epsilon(1,1)) / (1.0 - nu*nu);
   sigma(1,1) = (epsilon(1,1) + nu*epsilon(0,0)) / (1.0 - nu*nu);
   sigma(0,1) = epsilon(0,1) / (1.0 + nu);
   sigma(1,0) = sigma(0,1);
  }

  /// Get the principal stresses from the stress tensor
  void get_principal_stresses(const DenseMatrix<double>& sigma,
			      Vector<double>& eigenvals,
			      DenseMatrix<double>& eigenvecs)const
  {
   // Ensure that our eigenvectors are the right size
   eigenvals.resize(2);
   eigenvecs.resize(2);
  
   // Store the axial and shear stresses
   double s00 = sigma(0,0);
   double s01 = sigma(0,1);
   double s11 = sigma(1,1);

   // Calculate the principal stress magnitudes
   eigenvals[0] =
    0.5 * ( (s00 + s11) + sqrt((s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01)) );
   eigenvals[1] =
    0.5 * ( (s00 + s11) - sqrt((s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01)) );

   // Handle the shear free case
   if(s01==0.0)
    {
     eigenvecs(0,0)=1.0;
     eigenvecs(1,0)=0.0;
     eigenvecs(0,1)=0.0;
     eigenvecs(1,1)=1.0;
    }
  
   else
    {
     // TODO: better (more general) sign choice for streamlines

     // For max eval we choose y-positive evecs (suited to swelling sheet
     // problem)
     double sign = (eigenvals[0]-s00<0.0) ? -1.0 : 1.0;
     // Calculate the normalised principal stress direction for eigenvals[0]
     eigenvecs(0,0) =
      sign * (s01 / sqrt(s01*s01 + (eigenvals[0]-s00)*(eigenvals[0]-s00)));
     eigenvecs(1,0) =
      sign * ((eigenvals[0]-s00) / sqrt(s01*s01 + (eigenvals[0]-s00)*(eigenvals[0]-s00)));

     // For min eval we choose x-positive evecs (suited to swelling sheet
     // problem)
     sign = (s01<0.0) ? -1.0 : 1.0;
     // Calculate the normalised principal stress direction for eigenvals[1]
     eigenvecs(0,1) =
      sign * (s01 / sqrt(s01*s01 + (eigenvals[1]-s00)*(eigenvals[1]-s00)));
     eigenvecs(1,1) =
      sign * ((eigenvals[1]-s00) / sqrt(s01*s01 + (eigenvals[1]-s00)*(eigenvals[1]-s00)));
    }
  }

  /// Access function: Pointer to pressure function
  ScalarFctPt& pressure_fct_pt()
  {return Pressure_fct_pt;}

  /// Access function: Pointer to pressure function. Const version
  ScalarFctPt pressure_fct_pt() const
  {return Pressure_fct_pt;}

  /// Access function: Pointer to in plane forcing function
  VectorFctPt& in_plane_forcing_fct_pt()
  {return In_plane_forcing_fct_pt;}

  /// Access function: Pointer to in plane forcing function. Const version
  VectorFctPt in_plane_forcing_fct_pt() const
  {return In_plane_forcing_fct_pt;}

  /// Access function: Pointer to swelling function
  ScalarFctPt& swelling_fct_pt()
  {return Swelling_fct_pt;}

  /// Access function: Pointer to swelling function. Const version
  ScalarFctPt swelling_fct_pt() const
  {return Swelling_fct_pt;}

  /// Access function: Pointer to error metric function
  ErrorMetricFctPt& error_metric_fct_pt()
  {return Error_metric_fct_pt;}

  /// Access function: Pointer to multiple error metric function
  MultipleErrorMetricFctPt& multiple_error_metric_fct_pt() 
  {return Multiple_error_metric_fct_pt;}

  /// Access function: Pointer to error metric function function
  ErrorMetricFctPt error_metric_fct_pt() const {return Error_metric_fct_pt;}

  /// Access function: Pointer to multiple error metric function
  MultipleErrorMetricFctPt multiple_error_metric_fct_pt() const 
  {return Multiple_error_metric_fct_pt;}

  ///Access function to the Poisson ratio.
  const double*& nu_pt() {return Nu_pt;}

  ///Access function to the Poisson ratio (const version)
  const double& get_nu() const {return *Nu_pt;}

  ///Access function to the dampening coefficient.
  const double*& mu_pt() {return Mu_pt;}

  ///Access function to the dampening coefficient (const version)
  const double& get_mu() const {return *Mu_pt;}

  // [zdec] IN THE BIN
  /*
 /// Get the kth dof type at internal point l at time t(=0)
 virtual double get_w_bubble_dof(const unsigned& l,
 const unsigned& k,
 const unsigned& t = 0) const =0;

 /// Get the kth equation at internal point l
 virtual int local_w_bubble_equation(const unsigned& l, const unsigned& k)
 const =0;
  */
 
  /// Get pressure term at (Eulerian) position x. This function is
  /// virtual to allow overloading in multi-physics problems where
  /// the strength of the pressure function might be determined by
  /// another system of equations.
  inline virtual void get_pressure_foeppl_von_karman(const unsigned& ipt,
						     const Vector<double>& x,
						     double& pressure) const
  {
   //If no pressure function has been set, return zero
   if(Pressure_fct_pt==0)
    {
     pressure = 0.0;
    }
   else
    {
     // Get pressure strength
     (*Pressure_fct_pt)(x,pressure);
    }
  }

  /// Get pressure term at (Eulerian) position x. This function is
  /// virtual to allow overloading in multi-physics problems where
  /// the strength of the pressure function might be determined by
  /// another system of equations.
  inline virtual void get_in_plane_forcing_foeppl_von_karman(const unsigned& ipt,
							     const Vector<double>& x,
							     Vector<double>& pressure) const
  {
   //In plane is same as DIM of problem (2)
   pressure.resize(this->dim());
   //If no pressure function has been set, return zero
   if(In_plane_forcing_fct_pt==0)
    {
     pressure[0] = 0.0;
     pressure[1] = 0.0; 
    }
   else
    {
     // Get pressure strength
     (*In_plane_forcing_fct_pt)(x,pressure);
    }
  }

  /// Get swelling at (Eulerian) position x. This function is
  /// virtual to allow overloading.
  inline virtual void get_swelling_foeppl_von_karman(const Vector<double>& x,
						     double& swelling) const
  {
   //If no swelling function has been set, return zero
   if(Swelling_fct_pt==0)
    {
     swelling = 0.0;
    }
   else
    {
     // Get swelling magnitude
     (*Swelling_fct_pt)(x,swelling);
    }
  }
 
  /// \short Calculate the elastic energy of the element and return it as a
  /// double.
  virtual Vector<double> element_elastic_and_kinetic_energy();

  /// Add the element's contribution to its residual vector (wrapper)
  void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_foeppl_von_karman(
							   residuals,GeneralisedElement::Dummy_matrix,0);
  }


  /// Add the element's contribution to its residual vector and
  /// element Jacobian matrix (wrapper)
  void fill_in_contribution_to_jacobian(Vector<double> &residuals,
					DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_foeppl_von_karman(residuals,jacobian,1);
  }

  /// Add the element's contribution to its residual vector and
  /// element Jacobian matrix (wrapper)
  void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
							DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {
   //Call fill in Jacobian 
   fill_in_contribution_to_jacobian(residuals,jacobian);
   // There is no mass matrix: we will just want J w = 0

   // -- COPIED FROM DISPLACMENT FVK EQUATIONS --
   // Dummy diagonal (won't result in global unit matrix but
   // doesn't matter for zero eigenvalue/eigenvector
   unsigned ndof=mass_matrix.nrow();
   for (unsigned i=0;i<ndof;i++)
    {
     mass_matrix(i,i)+=1.0;
    }
  }
 
  /// Return FE representation of unknown values u(s)
  /// at local coordinate s
  virtual inline Vector<double> interpolated_u_foeppl_von_karman(const Vector<double> &s) const
  {
   //Find the number of unknown fields (should be 3 for FvK)
   //   const unsigned n_field = nfield();
   const unsigned w_index = 2;
   //Find the dimension of the element
   const unsigned dim = this->dim();
   //The number of first derivatives is the dimension of the element
   const unsigned n_deriv = dim;
   //The number of second derivatives is the triangle number of the dimension
   const unsigned n_2deriv = dim*(dim+1)/2;
   //Find out how many nodes there are
   const unsigned n_u_node = nnode_inplane();
   const unsigned n_w_node = nnode_outofplane(); 
   //Find out how many bubble nodes there are
   const unsigned n_b_type = nbubble_basis();
   //Find out how many nodes positional dofs there are
   const unsigned n_basis_type = nnodal_basis_type();
   const unsigned exactly_uno = 1; // [zdec] What am I supposed to do here?
   
   // [zdec] Generalise the retrieval of basis and test functions.
   // Set up memory for the shape and test functions
   Shape psi_u(n_u_node);
   Shape test_u(n_u_node);
   DShape dpsi_udxi(n_u_node,n_deriv);
   DShape dtest_udxi(n_u_node,n_deriv);
 
   //Local c1-shape funtion
   Shape psi_w(n_w_node,n_basis_type);
   Shape test_w(n_w_node,n_basis_type);
   Shape psi_b(n_b_type,exactly_uno);
   Shape test_b(n_b_type,exactly_uno);
   DShape dpsi_wdxi(n_w_node,n_basis_type,n_deriv);
   DShape dtest_wdxi(n_w_node,n_basis_type,n_deriv);
   DShape dpsi_b_dxi(n_b_type,exactly_uno,n_deriv);
   DShape dtest_b_dxi(n_b_type,exactly_uno,n_deriv);
   DShape d2psi_wdxi2(n_w_node,n_basis_type,n_2deriv);
   DShape d2test_wdxi2(n_w_node,n_basis_type,n_2deriv);
   DShape d2psi_b_dxi2(n_b_type,exactly_uno,n_2deriv);
   DShape d2test_b_dxi2(n_b_type,exactly_uno,n_2deriv);

   // [zdec] for the bin
   /*
   // Vectors of pointers to the appropriate nodal basis funtions so we can loop
   // over them. 
   Vector<Shape*> psi_field_pt(n_field, &psi_u);
   Vector<DShape*> dpsi_field_pt(n_field, &dpsi_udxi);
   Vector<DShape*> d2psi_field_pt(n_field, 0); // [zdec] Do we bother making this three long
   psi_field_pt[2]=&psi_w;
   dpsi_field_pt[2]=&dpsi_wdxi;
   d2psi_field_pt[2]=&d2psi_wdxi2;
   
   // Vectors of pointers to the appropriate internal basis funtions so we can
   // loop over them. Null for ux, uy [zdec] TODO generalise
   Vector<Shape*> psi_internal_field_pt(n_field, 0);
   Vector<DShape*> dpsi_internal_field_pt(n_field, 0);
   Vector<DShape*> d2psi_internal_field_pt(n_field, 0); // [zdec] Do we bother making this three long
   psi_internal_field_pt[2]=&psi_b;
   dpsi_internal_field_pt[2]=&dpsi_b_dxi;
   d2psi_internal_field_pt[2]=&d2psi_b_dxi2;
   */
 
   // Number of in-plane displacement fields is equal to #in-plane + #in-plane-derivs
   const unsigned n_u_fields = dim + dim*n_deriv;// DIM + DIM^2;
   // Number of out-of-plane displacement fields is equal 1 (value of deflection_
   // plus first (dim) and second deriv dim*(dim+1)/2
   const unsigned n_w_fields = 1 + n_deriv + n_2deriv;//1 + DIM + DIM(DIM+1)/2;

   //Initialise value of u
   Vector<double> interpolated_u(n_u_fields+n_w_fields,0.0);

   // The temporary field containers for our interpolation loop
   Vector<double> interpolated_f(n_field,0.0);
   DenseMatrix<double> interpolated_dfdxi(n_field,dim,0.0);
   DenseMatrix<double> interpolated_d2fdxi2(n_field,dim*(dim+1)/2,0.0);
   
   //Find values of c1-shape function
   d2shape_and_d2test_eulerian_foeppl_von_karman(s, psi_w, psi_b,
						 dpsi_wdxi, dpsi_b_dxi,
						 d2psi_wdxi2, d2psi_b_dxi2,
						 test_w, test_b,
						 dtest_wdxi, dtest_b_dxi,
						 d2test_wdxi2, d2test_b_dxi2);
  
   dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s, psi_u, dpsi_udxi,
						   test_u, dtest_udxi);
  
   //Interpolated unknown
   //=======================START OF INTERPOLATION=============================
   // Loop over the unknowns [u1,u2,w]
   for(unsigned i_field=0; i_field<n_field; i_field++)
    {
     //---------------NODAL DATA---------------
     // Nodal basis pointers for the current field
     // u1 and u2 share basis, w has a separate one
     Shape* basis_pt = psi_field_pt[i_field];
     DShape* dbasis_pt = dpsi_field_pt[i_field];
     DShape* d2basis_pt = d2psi_field_pt[i_field];
     // We only loop over the number of nodes that we need to for field i
     // (For Bell elements: all nodes for ux, uy; vertex nodes for w)
     unsigned n_node = nnode_for_field(i_field);
     // Get the local node numbers used in the interpolation of the field
     Vector<unsigned> local_nodes_for_field_i = nodes_for_field(i_field);

     // Loop over the nodes
     for(unsigned j_nodei=0; j_nodei<n_node; j_nodei++)
      {
       // Get the `j`-th node associated with field i
       // ()
       unsigned j_node = local_nodes_for_field_i[j_nodei];
       // Get the number of basis types for field i at node j
       unsigned n_type = ntype_for_field_at_node(i_field,j_node);
       // Loop over the types of basis function
       for(unsigned k_type=0; k_type<n_type; k_type++)
	{
	 double nodal_value =
	  nodal_value_for_field_at_node_of_type(i_field, j_node, k_type);
	  
	 // Add nodal contribution to unknown
	 interpolated_f[i_field] +=
	  nodal_value * (*basis_pt)(j_node,k_type);

	 // Add nodal contribution to first derivatives
	 for(unsigned l_deriv=0; l_deriv<n_deriv; l_deriv++)
	  {
	   interpolated_dfdxi(i_field,l_deriv) +=
	    nodal_value * (*dbasis_pt)(j_node, k_type, l_deriv);
	  } // End of loop over first derivatives (l_deriv)

	 // Add nodal contributions to second derivatives (for w only)
	 if(i_field==w_index)
	  {
	   for(unsigned l_2deriv=0; l_2deriv<n_2deriv; l_2deriv++)
	    {
	     interpolated_d2fdxi2(i_field,l_2deriv) +=
	      nodal_value * (*d2basis_pt)(j_node,k_type,l_2deriv);
	    } // End of loop over second derivatives l_2deriv
	  }
	} // End of loop over basis types (k_type)
      } // End of loop over nodes used by field i (j_nodei)
     //-----------END OF NODAL DATA------------
      
     //-------------INTERNAL DATA--------------
     Shape* internal_basis_pt = psi_internal_field_pt[i_field];
     DShape* internal_dbasis_pt = dpsi_internal_field_pt[i_field];
     DShape* internal_d2basis_pt = d2psi_internal_field_pt[i_field];

     // Begin loop over internal data [zdec] TODO generalise
     unsigned n_internal_types = ninternal_types_for_field(i_field);
     for(unsigned k_type=0; k_type<n_internal_types; k_type++)
      {
       double internal_value =
	internal_value_for_field_of_type(i_field, k_type);
	  
       // Add nodal contribution to unknown
       interpolated_f[i_field] +=
	internal_value * (*internal_basis_pt)(k_type,0);

       // Add nodal contribution to first derivatives
       for(unsigned l_deriv=0; l_deriv<n_deriv; l_deriv++)
	{
	 interpolated_dfdxi(i_field,l_deriv) +=
	  internal_value * (*internal_dbasis_pt)(k_type, 0, l_deriv);
	} // End of loop over first derivatives (l_deriv)

       // Add nodal contributions to second derivatives (for w only)
       if(i_field==w_index)
	{
	 for(unsigned l_2deriv=0; l_2deriv<n_2deriv; l_2deriv++)
	  {
	   interpolated_d2fdxi2(i_field,l_2deriv) +=
	    internal_value * (*internal_d2basis_pt)(k_type, 0,l_2deriv);
	  } // End of loop over second derivatives l_2deriv
	}
      } // End of loop over basis types (k_type)
     //----------END OF INTERNAL DATA----------
    } // End of loop over fields (i_field)
   //========================END OF INTERPOLATION==============================
   
   // Copy our interpolated fields into u in the order we want and return it.
   // ([zdec] w first?)
   interpolated_u[0] = interpolated_f[2];         // w
   interpolated_u[1] = interpolated_dfdxi(2,0);   // dwdx
   interpolated_u[1] = interpolated_dfdxi(2,1);   // dwdy
   interpolated_u[1] = interpolated_d2fdxi2(2,0); // d2wdx2
   interpolated_u[1] = interpolated_d2fdxi2(2,1); // d2wdxdy
   interpolated_u[1] = interpolated_d2fdxi2(2,2); // d2wdy2
   interpolated_u[1] = interpolated_f[0];         // ux
   interpolated_u[1] = interpolated_f[1];         // uy
   interpolated_u[1] = interpolated_dfdxi(0,0);   // duxdx
   interpolated_u[1] = interpolated_dfdxi(0,1);   // duxdy
   interpolated_u[1] = interpolated_dfdxi(1,0);   // duydx
   interpolated_u[1] = interpolated_dfdxi(1,1);   // duydy
   return(interpolated_u);
  }

  /// \short Self-test: Return 0 for OK
  unsigned self_test();

 protected:
  /// Pure virtual interface to the basis for the in--plane displacement 
  virtual void shape_u(const Vector<double> &s,  Shape &psi) const=0;

  /// \short Shape/test functions and derivs w.r.t. to global coords at
  /// local coord. s; return  Jacobian of mapping
  virtual double d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,
							       Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
							       DShape &d2psi_dx2,DShape& d2psi_b_dx2,
							       Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx,
							       DShape &d2test_dx2,DShape& d2test_b_dx2) const=0;

  /// \short Shape/test functions and derivs w.r.t. to global coords at
  /// local coord. s; return  Jacobian of mapping
  virtual double dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s,
							     Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
							     Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx) const=0;

  /// \short Shape/test functions at local coordinate s
  virtual void shape_and_test_foeppl_von_karman(const Vector<double> &s,
						Shape &psi, Shape& psi_b, Shape &test, Shape& test_b) const=0;
 
  /// \short in--plane Shape/test functions at local coordinate s
  virtual double dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s,
								 Shape &psi,DShape &dpsidx, Shape &test,DShape &dtestdx) const=0;

  /// \short Compute element residual Vector only (if flag=and/or element
  /// Jacobian matrix
  virtual void fill_in_generic_residual_contribution_foeppl_von_karman(
								       Vector<double> &residuals, DenseMatrix<double> &jacobian,
								       const unsigned& flag);

  /// Pointer to pressure function:
  ScalarFctPt Pressure_fct_pt;

  /// Pointer to in plane forcing function (i.e. the shear force applied to
  /// the face)
  VectorFctPt In_plane_forcing_fct_pt;
 
  /// Pointer to swelling function:
  ScalarFctPt Swelling_fct_pt; 

  /// Pointer to Poisson ratio, which this element cannot modify
  const double* Nu_pt;

  /// Pointer to the dampening coefficient, which this element cannot modify
  const double* Mu_pt;

  /// Pointer to global eta
  const double *Eta_pt;

  /// Default value for physical constant: Poisson ratio. 
  static const double Default_Nu_Value;

  /// Default value for constant: dampening coefficient. 
  static const double Default_Mu_Value;

  /// Default eta value so that we use 'natural' nondim and have no h dependence. 
  static const double Default_Eta_Value;

  /// Pointer to error metric
  ErrorMetricFctPt Error_metric_fct_pt;

  /// Pointer to error metric when we want multiple errors
  MultipleErrorMetricFctPt Multiple_error_metric_fct_pt;


 protected:

  /// Pointer to precomputed matrix that associates shape functions to monomials
  DenseMatrix<double> *Association_matrix_pt;
 };


} //end namespace oomph
#endif

