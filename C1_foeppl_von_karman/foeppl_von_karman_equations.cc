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
// Non--inline functions for the foeppl_von_karman equations
#include "foeppl_von_karman.h"

namespace oomph
{
 /// \short Default to incompressible sheet Poisson ratio
 const double FoepplVonKarmanEquations::Default_Nu_Value=0.5;
 
 /// \short Default to no damping
 const double FoepplVonKarmanEquations::Default_Mu_Value=0.0;
 
 /// \short Default to 'natural' FvK nondimensionalisation which is h 
 /// independent (NB NOT the same as in JFM paper)
 const double FoepplVonKarmanEquations::Default_Eta_Value=1;
 
 
 //======================================================================
 /// Fill in generic residual and jacobian contribution 
 //======================================================================
 void  FoepplVonKarmanEquations::
 fill_in_generic_residual_contribution_foeppl_von_karman(Vector<double> &residuals,
							 DenseMatrix<double> &jacobian,
							 const unsigned& flag)
 {
  // The index of in-plane and out-of-plane unknowns [zdec] is this even needed?
  const unsigned u_index = 0;
  const unsigned w_index = 2;
  const unsigned n_u_fields = 2;
  const unsigned n_w_fields = 1;
  
  // Find the dimension of the element [zdec] will this ever not be 2?
  const unsigned dim = this->dim();
  // The number of first derivatives is the dimension of the element
  const unsigned n_deriv = dim;
  // The number of second derivatives is the triangle number of the dimension
  const unsigned n_2deriv = dim*(dim+1)/2;
  
  // Find out how many nodes there are for each field
  const unsigned n_u_node = nu_node();
  const unsigned n_w_node = nw_node();
  // Get the vector of nodes used for each field
  const Vector<unsigned> u_nodes = get_u_nodes();
  const Vector<unsigned> w_nodes = get_w_nodes();
  // Find out how many basis types there are at each node
  // [zdec] this assumes that all nodes have the same number of types and that
  //        node 0 is used by both bases
  const unsigned n_u_nodal_type = nu_type_at_node(0);
  const unsigned n_w_nodal_type = nw_type_at_node(0);
  // Find out how many basis types there are internally
  const unsigned n_u_internal_type = nu_internal_type();
  const unsigned n_w_internal_type = nw_internal_type();
  
  //Get the Poisson ratio of the plate
  const double nu=get_nu();
  
  //Get the virtual dampening coefficient
  const double mu=get_mu();
  
  // [zdec] Do we bother separating psi (basis) and test functions?
  // Local in-plane nodal basis and test functions
  Shape psi_n_u(n_u_node, n_u_nodal_type);
  Shape test_n_u(n_u_node, n_u_nodal_type);
  DShape dpsi_n_udxi(n_u_node, n_u_nodal_type, n_deriv);
  DShape dtest_n_udxi(n_u_node, n_u_nodal_type, n_deriv); 
  /* [zdec] Do we include the internal basis for generality
  // Local in-plane internal basis and test functions
  Shape psi_i_u(n_u_internal_type);
  Shape test_i_u(n_u_internal_type);
  DShape dpsi_i_udxi(n_u_internal_type, n_deriv);
  DShape dtest_i_udxi(n_u_internal_type, n_deriv);
  */
  
  // Local out-of-plane nodal basis and test funtions
  Shape psi_n_w(n_w_node, n_w_nodal_type);
  Shape test_n_w(n_w_node, n_w_nodal_type);
  DShape dpsi_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
  DShape dtest_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
  DShape d2psi_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
  DShape d2test_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
   // Local out-of-plane internal basis and test functions
  Shape psi_i_w(n_w_internal_type);
  Shape test_i_w(n_w_internal_type);
  DShape dpsi_i_wdxi(n_w_internal_type, n_deriv);
  DShape dtest_i_wdxi(n_w_internal_type, n_deriv);
  DShape d2psi_i_wdxi2(n_w_internal_type, n_2deriv);
  DShape d2test_i_wdxi2(n_w_internal_type, n_2deriv);
  
  // Set the value of n_intpt
  const unsigned n_intpt = this->integral_pt()->nweight();
  // Integers to store the local equation and unknown numbers
  int local_eqn=0;
  int local_unknown=0;
  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);

    //Allocate and initialise to zero
    Vector<double> interp_x(dim,0.0);
    Vector<double> s(dim);
    s[0] = this->integral_pt()->knot(ipt,0);
    s[1] = this->integral_pt()->knot(ipt,1);
    interpolated_x(s,interp_x);
    //Call the derivatives of the shape and test functions for the unknown
    // [zdec] Make this general for each unknown.
    double J =
     d2shape_and_d2test_eulerian_foeppl_von_karman(s,
						   psi_n_w,
						   psi_i_w,
						   dpsi_n_wdxi,
						   dpsi_i_wdxi,
						   d2psi_n_wdxi2,
						   d2psi_i_wdxi2,
						   test_n_w,
						   test_i_w,
						   dtest_n_wdxi,
						   dtest_i_wdxi,
						   d2test_n_wdxi2,
						   d2test_i_wdxi2);
    
    dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,
						    psi_n_u,
						    dpsi_n_udxi,
						    test_n_u,
						    dtest_n_udxi);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;

    //=======================START OF INTERPOLATION=============================
    //Create space for in-plane interpolated unknowns
    Vector<double> interpolated_u(n_u_fields);
    Vector<double> interpolated_dudt(n_u_fields);
    DenseMatrix<double> interpolated_dudxi(n_u_fields, n_deriv);
    // Create space for out-of-plane interpolated unknowns
    Vector<double> interpolated_w(n_w_fields);
    Vector<double> interpolated_dwdt(n_w_fields);
    DenseMatrix<double> interpolated_dwdxi(n_w_fields, n_deriv);
    DenseMatrix<double> interpolated_d2wdxi2(n_w_fields, n_2deriv);
    
    //---Interpolate the in-plane unknowns-----------------------------------
    // Loop over nodes used by in-plane fields
    for(unsigned j_node=0; j_node<n_u_node; j_node++)
     {
      // Turn on damping if this field requires it AND we are not doing a
      // steady solve
      TimeStepper* timestepper_pt = this->node_pt(j_node)->time_stepper_pt();
      bool damping = u_is_damped() && !(timestepper_pt->is_steady());
      unsigned n_time = 0;
      if(damping)
       {
	n_time=timestepper_pt->ntstorage();
       }
      // Get the j-th node used by in-plane fields
      unsigned j_node_local = u_nodes[j_node];
      // Get the number of types used at the j-th node by in-plane fields
      unsigned n_type = nw_type_at_node(j_node_local);
      // Loop over types
      for(unsigned k_type=0; k_type<n_type; k_type++)
       {
	// Loop over in-plane unknowns
	for(unsigned alpha=0; alpha<n_u_fields; alpha++)
	 {
	  // Get the nodal value of type k
	  double u_value =
	   get_u_alpha_value_at_node_of_type(alpha, j_node, k_type);
	  // Add nodal contribution of type k to the interpolated displacement
	  interpolated_u[alpha] += u_value * psi_n_u(j_node, k_type);
	  // Loop over first derivatives
	  for(unsigned l_deriv=0; l_deriv<n_deriv; l_deriv++)
	   {
	    // Add the nodal contribution of type k to the derivative of the displacement
	    interpolated_dudxi(alpha,l_deriv)
	     += u_value * dpsi_n_udxi(j_node, k_type, l_deriv);
	   }
	  // (if damped n_time>0) Loop over the history values
	 }
       }
     }
    
    // Interpolate the out-of-plane unknowns
    for(unsigned j; j<n_w_node; j++)
     {
      // Get the 'j'-th node used by w
      unsigned j_local = w_nodes[j];
     }
    
    //========================END OF INTERPOLATION==============================    

    /*
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
	unsigned j_node = local_nodes_for_field_i[j_nodei];
	// Get the number of basis types for field i at node j
	unsigned n_type = ntype_for_field_at_node(i_field,j_node);
	// Turn on damping if this field requires it AND we are not doing a
	// steady solve
	TimeStepper* timestepper_pt = this->node_pt(j_node)->time_stepper_pt();
	bool damping = field_is_damped(i_field) && !(timestepper_pt->is_steady());
	unsigned n_time = 0;
	if(damping)
	 {
	  n_time=timestepper_pt->ntstorage();
	 }
	
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
	  
	  // Add nodal contribution to time derivative (n_time=0 if not damped)
	  double nodal_dudt_value=0.0;
	  for(unsigned t=0; t<n_time; t++)
	   {
	    double nodal_history_value =
	     nodal_value_for_field_at_node_of_type(t, i_field, j_node, k_type);
	    nodal_dudt_value += nodal_history_value * timestepper_pt->weight(1,t);
	   } // End of loop over history values (t)
	  interpolated_dfdt[i_field] += nodal_dudt_value * (*basis_pt)(j_node,k_type);
	 
	 } // End of loop over basis types (k_type)
       } // End of loop over nodes used by field i (j_nodei)
      //-----------END OF NODAL DATA------------
      
      //-------------INTERNAL DATA--------------
      Shape* internal_basis_pt = psi_internal_field_pt[i_field];
      DShape* internal_dbasis_pt = dpsi_internal_field_pt[i_field];
      DShape* internal_d2basis_pt = d2psi_internal_field_pt[i_field];
      // Turn on damping if this field requires it AND we are not doing a
      // steady solve
      TimeStepper* timestepper_pt =
       this->internal_data_for_field_pt(i_field)->time_stepper_pt();
      bool damping =
       field_is_damped(i_field) && !(timestepper_pt->is_steady());
      unsigned n_time = 0;
      if(damping)
       {
	n_time=timestepper_pt->ntstorage();
       }
      
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
	
	// Add internal contribution to time derivative (n_time=0 if not damped)
	double internal_dudt_value=0.0;
	for(unsigned t=0; t<n_time; t++)
	 {
	  double internal_history_value =
	   internal_value_for_field_of_type(t, i_field, k_type);
	  internal_dudt_value += internal_history_value * timestepper_pt->weight(1,t);
	 } // End of loop over history values (t)
	interpolated_dfdt[i_field] += internal_dudt_value * (*internal_basis_pt)(k_type,0);
	 
       } // End of loop over basis types (k_type)
      //----------END OF INTERNAL DATA----------
     } // End of loop over fields (i_field)
    */
    
    
    //======================[zdec] residuals and jacobians======================
    //Get the swelling function
    //--------------------
    double c_swell(0.0);
    get_swelling_foeppl_von_karman(interp_x,c_swell);
    
    // Get the stress
    DenseMatrix<double> sigma(2,2,0.0);
    get_sigma(sigma, interpolated_duidxj, interpolated_dwdxi, c_swell);
    
    //Get pressure function
    //-------------------
    double pressure(0.0);
    Vector<double> pressure_gradient(2,0.0);
    get_pressure_foeppl_von_karman(ipt,interp_x,pressure);
    get_in_plane_forcing_foeppl_von_karman(ipt,interp_x,pressure_gradient);

    // Loop over the nodal test functions
    for(unsigned l=0;l<n_w_node;l++)
     {
      for(unsigned k=0;k<n_basis_type;k++)
       {
	//Get the local equation
	local_eqn = this->nodal_local_eqn(l,k+2);
	//IF it's not a boundary condition
	if(local_eqn >= 0)
	 {
	  // Add virtual time derivative term for dampened buckling
	  residuals[local_eqn] += mu*interpolated_dwdt[0]*test_w(l,k)*W;
	  // Add body force/pressure term here
	  residuals[local_eqn] -= pressure*test_w(l,k)*W;
	  for(unsigned alpha=0;alpha<2;++alpha)
	   {
	    for(unsigned beta=0; beta<2;++beta)
	     {
	      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
	      residuals[local_eqn] += (1-nu)
	       *interpolated_d2wdxi2(0,alpha+beta)*d2test_wdxi2(l,k,alpha+beta)*W;
	      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
	      residuals[local_eqn] += (nu)
	       *interpolated_d2wdxi2(0,beta+beta)*d2test_wdxi2(l,k,alpha+alpha)*W;
	      // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
	      residuals[local_eqn] += eta()*sigma(alpha,beta)
	       *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)*W;
	     }
	   }
	  // Calculate the jacobian
	  //-----------------------
	  if(flag)
	   {
	    //Loop over the in--plane unknowns again
	    for(unsigned ll=0;ll<n_u_node;ll++)
	     {
	      // Loop over displacements
	      for(unsigned gamma=0;gamma<2;++gamma)
	       {
		local_unknown = this->nodal_local_eqn(ll,gamma);
		// If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Loop over displacements
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    // Loop over dimensions 
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Nonlinear Terms 
		      if(alpha==beta)
		       {
			jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)* 
			 interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(alpha==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
			 *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(beta==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
			 *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)/(1-nu*nu)*W;
		       }
		     }
		   }
		 }
	       }
	     }

	    //Loop over the test functions again
	    for(unsigned l2=0;l2<n_w_node;l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_basis_type;k2++)
	       {
		local_unknown = this->nodal_local_eqn(l2,k2+2);
		// If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi_w(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test_w(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // Linear Terms
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_wdxi2(l2,k2,alpha+beta)*d2test_wdxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_wdxi2(l2,k2,beta+beta)*d2test_wdxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
		       *dpsi_wdxi(l2,k2,alpha)*dtest_wdxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_wdxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_wdxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_wdxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_wdxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)*W;//DGR: eta^2 surely
		     }
		   }
		 }
	       }
	     }
	    //Loop over the internal test functions
	    for(unsigned l2=0;l2<nbubble_basis();l2++)
	     {
	      for(unsigned k2=0;k2<exactly_uno;k2++)
	       {
		local_unknown = local_internal_equation(w_index,
							l2);//local_w_bubble_equation(l2,k2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi_b(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test_w(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Linear Terms 
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_wdxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_wdxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_wdxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_wdxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_wdxi(l,k,beta)*W;//DGR: eta^2 surely
		     }
		   }
		 }
	       }
	     }
	   } // End of flag
	 }
       }
     }

    // Loop over the internal test functions
    for(unsigned l=0;l<nbubble_basis();l++)
     {
      for(unsigned k=0;k<exactly_uno;k++)
       {
	//Get the local equation
	local_eqn = local_internal_equation(w_index,
					    l);//local_w_bubble_equation(l,k);
	//IF it's not a boundary condition
	if(local_eqn >= 0)
	 {
	  // Add virtual time derivative term for dampened buckling
	  residuals[local_eqn] += mu*interpolated_dwdt[0]*test_b(l,k)*W;     
	  // Add body force/pressure term here
	  residuals[local_eqn] -= pressure*test_b(l,k)*W;

	  for(unsigned alpha=0;alpha<2;++alpha)
	   {
	    for(unsigned beta=0; beta<2;++beta)
	     {
	      // The Linear Terms 
	      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
	      residuals[local_eqn] += (1-nu)
	       *interpolated_d2wdxi2(0,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
	      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
	      residuals[local_eqn] += (nu)
	       *interpolated_d2wdxi2(0,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
	      // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
	      residuals[local_eqn] += eta()*sigma(alpha,beta)
	       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
	     }
	   }
	  // Calculate the jacobian
	  //-----------------------
	  if(flag)
	   {
	    //Loop over the in--plane unknowns again
	    for(unsigned ll=0;ll<n_u_node;ll++)
	     {
	      // Loop over displacements
	      for(unsigned gamma=0;gamma<2;++gamma)
	       {
		local_unknown = this->nodal_local_eqn(ll,gamma);
		// If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Loop over displacements
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    // Loop over dimensions 
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Nonlinear Terms 
		      if(alpha==beta)
		       {
			jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)*
			 interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(alpha==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
			 *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(beta==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
			 *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		     }
		   }
		 }
	       }

	     }
	    //Loop over the test functions again
	    for(unsigned l2=0;l2<n_w_node;l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_basis_type;k2++)
	       {
		local_unknown = this->nodal_local_eqn(l2,k2+2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi_w(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test_b(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Linear Terms 
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_wdxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_wdxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
		       *dpsi_wdxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_wdxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_wdxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_wdxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		     }
		   }
		 }
	       }
	     }
	    //Loop over the test functions again
	    for(unsigned l2=0;l2<nbubble_basis();l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<exactly_uno;k2++)
	       {
		local_unknown = local_internal_equation(w_index,
							l2);//local_w_bubble_equation(l2,k2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi_b(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test_b(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Linear Terms 
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		     }
		   }
		 }
	       }
	     }
	   } // End of flag
	 }
       }
     } //end loop over nodes
    for(unsigned l=0;l<n_u_node;++l)
     {
      // Now loop over displacement equations
      for(unsigned alpha=0; alpha<2 ; ++alpha)
       {
	//Get the local equation
	local_eqn = this->nodal_local_eqn(l,alpha);
	//IF it's not a boundary condition
	if(local_eqn >= 0)
	 {
	  // Now add on the stress terms
	  for(unsigned beta=0; beta<2;++beta)
	   {
	    // Variations in u_alpha
	    residuals[local_eqn] += sigma(alpha,beta)*dtest_udxi(l,beta)*W;
	   }
	  // Add the forcing term
	  residuals[local_eqn] += pressure_gradient[alpha]*psi_u(l)*W;

	  // Now loop over Jacobian
	  if(flag)
	   {
	    for(unsigned ll=0;ll<n_u_node;++ll)
	     {
	      // Now loop over displacement equations
	      for(unsigned gamma=0; gamma<2 ; ++gamma)
	       {
		//Get the local equation
		local_unknown = this->nodal_local_eqn(ll,gamma);
		// If the unknown is a degree of freedom 
		if(local_unknown>=0)
		 {
		  // Now add on the stress terms
		  for(unsigned beta=0; beta<2;++beta)
		   {
		    // The Linear Terms 
		    if(alpha==beta)
		     {
		      jacobian(local_eqn,local_unknown) += nu*dpsi_udxi(ll,gamma)
		       *dtest_udxi(l,beta)/(1-nu*nu)*W;
		     }
		    if(alpha==gamma)
		     {
		      jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,beta)
		       *dtest_udxi(l,beta)/(1-nu*nu)*W;
		     }
		    if(beta==gamma)
		     {
		      jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,alpha)
		       *dtest_udxi(l,beta)/(1-nu*nu)*W;
		     }
		   }// End loop beta
		 }// End if local unknown
	       }// End loop gamma dof
	     }// End loop nodal dofs
       
	    //Loop over the w nodal dofs
	    for(unsigned l2=0;l2<n_w_node;l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_basis_type;k2++)
	       {
		//If at a non-zero degree of freedom add in the entry
		//Get the local equation
		local_unknown = this->nodal_local_eqn(l2,k2+2);
		// If the unknown is a degree of freedom 
		if(local_unknown>=0)
		 {
		  // Now add on the stress terms
		  for(unsigned beta=0; beta<2;++beta)
		   {
		    // The NonLinear Terms 
		    jacobian(local_eqn,local_unknown) += nu*dpsi_wdxi(l2,k2,beta)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_wdxi(l2,k2,alpha)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_wdxi(l2,k2,beta)
		     *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		   }// End loop beta
		 }// End if local unknown
	       }// End loop position type dof
	     }// End loop w nodal dofs
       
	    //Loop over the w internal dofs
	    for(unsigned l2=0;l2<nbubble_basis();l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<exactly_uno;k2++)
	       {
		//If at a non-zero degree of freedom add in the entry
		//Get the local equation
		local_unknown = local_internal_equation(w_index,
							l2);// local_w_bubble_equation(l2,k2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Now add on the stress terms
		  for(unsigned beta=0; beta<2;++beta)
		   {
		    // The NonLinear Terms 
		    jacobian(local_eqn,local_unknown) += nu*dpsi_b_dxi(l2,k2,beta)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_b_dxi(l2,k2,alpha)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_b_dxi(l2,k2,beta)
		     *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		   } // End loop beta
		 } // End if local_unknown
	       } // End loop position type
	     } // End loop internal dofs
	   } // End if flag
	 } // End if local eqn
       } // End loop over alpha
     } // End loop over nodes
   } // End of loop over integration points
 }



 // Define the energy calculation function here
 Vector<double> FoepplVonKarmanEquations::element_elastic_and_kinetic_energy()
 {
  double element_elastic_energy = 0.0;
  double element_kinetic_energy = 0.0;

  //Find the number of unknown fields (should be 3 for FvK)
  const unsigned n_field = nfield();
  const unsigned u_index = 0; // [zdec] fix these hard coded indices
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

  //Get the Poisson ratio of the plate
  const double nu=get_nu();
 
  //Get the virtual dampening coefficient
  const double mu=get_mu();

  // [zdec] Generalise the retrieval of basis and test functions.
  // Set up memory for the shape and test functions
  Shape psi_u(n_u_node), test_u(n_u_node);
  DShape dpsi_udxi(n_u_node,n_deriv), dtest_udxi(n_u_node,n_deriv);
 
  //Local c1-shape funtion
  Shape psi_w(n_w_node,n_basis_type),test_w(n_w_node,n_basis_type),
   psi_b(n_b_type,exactly_uno),test_b(n_b_type,exactly_uno);
 
  DShape dpsi_wdxi(n_w_node,n_basis_type,n_deriv),dtest_wdxi(n_w_node,n_basis_type,n_deriv),
   dpsi_b_dxi(n_b_type,exactly_uno,n_deriv),dtest_b_dxi(n_b_type,exactly_uno,n_deriv),
   d2psi_wdxi2(n_w_node,n_basis_type,n_2deriv), d2test_wdxi2(n_w_node,n_basis_type,n_2deriv),
   d2psi_b_dxi2(n_b_type,exactly_uno,n_2deriv), d2test_b_dxi2(n_b_type,exactly_uno,n_2deriv);

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

  //Set the value of n_intpt
  const unsigned n_intpt = this->integral_pt()->nweight();

  
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);

    // [zdec] The interpolated field values
    Vector<double> interpolated_f(n_field,0.0);
    Vector<double> interpolated_dfdt(n_field,0.0); // Do we size this based on damped fields?
    DenseMatrix<double> interpolated_dfdxi(n_field,dim,0.0);
    DenseMatrix<double> interpolated_d2fdxi2(n_field,dim*(dim+1)/2,0.0);
 
    //Allocate and initialise to zero
    Vector<double> interp_x(dim,0.0);
    Vector<double> s(dim);
    s[0] = this->integral_pt()->knot(ipt,0);
    s[1] = this->integral_pt()->knot(ipt,1);
    interpolated_x(s,interp_x);
    //Call the derivatives of the shape and test functions for the unknown
    double J =
     d2shape_and_d2test_eulerian_foeppl_von_karman(s, psi_w, psi_b,
						   dpsi_wdxi, dpsi_b_dxi,
						   d2psi_wdxi2, d2psi_b_dxi2,
						   test_w, test_b,
						   dtest_wdxi, dtest_b_dxi,
						   d2test_wdxi2, d2test_b_dxi2);

    dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s, psi_u, dpsi_udxi,
						    test_u, dtest_udxi);
    
    //Premultiply the weights and the Jacobian
    double W = w*J;

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
	// Turn on damping if this field requires it AND we are not doing a
	// steady solve
	TimeStepper* timestepper_pt = this->node_pt(j_node)->time_stepper_pt();
	bool damping = field_is_damped(i_field) && !(timestepper_pt->is_steady());
	unsigned n_time = 0;
	if(damping)
	 {
	  n_time=timestepper_pt->ntstorage();
	 }
	
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
	  
	  // Add nodal contribution to time derivative (n_time=0 if not damped)
	  double nodal_dudt_value=0.0;
	  for(unsigned t=0; t<n_time; t++)
	   {
	    double nodal_history_value =
	     nodal_value_for_field_at_node_of_type(t, i_field, j_node, k_type);
	    nodal_dudt_value += nodal_history_value * timestepper_pt->weight(1,t);
	   } // End of loop over history values (t)
	  interpolated_dfdt[i_field] += nodal_dudt_value * (*basis_pt)(j_node,k_type);
	 
	 } // End of loop over basis types (k_type)
       } // End of loop over nodes used by field i (j_nodei)
      //-----------END OF NODAL DATA------------
      
      //-------------INTERNAL DATA--------------
      Shape* internal_basis_pt = psi_internal_field_pt[i_field];
      DShape* internal_dbasis_pt = dpsi_internal_field_pt[i_field];
      DShape* internal_d2basis_pt = d2psi_internal_field_pt[i_field];
      // Turn on damping if this field requires it AND we are not doing a
      // steady solve
      TimeStepper* timestepper_pt =
       this->internal_data_for_field_pt(i_field)->time_stepper_pt();
      bool damping =
       field_is_damped(i_field) && !(timestepper_pt->is_steady());
      unsigned n_time = 0;
      if(damping)
       {
	n_time=timestepper_pt->ntstorage();
       }

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
	
	// Add internal contribution to time derivative (n_time=0 if not damped)
	double internal_dudt_value=0.0;
	for(unsigned t=0; t<n_time; t++)
	 {
	  double internal_history_value =
	   internal_value_for_field_of_type(t, i_field, k_type);
	  internal_dudt_value += internal_history_value * timestepper_pt->weight(1,t);
	 } // End of loop over history values (t)
	interpolated_dfdt[i_field] += internal_dudt_value * (*internal_basis_pt)(k_type,0);
	 
       } // End of loop over basis types (k_type)
      //----------END OF INTERNAL DATA----------
     } // End of loop over fields (i_field)
    
    //========================END OF INTERPOLATION==============================
    //=============[zdec] move interpolated variables into u and w==============
    // Create containers for out-of-plane unknowns
    Vector<double> interpolated_w(1);
    Vector<double> interpolated_dwdt(1);
    DenseMatrix<double> interpolated_dwdxi(1,n_deriv);
    DenseMatrix<double> interpolated_d2wdxi2(1,n_2deriv);
    // Copy out-of-plane field into w container
    interpolated_w[0]    = interpolated_f[w_index];
    interpolated_dwdt[0] = interpolated_dfdt[w_index];
    for(unsigned i_deriv=0; i_deriv<n_deriv; i_deriv++)
     {
      interpolated_dwdxi(0,i_deriv) = interpolated_dfdxi(w_index,i_deriv);
     }
    for(unsigned i_2deriv=0; i_2deriv<n_2deriv; i_2deriv++)
     {
      interpolated_d2wdxi2(0,i_2deriv) = interpolated_d2fdxi2(w_index,i_2deriv);
     }
   
    //Create space for in-plane unknowns
    Vector<double> interpolated_ui(2);
    DenseMatrix<double> interpolated_duidxj(2,n_deriv);
    for(unsigned alpha=u_index; alpha<u_index+2; alpha++)
     {
      interpolated_ui[alpha] = interpolated_f[alpha];
      for(unsigned i_deriv=0; i_deriv<n_deriv; i_deriv++)
       {
	interpolated_duidxj(alpha,i_deriv) =
	 interpolated_dfdxi(u_index+alpha,i_deriv);
	oomph_info << alpha << " "
		   << i_deriv << " "
		   << interpolated_duidxj(alpha,i_deriv) << std::endl;
       }
     }
    //=======================end of interpolation transfer======================


    //----------------------------------------------------------------------
    // Add the elastic energy at ipt.

    // (Nondimensional) Lame parameters
    double lame2 = 1.0 / (2.0 * (1.0+nu));
    double lame1 = 2.0*nu*lame2 / (1.0 - 2.0*nu);

    //Get the swelling function
    double c_swell(0.0);
    get_swelling_foeppl_von_karman(interp_x,c_swell);
    
    // Truncated Green Lagrange strain tensor
    DenseMatrix<double> epsilon(dim,dim,0.0);
    for(unsigned alpha=0;alpha<dim;++alpha)
     {
      epsilon(alpha,alpha) -= c_swell;
      for(unsigned beta=0;beta<dim;++beta)
       {
	// Truncated Green Lagrange strain tensor
	epsilon(alpha,beta) += 0.5* interpolated_duidxj(alpha,beta)
	 + 0.5*interpolated_duidxj(beta,alpha)
	 + 0.5*interpolated_dwdxi(0,alpha)*interpolated_dwdxi(0,beta);
       }
     }

    // Add the plane strain energy contribution at ipt.
    for(unsigned alpha=0;alpha<dim;++alpha)
     {
      element_elastic_energy +=
       0.5*lame1*epsilon(alpha,alpha)*epsilon(alpha,alpha) * W;
      for(unsigned beta=0;beta<dim;++beta)
       {
	element_elastic_energy +=
	 eta()*lame2*epsilon(alpha,beta)*epsilon(alpha,beta) * W;
       }
     }
    
    // Add the bending energy contribution at ipt.
    element_elastic_energy +=
     0.5 * (
	    (interpolated_d2wdxi2(0,0)+interpolated_d2wdxi2(0,2))
	    *(interpolated_d2wdxi2(0,0)+interpolated_d2wdxi2(0,2))
	    + 2.0*(1.0-nu) * (
			      interpolated_d2wdxi2(0,1)*interpolated_d2wdxi2(0,1)
			      - interpolated_d2wdxi2(0,0)*interpolated_d2wdxi2(0,2)
			      )
	    ) * W;

    //----------------------------------------------------------------------
    // Add the kinetic energy at ipt.
    element_kinetic_energy += 0.5*mu*interpolated_dwdt[0]*interpolated_dwdt[0] * W;

   }
  Vector<double> energies(2,0.0);
  // First entry is the elastic energy and second is kinetic.
  energies[0] = element_elastic_energy;
  energies[1] = element_kinetic_energy;
  
  return energies;
 }




 //======================================================================
 /// Self-test:  Return 0 for OK
 //======================================================================
 unsigned  FoepplVonKarmanEquations::self_test()
 {

  bool passed=true;

  // Check lower-level stuff
  if (FiniteElement::self_test()!=0)
   {
    passed=false;
   }

  // Return verdict
  if (passed)
   {
    return 0;
   }
  else
   {
    return 1;
   }

 }

 //======================================================================
 /// Output function:
 ///
 ///   x,y,u,(epsilon),(sigma),(sigma_evals,sigma_evecs)
 ///
 /// nplot points in each coordinate direction
 //======================================================================
 void  FoepplVonKarmanEquations::output(std::ostream &outfile,
					const unsigned &nplot)
 {
  unsigned dim = this->dim();
 
  //Vector of local coordinates
  Vector<double> s(dim),x(dim);

  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);

  // Storage for variables
  double c_swell(0.0);
  Vector<double> u;
  DenseMatrix<double> interpolated_dwdxj(1,dim,0.0);
  DenseMatrix<double> interpolated_duidxj(2,dim,0.0);
  DenseMatrix<double> epsilon(2,2,0.0);
  DenseMatrix<double> sigma(2,2,0.0);
  DenseMatrix<double> sigma_eigenvecs(2,2,0.0);
  Vector<double> sigma_eigenvals(2,0.0);
  unsigned num_plot_points=this->nplot_points(nplot);
  Vector<double> r(3);

  // Loop over plot points
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    // Get local and global coordinates of plot point
    this->get_s_plot(iplot,nplot,s);
    interpolated_x(s,x);

    // Get interpolated unknowns
    u = interpolated_u_foeppl_von_karman(s);

    // Get degree of swelling for use in the strain tensor
    this->get_swelling_foeppl_von_karman(x,c_swell);

    // TODO: make the output customisable with flags e.g. output_axial_strain,
    //       output_principal_stress
    // TODO: make the indexing below more general and not hard coded.
    // Copy gradients from u into interpolated gradient matrices...
    interpolated_dwdxj(0,0) = u[1]; //dwdx1
    interpolated_dwdxj(0,1) = u[2]; //dwdx2
    interpolated_duidxj(0,0)= u[8]; //du1dx1
    interpolated_duidxj(0,1)= u[9]; //du1dx2
    interpolated_duidxj(1,0)= u[10]; //du2dx1
    interpolated_duidxj(1,1)= u[11]; //du2dx2
    // ...which are used to retrieve the strain tensor epsilon.
    get_epsilon(epsilon, interpolated_duidxj,interpolated_dwdxj, c_swell);

    // Use epsilon to find the stress sigma.
    get_sigma_from_epsilon(sigma, epsilon);

    // Get the principal values and the corresponding directions of stress.
    get_principal_stresses(sigma, sigma_eigenvals, sigma_eigenvecs);
   

    for(unsigned i=0;i<this->dim();i++)
     {
      outfile << x[i] << " ";
     }

    // Loop for variables
    for(Vector<double>::iterator it=u.begin();it!=u.end();++it)
     {
      outfile << *it << " " ;
     }
   
    // Output axial strains
    outfile << epsilon(0,0) << " " << epsilon(0,1) << " " << epsilon(1,1) << " ";

    // Output axial stress
    outfile << sigma(0,0) << " " << sigma(0,1) << " " << sigma(1,1) << " ";
   
    // Output principal stresses
    outfile << sigma_eigenvals[0] << " " << sigma_eigenvals[1] << " ";
   
    // Output principal stress directions
    outfile << sigma_eigenvecs(0,0) << " " << sigma_eigenvecs(1,0) << " "
	    << sigma_eigenvecs(0,1) << " " << sigma_eigenvecs(1,1) << " ";

    // End output line
    outfile << std::endl;
   }
  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);
 }


 // TODO: update further output functions for full output
 //======================================================================
 /// C-style output function: 
 ///
 ///   x,y,u   or    x,y,z,u
 ///
 /// nplot points in each coordinate direction
 //======================================================================
 void  FoepplVonKarmanEquations::output(FILE* file_pt,
					const unsigned &nplot)
 {
  //Vector of local coordinates
  Vector<double> s(this->dim()), x(this->dim());;

  // Tecplot header info
  fprintf(file_pt,"%s",this->tecplot_zone_string(nplot).c_str());

  // Loop over plot points
  Vector<double> u(this->required_nvalue(0),0.0);
  unsigned num_plot_points=this->nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    // Get local coordinates of plot point
    this->get_s_plot(iplot,nplot,s);

    // Get x position as Vector
    interpolated_x(s,x);

    for(unsigned i=0;i<this->dim();i++)
     {
      fprintf(file_pt,"%g ",x[i]);
     }
    u = interpolated_u_foeppl_von_karman(s);
    fprintf(file_pt,"%g \n",u[0]);
   }

  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(file_pt,nplot);
 }


 // [zdec] This is also called to output the pressure.
 // [TODO] Needs renaming.
 //======================================================================
 /// Output exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
 //======================================================================
 void FoepplVonKarmanEquations::output_fct(std::ostream &outfile,
					   const unsigned &nplot,
					   FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
 {
  //Vector of local coordinates
  Vector<double> s(this->dim());

  // Vector for coordintes
  Vector<double> x(this->dim());

  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);

  // Exact solution Vector (here a scalar)
  //
  // [zdec] Why was this->required_nvalue(0) called here? Caused range checking
  // errors due to required_nvalue(0) being 8 but exact_soln_pt resizing to a
  // vector of length one. Replaced with 1 for now as the above comment states
  // we already know this function to be scalar.
  // 
  // Vector<double> exact_soln(this->required_nvalue(0),0.0);
  unsigned nvalue=1;
  Vector<double> exact_soln(nvalue,0.0);
 
  // Loop over plot points
  unsigned num_plot_points=this->nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {

    // Get local coordinates of plot point
    this->get_s_plot(iplot,nplot,s);

    // Get x position as Vector
    interpolated_x(s,x);

    // Get exact solution at this point
    (*exact_soln_pt)(x,exact_soln);

    //Output x,y,...,u_exact
    for(unsigned i=0;i<this->dim();i++)
     {
      outfile << x[i] << " ";
     }
    // Loop over variables
    for(unsigned j=0;j<nvalue;j++)
     {
      outfile << exact_soln[j] << " ";
     }
    outfile <<  std::endl;
   }

  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);
 }


 //======================================================================
 /// Validate against exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
 /// HERE THIS MAY BE SUPERFLUOUS NOW
 //======================================================================
 void FoepplVonKarmanEquations::compute_error_in_deflection(std::ostream &outfile,
							    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
							    double& error, double& norm)
 {
  // Initialise
  error=0.0;
  norm=0.0;
  //Find out how many nodes there are
  // const unsigned n_u_node = this->nnode();
  const unsigned n_w_node = nnode_outofplane();
  //Find out how many bubble nodes there are
  const unsigned n_b_type = nbubble_basis();
  //Find out how many nodes positional dofs there are
  unsigned n_basis_type = nnodal_basis_type();
  // Find the internal dofs
  const unsigned exactly_uno = nbubble_basis_type();

  //Vector of local coordinates
  Vector<double> s(this->dim());

  // Vector for coordintes
  Vector<double> x(this->dim());

  //Set the value of n_intpt
  unsigned n_intpt = this->integral_pt()->nweight();

  // Exact solution Vector (here a scalar)
  Vector<double> exact_soln(this->required_nvalue(0),0.0);

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {

    //Assign values of s
    for(unsigned i=0;i<this->dim();i++)
     {
      s[i] = this->integral_pt()->knot(ipt,i);
     }

    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);

    // Get jacobian of mapping
    double J;
    {
     //Local c1-shape funtion
     Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
      psi_b(n_b_type,exactly_uno),test_b(n_b_type,exactly_uno);
   

     DShape dpsi_dxi(n_w_node,n_basis_type,this->dim()),dtest_dxi(n_w_node,n_basis_type,this->dim()),
      dpsi_b_dxi(n_b_type,exactly_uno,this->dim()),dtest_b_dxi(n_b_type,exactly_uno,this->dim()),
      d2psi_dxi2(n_w_node,n_basis_type,3), d2test_dxi2(n_w_node,n_basis_type,3),
      d2psi_b_dxi2(n_b_type,exactly_uno,3), d2test_b_dxi2(n_b_type,exactly_uno,3);

     J=this-> d2shape_and_d2test_eulerian_foeppl_von_karman(s,
							    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
							    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);
    }
    //Premultiply the weights and the Jacobian
    double W = w*J;
    // double Wlin = w*Jlin;

    // Get x position as Vector
    interpolated_x(s,x);

    // Get FE function value
    Vector<double> u_fe(this->required_nvalue(0),0.0);
    u_fe = interpolated_u_foeppl_von_karman(s);

    // Get exact solution at this point
    (*exact_soln_pt)(x,exact_soln);

    //Output x,y,...,error
    for(unsigned i=0;i<this->dim();i++)
     {
      outfile << x[i] << " ";
     }
    for(unsigned ii=0;ii<this->required_nvalue(0);ii++)
     {
      outfile << exact_soln[ii] << " " << exact_soln[ii]-u_fe[ii] << " ";
     }
    outfile << std::endl;

    // Loop over variables
    double tmp1 = 0.0, tmp2 =0.0;
    for(unsigned ii=0;ii<1;ii++)
     {
      // Add to error and norm
      tmp1 = (exact_soln[ii]*exact_soln[ii]*W);
      tmp2 = ((exact_soln[ii]-u_fe[ii])*(exact_soln[ii]-u_fe[ii])*W);
      norm += tmp1;
      error += tmp2;
     }
   } //End of loop over integration pts
 }
 //HERE OVERLOAD COMPUTE_ERROR SO WE CAN GET A VECTOR OF NORMS


 //======================================================================
 /// Validate against exact solution
 ///
 /// Solution is provided via function pointer.
 /// Plot error at a given number of plot points.
 ///
 //======================================================================
 void FoepplVonKarmanEquations::compute_error(std::ostream &outfile,
					      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
					      double& error, double& norm)
 {
  // Initialise
  error=0.0;
  norm=0.0;

  //Find out how many nodes there are
  // const unsigned n_u_node = this->nnode();
  const unsigned n_w_node = nnode_outofplane(); 
  //Find out how many bubble nodes there are
  const unsigned n_b_type = nbubble_basis();
  //Find out how many nodes positional dofs there are
  unsigned n_basis_type = nnodal_basis_type();
  // Find the internal dofs
  const unsigned exactly_uno = nbubble_basis_type();

  //Vector of local coordinates
  Vector<double> s(this->dim());

  // Vector for coordintes
  Vector<double> x(this->dim());

  //Set the value of n_intpt
  unsigned n_intpt = this->integral_pt()->nweight();

  // Tecplot
  //outfile << "ZONE" << std::endl;

  // Exact solution Vector (here a scalar)
  Vector<double> exact_soln(this->required_nvalue(0),0.0);

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {

    //Assign values of s
    for(unsigned i=0;i<this->dim();i++)
     {
      s[i] = this->integral_pt()->knot(ipt,i);
     }

    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);

    // Get jacobian of mapping
    double J;
    // double Jlin = this->J_eulerian1(s);// Nope
    {
     //Local c1-shape funtion
     Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
      psi_b(n_b_type,exactly_uno),test_b(n_b_type,exactly_uno);
   

     DShape dpsi_dxi(n_w_node,n_basis_type,this->dim()),dtest_dxi(n_w_node,n_basis_type,this->dim()),
      dpsi_b_dxi(n_b_type,exactly_uno,this->dim()),dtest_b_dxi(n_b_type,exactly_uno,this->dim()),
      d2psi_dxi2(n_w_node,n_basis_type,3), d2test_dxi2(n_w_node,n_basis_type,3),
      d2psi_b_dxi2(n_b_type,exactly_uno,3), d2test_b_dxi2(n_b_type,exactly_uno,3);


     J=this-> d2shape_and_d2test_eulerian_foeppl_von_karman(s,
							    psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
							    test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);
    }
    //Premultiply the weights and the Jacobian
    double W = w*J;
    // double Wlin = w*Jlin;

    // Get x position as Vector
    interpolated_x(s,x);

    // Get FE function value
    Vector<double> u_fe(this->required_nvalue(0),0.0);
    u_fe = interpolated_u_foeppl_von_karman(s);

    // Get exact solution at this point
    (*exact_soln_pt)(x,exact_soln);

    // Loop over variables
    double tmp1 = 0.0, tmp2 =0.0;
    if(Error_metric_fct_pt==0)
     {
      //Output x,y,...,error
      for(unsigned i=0;i<this->dim();i++)
       {
	outfile << x[i] << " ";
       }
      for(unsigned ii=0;ii<this->required_nvalue(0);ii++)
       {
	outfile << exact_soln[ii] << " " << exact_soln[ii]-u_fe[ii] << " ";
       }
      outfile << std::endl;

      for(unsigned ii=0;ii<1;ii++)
       {
	// Add to error and norm
	tmp1 = (exact_soln[ii]*exact_soln[ii]);
	tmp2 = ((exact_soln[ii]-u_fe[ii])*(exact_soln[ii]-u_fe[ii]));
       }
     }
    else
     {
      // Get the metric
      (*Error_metric_fct_pt)(x,u_fe,exact_soln,tmp2,tmp1);
     }
    norm += tmp1*W;
    error += tmp2*W;
   } //End of loop over integration pts
 }
}
