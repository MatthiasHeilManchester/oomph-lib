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

#ifndef OOMPH_BIHARMONIC_CURVED_BELL_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_CURVED_BELL_ELEMENTS_HEADER

#include "foeppl_von_karman_equations.h"
#include "../C1_basis/Bell_element_basis.h"
#include "../C1_basis/C1_curved_elements.h"
#include "../C1_basis/my_geom_object.h"
#include "../C1_basis/SubparametricTElement.h"

namespace oomph
{
 //=============================================================================
 /// FoepplVonKarmanC1CurvedBellElement elements are a subparametric scheme
 /// with  linear Lagrange interpolation for approximating the geometry and
 /// the C1-functions for approximating variables.
 //=============================================================================
 /* These elements use TElement<2,NNODE_1D> with Lagrangian bases to interpolate
 // the in-plane unknowns u_x, u_y although with sub-parametric simplex shape
 // interpolation.
 //
 // The Bell Hermite basis is used to interpolate the out-of-plane unknown w for
 // straight sided elements and is upgraded to Bernadou basis if the element is
 // on a curved boundary. This upgrade similarly corrects the shape
 // interpolation.
 //
 // As a result, all nodes contain the two in-plane Lagrangian dofs, vertex
 // nodes (i.e. nodes 0,1,2) each contain an additional six out-of-plane Hermite
 // dofs.
 //
 //  e.g FeopplVonKarmanC1CurvedBellElement<4>
 //                        | 0 | 1 | 2 | 3  | 4  |  5   |  6   |  7   |
 //          o <---------- |u_x|u_y| w |dwdx|dwdy|d2wdx2|d2wdxy|d2wdy2|
 //         / \
 //        x   x <-------- |u_x|u_y|
 //       /     \
 //      x   x   x  
 //     /         \
 //    o---x---x---o
 //                    
 */
 template <unsigned NNODE_1D>
 class FoepplVonKarmanC1CurvedBellElement :
  public virtual CurvableBellElement<NNODE_1D>,
  public virtual FoepplVonKarmanEquations
 {
 public:
  /// \short Function pointer to pressure function fct(x,f(x)) --
  /// x is a Vector!
  typedef void (*DisplacementFctPt)(const Vector<double>& x, double& f);

  /// \short function to pin all deflection dofs
  void pin_all_deflection_dofs() const;

  /// \short function to pin particular out--of--plane displacement dofs
  void fix_out_of_plane_displacement_dof(const unsigned& dof_number,
					 const unsigned& b,
					 const DisplacementFctPt& w);

  /// \short function to pin particular in--plane displacement dofs
  void fix_in_plane_displacement_dof(const unsigned& dof_number,
				     const unsigned& b,
				     const DisplacementFctPt& u);

  /// \short Function pointer to basis vectors function which sets  basis vectors
  /// b1 and b2 (which are in general functions of x)
  typedef void (*BasisVectorsFctPt) (const Vector<double>& x,
				     Vector<double>& b1,
				     Vector<double>& b2,
				     DenseMatrix<double>& Db1,
				     DenseMatrix<double>& Db2);

  /// \short Get the zeta coordinate
  inline void interpolated_zeta(const Vector<double> &s
				,Vector<double> &zeta) const // Higher
  {
   // If there is a macro element use it
   if(this->Macro_elem_pt!=0) {this->get_x_from_macro_element(s,zeta);}
   //Otherwise interpolate zeta_nodal using the shape functions
   else {interpolated_x(s,zeta);}
  }

  /// Lagrange interpolated shape for in--plane displacements
  void shape_u(const Vector<double> &s, Shape &psi) const;
  
  /// Lagrange interpolated d_shape for in--plane displacements
  void dshape_u_local(const Vector<double> &s,
		      Shape &psi, DShape &dpsids) const;

  /// Add the element's contribution to its residual vector (wrapper) with cached
  /// association matrix
  void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   // Store the expensive-to-construct matrix
   this->store_association_matrix();
   //Call the generic routine with the flag set to 1
   FoepplVonKarmanEquations::fill_in_contribution_to_residuals(residuals);
   // Remove the expensive-to-construct matrix
   this->delete_association_matrix();
  }
  
  /// Add the element's contribution to its residual vector and
  /// element Jacobian matrix (wrapper) with caching of association matrix
  void fill_in_contribution_to_jacobian(Vector<double> &residuals,
					DenseMatrix<double> &jacobian)
  {
   // Store the expensive-to-construct matrix
   this->store_association_matrix();
   //Call the generic routine with the flag set to 1
   FoepplVonKarmanEquations::fill_in_contribution_to_jacobian(residuals,jacobian);
   // Remove the expensive-to-construct matrix
   this->delete_association_matrix();
  }
  
  // ===================== New access functions [zdec] ==========================
  /// Number of fields (unknowns) interpolated by the element
  virtual unsigned nfield() const
  {
   return Nfield;
  }
  
  virtual bool field_is_bell_interpolated(const unsigned& i_field) const
  {
   if(i_field==2)
    {return true;}
   else
    {return false;}
  }
  
  // ------------------------ Nodal data ---------------------------------------
  /// Get the number of nodes that field i is interpolated over
  /// For FvK in-plane use all nodes, out-of-plane uses vertex nodes
  virtual unsigned nnode_for_field(const unsigned& i_field) const
  {
   return Nnode_for_field[i_field];
  }

  /// Get the nodes associated with interpolating field i
  virtual Vector<unsigned> nodes_for_field(const unsigned& i_field) const
  {
   return Node_index_for_field[i_field];
  }
  
  /// Get the number of basis type for field i at node j
  virtual unsigned ntype_for_field_at_node(const unsigned& i_field,
					   const unsigned& j_node) const
  {
   // The in-plane fields always use 1 unknown at each node
   // The out-of-plane fields use 6 unknowns at the vertex nodes
   return Ntype_for_field_at_node[i_field][j_node];
  }

  /// Get the dof index for field i type k at node j
  virtual unsigned nodal_value_index(const unsigned& i_field,
				     const unsigned& j_node,
				     const unsigned& k_type)
  {
   return Nodal_value_index[i_field][j_node][k_type];
  }
  
  /// Get the dof of the field i at node j of type k
  virtual double nodal_value_for_field_at_node_of_type(const unsigned& i_field,
						       const unsigned& j_nodei,
						       const unsigned& k_type) const
  {
   // As in-plane unknowns (i=0,1) are present at all nodes, they occupy the
   // first two positions (0,1 resp.). For the out-of-plane w (i=2), the k_th
   // type unknown will be at k+2.
   # ifdef PARANOID
   // [zdec] TODO check that 
   #endif

   // Get the local node number for the j_nodei-th node used by field i_field
   unsigned local_node_index = Node_index_for_field[i_field][j_nodei];
   // Get the index of the kth type for i_field at j_nodei
   unsigned local_value_index = Nodal_value_index[i_field][j_nodei][k_type];

   return raw_nodal_value(local_node_index,local_value_index);
  }
  
  /// Get the dof of the field i at node j of type k at time t
  virtual double nodal_value_for_field_at_node_of_type(const unsigned& t,
						       const unsigned& i_field,
						       const unsigned& j_nodei,
						       const unsigned& k_type) const
     {
   // As in-plane unknowns (i=0,1) are present at all nodes, they occupy the
   // first two positions (0,1 resp.). For the out-of-plane w (i=2), the k_th
   // type unknown will be at k+2.
   # ifdef PARANOID
   // [zdec] TODO check that 
   #endif

   // Get the local node number for the j_nodei-th node used by field i_field
   unsigned local_node_index = Node_index_for_field[i_field][j_nodei];
   // Get the index of the kth type for i_field at j_nodei
   unsigned local_value_index = Nodal_value_index[i_field][j_nodei][k_type];

   return raw_nodal_value(t,local_node_index,local_value_index);
  }
  
  
  // ----------------------- Internal data -------------------------------------
  // Each field has its own internal data which is resized according to the
  // number of dof/basis types required
  
  /// Get the number of internal data for field i
  virtual unsigned ninternal_types_for_field(const unsigned& i_field) const
  {
   return Ninternal_types_for_field[i_field];
  }

  /// Get the index of the internal data for field i
  virtual unsigned index_of_internal_data_for_field(const unsigned& i_field) const
  {
   return Index_of_internal_data_for_field[i_field];
  }
  
  /// Return the pointer to the internal data for field i
  virtual Data* internal_data_for_field_pt(const unsigned& i_field) const
  {
   return internal_data_pt(Index_of_internal_data_for_field[i_field]);
  }
  
  /// Return the value at the internal data for field i type k
  virtual double internal_value_for_field_of_type(const unsigned& i_field,
						    const unsigned& k_type) const
  {
   #ifdef PARANOID
   if(!this->element_is_curved())
    {
     throw OomphLibError("You should not be looking for internal values in an\
 element that isnt curved!",
			 OOMPH_CURRENT_FUNCTION,
			 OOMPH_EXCEPTION_LOCATION);
    }
   else if(i_field!=2)
    {
     throw OomphLibError("The only unknown with internal values is w [index 2]",
			 OOMPH_CURRENT_FUNCTION,
			 OOMPH_EXCEPTION_LOCATION);
    }
   #endif
   return internal_data_pt(i_field)->value(k_type);
  }
  
  /// Return the value at the internal data for field i type k
  virtual double internal_value_for_field_of_type(const unsigned& t,
						  const unsigned& i_field,
						  const unsigned& k_type) const
  {
   #ifdef PARANOID
   if(!this->element_is_curved())
    {
     throw OomphLibError("You should not be looking for internal values in an\
 element that isnt curved!",
			 OOMPH_CURRENT_FUNCTION,
			 OOMPH_EXCEPTION_LOCATION);
    }
   else if(i_field!=2)
    {
     throw OomphLibError("The only unknown with internal values is w [index 2]",
			 OOMPH_CURRENT_FUNCTION,
			 OOMPH_EXCEPTION_LOCATION);
    }
   #endif
   return internal_data_pt(i_field)->value(t,k_type);
  }
  
  /// Get the jth bubble dof at the lth internal point. Deliberately broken 
  /// for case when there is no curved edge.
  int local_internal_equation(const unsigned& i_field,
			      const unsigned& k_type) const
  {
   // Deliberately break this function for the below cases
   // If there is no curved edge then we cannot return anything meaningful
   // [zdec] Make this PARANOID?
   if(!this->element_is_curved())
    {
     throw OomphLibError(
			 "There are no time-dependent internal 'bubble' dofs for this element.",
			 OOMPH_CURRENT_FUNCTION,  OOMPH_EXCEPTION_LOCATION);
     // Return dummy value -2
     return -2;
    }
   // Now give the lth internal equation number
   else
    { return this->internal_local_eqn(index_of_internal_data_for_field(i_field),
				      k_type); }
  }

  // ========================= End new [zdec] ==================================

  // Is this element curved?
  bool element_is_curved() const
  {return CurvableBellElement<NNODE_1D>::element_is_curved();}
  
  // Get the number of nodal basis types, wrapper
  unsigned nnode_inplane() const
  {return this->nnode();};
  
  // Get the number of internal basis types, wrapper
  unsigned nnode_outofplane() const
  {return this->nvertex_node();};
  
  // Get the number of nodal basis types, wrapper
  unsigned nnodal_basis_type() const
  {return CurvableBellElement<NNODE_1D>::nnodal_basis_type();};
  
  // Get the number of internal basis types, wrapper
  unsigned nbubble_basis_type() const
  {return CurvableBellElement<NNODE_1D>::nbubble_basis_type();};

  // Get the number of internal bases, wrapper
  unsigned nbubble_basis() const
  {return CurvableBellElement<NNODE_1D>::nbubble_basis();};
 
 protected:
  /// Get rotation matrices that change the degrees of freedom to the basis set
  /// by Rotated_basis_fct_pt
  inline void rotation_matrix_at_node(const unsigned& inode, DenseDoubleMatrix&
				      rotation_matrix) const;
 
  /// \short Transform the shape functions so that they correspond to
  /// the new rotated dofs
  inline void rotate_shape(Shape& shape) const;

  /// \short Transform the shape functions and first derivatives so that they 
  /// correspond to the new rotated dofs
  inline void rotate_shape(Shape& shape, DShape& dshape) const;
  
  /// \short Transform the shape functions, first and second   derivatives so 
  /// that they correspond to the new rotated dofs
  inline void rotate_shape(Shape& shape, DShape& dshape,
			   DShape& d2shape) const;
  
  /// \short Get the jth bubble dof at the lth internal point at time t(=0).
  inline double get_w_bubble_dof(const unsigned& l,
				 const unsigned& j,
				 const unsigned& t = 0) const;
  
  /// \short Get the jth bubble dof at the lth internal point
  int local_w_bubble_equation(const unsigned& l, const unsigned& j) const;
  
 public:
  ///\short  Constructor: Call constructors for C1CurvedBellElement and
  /// Biharmonic equations
  FoepplVonKarmanC1CurvedBellElement() :
   FoepplVonKarmanEquations(), 
   Ninternal_types_for_field(Nfield,0),
   Index_of_internal_data_for_field(Nfield),
   Rotated_basis_fct_pt(0),
   Nnodes_to_rotate(0)
  {
   // Use the higher order integration scheme
   TGauss<2,4>* new_integral_pt = new TGauss<2,4>;
   this->set_integration_scheme(new_integral_pt); 

   // Add the (zero) internal data for each field
   for(unsigned i_field=0; i_field<Nfield; i_field++)
    {
     Index_of_internal_data_for_field[i_field] = this->add_internal_data(new Data(0));
    }

   Association_matrix_pt=0;
  }

  ///\short  Destructor: clean up alloacations
  ~FoepplVonKarmanC1CurvedBellElement()
  {
  }

  /// Broken copy constructor
  FoepplVonKarmanC1CurvedBellElement(const
				     FoepplVonKarmanC1CurvedBellElement<NNODE_1D>& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanC1CurvedBellElement");
  }

  /// Broken assignment operator
  void operator=(const FoepplVonKarmanC1CurvedBellElement<NNODE_1D>&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanC1CurvedBellElement");
  }

  /// \short Set up the rotated degrees of freedom
  inline void set_up_rotated_dofs(const unsigned& nnodes_to_rotate,
				  const Vector<unsigned>& nodes_to_rotate,
				  const BasisVectorsFctPt& basis_vectors_fct_pt);

  /// Upgrade the Bell element to a curved Bernadou element
  virtual void upgrade_element_to_curved(const MyC1CurvedElements::Edge& curved_edge, const double& s_ubar,
					 const double& s_obar,  CurvilineGeomObject* parametric_edge, 
					 const unsigned& boundary_order)
  {
   CurvableBellElement<NNODE_1D>::
    upgrade_element_to_curved(curved_edge, s_ubar, s_obar,
			      parametric_edge, boundary_order);
   Ninternal_types_for_field[2]+=
    this->bernadou_element_basis_pt()->n_internal_dofs();
  }
  
  /// \short  Required  # of `values' (pinned or dofs)
  /// at node n
  inline unsigned required_nvalue(const unsigned &n) const
  {
   return Initial_Nvalue[n];
  }
 
  /// \short Output function:
  ///  x,y,u   or    x,y,z,u
  void output(std::ostream &outfile)
  {
   FoepplVonKarmanEquations::output(outfile);
  }

  ///  \short Output function:
  ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
  void output(std::ostream &outfile, const unsigned &n_plot)
  {FoepplVonKarmanEquations::output(outfile,n_plot);}

  ///  \short Output function:
  ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
  void output_interpolated_exact_soln(std::ostream &outfile,
				      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt, const unsigned& n_plot);

  /// \short C-style output function:
  ///  x,y,u   or    x,y,z,u
  void output(FILE* file_pt)
  {FoepplVonKarmanEquations::output(file_pt);}

  ///  \short C-style output function:
  ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
  void output(FILE* file_pt, const unsigned &n_plot)
  {FoepplVonKarmanEquations::output(file_pt,n_plot);}


  /// \short Output function for an exact solution:
  ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
  void output_fct(std::ostream &outfile, const unsigned &n_plot,
		  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
   FoepplVonKarmanEquations::output_fct(outfile,n_plot,
					exact_soln_pt);
  }

  /// \short Output function for a time-dependent exact solution.
  ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
  /// (Calls the steady version)
  void output_fct(std::ostream &outfile, const unsigned &n_plot,
		  const double& time,
		  FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
   FoepplVonKarmanEquations::output_fct(outfile,n_plot,time,exact_soln_pt);
  }


 protected:

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  void shape_and_test_foeppl_von_karman(const Vector<double> &s, Shape &psi,
					Shape& psi_b, Shape& test, Shape& test_b) const;

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  inline double dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s,
							    Shape &psi, Shape &psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
							    Shape &test, Shape &test_b, DShape &dtest_dx, DShape &dtest_b_dx) const;

  /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
  inline double d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,
							      Shape &psi, Shape &psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
							      DShape &d2psi_dx2, DShape &d2psi_b_dx2,
							      Shape &test, Shape &test_b, DShape &dtest_dx, DShape &dtest_b_dx,
							      DShape &d2test_dx2, DShape &d2test_b_dx2) const;

  /// Get dshape and dtest for C0 (in--plane) displacements
  inline double dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s, 
								Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) const;

  // Private Data Members
 private:
  //============================================================================
  // [zdec] New private data storing the information about the fields nodes and
  // types
  
  /// The number of fields (unknowns) over the element [ux, uy, w]
  static const unsigned Nfield;
  
  /// The number of nodes
  static const unsigned Nnode;
  
  /// The number of nodes in each field = {nnode, nnode, nvertex}
  static const Vector<unsigned> Nnode_for_field;
  
  /// The index of the nodes used for each field [i_field][j_node]
  static const Vector< Vector<unsigned> > Node_index_for_field;
  
  /// The number of basis types for each field at each node [i_field][j_node]
  static const Vector< Vector<unsigned> > Ntype_for_field_at_node;

  /// The index of the basis types for each field
  static const Vector< Vector< Vector<unsigned> > > Nodal_value_index;
  
  /// The number of internal dofs/bases for each field
  Vector<unsigned> Ninternal_types_for_field;
  
  //============================================================================
  
  /// Enum to store which edge is curved set to none when element has no curved
  /// edges
  MyC1CurvedElements::Edge Curved_edge;

  /// Indices at which the added internal data is stored
  Vector<unsigned> Index_of_internal_data_for_field;

  /// Pointer to Stored Association matrix
  DenseMatrix<double>* Association_matrix_pt;

  /// Static int that holds the number of variables at
  /// nodes: always the same
  static const unsigned Initial_Nvalue[];

  /// A Pointer to the function that sets up the rotated basis at point x
  BasisVectorsFctPt Rotated_basis_fct_pt;

  /// Which nodes are we rotating
  Vector<unsigned> Nodes_to_rotate;

  /// Number of nodes to rotate
  unsigned Nnodes_to_rotate;


 };




 ////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////


 //==============================================================================
 /// Face geometry for the FoepplVonKarmanC1CurvedBellElement elements: The
 /// spatial dimension of the face elements is one lower than that of the bulk
 /// element but they have the same number of points along their 1D edges.
 //==============================================================================
 template < unsigned NNODE_1D>
 class FaceGeometry<FoepplVonKarmanC1CurvedBellElement<NNODE_1D> >:
  public virtual TElement<1,NNODE_1D>
 {

 public:

  /// \short Constructor: Call the constructor for the
  /// appropriate lower-dimensional TElement
  FaceGeometry() : TElement<1,NNODE_1D>() {}

 };



 /////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////

 //Inline functions:

 // //==============================================================================
 // /// Get the mapped position in the element. For straight sided elements this is
 // /// and affine mapping.
 // //==============================================================================
 // template < unsigned NNODE_1D>
 // void
 // FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::

 //==============================================================================
 /// Get the jth bubble dof at the lth internal point. Deliberately broken for
 /// the case where there is no curved edge
 //==============================================================================
 // [zdec] IN THE BIN
 /*
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::get_w_bubble_dof
 (const unsigned& l, const unsigned& j, const unsigned& t) const
 {return CurvableBellElement<NNODE_1D>::get_bubble_dof(l,j,t);} 
 */
 
 //==============================================================================
 /// Get the jth bubble dof at the lth internal point. Deliberately broken for
 /// case when there is no curved edge.
 //==============================================================================
 // [zdec] IN THE BIN
 /*
 template < unsigned NNODE_1D>
 int FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::local_w_bubble_equation(const
									   unsigned& l, const unsigned& j) const
 {return CurvableBellElement<NNODE_1D>::local_bubble_equation(l,j);}
 */
 
 //==============================================================================
 /// Set up the rotated degrees of freedom: includes a check for the number of
 /// rotation nodes being greater than three.
 //==============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::set_up_rotated_dofs(const unsigned&
									nnodes_to_rotate, const Vector<unsigned>& nodes_to_rotate, const
									BasisVectorsFctPt& basis_vectors_fct_pt)
 {
  // Change the member Nnode_to_rotate
  Nnodes_to_rotate=nnodes_to_rotate;
#ifdef PARANOID
  // Check that the number of nodes is smaller than 3
  if( nnodes_to_rotate > 3 )
   {
    throw OomphLibError(
			"There are only three nodes per element, so we cannot rotate more than\
 three ", OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
   }
#endif

  Nodes_to_rotate = nodes_to_rotate;

  // Point to the basis vectors function
  Rotated_basis_fct_pt = basis_vectors_fct_pt;
 }

 //==============================================================================
 /// Rotate the shape functions according to
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //==============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 rotation_matrix_at_node (const unsigned& inode, DenseDoubleMatrix&
			  rotation_matrix) const
 {
  // Initialise x normal and tangent
  Vector<double> x(2,0.0);

  // Get the node pointer
  Node* nod_pt=this->node_pt(inode);

  // Get the position of the vertex
  x[0]=nod_pt->x(0);
  x[1]=nod_pt->x(1);

  // Initialise the two basis vectors
  Vector<Vector<double> > bi(2,Vector<double>(2,0.0));
  Vector<DenseMatrix<double> > Dbi(2,DenseMatrix<double>(2,2,0.0));

  // Now find the two new basis vectors
  // Get the normal - overload if not continuous
  (*Rotated_basis_fct_pt)(x,bi[0],bi[1], Dbi[0], Dbi[1]);

  // Rotation matrix, B
  DenseMatrix<double> b1(2,2,0.0),b22(3,3,0.0), b21(3,2,0.0);

  // Fill in the submatrices
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
    for(unsigned beta=0; beta<2; ++beta)
     {
      // Fill in b1 - the Jacobian
      // Fill in the rotation of the first derivatives
      b1(alpha,beta) = bi[beta][alpha];
 
      // Avoid double counting the cross derivative
      if(alpha<=beta)
       {
	// Define row index
	const unsigned row = alpha + beta;
	for(unsigned gamma=0; gamma<2; ++gamma)
	 {
	  // Fill in b21 - the non affine part of the Jacobian derivative
	  // Define column index
	  unsigned col_b21 = gamma;
	  // Fill in the non-affine part of the rotation of the second derivatives
	  // if( beta>= alpha) ?
	  b21(row,col_b21) += Dbi[gamma](alpha,beta);
	  for(unsigned delta=0; delta<2; ++delta)
	   {
      
	    // Fill in b22 - the Affine part of the Jacobian derivative
	    // Redefine column index for the next submatrix
	    unsigned col_b22 = gamma + delta;
	    // Fill in the affine part of the rotation of the second derivatives
	    // if( beta>= alpha) ?
	    b22(row,col_b22) += bi[gamma][alpha] * bi[delta][beta];
	   }
	 }
       }
     }
   }


  // Fill in the submatrices to the full (6x6) matrix - we need to right multiply
  // this matrix so we need the transpose of the Jacobian
  // W dof remains the same
  rotation_matrix(0,0)=1.0;
  // Fill in b1
  for(unsigned i=0; i<2; ++i){
   for(unsigned j=0; j<2; ++j)
    { rotation_matrix(1+j,1+i) = b1(i,j); }
  }
  // Fill in b21
  for(unsigned i=0; i<3; ++i){
   for(unsigned j=0; j<2; ++j)
    { rotation_matrix(1+j,3+i) = b21(i,j); }
  }
  // Fill in b22
  for(unsigned i=0; i<3; ++i){
   for(unsigned j=0; j<3; ++j)
    { rotation_matrix(3+j,3+i) = b22(i,j); }
  }
 }


 //======================================================================
 /// Define the shape functions and test functions and derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //======================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::shape_and_test_foeppl_von_karman(
										     const Vector<double> &s, Shape &psi, Shape& psi_b,  Shape &test, Shape& test_b
										     ) const
 {
  throw OomphLibError(
		      "This still needs testing for curved elements.",
		      "void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::\
shape_and_test_foeppl_von_karman(...)", OOMPH_EXCEPTION_LOCATION); // HERE

  this->basis(s,psi,psi_b); 
 
  // Rotate the degrees of freedom
  rotate_shape(psi);
 
  // Galerkin
  // (Shallow) copy the basis functions
  test = psi;
  test_b = psi_b;

 }


 //======================================================================
 /// Define the shape functions and test functions and derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //======================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 dshape_and_dtest_eulerian_foeppl_von_karman(const Vector<double> &s, Shape &psi,
					     Shape& psi_b, DShape &dpsidx, DShape& dpsi_b_dx,  Shape &test, Shape& test_b,
					     DShape &dtestdx,DShape &dtest_b_dx) const
 {
  // Throw if called 
  throw OomphLibError(
		      "This still needs testing for curved elements.",
		      "void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::\
dshape_and_dtest_foeppl_von_karman(...)",
		      OOMPH_EXCEPTION_LOCATION);// HERE

  // Get the basis 
  double J=this->d_basis_eulerian(s,psi,psi_b,dpsidx,dpsi_b_dx);
 
  // Rotate the degrees of freedom
  rotate_shape(psi, dpsidx);
  // Galerkin
  // (Shallow) copy the basis functions
  test = psi;
  dtestdx= dpsidx;
  test_b = psi_b;
  dtest_b_dx= dpsi_b_dx;

  return J;
 }

 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 d2shape_and_d2test_eulerian_foeppl_von_karman(const Vector<double> &s,  Shape &psi,
					       Shape &psi_b, DShape &dpsidx, DShape &dpsi_bdx,  DShape &d2psidx,
					       DShape &d2psi_bdx,
					       Shape &test, Shape &test_b, DShape &dtestdx, DShape &dtest_bdx,
					       DShape &d2testdx, DShape &d2test_bdx) const
 {
  //Call the geometrical shape functions and derivatives
  double J=this->d2_basis_eulerian(s,psi,psi_b,dpsidx,dpsi_bdx,d2psidx,d2psi_bdx);
  // Rotate the dofs
  rotate_shape(psi, dpsidx, d2psidx);
 
  // Galerkin
  //Set the test functions equal to the shape functions (this is a shallow copy)
  test = psi;
  dtestdx= dpsidx;
  d2testdx = d2psidx;
  test_b = psi_b;
  dtest_bdx= dpsi_bdx;
  d2test_bdx = d2psi_bdx;

  //Return the jacobian
  return J;
 }

 //======================================================================
 /// Rotate the shape functions according to the specified basis on the 
 /// boundary. This function does a DenseDoubleMatrix solve to determine 
 /// new basis - which could be speeded up by caching the matrices higher
 /// up and performing the LU decomposition only once
 //======================================================================
 template < unsigned NNODE_1D>
 inline void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 rotate_shape(Shape& psi) const
 {
  const unsigned n_dof_types  = nnodal_basis_type();
  // Loop over the nodes with rotated dofs
  for(unsigned i=0; i<Nnodes_to_rotate; ++i)
   {
    // Get the nodes
    unsigned inode = Nodes_to_rotate[i];

    // Construct the vectors to hold the shape functions
    Vector<double> psi_vector(n_dof_types);

    // Get the rotation matrix
    DenseDoubleMatrix  rotation_matrix(n_dof_types,n_dof_types,0.0);
    this->rotation_matrix_at_node(inode,rotation_matrix);

    // Copy to the vectors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy to the vectors
      for(unsigned k=0;k<n_dof_types;++k)
       {
	// Copy over shape functions
	// psi_vector[l]=psi(inode,l);
	psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
       }
     }

    // Copy back to shape the rotated vetcors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
     }
   }
 }

 //======================================================================
 /// Rotate the shape functions according to the specified basis on the 
 /// boundary. This function does a DenseDoubleMatrix solve to determine 
 /// new basis - which could be speeded up by caching the matrices higher
 /// up and performing the LU decomposition only once
 //======================================================================
 template < unsigned NNODE_1D>
 inline void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 rotate_shape(Shape& psi, DShape& dpsidx) const
 {
  const unsigned n_dof_types  = nnodal_basis_type();
  const unsigned n_dim  = this->dim();
  // Loop over the nodes with rotated dofs
  for(unsigned i=0; i<Nnodes_to_rotate; ++i)
   {
    // Get the nodes
    unsigned inode = Nodes_to_rotate[i];

    // Construct the vectors to hold the shape functions
    Vector<double> psi_vector(n_dof_types);
    Vector<Vector<double> > dpsi_vector_dxi(n_dim,Vector<double>(n_dof_types));

    // Get the rotation matrix
    DenseDoubleMatrix  rotation_matrix(n_dof_types,n_dof_types,0.0);
    this->rotation_matrix_at_node(inode,rotation_matrix);

    // Copy to the vectors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy to the vectors
      for(unsigned k=0;k<n_dof_types;++k)
       {
	// Copy over shape functions
	// psi_vector[l]=psi(inode,l);
	psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
	// Copy over first derivatives
	for(unsigned i=0;i<n_dim; ++i)
	 {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
       }
     }

    // Copy back to shape the rotated vetcors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<n_dim; ++i)
       {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
     }
   }
 }

 //======================================================================
 /// Rotate the shape functions according to the specified basis on the 
 /// boundary. This function does a DenseDoubleMatrix solve to determine 
 /// new basis - which could be speeded up by caching the matrices higher
 /// up and performing the LU decomposition only once
 //======================================================================
 template < unsigned NNODE_1D>
 inline void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 rotate_shape(Shape& psi, DShape& dpsidx, DShape& d2psidx) const
 {
  const unsigned n_dof_types  = nnodal_basis_type();
  const unsigned n_dim  = this->dim();
  // n_dimth triangle number
  const unsigned n_2ndderiv = ((n_dim+1)*(n_dim))/2; // Guaranteed to be positive
  // Loop over the nodes with rotated dofs
  for(unsigned i=0; i<Nnodes_to_rotate; ++i)
   {
    // Get the nodes
    unsigned inode = Nodes_to_rotate[i];

    // Construct the vectors to hold the shape functions
    Vector<double> psi_vector(n_dof_types);
    Vector<Vector<double> > dpsi_vector_dxi(n_dim,Vector<double>(n_dof_types));
    Vector<Vector<double> > d2psi_vector_dxidxj(n_2ndderiv,Vector<double>(n_dof_types));

    // Get the rotation matrix
    DenseDoubleMatrix  rotation_matrix(n_dof_types,n_dof_types,0.0);
    this->rotation_matrix_at_node(inode,rotation_matrix);

    // Copy to the vectors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy to the vectors
      for(unsigned k=0;k<n_dof_types;++k)
       {
	// Copy over shape functions
	// psi_vector[l]=psi(inode,l);
	psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
	// Copy over first derivatives
	for(unsigned i=0;i<n_dim; ++i)
	 {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
	for(unsigned i=0;i<n_2ndderiv; ++i)
	 {d2psi_vector_dxidxj[i][l]+=d2psidx(inode,k,i)*rotation_matrix(l,k);}
       }
     }

    // Copy back to shape the rotated vetcors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<n_dim; ++i)
       {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
      // Copy over second derivatives
      for(unsigned i=0;i<n_2ndderiv; ++i)
       {d2psidx(inode,l,i)=d2psi_vector_dxidxj[i][l];}
     }
   }
 }

 //============================================================================
 /// Define the shape functions and test functions and derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //=============================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 dshape_u_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s,Shape &psi,
						 DShape &dpsidx,  Shape &test,  DShape &dtestdx) const
 {
  // Initialise
  const unsigned n_node = this->nnode();
  DShape dpsids(n_node,this->dim());
  dshape_u_local(s,psi,dpsids);
  double J;

  // Get the Jacobian
  DenseMatrix<double> jacobian(this->dim(),this->dim()),
   inverse_jacobian(this->dim(),this->dim());
  J = CurvableBellElement<NNODE_1D>::local_to_eulerian_mapping(s,jacobian,inverse_jacobian);
  // If the element has been upgraded
  // Now find the global derivatives
  for(unsigned l=0;l<n_node;++l)
   { 
    for(unsigned i=0; i<this->dim(); ++i)
     {
      // Initialise to zero
      dpsidx(l,i)=0.0;
      for(unsigned j=0; j<this->dim(); ++j)
       {
	// Convert to local coordinates
	dpsidx(l,i)+=inverse_jacobian(i,j)*dpsids(l,j);
       }
     }
   }
  test = psi;
  dtestdx = dpsidx;
  return J;
 }

 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::pin_all_deflection_dofs() const
 {
  const unsigned w_index = 2;
    
  // Get nodes
  const unsigned n_node = this->nnode_for_field(w_index);
  const Vector<unsigned> nodes = nodes_for_field(w_index);

  // Curved Bell elements only have deflection dofs at vertices
  for(unsigned j_nodei=0; j_nodei<n_node; j_nodei++)
   {
    // Get the j_nodei-th node used by i_field
    unsigned j_node = nodes[j_nodei];
    Node* nod_pt=this->node_pt(j_node);

    // Get the number of types at the current node
    unsigned n_type = ntype_for_field_at_node(w_index, j_nodei);
    
    // Check if it is on the boundary
    for(unsigned k_type=0; k_type<n_type; k_type++)
     {
      // Pin and set the value
      nod_pt->pin(2+k_type);
      nod_pt->set_value(2+k_type,0.0);
     }
   }

  // Now fix internal dofs
  unsigned n_internal = ninternal_types_for_field(w_index);
  for(unsigned k_type=0; k_type<n_internal; k_type++)
   {
    // Get node
    // Pin and set the value
    this->internal_data_for_field_pt(w_index)->pin(k_type);
   }
 }


 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 fix_out_of_plane_displacement_dof(const unsigned& dof_number, const unsigned&
				   b,const DisplacementFctPt& specified_deflection_fct_pt)
 {
  const unsigned n_vertices = nnode_outofplane(), n_dof_types = nnodal_basis_type();
  // Check that the dof number is a sensible value
  if(dof_number >= n_dof_types)
   {
    throw OomphLibError("Foppl von Karman elements only have 6 Hermite deflection degrees\
of freedom at internal points. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}",
                        "FoepplVonKarmanC1CurvedBellElement:fix_out_of_plane_dof()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<n_vertices; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(b);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0]=nod_pt->x(0);
      x[1]=nod_pt->x(1);
      // Get value
      double value;
      specified_deflection_fct_pt(x,value);
      // Pin and set the value
      nod_pt->pin(this->w_nodal_index_foeppl_von_karman()+dof_number);
      nod_pt->set_value(this->w_nodal_index_foeppl_von_karman()+dof_number,value);
     }
   }
 }

 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvedBellElement<NNODE_1D>::
 fix_in_plane_displacement_dof(const unsigned& dof_number, const unsigned& b, const
			       DisplacementFctPt& specified_displacement_fct_pt)
 {
  // Initialise constants that we use in this function
  const unsigned n_dof_types = 2;
  const unsigned dim = this->dim();
  const unsigned n_node = this->nnode();
  // Check that the dof number is a sensible value
  if(dof_number >= n_dof_types)
   {
    throw OomphLibError("Foppl von Karman elements only have 2 in-plane displacement degrees\
of freedom at internal points. They are {ux, uy}",
                        "FoepplVonKarmanC1CurvedBellElement:fix_out_of_plane_dof()",
                        OOMPH_EXCEPTION_LOCATION);
   }
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<n_node; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(b);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      // Since the element isn't necessarily isoparametric the nodes 'position'
      // is not necessarily correct?
      Vector<double> x(dim),s(dim);
      this->local_coordinate_of_node(n,s);
      interpolated_x(s, x); 
      // Fill in value
      double value;
      specified_displacement_fct_pt(x,value);
      // Pin and set the value
      nod_pt->pin(this->u_nodal_index_foeppl_von_karman()+dof_number);
      nod_pt->set_value(this->u_nodal_index_foeppl_von_karman()+dof_number,value);
     }
   }
 }
}


#endif
