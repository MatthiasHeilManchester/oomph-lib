// Header file for the Biharmonic Bell elements
#ifndef OOMPH_BIHARMONIC_CURVED_BELL_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_CURVED_BELL_ELEMENTS_HEADER

#include "kirchhoff_plate_bending_elements.h"
#include "../C1_basis/Bell_element_basis.h"
#include "../C1_basis/C1_curved_elements.h"
#include "../C1_basis/my_geom_object.h"

namespace oomph
{
//===============================================================================
/// KirchhoffPlateBendingC1CurvedBellElement elements are a subparametric scheme
/// with  linear Lagrange interpolation for approximating the geometry and
/// the C1-functions for approximating variables.
//==============================================================================

class KirchhoffPlateBendingC1CurvedBellElement : 
  public virtual CurvableBellElement<2>, 
  public virtual KirchhoffPlateBendingEquations
{
public:
 ///\short  Constructor: Call constructors for C1CurvedBellElement and
 /// Biharmonic equations
 KirchhoffPlateBendingC1CurvedBellElement() : CurvableBellElement(), 
  KirchhoffPlateBendingEquations(), Rotated_basis_fct_pt(0), Nnodes_to_rotate(0)
  {  }

 // Destructor
 ~KirchhoffPlateBendingC1CurvedBellElement() { }

 /// Broken copy constructor
 KirchhoffPlateBendingC1CurvedBellElement(const
  KirchhoffPlateBendingC1CurvedBellElement& dummy)
  { BrokenCopy::broken_copy(OOMPH_CURRENT_FUNCTION);}

 /// Broken assignment operator
 void operator=(const KirchhoffPlateBendingC1CurvedBellElement&)
  { BrokenCopy::broken_assign(OOMPH_CURRENT_FUNCTION);}

 /// \short Function pointer to basis vectors function which sets  basis vectors
 /// b1 and b2 (which are in general functions of x)
 typedef void (*BasisVectorsFctPt) (const Vector<double>& x, Vector<double>& b1,
  Vector<double>& b2, DenseMatrix<double>& Db1, DenseMatrix<double>& Db2);

 /// Add the element's contribution to its residual vector (wrapper) with cached
 /// association matrix
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   // Store the expensive-to-construct matrix
   this->store_association_matrix();
   //Call the generic routine with the flag set to 1
   KirchhoffPlateBendingEquations::fill_in_contribution_to_residuals(residuals);
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
   KirchhoffPlateBendingEquations::fill_in_contribution_to_jacobian(residuals,jacobian);
   // Remove the expensive-to-construct matrix
   this->delete_association_matrix();
  }

  // Get the number of nodal basis types, wrapper
  unsigned nnodal_basis_type() const {return CurvableBellElement::nnodal_basis_type();};

  // Get the number of internal basis types, wrapper
  unsigned nbubble_basis_type() const {return CurvableBellElement::nbubble_basis_type();};

  // Get the number of internal bases, wrapper
  unsigned nbubble_basis() const {return CurvableBellElement::nbubble_basis();};

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
 
 /// \short Get the jth bubble dof at the lth internal point.
 inline double get_w_bubble_dof(const unsigned& l, const unsigned& j) const
  {return CurvableBellElement::get_bubble_dof(l,j);} 

 /// \short Get the jth bubble dof at the lth internal point
 inline int local_w_bubble_equation(const unsigned& l, const unsigned& j) const
  {return CurvableBellElement::local_bubble_equation(l,j);}

public:
 /// \short Set up the rotated degrees of freedom
 inline void set_up_rotated_dofs(const unsigned& nnodes_to_rotate ,
  const Vector<unsigned>& nodes_to_rotate, const BasisVectorsFctPt&
  basis_vectors_fct_pt);

 /// \short  Required  # of `values' (pinned or dofs)
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const
  {return Initial_Nvalue;}

 /// \short Output function:
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {KirchhoffPlateBendingEquations::output(outfile);}

 ///  \short Output function:
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {KirchhoffPlateBendingEquations::output(outfile,n_plot);}

 ///  \short Output function:
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output_interpolated_exact_soln(std::ostream &outfile,
  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt, const unsigned& n_plot);

 /// \short C-style output function:
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {KirchhoffPlateBendingEquations::output(file_pt);}

 ///  \short C-style output function:
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {KirchhoffPlateBendingEquations::output(file_pt,n_plot);}


 /// \short Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
   KirchhoffPlateBendingEquations::output_fct(outfile,n_plot,
    exact_soln_pt);
  }

 /// \short Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
   KirchhoffPlateBendingEquations::output_fct(outfile,n_plot,time,
    exact_soln_pt);
  }

public:
 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline void shape_and_test_biharmonic(const Vector<double> &s, Shape &psi,
  Shape& psi_b, Shape& test, Shape& test_b) const;

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_and_dtest_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape &psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  Shape &test, Shape &test_b, DShape &dtest_dx, DShape &dtest_b_dx) const;

 /// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double d2shape_and_d2test_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape &psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  DShape &d2psi_dx2, DShape &d2psi_b_dx2,
  Shape &test, Shape &test_b, DShape &dtest_dx, DShape &dtest_b_dx,
  DShape &d2test_dx2, DShape &d2test_b_dx2) const;

// Private Data Members
private:

 /// \short Static int that holds the number of variables at
 /// nodes: always the same
 static const unsigned Initial_Nvalue;

//  /// \short unsigned that holds the index the 'bubble' dofs of the element are
//  // stored
//  unsigned Bubble_w_internal_index;

 /// A Pointer to the function that sets up the rotated basis at point x
 BasisVectorsFctPt Rotated_basis_fct_pt;

 /// Which nodes are we rotating
 Vector<unsigned> Nodes_to_rotate;

 // Number of nodes to rotate
 unsigned Nnodes_to_rotate;
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

//==============================================================================
/// Face geometry for the KirchhoffPlateBendingC1CurvedBellElement elements: The
/// spatial dimension of the face elements is one lower than that of the bulk
/// element but they have the same number of points along their 1D edges.
//==============================================================================
template<>
class FaceGeometry<KirchhoffPlateBendingC1CurvedBellElement > :
 public virtual TElement<1,2>
{

  public:

 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<1,2>() {}

};

/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//==============================================================================
/// Set up the rotated degrees of freedom: includes a check for the number of
/// rotation nodes being greater than three.
//==============================================================================
void KirchhoffPlateBendingC1CurvedBellElement::set_up_rotated_dofs(const unsigned&
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
void KirchhoffPlateBendingC1CurvedBellElement::
 rotation_matrix_at_node (const unsigned& inode, DenseDoubleMatrix&
 rotation_matrix) const
{
 // Initialise x normal and tangent
 Vector<double> x(2,0.0);

 // Get the node pointer
 Node* nod_pt=this->CurvableBellElement::node_pt(inode);

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
  void KirchhoffPlateBendingC1CurvedBellElement::shape_and_test_biharmonic(
   const Vector<double> &s, Shape &psi, Shape& psi_b,  Shape &test, Shape& test_b
   ) const
 {
  // Get dummy shape functions for the Bell call
  DShape dpsidx(3,6,2);
  DShape d2psidx(3,6,3);
   
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
  double KirchhoffPlateBendingC1CurvedBellElement::
  dshape_and_dtest_eulerian_biharmonic(const Vector<double> &s, Shape &psi,
  Shape& psi_b, DShape &dpsidx, DShape& dpsi_b_dx,  Shape &test, Shape& test_b,
  DShape &dtestdx,DShape &dtest_b_dx) const
 {
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
 
  double KirchhoffPlateBendingC1CurvedBellElement::
   d2shape_and_d2test_eulerian_biharmonic(const Vector<double> &s,  Shape &psi,
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
inline void KirchhoffPlateBendingC1CurvedBellElement::
 rotate_shape(Shape& psi) const
{
 // Loop over the nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];

   // Construct the vectors to hold the shape functions
   Vector<double> psi_vector(6,0);
   
   // Get the rotation matrix
   DenseDoubleMatrix  rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Copy to the vectors
   for(unsigned l=0;l<6;++l)
    for(unsigned k=0;k<6;++k)
    { psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); }

   // Copy back to shape the rotated vetcors
   for(unsigned l=0;l<6;++l)
    { psi(inode,l)=psi_vector[l]; }
  }
}

//======================================================================
/// Rotate the shape functions according to the specified basis on the 
/// boundary. This function does a DenseDoubleMatrix solve to determine 
/// new basis - which could be speeded up by caching the matrices higher
/// up and performing the LU decomposition only once
//======================================================================
inline void KirchhoffPlateBendingC1CurvedBellElement::
 rotate_shape(Shape& psi, DShape& dpsidx) const
{
 // Loop over the nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];

   // Construct the vectors to hold the shape functions
   Vector<double> psi_vector(6,0);
   Vector<Vector<double> > dpsi_vector_dxi(2,Vector<double>(6,0));

   // Get the rotation matrix
   DenseDoubleMatrix  rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Copy to the vectors
   for(unsigned l=0;l<6;++l)
    {
     // Copy to the vectors
     for(unsigned k=0;k<6;++k)
      {
      // Copy over shape functions
     // psi_vector[l]=psi(inode,l);
      psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
       {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
     }
    }

   // Copy back to shape the rotated vetcors
   for(unsigned l=0;l<6;++l)
    {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
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
inline void KirchhoffPlateBendingC1CurvedBellElement::
 rotate_shape(Shape& psi, DShape& dpsidx, DShape& d2psidx) const
{
 // Loop over the nodes with rotated dofs
 for(unsigned i=0; i<Nnodes_to_rotate; ++i)
  {
   // Get the nodes
   unsigned inode = Nodes_to_rotate[i];

   // Construct the vectors to hold the shape functions
   Vector<double> psi_vector(6,0);
   Vector<Vector<double> > dpsi_vector_dxi(2,Vector<double>(6,0));
   Vector<Vector<double> > d2psi_vector_dxidxj(3,Vector<double>(6,0));

   // Get the rotation matrix
   DenseDoubleMatrix  rotation_matrix(6,6,0.0);
   this->rotation_matrix_at_node(inode,rotation_matrix);

   // Copy to the vectors
   for(unsigned l=0;l<6;++l)
    {
     // Copy to the vectors
     for(unsigned k=0;k<6;++k)
      {
      // Copy over shape functions
      psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
       {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
      for(unsigned i=0;i<3; ++i)
       {d2psi_vector_dxidxj[i][l]+=d2psidx(inode,k,i)*rotation_matrix(l,k);}
      }
    }

   // Copy back to shape the rotated vetcors
   for(unsigned l=0;l<6;++l)
    {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<2; ++i)
        {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
      // Copy over second derivatives
      for(unsigned i=0;i<3; ++i)
        {d2psidx(inode,l,i)=d2psi_vector_dxidxj[i][l];}
    }
  }
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
}
#endif
