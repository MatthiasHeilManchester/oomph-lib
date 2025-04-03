// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#ifndef OOMPH_C1_HELPER_HEADER
#define OOMPH_C1_HELPER_HEADER

#include "unstructured_two_d_mesh_geometry_base.h"
#include "triangle_mesh.h"

namespace oomph
{

  //======================================================================
  /// Thin wrapper to curved line that defines curved boundaries
  /// in triangle meshes. Obtains higher derivatives from the underlying
  /// GeomObject. Used in C1 problems (e.g. plate bending) where we
  /// we need such information about the boundary parametrisation.
  //======================================================================
  class C1CurviLine 

  {
  public:
   
   /// Constructor: pointer to underlying TriangleMeshCurviline object
   C1CurviLine(TriangleMeshCurviLine* triangle_mesh_curviline_pt) :
    Triangle_mesh_curviline_pt(triangle_mesh_curviline_pt)
    {
    }
   
   /// Broken copy constructor
   C1CurviLine(const C1CurviLine& dummy) = delete;
   
   /// Broken assignment operator
   void operator=(const C1CurviLine&) = delete;
   
   /// (Empty) destructor
   virtual ~C1CurviLine() {}
   
   /// Position r as fct of zeta; forward to underlying geom object
   void position(const Vector<double>& zeta,
                 Vector<double>& r) const
    {
     Triangle_mesh_curviline_pt->geom_object_pt()->position(zeta,r); 
    }
   
   
   /// Derivative of position Vector w.r.t. to zeta:
   virtual void dposition(const Vector<double>& zeta,
                          Vector<double>& drdzeta) const
    {
     DenseMatrix<double> drdzeta_general(2,2);
     Triangle_mesh_curviline_pt->geom_object_pt()->dposition(zeta,drdzeta_general);
     drdzeta[0]=drdzeta_general(0,0);
     drdzeta[1]=drdzeta_general(0,1);
    }
   
   /// 2nd derivative of position Vector w.r.t. to coordinates:
   /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
   /// ddrdzeta(alpha,beta,i).
   /// Evaluated at current time.
   virtual void d2position(const Vector<double>& zeta,
                           Vector<double>& d2rdzeta) const
    {
     RankThreeTensor<double> d2rdzeta_general(1,1,2);
     Triangle_mesh_curviline_pt->geom_object_pt()->d2position(zeta,d2rdzeta_general);
     d2rdzeta[0]=d2rdzeta_general(0,0,0);
     d2rdzeta[1]=d2rdzeta_general(0,0,1);
    }
   
   
   /// The underlying GeomObject encodes r(zeta); here we invert this mapping
   /// with a tolerance to get zeta associated with the given point r_target.
   /// Default implemenation based on bisection and Newton's method.
   /// Overload with your own!
   virtual double get_zeta(const Vector<double>& r_target,
                           const double& tol=1.0e-8) const
    {
     Vector<double> r(2);
     Vector<double> zeta_vect(1);
     Vector<double> drdzeta(2);
     
     // Start Newton method in the middle
     double zeta_min=Triangle_mesh_curviline_pt->zeta_start();
     double zeta_max=Triangle_mesh_curviline_pt->zeta_end();
     double zeta=0.5*(zeta_max+zeta_min);

     // Do it
     unsigned max_iter=100;
     for (unsigned iter=0;iter<max_iter;iter++)
      {
       // Residual:
       zeta_vect[0]=zeta;
       position(zeta_vect,r);
       double res=sqrt( (r[0]-r_target[0])*(r[0]-r_target[0]) +
                        (r[1]-r_target[1])*(r[1]-r_target[1]) );
       if (res<tol)
        {
         if ( (zeta>zeta_min-tol) &&
              (zeta<zeta_max+tol) )
          {
           return zeta;
          }
         else
          {
           std::stringstream error_message;
           error_message << "Converged to zeta outside range: zeta = "
                         << zeta << std::endl
                         << "zeta_min/max = " << zeta_min << " "
                         << zeta_max
                         << std::endl;
           throw OomphLibError(error_message.str(),
                               OOMPH_CURRENT_FUNCTION,
                               OOMPH_EXCEPTION_LOCATION);
          }
        }

       
       // Get "Jacobian"
       dposition(zeta_vect,drdzeta);
       double dresdzeta=
        ( (r[0]-r_target[0])*drdzeta[0] +
          (r[1]-r_target[1])*drdzeta[1] )/res;

       // Newton correctino
       zeta-=res/dresdzeta;
      }

     // If we get here died:
     std::stringstream error_message;
     error_message
      << "newton method failed to converge after max_iter = "
      << max_iter << " iterations" << std::endl;
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     
     // dummy return
     return zeta;
    }
   

  private:

   /// Pointer to underlying triangle object
   TriangleMeshCurviLine* Triangle_mesh_curviline_pt;
   
  };



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

 

 //==============================================================================
 /// Base class for DuplicateNodeConstraintElement s (implemented in FvK and KS)
 //==============================================================================
class DuplicateNodeConstraintElement : public virtual GeneralisedElement
{

public:

 /// Constructor to ensure that the deformation is sufficiently smooth
 /// between different parts of a boundary (which may contain isolated
 /// kinks which make it C0 rather than C1). Pass:
 /// -- Pointers to nodes on the "left" and "right" boundary
 /// -- pointers to the C1Curvilines that describe the smooth
 ///    parts of the boundary on either side
 /// -- the boundary coordinates that identifies the corner point
 ///    relative the right and left boundary parametrisation
 ///    (specified via a one-sized vector). 
 DuplicateNodeConstraintElement(
  Node* const& left_node_pt,
  Node* const& right_node_pt,
  C1CurviLine* const& left_boundary_pt,
  C1CurviLine* const& right_boundary_pt,
  Vector<double> const& left_coord,
  Vector<double> const& right_coord) :
  Left_node_pt(left_node_pt),
  Right_node_pt(right_node_pt),
  Left_boundary_pt(left_boundary_pt),
  Right_boundary_pt(right_boundary_pt),
  Left_node_coord(left_coord),
  Right_node_coord(right_coord)
  {}
 

 /// Broken copy constructor
 DuplicateNodeConstraintElement(const C1CurviLine& dummy) = delete;
 
 /// Broken assignment operator
 void operator=(const DuplicateNodeConstraintElement&) = delete;
  
 /// This function must be called after all global boundary
 /// conditions have been applied to pin constraints that
 /// would enforce continuity that's already guaranteed.
 /// Specific implementation happens in derived classes
 /// for specific plate elements
 virtual void pin_redundant_constraints()=0;
 
protected: 

 /// Function to calculate Jacobian and Hessian of the coordinate mapping
 void get_jac_and_hess_of_coordinate_transform(
  DenseMatrix<double>& jac_of_transform,
  Vector<DenseMatrix<double>>& hess_of_transform)
  {
   //----------------------------------------------------------------------
   // We need the parametrisations either side of the vertex which define
   // the coordinates each node uses for its Hermite dofs.
   Vector<double> left_x(2, 0.0); // [zdec] debug
   Vector<double> right_x(2, 0.0); // [zdec] debug
   Vector<double> left_dxids(2, 0.0);
   Vector<double> left_d2xids2(2, 0.0);
   Vector<double> right_dxids(2, 0.0);
   Vector<double> right_d2xids2(2, 0.0);
   Left_boundary_pt->position(Left_node_coord, left_x); // [zdec] debug 
   Right_boundary_pt->position(Right_node_coord, right_x); // [zdec] debug
   Left_boundary_pt->dposition(Left_node_coord, left_dxids);
   Left_boundary_pt->d2position(Left_node_coord, left_d2xids2);
   Right_boundary_pt->dposition(Right_node_coord, right_dxids);
   Right_boundary_pt->d2position(Right_node_coord, right_d2xids2);

   // Get the speed of each parametrisation
   double left_mag =
    sqrt(left_dxids[0] * left_dxids[0] + left_dxids[1] * left_dxids[1]);
   double right_mag =
    sqrt(right_dxids[0] * right_dxids[0] + right_dxids[1] * right_dxids[1]);

   //----------------------------------------------------------------------
   // Normalise dxids to find the tangent vectors and their
   // derivatives either side of the vertex
   Vector<double> left_ti(2, 0.0);
   Vector<double> left_ni(2, 0.0);
   Vector<double> left_dtids(2, 0.0);
   Vector<double> left_dnids(2, 0.0);
   Vector<double> right_ti(2, 0.0);
   Vector<double> right_ni(2, 0.0);
   Vector<double> right_dtids(2, 0.0);
   Vector<double> right_dnids(2, 0.0);
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     // Fill in the tangents either side of the vertex
     left_ti[alpha] = left_dxids[alpha] / left_mag;
     right_ti[alpha] = right_dxids[alpha] / right_mag;
     // Fill in the derivatives of the (normalised) tangents either side of
     // the vertex
     left_dtids[alpha] =
      left_d2xids2[alpha] / std::pow(left_mag, 2) -
      (left_dxids[0] * left_d2xids2[0] + left_dxids[1] * left_d2xids2[1]) *
      left_dxids[alpha] / std::pow(left_mag, 4);
     right_dtids[alpha] = right_d2xids2[alpha] / std::pow(right_mag, 2) -
      (right_dxids[0] * right_d2xids2[0] +
       right_dxids[1] * right_d2xids2[1]) *
      right_dxids[alpha] / std::pow(right_mag, 4);
     // Use these to fill out the corresponding vectors for the normal
     // direction (nx,ny) = (ty,-tx)
    }
   // Use orthogonality to fill in normals and their derivatives
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     left_ni[alpha] = pow(-1, alpha) * left_ti[(alpha + 1) % 2];
     right_ni[alpha] = pow(-1, alpha) * right_ti[(alpha + 1) % 2];
     left_dnids[alpha] = pow(-1, alpha) * left_dtids[(alpha + 1) % 2];
     right_dnids[alpha] = pow(-1, alpha) * right_dtids[(alpha + 1) % 2];
    }

   //----------------------------------------------------------------------
   // We need to fill out the Jacobians and Hessians of the boundary
   // coordinates either side of the vertex
   DenseMatrix<double> left_jac(2, 2, 0.0);
   DenseMatrix<double> right_jac(2, 2, 0.0);
   Vector<DenseMatrix<double>> left_hess(2, DenseMatrix<double>(2, 2, 0.0));
   Vector<DenseMatrix<double>> right_hess(2, DenseMatrix<double>(2, 2, 0.0));
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     // Fill in Jacobians {{nx,tx},{ny,ty}}
     left_jac(alpha, 0) = left_ni[alpha];
     left_jac(alpha, 1) = left_ti[alpha];
     right_jac(alpha, 0) = right_ni[alpha];
     right_jac(alpha, 1) = right_ti[alpha];
     // Fill in Hessians
     // left_hess[alpha](0,0) = 0.0;
     left_hess[alpha](0, 1) = left_dnids[alpha];
     left_hess[alpha](1, 0) = left_dnids[alpha];
     left_hess[alpha](1, 1) = left_dtids[alpha];
     // right_hess[alpha](0,0) = 0.0;
     right_hess[alpha](0, 1) = right_dnids[alpha];
     right_hess[alpha](1, 0) = right_dnids[alpha];
     right_hess[alpha](1, 1) = right_dtids[alpha];
    }

   //----------------------------------------------------------------------
   // We need the inverse Jacobian and Hessian for the left parametrisation
   DenseMatrix<double> left_jac_inv(2, 2, 0.0);
   Vector<DenseMatrix<double>> left_hess_inv(2,
                                             DenseMatrix<double>(2, 2, 0.0));
   left_jac_inv(0, 0) = left_jac(1, 1);
   left_jac_inv(0, 1) = -left_jac(0, 1);
   left_jac_inv(1, 0) = -left_jac(1, 0);
   left_jac_inv(1, 1) = left_jac(0, 0);
   // Fill out inverse of Hessian
   // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     for (unsigned beta = 0; beta < 2; beta++)
      {
       for (unsigned gamma = 0; gamma < 2; gamma++)
        {
         for (unsigned alpha2 = 0; alpha2 < 2; alpha2++)
          {
           for (unsigned beta2 = 0; beta2 < 2; beta2++)
            {
             for (unsigned gamma2 = 0; gamma2 < 2; gamma2++)
              {
               left_hess_inv[alpha](beta, gamma) -=
                left_jac_inv(alpha, alpha2) *
                left_hess[alpha2](beta2, gamma2) *
                left_jac_inv(beta2, beta) * left_jac_inv(gamma2, gamma);
              }
            }
          }
        }
      }
    }

   //----------------------------------------------------------------------
   //----------------------------------------------------------------------
   // Use these to calculate the Jacobian of the left->right transform
   //     J = J_{left}^{-1}J_{right}
   // and the Hessian of the left->right transform
   //     H = H_{left}^{-1}J_{right}J_{right} + J_{left}^{-1}H_{right}
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     for (unsigned beta = 0; beta < 2; beta++)
      {
       for (unsigned gamma = 0; gamma < 2; gamma++)
        {
         // Add contribution to J
         jac_of_transform(alpha, beta) +=
          left_jac_inv(alpha, gamma) * right_jac(gamma, beta);
         for (unsigned mu = 0; mu < 2; mu++)
          {
           // Add second term contribution to H
           hess_of_transform[alpha](beta, gamma) +=
            left_jac_inv(alpha, mu) * right_hess[mu](beta, gamma);
           for (unsigned nu = 0; nu < 2; nu++)
            {
             // Add first term contribution to H
             hess_of_transform[alpha](beta, gamma) +=
              left_hess_inv[alpha](mu, nu) * right_jac(mu, beta) *
              right_jac(nu, gamma);
            }
          }
        }
      }
    }
  } // End get_jac_and_hess_of_coordinate_transform


 /// Store the index of the internal data keeping the Lagrange multipliers
 unsigned Index_of_lagrange_data;

 /// Store the index of the external data for the left node
 unsigned Index_of_left_data;

 /// Store the index of the external data for the right node
 unsigned Index_of_right_data;

 /// Pointer to the left node (before the vertex when anticlockwise)
 Node* Left_node_pt;

 /// Pointer to the right node (after the vertex when anticlockwise)
 Node* Right_node_pt;

 /// Pointer to the left node's boundary parametrisation
 C1CurviLine* Left_boundary_pt;

 /// Pointer to the right node's boundary parametrisation
 C1CurviLine* Right_boundary_pt;

 /// Coordinate of the left node on the left boundary
 Vector<double> Left_node_coord;

 /// Coordinate of the left node on the left boundary
 Vector<double> Right_node_coord;

 /// Tolerance for validating fully pinned constraints
 // [zdec] does this wnat to be the problem residual tolerance?
 double Constraint_tolerance = 1.0e-10;

 /// Tolerance for checking whether a dof has become decoupled from an
 /// equation.
 /// i.e. in the equation y=Ax, how small does A have to be before y no
 /// longer /numerically/ depends on x? This becomes relevant when derivative
 /// directions become orthogonal, we need to ensure they aren't considered
 /// linearly dependent. (We choose this to be slightly larger than machine
 /// precision and it shouldn't generally need to be touched)
 double Orthogonality_tolerance = 1.0e-15;
    
};


 
 
//==============================================================================
/// Class to specify boundary conditions for C1 (plate bending) problems
//==============================================================================
class BoundaryConditionForC1PlateBending
{

public:

 /// Constructor: Specify FD step for automatic evaluation of derivatives
 BoundaryConditionForC1PlateBending(const double& fd_step=1.0e-8) :
  FD_step(fd_step)
  {}

  /// Broken copy constructor
 BoundaryConditionForC1PlateBending(const C1CurviLine& dummy) = delete;
 
 /// Broken assignment operator
 void operator=(const BoundaryConditionForC1PlateBending&) = delete;
 
  /// Pure virtual function to specify value of the function (typically a
 /// displacement component) as a function of zeta, the 1D coordinate
 /// that parametrises the boundary
 virtual double f(const double& zeta) =0 ;


 /// Broken virtual function to specify the normal derivative of the function
 /// (typically a displacement component) w.r.t zeta, the 1D coordinate that
 /// parametrises the boundary. This is only needed for genuine C1 quantities
 /// (or for large-amplitude problems, e.g. Koiter Steigman, where all three
 /// displacement components need to be C1. Broken function will shout if it's
 /// evaluated.
 virtual double dfdn(const double& zeta)
  {
   std::stringstream error_message;
   error_message
    << "You'll have to implement the function that\n"
    << "specifies the normal derivative if you want to use it!\n"
    << "This broken virtual function was probably called \n"
    << "when you tried to implement clamping bcs in a plate problem.\n"
    << std::endl;
   throw OomphLibError(error_message.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);

  }
 
 /// Virtual function to specify the derivative of the function (typically a
 /// displacement component) w.r.t zeta, the 1D coordinate that parametrises the
 /// boundary. Defaults to fd evaluation. Overload with your own version if
 /// you prefer
 double dfdzeta(const double& zeta)
  {
   return (f(zeta+FD_step)-f(zeta))/FD_step;
  }
 

 /// Virtual function to specify the derivative of the normal derivatives 
 /// w.r.t zeta , the 1D coordinate that parametrises the boundary. Defaults to
 /// fd evaluation. Overload with your own version if you prefer.
 double d2fdndzeta(const double& zeta)
  {
   return (dfdn(zeta+FD_step)-dfdn(zeta))/FD_step;
  }

 /// Virtual function to specify the second derivative of the function
 /// w.r.t zeta , the 1D coordinate that parametrises the boundary. Defaults to
 /// fd evaluation. Overload with your own version if you prefer.
 double d2fdzeta2(const double& zeta)
  {
   return (f(zeta+FD_step)-2.0*f(zeta)+f(zeta-FD_step))/(FD_step*FD_step);
  }


private:

 /// FD step
 double FD_step;

};


 ////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////
 


 
 //===start of rotation helper class=========================================
 /// Helper class to contain all the rotation information in the element.
 //==========================================================================
 class RotatedBoundaryHelper
  {
  public:
   
   /// Constructor: just initialise the member data to their defaults (zeros)
   /// hierher Aidan: what's the parent element?
   RotatedBoundaryHelper(FiniteElement* const& parent_element_pt)
    : Parent_element_pt(parent_element_pt),
      Nnode(Parent_element_pt->nvertex_node()),
      Boundary_coordinate_of_node(3, 0.0),
      Nodal_boundary_parametrisation_pt(3, 0),
      Rotation_matrix_at_node(3, DenseMatrix<double>(6, 6, 0.0))
    {
    }
   
   /// Destructor
   ~RotatedBoundaryHelper() {}

   /// Return C1Curviline associated with the boundary
   /// parametrisatino that given local node is on
   /// hierher Aidan: please confirm
   C1CurviLine* nodal_boundary_parametrisation_pt(
    const unsigned& j_node)
    {
     return Nodal_boundary_parametrisation_pt[j_node];
    }
   
   /// Add a new boundary parametrisation to nodes all the nodes in the
   /// vector node_on_boundary
   void set_nodal_boundary_parametrisation(
    const Vector<unsigned>& node_on_boundary,
    const Vector<double>& boundary_coord_of_node,
    C1CurviLine* const& boundary_parametrisation_pt)
    {
     // Loop over all the nodes in node_on_boundary and add the boundary
     // pointer to their vector of boundaries
     unsigned n_node = node_on_boundary.size();
     for (unsigned j = 0; j < n_node; j++)
      {
       // The j-th node on the boundary
       unsigned j_node = node_on_boundary[j];
       
       // Set the boundary parametrisation data pointer for this node
       Nodal_boundary_parametrisation_pt[j_node] = boundary_parametrisation_pt;
       
       // Set the coordinate of node j on this boundary
       Boundary_coordinate_of_node[j_node] = boundary_coord_of_node[j];
       
       update_rotation_matrices();
      } // end of loop over nodes in node_on_boundary [j]
    } // end of set_nodal_boundary_parametrisation()
   

   
   /// Update all rotation matrices (checks if they are needed unless flag is
   /// true)
    void update_rotation_matrices()
    {
      // [zdec] hard coded the three vertex nodes
      unsigned n_vertex = 3;
      // Loop over each vertex
      for (unsigned j_node = 0; j_node < n_vertex; j_node++)
      {
        // If this node does not have a parametrisation (the pointer is still
        // null) skip over it, otherwise we go on to fill out the rotation
        // matrix
        if (!nodal_boundary_parametrisation_pt(j_node))
        {
          continue;
        }

        // Initialise the two basis vectors and their jacobians
        Vector<Vector<double>> bi(2, Vector<double>(2, 0.0));
        Vector<DenseMatrix<double>> dbidx(2, DenseMatrix<double>(2, 2, 0.0));

        // Our new coordinate system:
        //     (l, s)=(normal component, tangent component)
        // which we define in terms of basis vectors (rescaled)
        //     ni=dxi/dl / |n|           <-- Jacobian col 1
        //     ti=dxi/ds / |t|           <-- Jacobian col 2
        // and their derivatives
        //     dnidxj=d/dxj(dxi/dl / |n|) <-- Hessian `col' 1
        //     dtidxj=d/dxj(dxi/ds / |t|) <-- Hessian `col' 2

        // [zdec] we use i and j for brevity
        // but it should be alpha & beta
        // Need to write up how the transformation is done

        // Storage for our basis and derivatives
        Vector<double> ni(2, 0.0);
        Vector<double> ti(2, 0.0);
        Vector<double> dnids(2, 0.0);
        Vector<double> dtids(2, 0.0);

        // All tensors assumed evaluated on the boundary
        // Jacobian of inverse mapping
        DenseMatrix<double> jac_inv(2, 2, 0.0);
        // Hessian of mapping [zdec] (not needed because...)
        Vector<DenseMatrix<double>> hess(2, DenseMatrix<double>(2, 2, 0.0));
        // Hessian of inverse mapping [zdec] (...this can be found by
        // hand)
        Vector<DenseMatrix<double>> hess_inv(2, DenseMatrix<double>(2, 2, 0.0));

        // The basis is defined in terms of the boundary parametrisation
        Vector<double> boundary_coord = {Boundary_coordinate_of_node[j_node]};
        C1CurviLine* boundary_pt =
          Nodal_boundary_parametrisation_pt[j_node];
        Vector<double> x(2, 0.0);
        Vector<double> dxids(2, 0.0);
        Vector<double> d2xids2(2, 0.0);

        // Get position (debug) // hierher Aidan why debug?
        boundary_pt->position(boundary_coord, x);
        // Get tangent vector
        boundary_pt->dposition(boundary_coord, dxids);
        // Get second derivative
        boundary_pt->d2position(boundary_coord, d2xids2);

        double mag_t = sqrt(dxids[0] * dxids[0] + dxids[1] * dxids[1]);
        // ti is the normalised tangent vector
        ti[0] = dxids[0] / mag_t;
        ti[1] = dxids[1] / mag_t;
        // Derivative of (normalised) tangent
        dtids[0] = d2xids2[0] / std::pow(mag_t, 2) -
                   (dxids[0] * d2xids2[0] + dxids[1] * d2xids2[1]) * dxids[0] /
                     std::pow(mag_t, 4);
        dtids[1] = d2xids2[1] / std::pow(mag_t, 2) -
                   (dxids[0] * d2xids2[0] + dxids[1] * d2xids2[1]) * dxids[1] /
                     std::pow(mag_t, 4);
        // n = (t x e_z) implies
        ni[0] = ti[1];
        ni[1] = -ti[0];
        // Same for dnids
        dnids[0] = dtids[1];
        dnids[1] = -dtids[0];

        // Need inverse of mapping to calculate ds/dxi ----------------
        //   /  dx/dl  dx/ds  \ -1  ___  __1__ /  dy/ds -dx/ds \ .
        //   \  dy/dl  dy/ds  /     ---   det  \ -dy/dl  dx/dl /
        //
        //                          ___  /  dl/dx  dl/dy  \ .
        //                          ---  \  ds/dx  ds/dy  /
        //
        // Fill out inverse of Jacobian
        double det = (ni[0] * ti[1] - ni[1] * ti[0]);
        jac_inv(0, 0) = ti[1] / det;
        jac_inv(0, 1) = -ti[0] / det;
        jac_inv(1, 0) = -ni[1] / det;
        jac_inv(1, 1) = ni[0] / det;

        // Fill out the Hessian
        // (unneeded -- can calculate the inverse components by hand)
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // hess[alpha](0,0) = 0.0;
          hess[alpha](0, 1) = dnids[alpha];
          hess[alpha](1, 0) = dnids[alpha];
          hess[alpha](1, 1) = dtids[alpha];
        }

        // Fill out inverse of Hessian
        // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          for (unsigned beta = 0; beta < 2; beta++)
          {
            for (unsigned gamma = 0; gamma < 2; gamma++)
            {
              for (unsigned alpha2 = 0; alpha2 < 2; alpha2++)
              {
                for (unsigned beta2 = 0; beta2 < 2; beta2++)
                {
                  for (unsigned gamma2 = 0; gamma2 < 2; gamma2++)
                  {
                    hess_inv[alpha](beta, gamma) -=
                      jac_inv(alpha, alpha2) * hess[alpha2](beta2, gamma2) *
                      jac_inv(beta2, beta) * jac_inv(gamma2, gamma);
                  }
                }
              }
            }
          }
        }

        // Fill in the rotation matrix using the new basis
        fill_in_rotation_matrix_at_node_with_basis(j_node, jac_inv, hess_inv);

      } // end loop over vertices
    } // end of update_rotation_matrices()


    /// Access function to fill out rot_mat using rotation matrix
    void get_rotation_matrix_at_node(const unsigned& j_node,
                                     DenseMatrix<double>& rot_mat)
    {
      rot_mat = Rotation_matrix_at_node[j_node];
    }

  private:

    /// Helper function to fill in the rotation matrix for a given basis
    void fill_in_rotation_matrix_at_node_with_basis(
      const unsigned& j_node,
      const DenseMatrix<double>& jac_inv,
      const Vector<DenseMatrix<double>>& hess_inv)
    {
      // Rotation matrix, b constructed using submatrices b1, b12, b22
      DenseMatrix<double> b1(2, 2, 0.0), b22(3, 3, 0.0), b12(2, 3, 0.0);

      // Fill in the submatrices
      // Loop over the rotated first derivatives
      for (unsigned mu = 0; mu < 2; mu++)
      {
        // Loop over the unrotated first derivatives
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Fill in b1 - the Jacobian
          // Fill in the affine rotation of the first derivatives
          b1(mu, alpha) = jac_inv(mu, alpha);

          // Loop over unrotated second derivatives
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Avoid double counting the cross derivative
            if (alpha <= beta)
            {
              // Define column index
              const unsigned col = alpha + beta;

              // Fill in the non-affine part of the rotation of the first
              // derivatives
              b12(mu, col) += hess_inv[mu](alpha, beta);
              // [zdec] debug mixed derivative -- add extra
              if (alpha < beta)
              {
                // b12(mu, col) -= hess_inv[mu](alpha, beta);
              }
              // Loop over the rotated second derivatives
              for (unsigned nu = 0; nu < 2; nu++)
              {
                // // Avoid double counting the cross derivative
                // if (mu <= nu)
                {
                  // Fill in b22 - the Affine part of the Jacobian derivative
                  // Redefine row index for the next submatrix
                  unsigned row_b22 = mu + nu;
                  // Fill in the affine part of the rotation of the second
                  // derivatives [zdec] if( beta>= alpha) ?
                  b22(row_b22, col) += jac_inv(mu, alpha) * jac_inv(nu, beta);
                }
              }
            }
          }
        }
      }

      // Fill in the submatrices to the full (6x6) matrix
      Rotation_matrix_at_node[j_node](0, 0) = 1.0;
      // Fill in b1 --- the affine contribution to rotation of the
      // first derivatives
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned j = 0; j < 2; ++j)
        {
          Rotation_matrix_at_node[j_node](1 + i, 1 + j) = b1(i, j);
        }
      }
      // Fill in b21 --- the non-affine (second derivative dependent)
      // rotation of the first derivatives
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned j = 0; j < 3; ++j)
        {
          Rotation_matrix_at_node[j_node](1 + i, 3 + j) = b12(i, j);
        }
      }
      // Fill in b22 --- the rotation of the second derivatives
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned j = 0; j < 3; ++j)
        {
          Rotation_matrix_at_node[j_node](3 + i, 3 + j) = b22(i, j);
        }
      }
    } // end fill_in_rotation_matrix_at_node_with_basis

    /// Pointer to the `parent' finite element which this is a helper force
    FiniteElement* Parent_element_pt;

    /// The number of nodes (that we store rotation data for) in the fvk element
    /// that uses this helper
    unsigned Nnode;

    /// Vector containing boundary parametrised location for each node
    Vector<double> Boundary_coordinate_of_node;

    /// Vector containing boundary parametrisation at each node
    Vector<C1CurviLine*> Nodal_boundary_parametrisation_pt;

    /// Vector containing <rotation matrix at each node>
    Vector<DenseMatrix<double>> Rotation_matrix_at_node;
  };
 


 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==============================================================================
/// Namespace to deal update triangle meshes to deal with C1 elements
/// (beginning)
//==============================================================================
namespace C1Helper
{
 

 //=============================================================
 /// Enum to enumerate the possible edges that could be curved
 //=============================================================
 enum CurvedEdgeEnumeration
 {
  none = -1,
  zero = 0,
  one = 1,
  two = 2
 };

}


 
 //===========================================================================
 /// Template-free base class for curvable Bell Element
 //===========================================================================
class TemplateFreeCurvableBellElement : public virtual FiniteElement
 {

 public:
  
  /// Upgrade the element to be curved
  virtual void upgrade_element_to_curved(const C1Helper::CurvedEdgeEnumeration& curved_edge,
                                         const double& s_ubar,
                                         const double& s_obar,
                                         C1CurviLine* parametric_edge,
                                         const unsigned& boundary_order)=0;
  
  /// Access function to rotated boundary helper object
  virtual RotatedBoundaryHelper* rotated_boundary_helper_pt()=0;


  /// Clamp: i.e. pin the in-plane displacements and pin the out-of-plane
  /// displacement and its normal derivative. We also apply implied
  /// boundary conditions (e.g. specification of dw/dn also implies
  /// d^2w/dn/dzeta etc. boundary_values_pt[i] describes boundary conditions
  /// for the three displacement components... hierher Aidan: this is where
  /// it gets complicated: u_x or u_n etc?
  /// curviline_pt provides a pointer to the representation of the curvilinear
  /// boundary in the triangle mesh.
  virtual void fully_clamp_specified_boundary(
   const unsigned& b,
   const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
   TriangleMeshCurviLine* curviline_pt) = 0;
  
  
  /// Pin i.e. pin the in-plane and out of plane displacements only.
  /// We also apply implied boundary conditions (e.g. specification of w
  /// also implies dw/dt etc. boundary_values_pt[i] describes boundary conditions
  /// for the three displacement components... hierher Aidan: this is where
  /// it gets complicated: u_x or u_n etc?
  /// curviline_pt provides a pointer to the representation of the curvilinear
  /// boundary in the triangle mesh.
  virtual void pin_specified_boundary(
   const unsigned& b,
   const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
   TriangleMeshCurviLine* curviline_pt) = 0;


  /// Factory to create DuplicateNodeConstraintElement.
  /// which ensures that the deformation is sufficiently smooth
  /// between different parts of a boundary (which may contain isolated
  /// kinks which make it C0 rather than C1). Pass:
  /// -- Pointers to nodes on the "left" and "right" boundary
  /// -- pointers to the C1Curvilines that describe the smooth
  ///    parts of the boundary on either side
  /// -- the boundary coordinates that identifies the corner point
  ///    relative the right and left boundary parametrisation
  ///    (specified via a one-sized vector). 
  virtual DuplicateNodeConstraintElement* duplicate_constraint_element_factory(
   Node* const& left_node_pt,
   Node* const& right_node_pt,
   C1CurviLine* const& left_boundary_pt,
   C1CurviLine* const& right_boundary_pt,
   Vector<double> const& left_coord,
   Vector<double> const& right_coord)=0;
  
  
 };
 

//==============================================================================
/// Namespace to deal update triangle meshes to deal with C1 elements
/// (continued)
//==============================================================================
namespace C1Helper
{
 

 /// Boolean to suppress warning about polygonal boundaries
 extern bool Do_not_warn_about_polygonal_boundaries;

 /// Output stream to document split elements
 extern std::ofstream Split_elements_output_stream;
     
 /// Output stream to document duplicated nodes
 extern std::ofstream Duplicated_node_output_stream;
     
 /// Output stream to document elements whose boundaries have been
 /// upgraded to become curved
 extern std::ofstream Upgraded_to_curved_edge_element_stream;

 
//==============================================================================
/// Fct to duplicate nodes that span two boundaries of the triangle mesh pointed to by
/// bulk_mesh_pt. This makes sure that each node on the boundary is only 
/// associated with a single (smooth) boundary. Smooth boundaries  of the mesh/domain
/// are available  via the map c1_curviline_pt. The required continuity of
/// the solution is then imposed by a suitable DuplicateNodeConstraintElement
/// that is created and added to the Mesh pointed to by constraint_mesh_pt.
//==============================================================================
 extern void duplicate_corner_nodes(Mesh* bulk_mesh_pt, 
                                    std::map<unsigned,C1CurviLine*> c1_curviline_pt,
                                    Mesh* constraint_mesh_pt);
 

//==============================================================================
/// A function that upgrades straight sided elements to be curved. This involves
/// Setting up the parametric boundary, F(s) and the first derivative F'(s)
/// We also need to set the edge number of the upgraded element and the positions
/// of the nodes j and k (defined below) and set which edge (k) is to be exterior
///            @ k               
///           /(                 
///          /. \                
///         /._._)               
///      i @     @ j             
/// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
/// representation. That is we need to have a continuous 2nd derivative defined
/// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
/// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
/// as well.
/// The smooth mesh/domain boundaries are available via the map c1_curviline_pt.
/// hierher should we create/store/delete the c1_curviline_pt things in the
/// triangle mesh (either by default or when required). would probaly fix that
/// leak problem.
//=============================================================================
 extern void upgrade_edge_elements_to_curved_boundaries(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt);


 
//======================================================================
/// Function to set up rotated nodes on the boundaries of the triangle
/// mesh pointed to by bulk_mesh_pt. Necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
/// The smooth mesh/domain boundaries are available via the map c1_curviline_pt.
//======================================================================
 extern void rotate_edge_degrees_of_freedom(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt);

 
 
//======================================================================
/// Helper function to upgrade a triangle mesh to C1 computations.
/// -- split elements that have multiple edges on a boundary (so the number
///    of elements in the mesh changes!
/// -- duplicate corner nodes (so the number of nodes in the mesh changes
/// -- create DuplicateNodeConstraintElements that constrain the degrees of freedom
///    at the newly created nodes to ensure continuity of the displacements
///    These elements are added to the mesh pointed to by constraint_mesh_pt
/// -- rotate the degrees of freedom on the curvilinear boundaries so they
///    represent displacements and derivatives in the tangential and
///    normal direction // hierher Aidan true?
/// The smooth mesh/domain boundaries are available via the map c1_curviline_pt.
//======================================================================
 template <class ELEMENT>
 void upgrade_triangle_mesh_for_c1_plate_bending(
  TriangleMeshBase* bulk_mesh_pt,
  Mesh* constraint_mesh_pt)
 {
  // Get map to curvline boundaries of mesh
  std::map<unsigned, TriangleMeshCurviLine*> curviline_boundary_pt =
   bulk_mesh_pt->curviline_boundary_pt();

  unsigned nb=bulk_mesh_pt->nboundary();
  if (!C1Helper::Do_not_warn_about_polygonal_boundaries)
   {
    if (curviline_boundary_pt.size()!=nb)
     {
      std::stringstream warning_message;
      warning_message
       << "Black box helper function upgrade_triangle_mesh_for_c1_plate_bending()\n"
       << "will only rotate degrees of freedom on boundaries that are described\n"
       << "by TriangleMeshCurviLine. It seems that in your triangle mesh only\n"
       << curviline_boundary_pt.size() << " of " << nb << " boundaries \n"
       << "are of this type. This only matters if you want to apply clamping-type\n"
       << "boundary conditions along those boundaries. Continue at your own risk\n"
       << "and/or make this message disappear by setting\n\n"
       << "     C1Helper::Do_not_warn_about_polygonal_boundaries\n\n"
       << "to false. More intelligently, replace the polygonal boundaries\n"
       << "by TriangleMeshCurviLines. Alternatively, you can deal with this\n"
       << "yourself, of course. The black box function is only provided for\n"
       << "convenience."
       << std::endl;
      OomphLibWarning(warning_message.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
      
     }
   }

  
  // Map as "sparse vector" for C1 curvilines
  // Note: don't delete these. They're retained in the
  // Constraint elements. Leak... 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt;
  for (auto [ b, curviline_pt] :  curviline_boundary_pt)
   {
    c1_curviline_pt[b] = new C1CurviLine(curviline_pt);
   }
  
  // Split elements that have multiple edges on a boundary
  // Note: Sets up the boundary loopup scheme too.
  bulk_mesh_pt->
   template split_elements_with_multiple_boundary_edges<ELEMENT>
   (Split_elements_output_stream);
  
  // Duplicate corner nodes
  duplicate_corner_nodes(bulk_mesh_pt,
                         c1_curviline_pt,
                         constraint_mesh_pt);
  
  // Re-setup boundary cooordinates (bypass check of discrepancy
  // between polygon/smooth representation of boundary.
  double backup=ToleranceForVertexMismatchInPolygons::Tolerable_error;
  ToleranceForVertexMismatchInPolygons::Tolerable_error=DBL_MAX;
  for (unsigned b=0;b<nb;b++)
   {
    bulk_mesh_pt->template setup_boundary_coordinates<ELEMENT>(b);
   }
  ToleranceForVertexMismatchInPolygons::Tolerable_error=backup;

  
  // Upgrade
  upgrade_edge_elements_to_curved_boundaries(
   bulk_mesh_pt,
   c1_curviline_pt);
  
  // Rotate degrees of freedom (only really needed for clamped bcs)
  rotate_edge_degrees_of_freedom(
   bulk_mesh_pt,
   c1_curviline_pt);
  
 }

 

 } // end namespace C1Helper


} // namespace oomph

#endif
