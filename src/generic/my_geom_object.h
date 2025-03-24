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
#ifndef OOMPH_MY_GEOM_OBJECTS_HEADER // hierher C1_curviline.h
#define OOMPH_MY_GEOM_OBJECTS_HEADER

#include "unstructured_two_d_mesh_geometry_base.h"


namespace oomph
{


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

 // hierher break copy constructor etc.
 
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
   // hierher throw
   abort();
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
 
  //======================================================================
  /// hierher
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


   // hierher in all of these, why is zeta still a vector?
   
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
           oomph_info << "Return with zeta = " << zeta << std::endl;
           return zeta;
          }
         else
          {
           oomph_info << "Converged to zeta outside range: zeta = "
                      << zeta << std::endl
                      << "zeta_min/max = " << zeta_min << " " << zeta_max
                      << std::endl;
           abort(); // hierher
          }
        }
       
       oomph_info << "iter, zeta, res "
                  << iter << " "
                  << zeta << " "
                  << res << " "
                  << std::endl;
       
       // Get "Jacobian"
       dposition(zeta_vect,drdzeta);
       double dresdzeta=
        ( (r[0]-r_target[0])*drdzeta[0] +
          (r[1]-r_target[1])*drdzeta[1] )/res;

       // Newton correctino
       zeta-=res/dresdzeta;
      }

     // If we get here died:
     oomph_info << "newton method failed to converge after max_iter = "
                << max_iter << " iterations" << std::endl;
     abort(); // hierher

     // dummy return
     return zeta;
    }


   // hierher may retain but all functions should have been wrapped now.
   // /// Pointer to underlying triangle object
   // TriangleMeshCurviLine* triangle_mesh_curviline_pt() const
   //  {
   //   return Triangle_mesh_curviline_pt;
   //  }
   

  private:

   /// Pointer to underlying triangle object
   TriangleMeshCurviLine* Triangle_mesh_curviline_pt;
   
  };





////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



  //===start of rotation helper class=========================================
  /// Helper class to contain all the rotation information in the element.
  class RotatedBoundaryHelper
  {
  public:
    /// Constructor: just initialise the member data to their defaults (zeros)
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

    C1CurviLine* nodal_boundary_parametrisation_pt(
      const unsigned& j_node)
    {
      return Nodal_boundary_parametrisation_pt[j_node];
    }

   // hierher isn't this the same as in fkv?

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

        // Get position (debug) // hierher why debug?
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


        // // [zdec] debug
        // std::ofstream jac_and_hess;
        // jac_and_hess.open("jac_and_hess_new.csv", std::ios_base::app);
        // jac_and_hess << "Jacobian inverse:" << std::endl
        // 		   << bi[0][0] << " " << bi[0][1] << std::endl
        // 		   << bi[1][0] << " " << bi[1][1] << std::endl
        // 		   << "Hessian inverse [x]:" << std::endl
        // 		   << Dbi[0](0,0) << " " << Dbi[0](0,1) << std::endl
        // 		   << Dbi[0](1,0) << " " << Dbi[0](1,1) << std::endl
        // 		   << "Hessian inverse [y]:" << std::endl
        // 		   << Dbi[1](0,0) << " " << Dbi[1](0,1) << std::endl
        // 		   << Dbi[1](1,0) << " " << Dbi[1](1,1) << std::endl <<
        // std::endl;


        // // [zdec] debug
        // std::ofstream jac_and_hess;
        // jac_and_hess.open("jac_and_hess_new.csv", std::ios_base::app);
        // jac_and_hess << "Jacobian inverse:" << std::endl
        //              << jac_inv(0, 0) << " " << jac_inv(0, 1) << std::endl
        //              << jac_inv(1, 0) << " " << jac_inv(1, 1) << std::endl
        //              << "Hessian inverse [x]:" << std::endl
        //              << hess_inv[0](0, 0) << " " << hess_inv[0](0, 1)
        //              << std::endl
        //              << hess_inv[0](1, 0) << " " << hess_inv[0](1, 1)
        //              << std::endl
        //              << "Hessian inverse [y]:" << std::endl
        //              << hess_inv[1](0, 0) << " " << hess_inv[1](0, 1)
        //              << std::endl
        //              << hess_inv[1](1, 0) << " " << hess_inv[1](1, 1)
        //              << std::endl
        //              << std::endl;
        // jac_and_hess.close();

        // // [zdec] debug
        // std::ofstream debug_stream;
        // debug_stream.open("norm_and_tan.dat", std::ios_base::app);
        // debug_stream << x[0] << " " << x[1] << " " << ni[0] << " " << ni[1]
        //              << " " << ti[0] << " " << ti[1] << " " << dnids[0] << "
        //              "
        //              << dnids[1] << " " << dtids[0] << " " << dtids[1] << " "
        //              << d2xids2[0] << " " << d2xids2[1] << std::endl;
        // debug_stream.close();

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

} // namespace oomph

#endif
