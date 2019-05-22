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

//oomph-lib headers
#include "../generic/shape.h"
#include "C1_curved_elements.h"
#include "my_geom_object.h"
#include "../generic/prettyprint98.h"

#ifndef SUBPARAMETRIC_TELEMENTS
#define SUBPARAMETRIC_TELEMENTS
namespace oomph
{
//==============================================================================
// SubparametricTElements are subparametric finite elements that provide both
// basis functions to interpolate the unknowns and shape functions to 
// interpolate the geometry. In general these elements have both nodal basis
// functions and bubble basis functions (corresponding to internal dofs)
//==============================================================================

/// \short Class for subparametric element that has a seprate basis for the 
// unknowns and a shape for the interpolation of the coordinates.
class SubparametricTElement : public TElement<2,2>
{
public:
 /// Constructor
 SubparametricTElement() {}; 

 /// Destructor
 ~SubparametricTElement() {}; 

 /// Return total number of basis functions 
 virtual unsigned nbasis()const 
   {return nvertex_node()*nnodal_basis_type()+nbubble_basis()*nbubble_basis_type();} 

 /// Return number of basis function types 
 virtual unsigned nnodal_basis_type()const = 0 ; 

 /// Return number of bubble basis functions 
 virtual unsigned nbubble_basis()const = 0 ; 

 /// Return number of bubble basis functions 
 virtual unsigned nbubble_basis_type() const = 0 ; 
  
 /// Return number of internal dofs (may be different from basis functions in
 /// general if multiple unknowns are interpolated by SubparametricTElements)
 virtual unsigned ninternal_dofs() const = 0 ; 

 /// Get the itype th bubble dof at internal point ibdof
 virtual double get_bubble_dof(const unsigned ibdof, const unsigned
   itype) const = 0;

 /// Get the basis for the  unknowns
 virtual void basis(const Vector<double>& s, Shape& nodal_basis, Shape& bubble_basis) const = 0;

 /// Get the local derivative of the basis for the unknowns
 virtual void d_basis_local(const Vector<double>& s, Shape& nodal_basis, Shape& bubble_basis, DShape& dnodal_basis,
    DShape& dbubble_basis) const = 0;

 /// Get the local second derivative of the basis for the unknowns
 virtual void d2_basis_local(const Vector<double>& s, Shape& nodal_basis, 
   Shape& bubble_basis, DShape& dnodal_basis, DShape& dbubble_basis, DShape&
    d2nodal_basis, DShape& d2bubble_basis) const  = 0;

 /// Get the global (Eulerian) derivative of the basis for the unknowns
 virtual double d_basis_eulerian(const Vector<double>& s, Shape& nodal_basis, Shape& bubble_basis, DShape& dnodal_basis,
    DShape& dbubble_basis) const = 0;

 /// Get the global (Eulerian) second derivative of the basis for the unknowns
 virtual double d2_basis_eulerian(const Vector<double>& s, Shape& nodal_basis, 
   Shape& bubble_basis, DShape& dnodal_basis, DShape& dbubble_basis, DShape&
    d2nodal_basis, DShape& d2bubble_basis) const  = 0;

}; //End of Subparametric TElement class

/// \short Curvable Bell element. By default this element is the standard 
/// simplex Bell element with 6dofs per node. The element can be upgraded
/// to provide accurate representation of boundary conditions, which adds
/// additional (bubble) unknowns to the interior.
class CurvableBellElement : public SubparametricTElement
{
public:
 /// Constructor
 CurvableBellElement() : Curved_edge(MyC1CurvedElements::none)
  { 
   // Add the (zero) bubble dofs
   Internal_data_index = this->add_internal_data(new Data(0));
   Bernadou_element_basis_pt=0;
   Association_matrix_pt=0;
  //  // Store the vertices (TMP) for Bell
  //  Verts = (Vector<Vector<double> >(nvertex_node(),Vector<double>(dim(),0.0)));
  //  // Fill in
  //  for(unsigned ivert=0;ivert<nvertex_node();++ivert)
  //    {
  //    Verts[ivert][0] =nodal_position(ivert,0);
  //    Verts[ivert][1] =nodal_position(ivert,1);
  //    }
  // std::cout<<"{"<<Verts[0][0]<<","<<Verts[0][1]<<"} ";
  // std::cout<<"{"<<Verts[1][0]<<","<<Verts[1][1]<<"} ";
  // std::cout<<"{"<<Verts[1][0]<<","<<Verts[1][1]<<"}\n";
  };

 ///Destructor 
 ~CurvableBellElement()
  {
   // Clean Up
   delete Association_matrix_pt;
   delete Bernadou_element_basis_pt;
  }

 /// \short Alias for enum to enumerate the possible edges that could be curved
 typedef typename MyC1CurvedElements::Edge Edge; 

 ///  Boolean function indicating whether element is curved or not
 bool element_is_curved() const {return Curved_edge != MyC1CurvedElements::none; }

 /// \short get the coordinate i
 void interpolated_x (const Vector<double>& s, Vector<double>& x) const 
  {
   // Wrapper: call the BernadouElement curved basis if curved 
   if(element_is_curved()) 
    {Bernadou_element_basis_pt->coordinate_x(s,x);}
   else 
    {TElement<2,2>::interpolated_x(s,x);}
  }

 /// Overloaded shape. Thin wrapper which breaks shape for upgraded elements,
 /// which have a mapping but no shape function defined.
 virtual void shape(const Vector<double>& s, Shape& shape)
  {
   if(element_is_curved())
     { 
      throw OomphLibError("No shape defined for these elements: use interpolated_x \
to access interpolated eulerian coordinate",
       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
     }
   // Default to TElement shape
   else
     {TElement<2,2>::shape(s,shape);}
  }

 /// Overloaded shape. Thin wrapper which breaks shape for upgraded elements,
 /// which have a mapping but no shape function defined.
 virtual void dshape_local(const Vector<double>& s, Shape& shape, DShape& dshape)
  {
   if(element_is_curved())
     { 
      throw OomphLibError("No shape defined for these elements: use interpolated_x \
to access interpolated eulerian coordinate",
       OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
     }
   else
     {TElement<2,2>::dshape_local(s,shape,dshape);}
  }

 /// Local_to_eulerian mapping with local coordinate argument: when upgraded this 
 /// uses the Bernadou implementation of the Jacobian 
 virtual double local_to_eulerian_mapping(const Vector<double>& s, DenseMatrix<double>& jacobian,
   DenseMatrix<double>& inverse_jacobian)
  {
   if(element_is_curved())
    { 
     // Fill in the Jacobian
     Bernadou_element_basis_pt->get_jacobian(s,jacobian); 
     // Invert the Jacobian and return the determinant
     return invert_jacobian<2>(jacobian,inverse_jacobian);
    }
   else
     {
      // Assemble dshape to get the mapping
      Shape psi(nvertex_node());
      DShape dpsi(nvertex_node(),dim());
      TElement<2,2>::dshape_local(s,psi,dpsi);
      return TElement<2,2>::local_to_eulerian_mapping(dpsi,jacobian,inverse_jacobian);
     }
  }

 /// Get the basis for the  unknowns
 virtual void basis(const Vector<double>& s, Shape& nodal_basis, Shape& bubble_basis) const
   {
    if(element_is_curved())
      { 
       Bernadou_element_basis_pt->shape(s,nodal_basis,bubble_basis);
      }
    else
      { //HERE
      Vector<Vector<double> > Verts = (Vector<Vector<double> >(nvertex_node(),Vector<double>(dim(),0.0)));
      for(unsigned ivert=0;ivert<nvertex_node();++ivert)
         for(unsigned icoord=0;icoord<dim();++icoord)
           Verts[ivert][icoord] =nodal_position(ivert,icoord);
      DShape dummydshape(nvertex_node(),nnodal_basis_type(),dim());
      DShape dummyd2shape(nvertex_node(),nnodal_basis_type(),dim()*dim()-1);
      Bell_basis.d2_basis_eulerian(s,Verts,nodal_basis,dummydshape,dummyd2shape);
      }
   };

 /// Get the local derivative of the basis for the unknowns
 virtual void d_basis_local(const Vector<double>& s, Shape& nodal_basis, Shape& bubble_basis, DShape& dnodal_basis,
    DShape& dbubble_basis) const
  {
   throw OomphLibError("Needs implementing. BLAME DAVID.",OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

 /// Get the local second derivative of the basis for the unknowns
 virtual void d2_basis_local(const Vector<double>& s, Shape& nodal_basis, 
   Shape& bubble_basis, DShape& dnodal_basis, DShape& dbubble_basis, DShape&
    d2nodal_basis, DShape& d2bubble_basis) const
  {
   throw OomphLibError("Needs implementing. BLAME DAVID.",OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

 /// Get the global (Eulerian) derivative of the basis for the unknowns
 virtual double d_basis_eulerian(const Vector<double>& s, Shape& nodal_basis, Shape& bubble_basis, DShape& dnodal_basis,
    DShape& dbubble_basis) const
    {
    if(element_is_curved())
      { 
       return Bernadou_element_basis_pt->d_shape_dx(s,nodal_basis, bubble_basis,
        dnodal_basis, dbubble_basis);
      }
    else
      { //HERE
      Vector<Vector<double> > Verts = (Vector<Vector<double> >(nvertex_node(),Vector<double>(dim(),0.0)));
      for(unsigned ivert=0;ivert<nvertex_node();++ivert)
         for(unsigned icoord=0;icoord<dim();++icoord)
           Verts[ivert][icoord] =nodal_position(ivert,icoord);
      DShape dummyd2shape(nvertex_node(),nnodal_basis_type(),dim()*dim()-1);
      Bell_basis.d2_basis_eulerian(s,Verts,nodal_basis,dnodal_basis,dummyd2shape);
      return TElement<2,2>::J_eulerian(s);
      }
     }

 /// Get the global (Eulerian) second derivative of the basis for the unknowns
 virtual double d2_basis_eulerian(const Vector<double>& s, Shape& nodal_basis, 
   Shape& bubble_basis, DShape& dnodal_basis, DShape& dbubble_basis, DShape&
    d2nodal_basis, DShape& d2bubble_basis) const  
    {
    if(element_is_curved())
      {
       // IF the association matrix_pt is equal to null pointer we don't have one
       // so use the `slow' version
       if(Association_matrix_pt==0) 
        {
         return Bernadou_element_basis_pt->d2_shape_dx2(s,nodal_basis,bubble_basis, dnodal_basis,
         dbubble_basis,d2nodal_basis,d2bubble_basis);
        }
       // Else use the cached association matrix to compute
       else 
        {
         return Bernadou_element_basis_pt->d2_shape_dx2(s,nodal_basis,bubble_basis, dnodal_basis,
         dbubble_basis,d2nodal_basis,d2bubble_basis,*Association_matrix_pt);
        }
      }
     // Use the Bell basis functions if not upgraded
    else
      { 
      Vector<Vector<double> > Verts = (Vector<Vector<double> >(nvertex_node(),Vector<double>(dim(),0.0)));
      for(unsigned ivert=0;ivert<nvertex_node();++ivert)
         for(unsigned icoord=0;icoord<dim();++icoord)
           Verts[ivert][icoord] =nodal_position(ivert,icoord);
      Bell_basis.d2_basis_eulerian(s,Verts,nodal_basis,dnodal_basis,d2nodal_basis);
      return TElement<2,2>::J_eulerian(s);
      }
     }

 /// \short Precompute the association matrix for the curved basis. This is 
 /// important for optimisation, as its contruction is expensive. It only 
 /// needs constructing once per get_residuals or get_jacobian call.
 void store_association_matrix()
   {
    // If the element has been upgraded
    if(element_is_curved())
     { 
      // Create the matrix
      Association_matrix_pt =  new DenseMatrix<double>(Bernadou_element_basis_pt->n_basis_functions(),
        Bernadou_element_basis_pt->n_basic_basis_functions(),0.0);
      // Fill in the matrix
      Bernadou_element_basis_pt->fill_in_full_association_matrix(*Association_matrix_pt);
     }
    else
     { /* No association matrix for the Bell Elements. Throw here*/}
   }; 

 /// \short Delete the association matrix for the curved basis. This is 
 /// important for optimisation, as we don't want to be storing a large matrix 
 /// for every element. 
 void delete_association_matrix()
   {
    // If the element has been upgraded
    if(element_is_curved())
     { 
      delete Association_matrix_pt;
      Association_matrix_pt=0;
     }
    else
     { /* No association matrix for the Bell Elements. Throw here ?*/}
   }; 

 /// \short Return number of bubble basis functions. This will be zero
 /// for un-upgraded elements. Upgrading introduces nbubble_dof additional basis 
 /// functions depending on the basis 
 unsigned nbubble_basis() const 
  {return (element_is_curved() ? Bernadou_element_basis_pt->n_internal_dofs() : 0);}; 

 // Return number of bubble dof types
 unsigned nbubble_basis_type() const {return 1;}; 

 /// Return number of nodal dof types 
 unsigned nnodal_basis_type()const { return 6; }; 
  
 /// \short Return number of nodal dof types. The default implementation in these
 /// elements is to use nbasis_type but in derived classes these two values may
 /// not be the same
 virtual unsigned ndof_type()const { return nnodal_basis_type(); }; 

 /// \short Return number of internal dofs. This will be zero
 /// for un-upgraded elements. The default implementation in these elements
 /// adds nbubble_basis additional basis functions
 virtual unsigned ninternal_dofs() const { return nbubble_basis(); }; 

 /// Get the bubble dof
 virtual double get_bubble_dof(const unsigned ibdof, const unsigned
   itype) const 
  {
  // Deliberately break this function for the below cases
  // If there is no curved edge then we cannot return anything meaningful
  if(Curved_edge==MyC1CurvedElements::none)
  {
  throw OomphLibError("There are no time-dependent internal 'bubble' dofs for \
this element.",OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  // Return dummy value 0.0 to shut up compiler
  return 0;
  }
  // For these elements we only have a single dof at each internal point
  else if(itype!=0)
  {
  throw OomphLibError(
   "There is only a single degree of freedom at the internal points in this \
element.", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  // Return dummy value 0.0 to shut up compiler
  return 0;
  }
  // Now give the lth internal degree of freedom
  else
   {return this->internal_data_pt(Internal_data_index)->value(ibdof);}
  }
 
 /// Index of the internal data
 virtual unsigned index_of_internal_data() const 
  {return Internal_data_index;}
  
 /// Add the a curved element pointer of type BERNADOU_BASIS
 template<typename BERNADOU_BASIS>
 void add_new_curved_basis()
  {
   if(Bernadou_element_basis_pt ==0)
    {Bernadou_element_basis_pt =  new BERNADOU_BASIS;}
   else
    {
    throw OomphLibError(
    "There is already a curved basis for this curvable Bell element.",
     OOMPH_CURRENT_FUNCTION,  OOMPH_EXCEPTION_LOCATION);
    }
  } 

 /// Upgrade the element to be curved
 virtual void upgrade_element_to_curved(const Edge& curved_edge, const double& s_ubar,
  const double& s_obar,  CurvilineGeomObject* parametric_edge, const unsigned& boundary_order)
  {
   using namespace MyC1CurvedElements;
  #ifdef PARANOID
  // Check that we haven't upgraded this element already
  if(Curved_edge != none)
   {
    throw OomphLibError(
    "Cannot upgrade more than a single edge to be curved in C1 Curved Bell \
Elements.",OOMPH_CURRENT_FUNCTION,  OOMPH_EXCEPTION_LOCATION);
   }
  #endif
  // Add the curved edge
  Curved_edge = curved_edge;

  Integral* new_integral_pt; 
  // Switch for the boundary order
  switch(boundary_order)
   {
   case 3:
    add_new_curved_basis<BernadouElementBasis<3> >();
    new_integral_pt = new TGauss<2,13>;
   break;
   case 5:
    add_new_curved_basis<BernadouElementBasis<5> >();
    new_integral_pt = new TGauss<2,16>;
   break;
   default:
    throw OomphLibError(
    "Currently only BernadouElementBasis<3> and BernadouElementBasis<5> are implemented."
     ,OOMPH_CURRENT_FUNCTION,  OOMPH_EXCEPTION_LOCATION);
   } 
   this->set_integration_scheme(new_integral_pt); 
  // Set the number of internal dofs to nbubble
  // Number_of_internal_dof_types = 1;
  // Number_of_internal_dofs = basis_pt->ninternal_dofs;
  // Nbubble_basis =  nbubble;
  unsigned n_bubble = Bernadou_element_basis_pt->n_internal_dofs();
  Internal_data_index = this->add_internal_data(new Data(n_bubble));
 
  // Vector to store nodes
  Vector<Vector<double> > vertices(Vector<Vector<double> >(nvertex_node(),Vector<double>(dim(),0.0)));
  // Now switch to upgrade
  // The shape functions are designed such that the curved edge is always edge
  // two. So this is where we set that up. 
  // Throw an error if an edge is upgraded to none
   if(Curved_edge== none)
    {
     throw OomphLibError( "Cannot upgrade edge 'none'. Curved elements must have\
 one side defined by a parametric function.", OOMPH_CURRENT_FUNCTION,  
      OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
    // Loop over coordinate components then vertices and permute vertices
    for(unsigned i=0;i<2;++i)
     {
     for(unsigned ivert=0;ivert<nvertex_node();++ivert)
      {vertices[ivert][i]=this->node_pt((ivert+unsigned(curved_edge)+1)%3)->x(i);}
     }
    }
  // Add the vertices to make the shape functions fully functional
  Bernadou_element_basis_pt->upgrade_element(vertices, s_ubar,s_obar,
   curved_edge,*parametric_edge);
 }

private:
 /// Enum to store which edge is curved set to none when element has no curved
 /// edges
 MyC1CurvedElements::Edge Curved_edge;

 /// Index at which the added internal data is storedj
 unsigned Internal_data_index;

 /// Pointer to Bernadou Element Basis
 MyC1CurvedElements::BernadouElementBasisBase* Bernadou_element_basis_pt;

 /// Basis functions
 MyShape::BellElementBasis Bell_basis;

 /// Pointer to Stored Association matrix
 DenseMatrix<double>* Association_matrix_pt;

 /// HERE these neeedn't be here
// Vector<Vector<double> > Verts;

}; //End of Subparametric TElement class


} //End of namespace extension
#endif
