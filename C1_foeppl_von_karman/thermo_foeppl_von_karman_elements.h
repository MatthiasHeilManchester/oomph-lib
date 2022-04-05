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
#ifndef OOMPH_THERMO_FVK_ELEMENTS_HEADER
#define OOMPH_THERMO_FVK_ELEMENTS_HEADER

#include<sstream>

//OOMPH-LIB headers
#include "foeppl_von_karman_elements.h"


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
class ThermoFoepplVonKarmanEquations : public virtual FoepplVonKarmanEquations
{

public:

 /// \short Function pointer to temperature function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*TemperatureFctPt)(const Vector<double>& x, double& f);

 /// Constructor (must initialise the Pressure_fct_pt to null)
 ThermoFoepplVonKarmanEquations() :
  FoepplVonKarmanEquations()
  {
   Alpha_pt = &Default_Alpha_Value;
  }

 /// Broken copy constructor
 ThermoFoepplVonKarmanEquations(const ThermoFoepplVonKarmanEquations& dummy)
  {
   BrokenCopy::broken_copy("ThermoFoepplVonKarmanEquations");
  }

 /// Broken assignment operator
 void operator=(const ThermoFoepplVonKarmanEquations&)
  {
   BrokenCopy::broken_assign("ThermoFoepplVonKarmanEquations");
  }

 /// Fill in the stress tensor
 void get_sigma(DenseMatrix<double>& sigma,
		const DenseMatrix<double>& grad_u,
		const DenseMatrix<double>& grad_w,
		const double& alphaDeltaT)const
  {
   // Poisson ratio
   double nu(get_nu());
   // Truncated Green Lagrange strain tensor
   DenseMatrix<double> epsilon(this->dim(),this->dim(),0.0);
   for(unsigned alpha=0;alpha<this->dim();++alpha)
    {
     for(unsigned beta=0;beta<this->dim();++beta)
      {
       // Truncated Green Lagrange strain tensor
       epsilon(alpha,beta) += 0.5* grad_u(alpha,beta) + 0.5*grad_u(beta,alpha)
                           +0.5*grad_w(0,alpha)*grad_w(0,beta);
      }
    }
   
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
     // Thermal stresses
     sigma(alpha,alpha) += -alphaDeltaT;
    }
 
  }
 /// Access function: Pointer to temperature function
 TemperatureFctPt& temperature_fct_pt() {return Temperature_fct_pt;}

 /// Access function: Pointer to temperature function. Const version
 TemperatureFctPt temperature_fct_pt() const {return Temperature_fct_pt;}

 /// Access function for alpha
 const double &get_alpha() const {return *Alpha_pt;}

 /// Access function for a pointer to eta
 const double* &alpha_pt() {return Alpha_pt;}

 /// Get temperature at (Eulerian) position x. This function is
 /// virtual to allow overloading.
 inline virtual void get_temperature_foeppl_von_karman(const unsigned& ipt,
                                        const Vector<double>& x,
                                        double& temperature) const
  {
   //If no pressure function has been set, return zero
   if(Temperature_fct_pt==0)
    {
     temperature = 0.0;
    }
   else
    {
     // Get pressure strength
     (*Temperature_fct_pt)(x,temperature);
    }
  }

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_thermo_foeppl_von_karman(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_thermo_foeppl_von_karman(residuals,jacobian,1);
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


protected:
 
 /// \short Compute element residual Vector only (if flag=and/or element
 /// Jacobian matrix
 virtual void fill_in_generic_residual_contribution_thermo_foeppl_von_karman(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);
 
 
 /// Pointer to temperature function:
 TemperatureFctPt Temperature_fct_pt;

 /// Pointer to alpha, thermal expansion coefficient
 const double *Alpha_pt;

 /// Default value for material parameter alpha
 static const double Default_Alpha_Value;
 
};


} //end namespace oomph
#endif
