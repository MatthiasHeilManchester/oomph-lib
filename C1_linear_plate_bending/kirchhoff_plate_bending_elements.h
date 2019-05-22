// Header file for the Biharmonic Bell elements
#ifndef OOMPH_MYBIHARMONIC_ELEMENTS_HEADER
#define OOMPH_MYBIHARMONIC_ELEMENTS_HEADER


#include<sstream>

//OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"
#include "../C1_basis/SubparametricTElement.h"

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
class KirchhoffPlateBendingEquations : public virtual FiniteElement 
{
public:
 /// \short Function pointer to pressure function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*PressureFctPt)(const Vector<double>& x, double& f);

 /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
 /// x is a Vector!
 typedef void (*PressureFctGradientPt)(const Vector<double>& x,
                                            Vector<double>& gradient);

 /// Constructor (must initialise the Pressure_fct_pt to null)
 KirchhoffPlateBendingEquations() : Pressure_fct_pt(0), Pressure_fct_gradient_pt(0)
    {}

 /// Broken copy constructor
 KirchhoffPlateBendingEquations(const KirchhoffPlateBendingEquations& dummy)
  {
   BrokenCopy::broken_copy("KirchhoffPlateBendingEquations");
  }

 /// Broken assignment operator
 void operator=(const KirchhoffPlateBendingEquations&)
  {
   BrokenCopy::broken_assign("KirchhoffPlateBendingEquations");
  }

 /// \short Return the index at which the unknown value
 /// is stored.
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 virtual inline unsigned u_index_biharmonic(const unsigned& l) const 
  {return  l;}

 /// Output with default number of plot points
 void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   KirchhoffPlateBendingEquations::output(outfile,n_plot);
  }

 /// \short Output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   KirchhoffPlateBendingEquations::output(file_pt,n_plot);
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


 /// Dummy, time dependent error checker
 void compute_error(std::ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time, double& error, double& norm)
  {
   throw OomphLibError(
    "There is no time-dependent compute_error() for these elements",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

 /// Access function: Pointer to pressure function
 PressureFctPt& pressure_fct_pt() {return Pressure_fct_pt;}

 /// Access function: Pointer to pressure function. Const version
 PressureFctPt pressure_fct_pt() const {return Pressure_fct_pt;}

 /// Access function: Pointer to gradient of pressure function
 PressureFctGradientPt& pressure_fct_gradient_pt()
  {return Pressure_fct_gradient_pt;}

 /// Access function: Pointer to gradient pressure function. Const version
 PressureFctGradientPt pressure_fct_gradient_pt() const
  {return Pressure_fct_gradient_pt;}

 ///Access function to the Poisson ratio.
 const double*& nu_pt() {return Nu_pt;}

 ///Access function to the Poisson ratio (const version)
 const double& get_nu() const {return *Nu_pt;}

 // Get the kth dof type at internal point l
 virtual double get_w_bubble_dof(const unsigned& l, const unsigned& k) const =0;

 // Get the kth equation at internal point l
 virtual int local_w_bubble_equation(const unsigned& l, const unsigned& k)
   const =0;

 // Get the number of basis functions, pure virtual
 virtual unsigned nnodal_basis_type() const =  0;

 // Get the number of internal basis functions, pure virtual
 virtual unsigned nbubble_basis() const = 0;

 // Get the number of internal basis functions, pure virtual
 virtual unsigned nbubble_basis_type() const = 0;

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline virtual void get_pressure_biharmonic(const unsigned& ipt,
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

 /// Add the element's contribution to its residual vector (wrapper)
 virtual void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_biharmonic(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_biharmonic(residuals,jacobian,1);
  }

 /// \short Return FE representation of unknown values u(s)
 /// at local coordinate s
 inline void interpolated_u_biharmonic(const Vector<double> &s, Vector<double>&
interpolated_u) const
  {
   //Find number of position dofs
   const unsigned n_basis_type = this->nnodal_basis_type();
   // Find the internal dofs
   const unsigned n_b_position_type = this->nbubble_basis_type();
   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   //Find out how many internal points there are
   const unsigned n_b_node = this->nbubble_basis();

   //Get the index at which the unknown is stored
   //const unsigned u_length= 6; //Hierher

   //Local c1-shape funtion
   Shape psi(n_node,n_basis_type),test(n_node,n_basis_type),
    psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);

   DShape dpsi_dxi(n_node,n_basis_type,dim()),dtest_dxi(n_node,n_basis_type,dim()),
    dpsi_b_dxi(n_b_node,n_b_position_type,dim()),dtest_b_dxi(n_b_node,n_b_position_type,dim()),
    d2psi_dxi2(n_node,n_basis_type,3), d2test_dxi2(n_node,n_basis_type,3),
    d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

   //Initialise value of u
   //Find values of c1-shape function
   d2shape_and_d2test_eulerian_biharmonic(s,psi,psi_b,dpsi_dxi,dpsi_b_dxi,
    d2psi_dxi2,d2psi_b_dxi2,test,test_b,dtest_dxi,dtest_b_dxi,d2test_dxi2,
    d2test_b_dxi2);

   //Interpolated unknown
   for(unsigned l=0;l<n_node;l++)
   {
    for(unsigned k=0;k<n_basis_type;k++)
     {
      // Get the kth nodal value at node l for displacement i
      double u_value_kl = this->raw_nodal_value(l,u_index_biharmonic(k));
      interpolated_u[0] += u_value_kl*psi(l,k);
      interpolated_u[1] += u_value_kl*dpsi_dxi(l,k,0);
      interpolated_u[2] += u_value_kl*dpsi_dxi(l,k,1);
      interpolated_u[3] += u_value_kl*d2psi_dxi2(l,k,0);
      interpolated_u[4] += u_value_kl*d2psi_dxi2(l,k,1);
      interpolated_u[5] += u_value_kl*d2psi_dxi2(l,k,2);
     }
   }

   // Bubble dofs
   for(unsigned l=0;l<nbubble_basis();l++)
   {
    for(unsigned k=0;k<nbubble_basis_type();k++)
     {
      double u_value = get_w_bubble_dof(l,k);
      interpolated_u[0] += u_value * psi_b(l,k);
      interpolated_u[1] += u_value*dpsi_b_dxi(l,k,0);
      interpolated_u[2] += u_value*dpsi_b_dxi(l,k,1);
      interpolated_u[3] += u_value*d2psi_b_dxi2(l,k,0);
      interpolated_u[4] += u_value*d2psi_b_dxi2(l,k,1);
      interpolated_u[5] += u_value*d2psi_b_dxi2(l,k,2);
     }
   }
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();

protected:
 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  DShape &d2psi_dx2,DShape& d2psi_b_dx2,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx,
  DShape &d2test_dx2,DShape& d2test_b_dx2) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx) const=0;

 /// \short Shape/test functions at local coordinate s
 virtual void shape_and_test_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, Shape &test, Shape& test_b) const=0;

 /// \short Compute element residual Vector only (if flag=and/or element
 /// Jacobian matrix
 virtual void fill_in_generic_residual_contribution_biharmonic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);

 /// Pointer to pressure function:
 PressureFctPt Pressure_fct_pt;

 /// Pointer to gradient of pressure function
 PressureFctGradientPt Pressure_fct_gradient_pt;

 /// Pointer to Poisson ratio, which this element cannot modify
 const double* Nu_pt;
};
} //end namespace oomph
#endif
