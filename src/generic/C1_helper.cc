#include "C1_helper.h"

namespace oomph
{
 
//==============================================================================
/// Namespace to deal update triangle meshes to deal with C1 elements
//==============================================================================
namespace C1Helper
{
 

 
 /// Boolean to suppress warning about polygonal boundaries
 bool Do_not_warn_about_polygonal_boundaries=false;
  
 /// Output stream to document split elements
 std::ofstream Split_elements_output_stream;
     
 /// Output stream to document duplicated nodes
 std::ofstream Duplicated_node_output_stream;
     
 /// Output stream to document elements whose boundaries have been
 /// upgraded to become curved
 std::ofstream Upgraded_to_curved_edge_element_stream;

 
//==============================================================================
// hierher update
/// Duplicate nodes at corners in order to properly apply boundary
/// conditions from each edge. Also adds (8) Lagrange multiplier dofs to the
/// problem in order to constrain continuous interpolation here across its (8)
/// vertex dofs. (Note "corner" here refers to the meeting point of any two
/// sub-boundaries in the closed external boundary)
//==============================================================================
 void duplicate_corner_nodes(Mesh* bulk_mesh_pt, 
                             std::map<unsigned,C1CurviLine*> c1_curviline_pt,
                             Mesh* constraint_mesh_pt)
 {

  // hierher check if mesh is distributed!
  
  // Collection of nodes that occupy two boundaries together with the boundary IDs
  // (ordered: first < second)
  std::map<Node*,std::pair<unsigned,unsigned>> boundaries_of_boundary_node_pt;

  // Loop over the curvilinear parts of the outer boundary
  for (const auto& [i_bound, para] : c1_curviline_pt)
  {
   unsigned n_b_node = bulk_mesh_pt->nboundary_node(i_bound);
   for(unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
    {
     // Store the node we are checking
     Node* node_pt = bulk_mesh_pt->boundary_node_pt(i_bound,i_b_node);

     // Pointer to set that contains the boundaries we're on
     std::set<unsigned>* boundaries_pt=0;
     node_pt->get_boundaries_pt(boundaries_pt);
     if (boundaries_pt!=0)
      {
       if (boundaries_pt->size()==2)
        {
         unsigned b_min=UINT_MAX;
         unsigned b_max=0;
         for (unsigned b : (*boundaries_pt))
          {
           if (b<b_min) b_min=b;
           if (b>b_max) b_max=b;           
          }
         // Ordered!
         boundaries_of_boundary_node_pt[node_pt].first=b_min;
         boundaries_of_boundary_node_pt[node_pt].second=b_max;
        }
      }
    }
  }
  
  // Here are the nodes that need to be duplicated. We duplicate them
  // on the lower of its two boundaries (this is stored first)
  for (auto a : boundaries_of_boundary_node_pt)
   {
    Node* node_to_be_duplicated_pt=a.first;
    unsigned boundary_on_which_node_is_duplicated=a.second.first;
    unsigned boundary_on_which_node_is_left=a.second.second;

    // Doc?
    if (Duplicated_node_output_stream.is_open())
     {
      Duplicated_node_output_stream
       << node_to_be_duplicated_pt->x(0) << " "
       << node_to_be_duplicated_pt->x(1) << " "
       << boundary_on_which_node_is_duplicated << " "
       << boundary_on_which_node_is_left
       << std::endl;
     }
    
    // Find the boundary element that contains the node to be duplicated on
    // the boundary where the node is to be duplicated
    FiniteElement* el_where_node_is_to_be_duplicated_pt=0;
    unsigned n_b_el = bulk_mesh_pt->nboundary_element(boundary_on_which_node_is_duplicated);
    for (unsigned i_b_el = 0; i_b_el < n_b_el; i_b_el++)
    {
      // Get the element pointer
      FiniteElement* el_pt = bulk_mesh_pt->boundary_element_pt
       (boundary_on_which_node_is_duplicated, i_b_el);
      // If the corner node pt is in the element we have found the right
      // element
      if (el_pt->get_node_number(node_to_be_duplicated_pt) != -1)
       {
        el_where_node_is_to_be_duplicated_pt = el_pt;
        break;
      }
    }
    
    // Now we need to create a new node and substitute the element's
    // old corner node for this new one
    Node* new_node_pt = el_where_node_is_to_be_duplicated_pt->construct_boundary_node(
     el_where_node_is_to_be_duplicated_pt->get_node_number(node_to_be_duplicated_pt));
    
    // Copy the position and other info from the old node into the new node
    new_node_pt->x(0)=node_to_be_duplicated_pt->x(0);
    new_node_pt->x(1)=node_to_be_duplicated_pt->x(1);

    // Then we add this node to the mesh
    bulk_mesh_pt->add_node_pt(new_node_pt);

    // Then replace the old node for the new one on the boundary
    bulk_mesh_pt->remove_boundary_node(boundary_on_which_node_is_duplicated,node_to_be_duplicated_pt);
    bulk_mesh_pt->   add_boundary_node(boundary_on_which_node_is_duplicated,new_node_pt);

    // The final job is to constrain this duplication using the specialised
    // Lagrange multiplier elements which enforce equality of displacement and
    // its derivatives either side of this corner.
    C1CurviLine* left_parametrisation_pt  = c1_curviline_pt[boundary_on_which_node_is_left];
    C1CurviLine* right_parametrisation_pt = c1_curviline_pt[boundary_on_which_node_is_duplicated];

    // Get the coordinates on each node on their respective boundaries
    Vector<double> left_boundary_coordinate =
     {left_parametrisation_pt->get_zeta(node_to_be_duplicated_pt->position())};
    Vector<double> right_boundary_coordinate =
     {right_parametrisation_pt->get_zeta(new_node_pt->position())};

    // Create the constraining element using the first bulk element
    // in them mesh
    TemplateFreeCurvableBellElement* bulk_el_pt=dynamic_cast<TemplateFreeCurvableBellElement*>
     (bulk_mesh_pt->element_pt(0));
    
    DuplicateNodeConstraintElement* constraint_element_pt =
     bulk_el_pt->duplicate_constraint_element_factory(node_to_be_duplicated_pt,
                                                      new_node_pt,
                                                      left_parametrisation_pt,
                                                      right_parametrisation_pt,
                                                      left_boundary_coordinate,
                                                      right_boundary_coordinate);
    
    // Add the constraining element to the mesh
    constraint_mesh_pt->add_element_pt(constraint_element_pt);
   }   
 }



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
//=============================================================================
 void upgrade_edge_elements_to_curved_boundaries(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt) 
 {
  
  // Loop over the curvilinear parts of the outer boundary
  for (const auto& [ibound, c1_curve_pt] : c1_curviline_pt)
   {
    
    // Loop over the bulk elements adjacent to boundary ibound
    const unsigned n_els=bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<n_els; e++)
     {
      // Get pointer to bulk element adjacent to b
      FiniteElement* bulk_el_pt =
       bulk_mesh_pt->boundary_element_pt(ibound,e);

      // Initialise enum for the curved edge
      C1Helper::CurvedEdgeEnumeration edge=C1Helper::CurvedEdgeEnumeration::none;
      
      // Loop over all (three) vertex nodes of the element and
      // identify single node that is interior (i.e. not on any
      // of the outer boundaries
      unsigned index_of_interior_node = 3;
      unsigned nnode_not_on_any_outer_boundary = 0;
      const unsigned nnode = 3;
      Vector<Vector<double> > xn(nnode,Vector<double>(2,0.0));
      for(unsigned n=0;n<nnode;++n)
       {
        Node* nod_pt = bulk_el_pt->node_pt(n);
        xn[n][0]=nod_pt->x(0);
        xn[n][1]=nod_pt->x(1);
        
        // Check if it is on any of the outer boundaries
        bool node_is_on_some_outer_boundary=false;
        for (const auto& [b, dummy_c1_curve_pt] : c1_curviline_pt)
         {
          if (nod_pt->is_on_boundary(b))
           {
            node_is_on_some_outer_boundary=true;
            break;
           }
         }
        if (!node_is_on_some_outer_boundary)
         {
          index_of_interior_node = n;
          nnode_not_on_any_outer_boundary++;
         }
       }// end record boundary nodes
      
      // hierher shouldn't these be called zeta (everywhere; sigh)
      // boundary coordinate at the next (cyclic) node after interior
      const double s_ubar =
       c1_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);
      
      // boundary coordinate at the previous (cyclic) node before interior
      const double s_obar =
       c1_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);
      
      // Assign edge case
      edge = static_cast<C1Helper::CurvedEdgeEnumeration>(index_of_interior_node);
      
#ifdef PARANOID
      // Check nnode_on_neither_boundary
      if (nnode_not_on_any_outer_boundary == 0)
       {
        throw OomphLibError(
         "No interior nodes. One node per CurvedElement must be interior.",
         OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
       }
      else if (nnode_not_on_any_outer_boundary> 1)
       {
        throw OomphLibError(
         "Multiple interior nodes. Only one node per CurvedElement can be interior.",
         OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
       }
      
      // Check for inverted elements
      if (s_ubar>s_obar)
       {
        throw OomphLibError(
         "Decreasing parametric coordinate. Parametric coordinate must increase as the edge is traversed anti-clockwise.",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       } // end checks
#endif
      
      // Upgrade it
      TemplateFreeCurvableBellElement* curv_el_pt=
       dynamic_cast<TemplateFreeCurvableBellElement*>(bulk_el_pt);
#ifdef PARANOID
      if (curv_el_pt==0)
       {
        throw OomphLibError(
         "Cast to TemplateFreeCurvableBellElement failed",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
#endif

      // hierher shouldn't we hard code this to only allow specific options
      // this can't be any number, right?
      unsigned boundary_order=5;
      curv_el_pt->upgrade_element_to_curved(edge, s_ubar, s_obar,
                                            c1_curve_pt,
                                            boundary_order);


      // Doc?
      if (Upgraded_to_curved_edge_element_stream.is_open())
       {
        unsigned nplot=5;
        curv_el_pt->output(Upgraded_to_curved_edge_element_stream,nplot);
       }
      
      
     }

   } // end of loop over outer boundaries
  
 } // end_upgrade_elements


 
//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
//======================================================================
 void rotate_edge_degrees_of_freedom(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt) 
{
 
 // Loop over the bulk elements: Yes, really because we also need to deal with those that only
 // have a single node on the boundary!
 unsigned n_element = bulk_mesh_pt-> nelement();
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element 
   FiniteElement* el_pt = bulk_mesh_pt->finite_element_pt(e);
   
   // Loop over the curvilinear parts of the outer boundary
   for (const auto& [b, c1_curve_pt] : c1_curviline_pt)
    {
     // local node numbers of nodes on external boundaries
     Vector<unsigned> boundary_node;
     
     // Boundary coordinates of nodes on the external boundaries
     Vector<double> boundary_coordinate_of_node;
     
     // Loop over vertex nodes (they come first)
     const unsigned nnode=3;
     for (unsigned n=0; n<nnode;++n)
      {
       // If on external boundary b
       if (el_pt->node_pt(n)->is_on_boundary(b))
        {
         boundary_node.push_back(n);
         double coord = c1_curve_pt->get_zeta(el_pt->node_pt(n)->position());
         boundary_coordinate_of_node.push_back(coord);
        }
      }
     
     // If the element has nodes on the boundary, rotate the Hermite dofs
     if(!boundary_node.empty())
      {
       // Rotate the nodes by passing the index of the nodes and the
       // normal / tangent vectors to the element
       
       // Upgrade it
       TemplateFreeCurvableBellElement* curv_el_pt=
        dynamic_cast<TemplateFreeCurvableBellElement*>(el_pt);
#ifdef PARANOID
       if (curv_el_pt==0)
        {
         throw OomphLibError(
          "Cast to TemplateFreeCurvableBellElement failed",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
       
       curv_el_pt->
        rotated_boundary_helper_pt()->
        set_nodal_boundary_parametrisation(boundary_node,
                                           boundary_coordinate_of_node,
                                           c1_curve_pt);
      }
    }
  }
 
} // end rotate_edge_degrees_of_freedom

 

} // end namespace C1Helper

} // end namespace oomph
