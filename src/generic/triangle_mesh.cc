// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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

#include <algorithm>
#include "map_matrix.h"
#include "unstructured_two_d_mesh_geometry_base.h"
#include "triangle_mesh.h"

namespace oomph
{
#ifdef OOMPH_HAS_TRIANGLE_LIB

  //==============================================================
  /// Dump the triangulateio structure to a dump file and
  /// record boundary coordinates of boundary nodes
  //==============================================================
  void TriangleMeshBase::dump_triangulateio(std::ostream& dump_file)
  {
    TriangleHelper::dump_triangulateio(Triangulateio, dump_file);

#ifdef OOMPH_HAS_MPI
    // If the mesh is not distributed then process what follows
    if (!this->is_mesh_distributed())
    {
#endif // #ifdef OOMPH_HAS_MPI

      // Loop over all boundary nodes and dump out boundary coordinates
      // if they exist
      Vector<double> zeta(1);
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        if (Boundary_coordinate_exists[b])
        {
          dump_file << "1 # Boundary coordinate for boundary " << b
                    << " does exist\n";
          unsigned nnod = nboundary_node(b);
          dump_file << nnod << " # Number of dumped boundary nodes\n";
          for (unsigned j = 0; j < nnod; j++)
          {
            Node* nod_pt = boundary_node_pt(b, j);
            nod_pt->get_coordinates_on_boundary(b, zeta);
            dump_file << zeta[0] << std::endl;
          }
          dump_file << "-999 # Done boundary coords for boundary " << b << "\n";
        }
        else
        {
          dump_file << "0 # Boundary coordinate for boundary " << b
                    << " does not exist\n";
        }
      }

#ifdef OOMPH_HAS_MPI
    }
#endif // #ifdef OOMPH_HAS_MPI
  }


  //==============================================================
  /// Regenerate the mesh from a dumped triangulateio file
  /// and dumped boundary coordinates of boundary nodes
  //==============================================================
  void TriangleMeshBase::remesh_from_triangulateio(std::istream& restart_file)
  {
#ifdef PARANOID
    // Record number of boundaries
    unsigned nbound_old = nboundary();
#endif

    // Clear the existing triangulate io
    TriangleHelper::clear_triangulateio(Triangulateio);

    // Read the data into the file
    TriangleHelper::read_triangulateio(restart_file, Triangulateio);

    // Now remesh from the new data structure
    this->remesh_from_internal_triangulateio();

#ifdef OOMPH_HAS_MPI
    // If the mesh is not distributed then process what follows
    if (!this->is_mesh_distributed())
    {
#endif // #ifdef OOMPH_HAS_MPI

#ifdef PARANOID
      // Record number of boundary nodes after remesh
      unsigned nbound_new = nboundary();
      if (nbound_new != nbound_old)
      {
        std::ostringstream error_stream;
        error_stream
          << "Number of boundaries before remesh from triangulateio, "
          << nbound_new << ",\ndoesn't match number boundaries afterwards, "
          << nbound_old
          << ". Have you messed \naround with boundary nodes in the "
          << "derived mesh constructor (or after calling \nit)? If so,"
          << " the dump/restart won't work as written at the moment.";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif


      // Loop over all boundary nodes and read boundary coordinates
      // if they exist
      Vector<double> zeta(1);
      std::string input_string;
      unsigned nb = nboundary();
      for (unsigned b = 0; b < nb; b++)
      {
        // Read line up to termination sign
        getline(restart_file, input_string, '#');

        // Ignore rest of line
        restart_file.ignore(80, '\n');

        // Did boundary coordinate exist?
        const unsigned bound_coord_exists = atoi(input_string.c_str());
        if (bound_coord_exists == 1)
        {
          // Remember it!
          Boundary_coordinate_exists[b] = true;

          // Read line up to termination sign
          getline(restart_file, input_string, '#');

          // Ignore rest of line
          restart_file.ignore(80, '\n');

          // How many nodes did we dump?
          const unsigned nnod_dumped = atoi(input_string.c_str());

          // Does it match?
          unsigned nnod = nboundary_node(b);
          if (nnod != nnod_dumped)
          {
            std::ostringstream error_stream;
            error_stream << "Number of dumped boundary nodes " << nnod_dumped
                         << " doesn't match number of nodes on boundary " << b
                         << ": " << nnod << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }

          // Loop over all nodes
          for (unsigned j = 0; j < nnod; j++)
          {
            // Read line up to termination sign
            getline(restart_file, input_string);

            // Boundary coordinate
            zeta[0] = atof(input_string.c_str());

            // Set it
            Node* nod_pt = boundary_node_pt(b, j);
            nod_pt->set_coordinates_on_boundary(b, zeta);
          }

          // Read line up to termination sign
          getline(restart_file, input_string, '#');

          // Ignore rest of line
          restart_file.ignore(80, '\n');

          // Have we reached the end?
          const int check = atoi(input_string.c_str());
          if (check != -999)
          {
            std::ostringstream error_stream;
            error_stream << "Haven't read all nodes on boundary " << b
                         << std::endl;
            throw OomphLibError(error_stream.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        }
        else
        {
          oomph_info << "Restart: Boundary coordinate for boundary " << b
                     << " does not exist.\n";
        }
      }
#ifdef OOMPH_HAS_MPI
    } // if (!this->is_mesh_distributed())
#endif // #ifdef OOMPH_HAS_MPI
  }


  //==============================================================
  /// Write a Triangulateio_object file of the TriangulateIO object
  /// String s is add to assign a different value for
  /// input and/or output structure.
  /// The function give the same result of the "report" function
  /// included in the tricall.c, esternal_src.
  //==============================================================
  void TriangleMeshBase::write_triangulateio(TriangulateIO& triangle,
                                             std::string& s)
  {
    std::ofstream outfile;
    char filename[100];

    sprintf(filename, "Triangulateio_object_%s.dat", s.c_str());
    outfile.open(filename);
    outfile << "# Triangulateio object values:\n\n" << std::endl;

    // Write points coordinates
    if (triangle.numberofpoints != 0)
    {
      outfile << "# Triangulateio number of points is:"
              << triangle.numberofpoints << std::endl;
    }
    if (triangle.pointlist != NULL)
    {
      outfile << "# Vertex coordinates are:" << std::endl;
      for (int k = 0; k < triangle.numberofpoints * 2; k += 2)
      {
        outfile << (k * 0.5) + 1 << " " << triangle.pointlist[k] << " "
                << triangle.pointlist[k + 1] << std::endl;
      }
    }

    // Write points attribute list
    if (triangle.numberofpointattributes != 0)
    {
      outfile << "# Triangulateio number of points attributelist is:"
              << triangle.numberofpointattributes << std::endl;
    }
    if (triangle.pointattributelist != NULL)
    {
      outfile << "# Vertex attribute are:" << std::endl;
      for (int k = 0; k < triangle.numberofpointattributes; k++)
      {
        outfile << triangle.pointattributelist[k] << std::endl;
      }
    }

    // Write point markers list
    if (triangle.pointmarkerlist != NULL)
    {
      outfile << "# Vertex Markers are:" << std::endl;
      for (int k = 0; k < triangle.numberofpoints; k++)
      {
        outfile << triangle.pointmarkerlist[k] << std::endl;
      }
    }

    // Write the 1.node file used by the showme function
    std::ofstream nodefile;
    char nodename[100];

    sprintf(nodename, "file_%s.1.node", s.c_str());
    nodefile.open(nodename);
    nodefile << triangle.numberofpoints << " 2 "
             << triangle.numberofpointattributes << " 0" << std::endl;
    for (int j = 0; j < triangle.numberofpoints * 2; j += 2)
    {
      nodefile << (j / 2) + 1 << " " << triangle.pointlist[j] << " "
               << triangle.pointlist[j + 1] << std::endl;
    }
    nodefile.close();


    // Write segments edge elements
    if (triangle.numberofsegments != 0)
    {
      outfile << "# The number of segments is:" << triangle.numberofsegments
              << std::endl;
    }
    if (triangle.segmentlist != NULL)
    {
      outfile << "# Segments are:" << std::endl;
      for (int k = 0; k < triangle.numberofsegments * 2; k += 2)
      {
        outfile << triangle.segmentlist[k] << "  "
                << triangle.segmentlist[k + 1] << std::endl;
      }
    }

    // Write segments markers list
    if (triangle.segmentmarkerlist != NULL)
    {
      outfile << "# Segments Markers are:" << std::endl;
      for (int k = 0; k < triangle.numberofsegments; k++)
      {
        outfile << triangle.segmentmarkerlist[k] << std::endl;
      }
    }

    // Write regions
    if (triangle.numberofregions != 0)
    {
      outfile << "# The number of region is:" << triangle.numberofregions
              << std::endl;
    }

    // Write holes
    if (triangle.numberofholes != 0)
    {
      outfile << "# The number of holes is:" << triangle.numberofholes
              << std::endl;
    }
    if (triangle.holelist != NULL)
    {
      outfile << "#  Holes are:" << std::endl;
      for (int k = 0; k < triangle.numberofholes * 2; k += 2)
      {
        outfile << triangle.holelist[k] << "  " << triangle.holelist[k + 1]
                << std::endl;
      }
    }

    // Write triangles
    if (triangle.numberoftriangles != 0)
    {
      outfile << "# Triangulateio number of triangles:"
              << triangle.numberoftriangles << std::endl;
    }
    if (triangle.numberofcorners != 0)
    {
      outfile << "# Triangulateio number of corners:"
              << triangle.numberofcorners << std::endl;
    }
    if (triangle.numberoftriangleattributes != 0)
    {
      outfile << "# Triangulateio number of triangles attributes:"
              << triangle.numberoftriangleattributes << std::endl;
    }
    if (triangle.trianglelist != NULL)
    {
      outfile << "# Traingles are:" << std::endl;
      for (int k = 0; k < triangle.numberoftriangles * 3; k += 3)
      {
        outfile << triangle.trianglelist[k] << " "
                << triangle.trianglelist[k + 1] << " "
                << triangle.trianglelist[k + 2] << std::endl;
      }
    }

    if (triangle.trianglearealist != NULL)
    {
      outfile << "# Triangle's areas are:" << std::endl;
      for (int k = 0; k < triangle.numberoftriangles; k++)
      {
        outfile << triangle.trianglearealist[k] << std::endl;
      }
    }

    if (triangle.trianglelist != NULL)
    {
      // Write the 1.ele file used by the showme function
      std::ofstream elefile;
      char elename[100];

      sprintf(elename, "file_%s.1.ele", s.c_str());
      elefile.open(elename);
      elefile << triangle.numberoftriangles << " 3 0" << std::endl;
      for (int j = 0; j < triangle.numberoftriangles * 3; j += 3)
      {
        elefile << (j / 3) + 1 << " " << triangle.trianglelist[j] << " "
                << triangle.trianglelist[j + 1] << " "
                << triangle.trianglelist[j + 2] << std::endl;
      }
      elefile.close();
    }

    outfile.close();
  }

#endif

  //================================================================
 /// Setup lookup schemes which establish which elements are located
 /// next to which boundaries (Doc to outfile if it's open).
 //================================================================
 void TriangleMeshBase::setup_boundary_element_info(std::ostream& outfile)
 {
  
  // Should we document the output here
  bool doc = false;
  
  if (outfile) doc = true;
  
  // Number of boundaries
  unsigned nbound = nboundary();
  
  // Wipe/allocate storage for arrays
  Boundary_element_pt.clear();
  Face_index_at_boundary.clear();
  Boundary_element_pt.resize(nbound);
  Face_index_at_boundary.resize(nbound);
  
  
  // Loop over elements
  //-------------------
  unsigned nel = nelement();
  
  for (unsigned e = 0; e < nel; e++)
   {
    // Get pointer to element
    FiniteElement* fe_pt = finite_element_pt(e);
    
    if (doc)
     {
      outfile << "Element: " << e << " " << fe_pt << std::endl;
     }
    
    // Only include 2D elements! Some meshes contain interface elements too.
    if (fe_pt->dim() == 2)
     {
      // Loop over the element's nodes and find out which boundaries they're
      // on
      // ----------------------------------------------------------------------
      Vector<std::set<unsigned>*> boundaries_pt(3, 0);
      
      // We need only loop over the corner nodes
      for (unsigned i = 0; i < 3; i++)
       {
        fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
        if (doc)
         {
          outfile << "Node " << i << " in element: " << e
                  << " at "
                  << fe_pt->node_pt(i)->x(0) << " "
                  << fe_pt->node_pt(i)->x(1) << " "
                  << " is on boundaries: ";
          if (boundaries_pt[i]!=0)
           {
            for (unsigned b : *boundaries_pt[i])
             {
              outfile << b << " ";
             }
           }
          outfile << std::endl;
         }
       }
      
      
      // Edge 0 connects points 1 and 2
      //-----------------------------
      
      if (boundaries_pt[1] && boundaries_pt[2])
       {
        // Find the common boundaries of each edge
        std::set<unsigned> edge_boundary;
        std::set_intersection(boundaries_pt[1]->begin(),
                              boundaries_pt[1]->end(),
                              boundaries_pt[2]->begin(),
                              boundaries_pt[2]->end(),
                              std::insert_iterator<std::set<unsigned>>(
                               edge_boundary, edge_boundary.begin()));
        std::set<unsigned>::iterator it = edge_boundary.begin();
        
        // Edge does exist:
        if (edge_boundary.size() > 0)
         {
#ifdef PARANOID
          if (edge_boundary.size() > 1)
           {
            std::ostringstream error_stream;
            error_stream
             << "This shouldn't happen!";
            throw OomphLibError(
             error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
           }
#endif
          // The next element on this boundary (pointed to by it0) has face index 0 
          Boundary_element_pt[*it].push_back(fe_pt);
          Face_index_at_boundary[*it].push_back(0);
         }
       }
      
      // Edge 1 connects points 0 and 2
      //-----------------------------
      
      if (boundaries_pt[0] && boundaries_pt[2])
       {
        std::set<unsigned> edge_boundary;
        std::set_intersection(boundaries_pt[0]->begin(),
                              boundaries_pt[0]->end(),
                              boundaries_pt[2]->begin(),
                              boundaries_pt[2]->end(),
                              std::insert_iterator<std::set<unsigned>>(
                               edge_boundary, edge_boundary.begin()));
        std::set<unsigned>::iterator it = edge_boundary.begin();
        
        // Edge does exist:
        if (edge_boundary.size() > 0)
         {
#ifdef PARANOID
          if (edge_boundary.size() > 1)
           {
            std::ostringstream error_stream;
            error_stream
             << "This shouldn't happen!";
            throw OomphLibError(
             error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
           }
#endif
          // The next element on this boundary (pointed to by it) has face index 1 
          Boundary_element_pt[*it].push_back(fe_pt);
          Face_index_at_boundary[*it].push_back(1);           
         }
       }
      
      // Edge 2 connects points 0 and 1
      //-----------------------------
      
      if (boundaries_pt[0] && boundaries_pt[1])
       {
        std::set<unsigned> edge_boundary;
        std::set_intersection(boundaries_pt[0]->begin(),
                              boundaries_pt[0]->end(),
                              boundaries_pt[1]->begin(),
                              boundaries_pt[1]->end(),
                              std::insert_iterator<std::set<unsigned>>(
                               edge_boundary, edge_boundary.begin()));
        
        std::set<unsigned>::iterator it = edge_boundary.begin();
        
        // Edge does exist:
        if (edge_boundary.size() > 0)
         {
#ifdef PARANOID
          if (edge_boundary.size() > 1)
           {
            std::ostringstream error_stream;
            error_stream
             << "This shouldn't happen!";
            throw OomphLibError(
             error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
           }
#endif
          // The next element on this boundary (pointed to by it) has face index 2 
          Boundary_element_pt[*it].push_back(fe_pt);
          Face_index_at_boundary[*it].push_back(2);
          
          
         }
        
       }
      
     }
    
   }
  
  
  // Doc?
  //-----
  if (doc)
   {
    // Loop over boundaries
    for (unsigned i = 0; i < nbound; i++)
     {
      unsigned nel = Boundary_element_pt[i].size();
      outfile << "Boundary: " << i << " is adjacent to " << nel << " elements"
              << std::endl;
      
      // Loop over elements on given boundary
      for (unsigned e = 0; e < nel; e++)
       {
        FiniteElement* fe_pt = Boundary_element_pt[i][e];
        outfile << "Boundary element:" << fe_pt
                << " Face index of boundary is "
                << Face_index_at_boundary[i][e] << std::endl;
       }
     }
   }
  
  // Lookup scheme has now been setup yet
  Lookup_for_elements_next_boundary_is_setup = true;
    
 }
 
 
} // namespace oomph
