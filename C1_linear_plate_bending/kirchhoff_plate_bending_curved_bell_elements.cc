// Non--inline functions for BellBiharmonic elements
#include "kirchhoff_plate_bending.h"

namespace oomph
{

//======================================================================
/// Set the data for the number of Variables at each node
//======================================================================
 template<unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
 const unsigned KirchhoffPlateBendingC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::Initial_Nvalue = 6;


//====================================================================
// Force build of templates
//====================================================================
template class KirchhoffPlateBendingC1CurvedBellElement<2,2,3>;
template class KirchhoffPlateBendingC1CurvedBellElement<2,2,5>;

}
