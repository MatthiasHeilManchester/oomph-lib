# A few file definitions
#-----------------------

# Define the sources
sourcesLPB = kirchhoff_plate_bending_elements.cc kirchhoff_plate_bending_curved_bell_elements.cc

# Define the headers
headersLPB = kirchhoff_plate_bending_elements.h kirchhoff_plate_bending_curved_bell_elements.h

# Define name of library
libnameLPB = C1_linear_plate_bending

# Combine headers and sources
headersLPB_and_sourcesLPB = $(headersLPB) $(sourcesLPB)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libC1_linear_plate_bending.la

# Sources that the library depends on: 
#-------------------------------------
libC1_linear_plate_bending_la_SOURCES = $(headersLPB_and_sourcesLPB)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
C1_linear_plate_bendinginclude_HEADERS = $(libnameLPB).h 
# C1_basisincludedir = $(includedir) 
C1_linear_plate_bendingincludedir = $(includedir)

#Here's the subdirectory where the actual header files are placed
library_C1_linear_plate_bendingincludedir=$(includedir)/C1_linear_plate_bending

#These are the header files that are to be placed in subdirectory
library_C1_linear_plate_bendinginclude_HEADERS=$(headersLPB) 


# Required libraries -- [assuming that we want to link against stuff in generic
#-------------------     add other oomph-lib libraries if you need them....]
# Include path for library headers -- need to refer to include files
# in their respective source directories as they will not have been 
# installed yet! 
AM_CPPFLAGS += -I$(top_builddir)/src/generic # -g -DPARANOID -DRANGE_CHECK


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libnameLPB).h: dummy_$(libnameLPB).h

dummy_$(libnameLPB).h:  $(headersLPB)
	echo  $(libnameLPB) $(headersLPB) > all_$(libnameLPB).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libnameLPB).aux > $(libnameLPB).h
	rm all_$(libnameLPB).aux

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_C1_linear_plate_bendingincludedir) && rm -f $(headersLPB) )
	(echo "$(headersLPB)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_C1_linear_plate_bendingincludedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Tidy up
#--------
clean-local:
	rm -f $(libnameLPB).h


