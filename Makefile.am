# A few file definitions
#-----------------------

# Define the sources
sourcesC1 =  \
Bell_element_basis.cc C1_curved_elements.cc C1_curved_basic_basis.cc \
C1_curved_association_matrix.cc

sourcesFvK =  \
foeppl_von_karman_elements.cc foeppl_von_karman_bell_elements.cc \
 foeppl_von_karman_curved_bell_elements.cc

# Define the headers
headersC1 =  \
Bell_element_basis.h C1_curved_elements.h my_geom_object.h

headersFvK =  \
foeppl_von_karman.h foeppl_von_karman_bell_elements.h \
foeppl_von_karman_elements.h foeppl_von_karman_curved_bell_elements.h

# Define name of library
libnameC1 = C1_basis
libnameFvK = C1_foeppl_von_karman

# Combine headers and sources
headersC1_and_sourcesC1 = $(headersC1) $(sourcesC1)
headersFvK_and_sourcesFvK = $(headersFvK) $(sourcesFvK)

# Define the library that gets placed in lib directory
#-----------------------------------------------------
lib_LTLIBRARIES = libC1_basis.la libC1_foeppl_von_karman.la

# Sources that the library depends on: 
#-------------------------------------
libC1_basis_la_SOURCES = $(headersC1_and_sourcesC1)
libC1_foeppl_von_karman_la_SOURCES = $(headersFvK_and_sourcesFvK)

# The library's include headers:
#-------------------------------
# Headers that are to be included in the $(includedir) directory:
# This is the combined header which contains "#include<...>" commands
# for the real headers in the subdirectory below $(includedir)
C1_basisinclude_HEADERS = $(libnameC1).h 
C1_foeppl_von_karmaninclude_HEADERS = $(libnameFvK).h 
C1_basisincludedir = $(includedir) 
C1_foeppl_von_karmanincludedir = $(includedir)

#Here's the subdirectory where the actual header files are placed
library_C1_basisincludedir=$(includedir)/C1_basis
library_C1_foeppl_von_karmanincludedir=$(includedir)/C1_foeppl_von_karman

#These are the header files that are to be placed in subdirectory
library_C1_basisinclude_HEADERS=$(headersC1) 
library_C1_foeppl_von_karmaninclude_HEADERS=$(headersFvK) 


# Required libraries -- [assuming that we want to link against stuff in generic
#-------------------     add other oomph-lib libraries if you need them....]
# Include path for library headers -- need to refer to include files
# in their respective source directories as they will not have been 
# installed yet! 
AM_CPPFLAGS += -I$(top_builddir)/src/generic \
  -I$(top_builddir)/build/include/C1_basis/


# Combined header file
#---------------------
# Rule for building combined header (dummy target never gets made so
# the combined header is remade every time)
$(libnameC1).h: dummy_$(libnameC1).h
$(libnameFvK).h: dummy_$(libnameFvK).h

dummy_$(libnameC1).h:  $(headersC1)
	echo  $(libnameC1) $(headersC1) > all_$(libnameC1).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libnameC1).aux > $(libnameC1).h
	rm all_$(libnameC1).aux

dummy_$(libnameFvK).h:  $(headersFvK)
	echo  $(libnameFvK) $(headersFvK) > all_$(libnameFvK).aux
	$(AWK) -f $(top_builddir)/bin/headers.awk < \
	       all_$(libnameFvK).aux > $(libnameFvK).h
	rm all_$(libnameFvK).aux

# Extra hook for install: Optionally replace headers by symbolic links
#---------------------------------------------------------------------
if SYMBOLIC_LINKS_FOR_HEADERS
install-data-hook: 
	(cd $(library_includedir) && rm -f $(headers) )
	(echo "$(headers)" > include_files.list )
	($(top_builddir)/bin/change_headers_to_links.sh `pwd`)
	($(LN_S) `cat include_files.list.aux` $(library_includedir) )
	(rm -r include_files.list.aux include_files.list ) 
else
install-data-hook: 
endif


# Tidy up
#--------
clean-local:
	rm -f $(libname).h


