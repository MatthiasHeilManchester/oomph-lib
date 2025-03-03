# Initialisation script for self-tests wipes the master log validation.log
TESTS = validate.sh


# My own target for distributions that allow long filenames
# assuming there is some general purpose tar without the
# idiotic restrictions to the length of the filenames
# available on the machine
mydist: distdir
	tar cvfz $(distdir).tar.gz $(distdir)
	$(am__remove_distdir)

# Do we want to build the demo_drivers?
#--------------------------------------
if WANT_DEMO
   demo_dirs = demo_drivers
else
   demo_dirs =
endif



# Do we want to build the user sources?
#-------------------------------------
if WANT_USER_SRC
   user_src_dirs = user_src
else
   user_src_dirs =
endif



# Directories which should be built
#------------------------------------------
# NOTE: Order matters! 
#       (0) bin
#       (1) external_src, src and user_src first because libraries 
#           need to built before demo and/or user codes can be built.
#           Src order: External libraries cannot depend on 
#           internal ones and internal ones cannot depend on
#           user ones.
#       (2) . because test script in top level directory (=here) needs
#           to be run before the tests in the sub-directories
#       (3) demo driver codes.
sub_dirs = bin external_src external_distributions src $(user_src_dirs) . \
                   $(demo_dirs) 



# Doc directory (only make doc if doxygen is available and if it's wanted)
#-------------------------------------------------------------------------
# [NOTE: automake only allows single conditionals]
#-------------------------------------------------
if HAVE_DOXYGEN 
if WANT_DOC
   doc_dirs = doc
else
   doc_dirs =
endif
else
  doc_dirs =
endif

# The ultimate list of sub-directories (does not include user drivers!)
#----------------------------------------------------------------------
SUBDIRS = $(sub_dirs) $(doc_dirs) self_test

# The ultimate list of sub-directories that are to be distributed
#----------------------------------------------------------------
# includes the private stuff and the user drivers
#------------------------------------------------
# Note that the PRIVATE_DIR is substituted from the configure.ac file
# ---------------------------------------------------------------------
# because otherwise automake breaks
# ---------------------------------------------------------------------
DIST_SUBDIRS = user_drivers $(sub_dirs) $(doc_dirs) self_test $(PRIVATE_DIR)

# A few extract scripts to be included into the distribution
#-----------------------------------------------------------
EXTRA_DIST = \
validate.sh \
autogen.sh \
non_interactive_autogen.sh \
config \
depcomp config.guess config.sub install-sh .gitignore


# Customise the distribution
dist-hook:
# Kill stray symbolic links to private configure options
	(cd config/configure_options && find . -type l -exec rm {} \; )
	($(top_builddir)/bin/customise_dist.sh  `(cd $(top_builddir); pwd .)`  `(cd $(distdir); pwd .)` $(distdir) )
	(cd  $(distdir)/config/configure_options && cp -f default current)


clean-local: 
	rm -rf build validation.log
# We've taken the user drivers out of the SUBDIRS variable
# so they need to be cleaned by hand for normal clean.
# However, only do a clean if a Makefile exists!
	(cd user_drivers && if test -e Makefile; then ${MAKE} clean; fi )
# Kill stray symbolic links to private configure options
	(cd config/configure_options && find . -type l -exec rm {} \; )
# The private files are also not in the SUBDIRS variable so they also need
# to be cleaned by hand. We put a `-` in front of the command to say that
# it is allowed to fail.
	-(cd private && if test -e Makefile; then ${MAKE} clean; fi)

# Not necessary to force a descend in to user_drivers as they're
# cleaned automatically during distclean (because they're mentioned
# in DIST_SUBDIRS
distclean-local:
#We need to distclean the private directories because they are not
#mentioned in DIST_SUBDIRS
#	-(cd private && if test -e Makefile; then ${MAKE} distclean; fi)
	rm -f old_symbolic_links.txt
	rm -f config/configure.ac_scripts/private*.dir_list


# Sorry -- this is a total mess! I'm surrounding the actual output (the length
# in characters of the mpi run command with whitespace stripped out) by
# START_BLA and END_BLA because when this is called from the analyse
# script, make is actually called recursively and therefore additional
# output appears (which I can now strip out). One beer is on offer
# for whoever wants to do this properly. For now it works... 
do-we-have-mpi:
	@echo "$(MPI_RUN_COMMAND)" | tr -d '[[:space:]]' | wc -m > .tmp.wc; echo "START_BLA"`cat .tmp.wc`"END_BLA"; rm .tmp.wc



# Check if we have libtool installed (needed for updating
# distribution when new directories are added etc.)
check_if_we_have_have_libtool_installed_and_say_that_we_cant_reconfigure_if_not:
	@if $(have_libtool); then \
            echo "Have libtool so can reconfigure oomph-lib"; \
            echo "e.g. after adding new directories."; \
            exit 0; \
        else \
            echo " "; \
            echo " "; \
            echo "==================================================="; \
            echo " "; \
            echo " "; \
            echo " "; \
            echo "Don't have libtool so we can't reconfigure oomph-lib"; \
            echo "e.g. after adding new directories. Please install it"; \
            echo "and re-run autogen.sh!"; \
            echo " "; \
            echo " "; \
            echo " "; \
            echo "==================================================="; \
            echo " "; \
            echo " "; \
            exit 1; \
        fi;
