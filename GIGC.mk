# Make sure ESMADIR is defined
# ----------------------------
ifndef BASEDIR
export BASEDIR=$(PWD)/GIGC
endif 
ifndef ESMADIR
export ESMADIR=$(PWD)/GIGC/Shared
endif
ifndef ESMF_DIR
export ESMF_DIR=$(PWD)/GIGC/ESMF
endif
       RCDIR   :=$(BASEDIR)/Registry/

# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies

# ESMF-specific settings
# ----------------------------
export ESMF_COMPILER=intel
export ESMF_COMM=openmpi
export ESMF_INSTALL_PREFIX=$(ESMF_DIR)/$(ARCH)
export ESMF_INSTALL_LIBDIR=$(ESMF_DIR)/$(ARCH)/lib
export ESMF_INSTALL_MODDIR=$(ESMF_DIR)/$(ARCH)/mod
export ESMF_INSTALL_HEADERDIR=$(ESMF_DIR)/$(ARCH)/include
export ESMF_F90COMPILEOPTS=-align all -fPIC
export ESMF_CXXCOMPILEOPTS=-fPIC
export ESMF_OPENMP=OFF
