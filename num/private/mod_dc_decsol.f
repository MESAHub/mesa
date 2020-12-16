      module mod_dc_decsol
      use const_def, only: dp
      use utils_lib, only: mesa_error
      
      contains

#include "decomc.dek"
#include "decomr.dek"
#include "estrad.dek"
#include "decsol_done.dek"
#include "slvrad.dek"
#include "slvrod.dek"
      
      end module mod_dc_decsol
      
      
