#include "ittnotify.h"

void fortran_sde_start()
{
  __SSC_MARK(0x111);
}  

void fortran_sde_stop()
{
  __SSC_MARK(0x222);
}  

void fortran_itt_resume()
{
  __itt_resume();
}

void fortran_itt_pause()
{
  __itt_pause();
}
