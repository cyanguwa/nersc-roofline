MODULE ITT_SDE_FORTRAN
USE, INTRINSIC :: ISO_C_BINDING

INTERFACE
   
   SUBROUTINE FORTRAN_ITT_RESUME() &
      BIND(C, NAME='fortran_itt_resume')
   END SUBROUTINE FORTRAN_ITT_RESUME

   SUBROUTINE FORTRAN_ITT_PAUSE() &
      BIND(C, NAME='fortran_itt_pause')
   END SUBROUTINE FORTRAN_ITT_PAUSE

   SUBROUTINE FORTRAN_SDE_START() &
      BIND(C, NAME='fortran_sde_start')
   END SUBROUTINE FORTRAN_SDE_START

   SUBROUTINE FORTRAN_SDE_STOP() &
      BIND(C, NAME='fortran_sde_stop')
   END SUBROUTINE FORTRAN_SDE_STOP
END INTERFACE

contains

   subroutine start_collection()
     call fortran_sde_start()
     call fortran_itt_resume()
   end subroutine start_collection

   subroutine stop_collection() 
    call fortran_itt_pause()
    call fortran_sde_stop()
   end subroutine stop_collection

END MODULE
