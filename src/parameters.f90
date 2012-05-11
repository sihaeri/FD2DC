MODULE parameters

!--Number of variables to solve 
INTEGER,PARAMETER       :: nphi = 4
INTEGER,PARAMETER       :: max_char_len = 20

!--File units
INTEGER,PARAMETER       :: set_unit = 1 !--Problem setup file
INTEGER,PARAMETER       :: out_unit = set_unit + 1 !--OUtput information at each iteration
INTEGER,PARAMETER       :: grd_unit = out_unit + 1 !--grid file
INTEGER,PARAMETER       :: sres_unit = grd_unit + 1 !--read restart file
INTEGER,PARAMETER       :: eres_unit = sres_unit + 1 !--write restart file
INTEGER,PARAMETER       :: plt_unit = eres_unit + 1 !-- tect plot results
INTEGER,PARAMETER       :: ubar_unit = plt_unit + 1 
INTEGER,PARAMETER       :: vbar_unit = ubar_unit + 1 
INTEGER,PARAMETER       :: u2bar_unit = vbar_unit + 1 
INTEGER,PARAMETER       :: v2bar_unit = u2bar_unit + 1 
INTEGER,PARAMETER       :: uvbar_unit = v2bar_unit + 1 
INTEGER,PARAMETER       :: alloc_create = 1
INTEGER,PARAMETER       :: alloc_destroy = alloc_create+1
INTEGER,PARAMETER       :: max_len_tecline  = 32000

INTEGER,PARAMETER       :: force_predict = 1
INTEGER,PARAMETER       :: force_correct = force_predict + 1

INTEGER,PARAMETER       :: solver_sparsekit = 1
INTEGER,PARAMETER       :: solver_sip       = solver_sparsekit + 1
INTEGER,PARAMETER       :: solver_cg        = solver_sip + 1
INTEGER,PARAMETER       :: solver_hypre     = solver_cg + 1

INTEGER,PARAMETER       :: maxmove = 2

INTEGER,PARAMETER       :: do_collect_stat = 1
INTEGER,PARAMETER       :: end_collect_stat = do_collect_stat + 1

END MODULE parameters