#! /bin/ksh
#

always_oas(){
rout "${cblue}>> always_oas${cnormal}"
  liboas=""
  libpsmile="$oasdir/$platform/lib/libpsmile.MPI1.a $oasdir/$platform/lib/libmct.a $oasdir/$platform/lib/libmpeu.a $oasdir/$platform/lib/libscrip.a"
  incpsmile="-I$oasdir/$platform/build/lib/psmile.MPI1"
rout "${cblue}<< always_oas${cnormal}"
}

substitutions_oas(){
rout "${cblue}>> substitutions_oas${cnormal}"
    c_substitutions_oas
  comment "   cp new  mod_oasis_method.F90 to psmile/src"
    cp $rootdir/bldsva/intf_oas3/oasis3-mct/arch/JURECA/src/mod_oasis_method.F90 ${oasdir}/lib/psmile/src >> $log_file 2>> $err_file
  check
  comment "   cp new  mod_oasis_grid.F90 to psmile/src"
    cp $rootdir/bldsva/intf_oas3/oasis3-mct/arch/JURECA/src/mod_oasis_grid.F90 ${oasdir}/lib/psmile/src >> $log_file 2>> $err_file
  check
  comment "   sed prism_get_freq functionality to mod_prism.F90"
    sed -i "/oasis_get_debug/a   use mod_oasis_method ,only: prism_get_freq            => oasis_get_freq" ${oasdir}/lib/psmile/src/mod_prism.F90  >> $log_file 2>> $err_file   # critical anchor
  check
    # set prefix vor all mct files to run with other mct
    prefix="oas_"
  comment "   cd to ${oasdir}/lib/mct/mct"
    cd ${oasdir}/lib/mct/mct >> $log_file 2>> $err_file
  check
  comment "   mv all m_* to ${prefix}m_*"
    for i in m_* ; do mv $i ${prefix}${i}  >> $log_file 2>> $err_file ; check ; done
  comment "   mv mct_mod.F90 to ${prefix}mct_mod.F90"
    mv mct_mod.F90 ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
  comment "   rename all interfaces in mct_mod.F90 and Makefile"
    sed -i "s/\(${prefix}\)*mct_/${prefix}mct_/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/\(${prefix}\)*mct_/${prefix}mct_/g" Makefile >> $log_file 2>> $err_file
  check
  comment "   rename all sources in mct_mod.F90 and Makefile"
    sed -i "s/use m_/use ${prefix}m_/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/\(${prefix}\)*m_/${prefix}m_/g" Makefile >> $log_file 2>> $err_file
  check
  comment "   except those which are redirecting mpeu sources"
    sed -i "s/${prefix}m_List/m_List/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_string/m_string/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_die/m_die/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_MergeSort/m_MergeSort/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_inpak90/m_inpak90/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
    sed -i "s/${prefix}m_Permuter/m_Permuter/g" ${prefix}mct_mod.F90 >> $log_file 2>> $err_file
  check
  comment "   rename all module names"
    sed -i "s/module m_/module ${prefix}m_/g" * >> $log_file 2>> $err_file
  check
  comment "   rename all dependencies"
    sed -i "s/use m_MCTW/use ${prefix}m_MCTW/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Glob/use ${prefix}m_Glob/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_AttrV/use ${prefix}m_AttrV/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Accum/use ${prefix}m_Accum/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Gener/use ${prefix}m_Gener/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Spa/use ${prefix}m_Spa/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Nav/use ${prefix}m_Nav/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Con/use ${prefix}m_Con/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Ex/use ${prefix}m_Ex/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Rou/use ${prefix}m_Rou/g" * >> $log_file 2>> $err_file
  check
    sed -i "s/use m_Rear/use ${prefix}m_Rear/g" * >> $log_file 2>> $err_file
  check
  comment "   cd in psmile-source dir"
    cd ${oasdir}/lib/psmile/src >> $log_file 2>> $err_file
  check
  comment "   rename all mct references"
    sed -i "s/\(${prefix}\)*mct_/${prefix}mct_/g" * >> $log_file 2>> $err_file
  check
    
rout "${cblue}<< substitutions_oas${cnormal}"
}

configure_oas(){
rout "${cblue}>> configure_oas${cnormal}"
  file=${oasdir}/util/make_dir/make.oas3
  comment "   cp jureca oasis3-mct makefile to /util/make_dir/"
    cp $rootdir/bldsva/intf_oas3/oasis3-mct/arch/$platform/config/make.intel_jureca_oa3 $file >> $log_file 2>> $err_file
  check
  c_configure_oas
  comment "   sed new psmile includes to Makefile"
    sed -i 's@__inc__@-I$(LIBBUILD)/psmile.$(CHAN) -I$(LIBBUILD)/scrip  -I$(LIBBUILD)/mct'" -I$ncdfPath/include@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ldflg to oas Makefile"
    sed -i "s@__ldflg__@@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comF90 to oas Makefile"
    sed -i "s@__comF90__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  check
  comment "   sed comCC to oas Makefile"
    sed -i "s@__comCC__@$mpiPath/bin/mpicc@" $file >> $log_file 2>> $err_file
  check
  comment "   sed ld to oas Makefile"
    sed -i "s@__ld__@$mpiPath/bin/mpif90@" $file >> $log_file 2>> $err_file
  check
  comment "   sed libs to oas Makefile"
    sed -i "s@__lib__@-L$ncdfPath/lib/ -lnetcdff@" $file >> $log_file 2>> $err_file
  check
  comment "   sed precision to oas Makefile"
    sed -i "s@__precision__@-i4 -r8@" $file >> $log_file 2>> $err_file
  check

rout "${cblue}<< configure_oas${cnormal}"
}

make_oas(){
rout "${cblue}>> make_oas${cnormal}"
  c_make_oas
rout "${cblue}<< make_oas${cnormal}"
}


setup_oas(){
rout "${cblue}>> setupOas${cnormal}"
  comment "   copy cf_name_table to rundir"
    cp $rootdir/bldsva/data_oas3/cf_name_table.txt $rundir >> $log_file 2>> $err_file
  check
  comment "   copy oas Makefile to rundir"
    cp $namelist_oas $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed procs, gridsize & coupling freq into namcouple"
  if [[ $withPFL == "true" && $withCOS == "true" ]] then
    ncpl_exe1=$nproc_cos
    ncpl_exe2=$nproc_pfl
    ncpl_exe3=$nproc_clm
    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe3/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe3/$ncpl_exe3/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq2/$cplfreq2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    
    sed "s/ngpflx/$gx_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngpfly/$gy_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosx/$(($gx_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosy/$(($gy_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

  fi  

  if [[ $withPFL == "true" && $withCOS == "false" ]] then
    ncpl_exe1=$nproc_pfl
    ncpl_exe2=$nproc_clm

    sed "s/nproc_exe1/$nproc_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq2/$cplfreq2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngpflx/$gx_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngpfly/$gy_pfl/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  fi
  if [[ $withPFL == "false" && $withCOS == "true" ]] then
    ncpl_exe1=$nproc_cos
    ncpl_exe2=$nproc_clm
    sed "s/nproc_exe1/$nproc_cos/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe1/$ncpl_exe1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/nproc_exe2/$nproc_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ncpl_exe2/$ncpl_exe2/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/cplfreq1/$cplfreq1/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

    sed "s/ngcosx/$(($gx_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngcosy/$(($gy_cos-($nbndlines*2)))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmx/$gx_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
    sed "s/ngclmy/$gy_clm/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

  fi
  comment "   sed sim time into namcouple"
    sed "s/totalruntime/$(($runhours*3600))/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check
  comment "   sed startdate into namcouple"
    sed "s/yyyymmdd/${yyyy}${mm}${dd}/" -i $rundir/namcouple >> $log_file 2>> $err_file
  check

rout "${cblue}<< setupOas${cnormal}"
}

