
1. copy this directory, /work/00315/tg455890/stampede2/EPIC/Program
   Please make sure all executables in Program/exec are in place.

2. Link GDAS file
   2.1 Modify WRKDIR in scripts/retrieve_gdas.sh
   2.1 Execute
       $> retrieve_gdas.sh 2019043018
       you will get $WRKDIR/2019043018_memxxx

3. Prepare for ICS
   3.1 Modify script scripts/run_ics.sh
       "wrkroot",
       "gdasdir" (to match WRKDIR above in retrieve_gdas.sh ),
       "rootdir" directory for "Program" copied above
   3.2 Execute
       $> scripts/run_ics.sh YYYYMMDDHH NN [test_runs|test_mp|test_spp|test_mspp]
       where NN is the number of ensembles, 40
       "test_runs": single physics
       "test_mp":   mulitple physics
       "test_spp":  single physics spp run
       "test_mspp": multiple physics spp run
   3.3 You can modify the run-time job parameters in scripts/driver_chgres_regional.c768

4. Prepare for LBCs
   4.1 Modify script, scripts/run_lbc.sh
       "wrkroot"
       "rootdir"
       "gdasdir"
   4.2 Execute
       $> scripts/run_lbc.sh
   4.3 Modify run-time parameters in templates/chgres_driver.job

5. Run FV3
   5.1 Modify scripts/run_sar.sh
       "wrkroot"
	   "rootdir"

	   layout_x="14"
	   layout_y="10"
	   nquilt="28"
	   io_layout="1,1"

	   nodes=$((npes/24))
	   m=$((npes%24))

	   etc for optimal performance on Stampede
   5.2 Execute
       $> scripts/run_sar.sh YYYYMMDDHH NN test_runs

   5.3 run-time parameters in templates/run_fv3_on_Stampede.job
