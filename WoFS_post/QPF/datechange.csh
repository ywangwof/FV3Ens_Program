#!/bin/csh -f
set echo

#set base_dir = '/work/brian.matilla/python_realtime/'
set base_dir = '/home/brian.matilla/work_backup/python_realtime/'

#set pat_dir = '/work/brian.matilla/python_realtime_PAT/'
set pat_dir = '/home/brian.matilla/work_backup/python_realtime_PAT/'


set date = `date -u --date='0 minutes ago' +%Y%m%d` 
set prev_date = `date -u --date='1 day ago' +%Y%m%d`

#setenv yyyy `date -u +%Y`
#setenv mm `date -u +%m`
#setenv dd `date -u +%d `
#
#setenv ddy `date -d "$date 1 day ago" +%d`
#setenv ddt `date -d "$date 1 day" +%d`
#
##set date = $yyyy$mm$dd
#set prev_date = $yyyy$mm$ddy
#set next_date = $yyyy$mm$ddt

cd $base_dir

#exec sed -i -e "18s/$prev_date/$date/g" run_mrms_remap.csh -i -e "19s/$date/$next_date/g" run_mrms_remap.csh -i -e "s/$prev_date/$date/g" run_mrms_qpe_plot.csh 
exec sed -i -e "s/$prev_date/$date/g" run_newse_post_swath_rainfall.job -i -e "s/$prev_date/$date/g" run_newse_plot_swath_rainfall.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_post_ensemble.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_post_swath.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_post_environment.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_post_sounding.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_swath.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_sounding2.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_sounding3.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_sounding4.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_sounding5.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_sounding.job -i -e "s/$prev_date/$date/g" ${pat_dir}run_newse_plot_timestep.job -i -e "s/$prev_date/$date/g" run_newse_basemap.job -i -e "s/$prev_date/$date/g" /home/brian.matilla/work_backup/GLM/scripts_rt/create_GOES16_GLM_OBS_rlt.csh 

#exec sed -i -e "s/$prev_date/$date/g" run_mrms_remap_v2.csh -i -e "s/$prev_date/$date/g" run_mrms_qpe_plot_agg.csh -i -e "s/$prev_date/$date/g" run_mrms_agg_realtime.csh -i -e "s/$prev_date/$date/g" run_paint.csh -i -e "s/$prev_date/$date/g" run_verif_obj.csh "s/$prev_date/$date/g" run_mrms_createnans.csh
