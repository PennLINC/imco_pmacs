bsub -o /project/imco/baller/PBP_graphics/scripts/logfiles/outputlogjob_pbp_viz.out -e /project/imco/baller/scripts/PBP_graphics/logfiles/outputlogjob_pbp_viz.error -R "rusage[mem=128G]" < /project/imco/baller/scripts/PBP_graphics/run_visualization_scripts_from_command_line.sh