bsub -o /project/imco/baller/scripts/logfiles/outputlogjob.out -e /project/imco/baller/scripts/logfiles/outputlogjob.error -R "rusage[mem=128G]" < /project/imco/baller/scripts/bcp_wrapper.sh
