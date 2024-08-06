#!/bin/sh

function mapp() {
#  This is from, http://prll.sourceforge.net/shell_parallel.html
  if [[ -z $MAPP_NR_CPUS ]] ; then
#	local MAPP_NR_CPUS=$(grep "processor	:" < /proc/cpuinfo | wc -l)
#   max core for calculation; modified by Toshiyuki Sumikama
      local MAPP_NR_CPUS=8 # was previously 8
  fi
  local mapp_pid=$(exec bash -c 'echo $PPID')
  local mapp_funname=$1
  local -a mapp_params
  mapp_params=("$@")
#   mapp_nr_args; modified by Toshiyuki Sumikama
#  local mapp_nr_args=${#mapp_params[@]}
  local mapp_nr_args=`expr ${#mapp_params[@]} - 1`
  local mapp_current=0
  function mapp_trap() {
#  echo "MAPP PROGRESS: $((mapp_current*100/mapp_nr_args))%" 1>&2
	if [[ $mapp_current -lt $mapp_nr_args ]] ; then
	    let mapp_current+=1
	    (
		$mapp_funname "${mapp_params[$mapp_current]}"
		kill -USR1 $mapp_pid
	    ) &
	fi
  }

  trap mapp_trap SIGUSR1
  while [[ $mapp_current -lt $mapp_nr_args ]]; do
	wait
	if [[ $mapp_current -lt $mapp_nr_args && $? -lt 127 ]] ; then
	    sleep 1
	    local mapp_tmp_count=$mapp_current
	    wait
	    if [[ $mapp_tmp_count -eq $mapp_current ]] ; then
		echo "   MAPP_FORCE" 1>&2
		for i in $(seq 1 ${MAPP_NR_CPUS}) ; do
		    mapp_trap
		done
	    fi
	fi
  done
  for i in $(seq 1 ${MAPP_NR_CPUS}) ; do
	wait
  done
  trap - SIGUSR1
  unset -f mapp_trap
}

myfun()
{
(
    cd $HOME/nucnet-tools-code/my_examples/network # where comannds are executed

    # FOR EDITING ZONE.XML
    # Define the file path and the new values





    # echo $val
    FILE="$HOME/nucnet-tools-code/data_pub/zone.xml"
    # my_array=(1.4985e0 1.4985e1 1.4985e2 1.4985e3 1.4985e4 1.4985e5 1.4985e6 xxx) # current value at very end
    
    # prev=${my_array[$1 - 1]}
    # next=${my_array[$1]}
    
    # turns out i need to make a new zone.xml file entirely, because otherwise
    # ./run_single_zone might use the wrong one due to parallelism
    # # LOCK FILE SO CAN RUN IN PARALLEL
    # # Acquire an exclusive lock on the XML file
    # exec 200>$FILE.lock
    # flock -x 200
    new_value=$(bc <<< "1498500 + $1*10")
    xmlstarlet ed -u "/zone_data/zone/optional_properties/property[@name='tau_1']" -v "$new_value" $FILE > $HOME/nucnet-tools-code/ZZZ_my_data/zone$1.xml
    # # Release the lock
    # flock -u 200
    

    # RUNNING ZONE
    ./run_single_zone $HOME/nucnet-tools-code/data_pub/net_with_fission.xml $HOME/nucnet-tools-code/ZZZ_my_data/zone$1.xml $HOME/nucnet-tools-code/ZZZ_my_data/my_output_with_fission$1.xml "[z <= 90]" >> /dev/null
    ./run_single_zone $HOME/nucnet-tools-code/data_pub/my_net.xml $HOME/nucnet-tools-code/ZZZ_my_data/zone$1.xml $HOME/nucnet-tools-code/ZZZ_my_data/my_output$1.xml "[z <= 90]" >> /dev/null


    # GETTING TXT FILES
    $HOME/nucnet-tools-code/examples/analysis/print_abundances_vs_nucleon_number $HOME/nucnet-tools-code/ZZZ_my_data/my_output_with_fission$1.xml a "[last()]" > $HOME/nucnet-tools-code/ZZZ_my_data/ya_with_fission$1.txt 
    $HOME/nucnet-tools-code/examples/analysis/print_abundances_vs_nucleon_number $HOME/nucnet-tools-code/ZZZ_my_data/my_output$1.xml a "[last()]" > $HOME/nucnet-tools-code/ZZZ_my_data/ya$1.txt
    
    # FINDING MIN ERROR
    # cd $HOME/nucnet-tools-code/ZZZ_my_data
    # result=$(python3 Iterative_R_Process.py "ya")
    # echo $result
    echo $1' was finished.'
)
}

for eachrun in {0..100} # 0000..0010
do
    list=$list' '$eachrun
done

echo $list
# val=$1 # takes first argument, which is input from command line

# (mapp myfun $val $list)
mapp myfun $list $val

