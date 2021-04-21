#
# Script to run velocity analysis for camera images.
#

Cameras="LocA LocV01 LocV02 LocV03"
#Cameras="LocV01 LocV02 LocV03"
#Cameras="LocV01"
#Cameras="LocA"

for Camera in $Cameras
do

	if [ $Camera == 'LocA' ]
	then
		timestart=001
		COL_NUM1s="95" # 145 195 245 295 345"
		istep=5
	elif [[ $Camera == 'LocV01' || $Camera == 'LocV02' || $Camera == 'LocV03' ]]
	then
		timestart=001
		COL_NUM1s="5" # 15 25 35 45 55 65 75"
		if [ $Camera == 'LocV01' ]
		then
			istep=5
		elif [ $Camera == 'LocV02' ]
		then
			istep=5
		elif [ $Camera == 'LocV03' ]
		then
			istep=5
		fi
	fi

	for COL_NUM1 in $COL_NUM1s
	do
		COL_NUM2=$(($COL_NUM1+$istep))
		COL_NUM3=$(($COL_NUM2+$istep))
		echo $COL_NUM1, $COL_NUM2, $COL_NUM3, $Camera
		python pyplis_velocity_analysis.py --Camera=$Camera --COL_NUM1=$COL_NUM1 --COL_NUM2=$COL_NUM2 --COL_NUM3=$COL_NUM3 --timestart=$timestart
	done
done

## OLD STUFF, IGNORE FOR NOW

# #python pyplis_read_plot.py ./ LocA --experimenttype=SUNSZA40WQBLC --timestampBG=005 --InputFolderBG=./ --experimenttypeBG=SUNSZA40W_BG --PlotType=Velocity --timestart=001 --timeend=005 --timestep=1

# experimenttype=SUNSZA40WQBLC
# experimenttypeBG=SUNSZA40WQBLC_BG
# InputFolder=/home/aky/NILU/xnilu_wrk/users/aky/Projects/COMTESSA/Experiments/palm_tomo_Rena_HT_Velo/
# InputFolderBG=/home/aky/NILU/xnilu_wrk/users/aky/Projects/COMTESSA/Experiments/palm_tomo_Rena_HT_Velo/
# Locations="LocV01 LocV02 LocV03"
# Times="1 2 3 4 5 6 7 8 9 10" # 002"

# # for Location in $Locations
# # do
# # 	python pyplis_read_plot.py $InputFolder $Location --experimenttype=$experimenttype --timestampBG=005 --InputFolderBG=$InputFolderBG --experimenttypeBG=$experimenttypeBG  --PlotType=Velocity --timestart=001 --timeend=010 --timestep=1
# # done


# for Location in $Locations
# do
# 	for Time in $Times
# 	do
# 		python pyplis_read_plot.py $InputFolder $Location --experimenttype=$experimenttype --timestampBG=005 --InputFolderBG=$InputFolderBG --experimenttypeBG=$experimenttypeBG  --PlotType=Images --timestart=$Time  --timestep=1 --pngfolder=./tmpvelo/
# 	done
# done

# convert +append ./tmpvelo/SUNSZA40WQBLC_LocV01_???_appabs_V1.png  tmpvelo/LocV01_appabs.png
