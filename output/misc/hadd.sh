#!/bin/bash

for i in $(seq 0 15);
do
	if [ $i -lt 10 ]
	then
		idx=$i;

	elif [ $i -eq 10 ]
	then
		idx = a;

	elif [ $i -eq 11 ]
	then
		idx=b;

	elif [ $i -eq 12 ]
	then
		idx=c;

	elif [ $i -eq 13 ]
	then
		idx=d;

	elif [ $i -eq 14 ]
	then
		idx=e;

	elif [ $i -eq 15 ]
	then
		idx=f;

	fi
	
	for j in $(seq 10 15);
	do

		if [ $j -lt 10 ]
		then
			jdx=$j;
			hadd /data/ooconnor/background/reco2-$idx$j.root /data/sbnd/background/reco2-$idx$j*;
		fi
		

		if [ $j -eq 10 ]
		then
			jdx=a
			hadd /data/ooconnor/background/reco2-$idx$jdx.root /data/sbnd/background/reco2-$idx$jdx*;
		fi

		if [ $j -eq 11 ]
		then
			jdx=b
			hadd /data/ooconnor/background/reco2-$idx$jdx.root /data/sbnd/background/reco2-$idx$jdx*;
		fi

		if [ $j -eq 12 ]
		then
			jdx=c
			hadd /data/ooconnor/background/reco2-$idx$jdx.root /data/sbnd/background/reco2-$idx$jdx*;	
		fi

		if [ $j -eq 13 ]
		then
			jdx=d
			hadd /data/ooconnor/background/reco2-$idx$jdx.root /data/sbnd/background/reco2-$idx$jdx*;
		fi

		if [ $j -eq 14 ]
		then
			jdx=e
			hadd /data/ooconnor/background/reco2-$idx$jdx.root /data/sbnd/background/reco2-$idx$jdx*;
		fi

		if [ $j -eq 15 ]
		then
			jdx=f
			hadd /data/ooconnor/background/reco2-$idx$jdx.root /data/sbnd/background/reco2-$idx$jdx*;
		fi
	
	done
done

