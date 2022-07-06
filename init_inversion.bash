#!/bin/bash
# set up running directories
TsuData="TsunamiData"   # S-net or TsunamiData
curdir=`pwd`
if [ -d forward ]; then 
 rm -r forward
fi
if [ -d inversion ]; then
 rm -r inversion
fi
if [ -d reversal ]; then
 rm -r reversal
fi
if [ -d data_res ]; then
 rm -r data_res
fi

mkdir forward
mkdir inversion
mkdir reversal
mkdir data_res

# initialize running directories

# forward
ln -s $curdir/pcomcot/pcomcot ./forward
ln -s $curdir/bathy_pre/layer01.xyz ./forward
cp $curdir/data_pre/Stations.ctl ./forward

# reversal
ln -s $curdir/pcomcot_timereversal/pcomcot ./reversal
ln -s $curdir/bathy_pre/layer01.xyz ./reversal
cp $curdir/data_pre/Stations.ctl ./reversal
cp $curdir/data_pre/Station*.dat ./reversal
cp $curdir/pcomcot_timereversal/FaultParameters.ctl ./reversal
ln -s $curdir/aux/conv_init_waterelev ./reversal

# inversion
ln -s $curdir/aux/conv_init_waterelev ./inversion

# data_res
ln -s $curdir/aux/convert_data_res ./data_res

# pcomcot_fault
# cp $curdir/data_pre/Stations.ctl ./pcomcot_fault
#rm ./pcomcot_fault/layer01.xyz
#cp $curdir/bathy_pre/layer01.xyz ./pcomcot_fault
#cp $curdir/data_pre/$TsuData/data_obs.mat ./pcomcot_fault
#cp $curdir/data_pre/$TsuData/data_obs_ori.mat ./pcomcot_fault
