#!/bin/bash

ver=`cat ../../../../Program.cs|grep "string version ="|sed 's/.*= "//; s/".*//'`;
mkdir -p ~/temp-SELDLA
mv win mac linux ~/temp-SELDLA
cd ~/temp-SELDLA
chmod 644 linux/* mac/*
chmod 755 linux/SELDLA mac/SELDLA
mv win SELDLA-${ver}-win-x64
mv mac SELDLA-${ver}-mac-x64
mv linux SELDLA-${ver}-linux-x64
zip -r SELDLA-${ver}-win-x64.zip SELDLA-${ver}-win-x64 &
tar zvcf SELDLA-${ver}-linux-x64.tar.gz SELDLA-${ver}-linux-x64 &
tar zvcf SELDLA-${ver}-mac-x64.tar.gz SELDLA-${ver}-mac-x64 &
wait

scp SELDLA-${ver}-win-x64.zip SELDLA-${ver}-linux-x64.tar.gz SELDLA-${ver}-mac-x64.tar.gz sakura:www/software/SELDLA
