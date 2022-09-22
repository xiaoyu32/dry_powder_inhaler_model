#!/bin/bash

if test "`id -u`" -ne 0
	then 
	echo "ERROR: You need to run this script as root" 
	exit 
fi

echo "checking for old installations..."
if [[ -L /usr/bin/lpp ]]; then
	echo "ERROR: you need to uninstall old lpp versions first."
	echo "you can do by so executing"
	echo "sudo rm /usr/bin/lpp /usr/bin/pizza"
	exit 1
fi

echo "setting symbolic links..."
#ln -s $(pwd)/src/lpp.py /usr/bin/lpp
ln -s $(pwd)/src/lpp.py /opt/local/bin/lpp

echo "checking installation..."
echo "which lpp"
echo $(which lpp)
