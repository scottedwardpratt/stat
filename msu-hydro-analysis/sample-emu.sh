#!/bin/sh
## sample the emulator using our test program
if [ ! -e $1/Emulator.statefile ]; then 
		echo "# no snapshot file, run train-emu.sh first"
		exit -1
fi

test_Emu++ $1/Emulator.statefile
