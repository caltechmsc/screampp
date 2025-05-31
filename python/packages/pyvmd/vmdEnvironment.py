'''vcvicek setup for vmd in python
the default display is WIN unless VMDDISPLAYDEVICE is set (then it can be TEXT)'''

import os
# there is some funny bussiness with environmental variables
# just using putenv does not work, so I use this
# note that there might be problems with this on MAC
def putenv(var,value):
    os.environ[var] = value

vmddir = os.getenv('VASEKINSTALL') + '/vmd'
putenv('VMDDIR',vmddir)
if os.getenv('VMDDISPLAYDEVICE') is None:
    putenv('VMDDISPLAYDEVICE','WIN')  # TEXT or WIN
putenv('SURF_BIN',vmddir+'/surf_LINUXAMD64')
putenv('TACHYON_BIN',vmddir+'/tachyon_LINUXAMD64')
putenv('STRIDE_BIN',vmddir+'/stride_LINUXAMD64')

# bash version
# export VMDDIR=$VASEKINSTALL/vmd
# export VMDDISPLAYDEVICE=WIN  # TEXT or WIN
# export SURF_BIN=$VMDDIR/surf_LINUXAMD64
# export TACHYON_BIN=$VMDDIR/tachyon_LINUXAMD64
# export STRIDE_BIN=$VMDDIR/stride_LINUXAMD64

def setWin():
    putenv('VMDDISPLAYDEVICE','WIN')  # TEXT or WIN

def setText():
    putenv('VMDDISPLAYDEVICE','TEXT')  # TEXT or WIN

