"""python 2.7 does not have timing module anymore as in version 2.4, this is to emulate the old timing module for the scream code"""

import time

timestart = 0.0
timefinish = 0.0

def start():
    global timestart 
    timestart = time.clock()
def finish():
    global timefinish
    timefinish = time.clock()
def micro():
    global timestart
    global timefinish
    return (timefinish-timestart)*1e6

    

