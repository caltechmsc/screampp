# GENERAL
setenv biosoft /project/Biogroup/Software
set ARCH = `uname -m`
if ($ARCH == x86_64) then
  set path = ( /exec/python/pythonEPD-7.0-2-rh3-x86_64/bin $path )
  alias python /exec/python/pythonEPD-7.0-2-rh3-x86_64/bin/python
else
  set path = ( /exec/python/pythonEPD-7.0-2-rh3-x86/bin $path )
  alias python /exec/python/pythonEPD-7.0-2-rh3-x86/bin/python
endif

### Required SCREAM Variables
setenv SCREAM_NEW     /project/Biogroup/Software/scream3/
setenv SCREAM_NEW_CHG amber
setenv SCREAM_NEW_LIB ${SCREAM_NEW}/lib/amber/
setenv SCREAM_NEW_CNN ${SCREAM_NEW}/lib/cnn/
setenv SCREAM_NEW_RTF ${SCREAM_NEW}/lib/rtf/

### Aliases
alias  scream      ${SCREAM_NEW}/python/scream.py
alias  screamwrap  ${SCREAM_NEW}/python/scream_wrap.py
alias  screammulti ${SCREAM_NEW}/python/scream_multi.py

### Paths
if ($ARCH == x86_64) then
  setenv PYTHONPATH ${SCREAM_NEW}/build/lib.linux-x86_64-2.7:${SCREAM_NEW}/python/packages:${PYTHONPATH}
else
  setenv PYTHONPATH ${SCREAM_NEW}/build/lib.linux-i686-2.7:${SCREAM_NEW}/python/packages:${PYTHONPATH}
endif
set path = ( ${SCREAM_NEW}/python $path )

# Print Help
echo "SCREAM usage"
echo "scream      {scream.par}"
echo "screamwrap  {bgf} {rotlib} {residue list}"
echo "screammulti {bgf} {rotlib} {# to output} {resiude list}"
