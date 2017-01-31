# generate prototypes for Samba C code
# originally writen by A. Tridge, June 1996
# this modified version was taken from rsync 2.4.6

# removed some things that I do not need in SiPSMG ( ltz, 2009 )
# modified by Chensong Zhang for FASP ( 04/04/2010 )
# modified by Chensong Zhang for FASP+Doxygen ( 05/22/2010 )
# modified by Xiaozhe Hu for BSR format ( 10/26/2010 )
# modified by Chensong Zhang for the new format ( 03/03/2011 )
# modified by Xiaozhe Hu for CUDA ( 09/20/2011 )
# modified by Chensong Zhang to update output format ( 01/31/2017 )

BEGIN {
  inheader=0;
  print "/*! \\file  " name
  print " *"
  print " *  \\brief Function decoration for the FASP package"
  print " *"
  print " *---------------------------------------------------------------------------------"
  print " *  Copyright (C) 2008--2017 by the FASP team. All rights reserved.                "
  print " *  Released under the terms of the GNU Lesser General Public License 3.0 or later."
  print " *---------------------------------------------------------------------------------"
  print " *"
  print " *  \\warning DO NOT EDIT!!! This file is automatically generated!"
  print " */ \n"

  print "#include \"fasp.h\" "
  print "#include \"fasp_block.h\" "
  print "#include \"fasp4cuda.h\" "
  print "#ifdef __cplusplus"
  print "extern \"C\" { "
  print "#endif" 
}

{
  if (inheader) {
    if (match($0,"[)][ \t]*$")) {
      inheader = 0;
      printf "%s;\n\n",$0;
    } else {
      printf "%s\n",$0;
    }
    next;
  }
}

/\/*! \\file/ {
    printf "\n/*-------- In file: %s --------*/\n\n",$3;  
}

/^static|^extern/ || !/^[a-zA-Z]/ || /[;]/ {
  next;
}

!/^INT|^REAL|^FILE|^OFF_T|^size_t|^off_t|^pid_t|^unsigned|^mode_t|^DIR|^user|^int|^char|^uint|^struct|^BOOL|^void|^double|^time|^dCSRmat|^dCOOmat|^dvector|^iCSRmat|^ivector|^AMG_data|^ILU_data|^dSTRmat|^dBSRmat|^dCSRLmat|^precond|^cudvector|^cuivector|^cudCSRmat/ {
  next;
}

/[(].*[)][ \t]*$/ {
    printf "%s;\n\n",$0;
    next;
}

/[(]/ {
  inheader=1;
  printf "%s\n",$0;
  next;
}

END {
  print "#ifdef __cplusplus"
  print "} "
  print "#endif"
  print "/* End of " name " */"
}
