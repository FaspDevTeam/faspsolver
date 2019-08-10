#!/usr/bin/wish
########################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
# Top level tcl install script. Installs the libraries, tutorial and
# test programs
#
########################################################################

proc InstallError { x } {
	    puts stderr {*****************************************************************}
	    puts stderr "ERROR: $x"
	    puts stderr {*****************************************************************}
	    exit
}

proc InstallStatus { s } {
global status
set status $s
update idletasks
}

proc InstallHelper { s } {
global helper
set helper $s
update idletasks
}

proc Help {} {
    global prompthelp
    if {[info exists prompthelp(ok)]} {unset prompthelp(ok)}
    set f [toplevel .prompthelp -borderwidth 10 -width 30]
    set b [frame $f.buttons -bd 10]
    ScrolledText $f 72 24
    button $b.close -text close -command {set prompthelp(ok) 1}
    $f.text configure \
	    -font -*-Helvetica-bold-r-*--*-120-* \
	    -background white -foreground black -state normal
    label $f.label -relief flat\
	    -text "FASP install help" -border 0 \
	    -font  -*-Helvetica-bold-r-*--*-120-* 
    pack $b.close -side left
    pack $f.label -side top
    pack $f.buttons -side top
    pack $f.yscroll -side right -fill y
    pack $f.label  -side top -fill x -anchor e
    pack $f.text  -side top -fill both -expand true
    pack $f.xscroll -side bottom -fill x
    pack $f.buttons -side bottom
    if {[catch {open "FASP_GUI_help.txt" r} in1]} {
	InstallHelper "Can not read help file FASP_GUI_help.txt" 
    } else {
	while { [gets $in1 line] >=0 } {
	    $f.text insert end "$line\n"
	}
	close $in1
    }
    $f.text configure -state disabled
    grab $f
    tkwait variable prompthelp(ok)
    grab release $f
    if {$prompthelp(ok)} {destroy $f}
}

proc ScrolledText { parent width height } {
# This proc follows closely the Brent's book
    catch {frame $parent}
    text $parent.text
    $parent.text configure -setgrid true \
	    -xscrollcommand [list $parent.xscroll set] \
	    -yscrollcommand [list $parent.yscroll set] \
	    -wrap word -width $width -height $height
    scrollbar $parent.xscroll -orient horizontal \
	    -command [list $parent.text xview]
    scrollbar $parent.yscroll -orient vertical \
	    -command [list $parent.text yview]
}

proc ScrolledText1 { parent width height } {
# This proc follows closely the Brent's book
    catch {frame $parent}
    text $parent.text 
    $parent.text configure -setgrid true \
	    -xscrollcommand [list $parent.xscroll set] \
	    -yscrollcommand [list $parent.yscroll set] \
	    -wrap word -width $width -height $height \
            -borderwidth 5 -relief sunken
    scrollbar $parent.xscroll -orient horizontal \
	    -command [list $parent.text xview]
    scrollbar $parent.yscroll -orient vertical \
	    -command [list $parent.text yview]
    pack $parent.yscroll -side right -fill y
    pack $parent.xscroll -side bottom -fill x
}

proc ItemsSelect { w y } {
#    $w select set anchor [$w nearest $y]
    set ix [$w nearest $y]
    $w select set $ix
    $w see $ix
}

proc FindFile { startDir namePat } {
    global filel
# Find file: Unix chapter: Brent Welch http://www.beedub.com/book/
# here we glob for files only in the current dir pwd
    set pwd [pwd]
    if [catch {cd $startDir} err] {
	InstallError [list "Error in find file" $err]
    }
    foreach match [glob -nocomplain -- $namePat] {
	lappend filel $startDir/$match
    }
    foreach file [glob -nocomplain *] {
	if [file isdirectory $file] {
	    FindFile $startDir/$file $namePat
	}
    }
    cd $pwd
}

proc DistClean { } {
    global build_dir
    global build_type
    global status
    global filel
    set filel {}
    update idletasks
    set current_dir [pwd]
    if {[file isdirectory "$current_dir/BUILD_FASP"]} {
	catch "exec /bin/rm -rf $current_dir/BUILD_FASP" c3
    }  
    if {[info exists build_dir]} {unset build_dir}
    #### complicated delete of emacs backup files
    eval [list FindFile $current_dir *~]
    foreach n [list $filel] {
	catch "exec /bin/rm -rf $n" c4	     
    }
}

proc  Uninstall { } {
global status
global build_dir
    update idletasks
    set current_dir [pwd]
    set build_dir "$current_dir/BUILD_FASP"
    if {[file isdirectory $build_dir]} {
	    if {[file exists $build_dir/install_manifest.txt ]} {
	        catch "exec xargs /bin/rm < $build_dir/install_manifest.txt" c1
	    }
	    catch "exec /bin/rm -rf $build_dir/install_manifest.txt doc/htdocs" c2
    } 
    DistClean
    InstallHelper "DONE Uninstalling FASP."
}

proc Config { f } {
    global build_dir
    global build_type
    global c_options
    global compiler_flags
    global extra_lib_path
    global option_lbl
    global maxlen
    global status
## Create Build dir (name should be set here). 
    DistClean
    set build_dir BUILD_FASP
    if {[catch {exec /bin/mkdir -p "$build_dir"}]} {
	InstallError "Cannot create $build_dir. Exiting "
    }
#
    InstallHelper "Configuring..."

    set cmake_flags [list "-DCMAKE_BUILD_TYPE=$build_type"]
##   
    if {$c_options(verbose)} {
	set  cmake_flags [lappend cmake_flags \
			       "-DCMAKE_VERBOSE_MAKEFILE=ON" \
			       "-DCMAKE_RULE_MESSAGES=ON"]
    } else {
	set cmake_flags [lappend cmake_flags \
			     "-DCMAKE_VERBOSE_MAKEFILE=OFF" \
			     "-DCMAKE_RULE_MESSAGES=ON"]
    }
    if {$c_options(shared)} {
	set cmake_flags [lappend cmake_flags "-DSHARED=$c_options(shared)"]
    }
    if {$c_options(openmp)} {
	set cmake_flags [lappend cmake_flags "-DUSE_OPENMP=$c_options(openmp)"]
    }
    if {$c_options(doxywizard)} {
	set cmake_flags [lappend cmake_flags "-DDOXYWIZARD=$c_options(doxywizard)"]
    }
    if {$c_options(umfpack)} {
	set cmake_flags [lappend cmake_flags "-DUSE_UMFPACK=$c_options(umfpack)"]
	if {![regexp {^[ |\t]*$} [join $extra_lib_path(suitesparse_dir)]]} {
	    set cmake_flags \
	    	[concat $cmake_flags  " -DSUITESPARSE_DIR=\\\"$extra_lib_path(suitesparse_dir)\\\""]
    	}
    }	
##
    if {![regexp {^[ |\t]*$} [join $compiler_flags(c)]]} {
	set cmake_flags \
	    [concat $cmake_flags " -DADD_CFLAGS=\\\"$compiler_flags(c)\\\""]
    } 
    if {![regexp {^[ |\t]*$} [join $compiler_flags(cxx)]]} {
	set cmake_flags \
	    [concat $cmake_flags " -DADD_CXXFLAGS=\\\"$compiler_flags(cxx)\\\""]
    }
    if {![regexp {^[ |\t]*$} [join $compiler_flags(fortran)]]} {
	set cmake_flags \
	    [concat $cmake_flags  " -DADD_FFLAGS=\\\"$compiler_flags(fortran)\\\""]
    }
    
#### SHELL command 

##   This needs to be changed for Windows
    set xcommand(cmake)  "|/bin/sh -c  \"cd $build_dir && cmake .. $cmake_flags\" "
#    puts [list $xcommand(cmake) [pwd] ]
    set in  [eval {open} [list $xcommand(cmake) r]]
    $f configure -state normal
    while {[gets $in line] >=0}  {
	$f insert end "$line\n"
	update idletasks
	$f see end
    }
    $f configure -state disabled 
    InstallHelper "DONE Configuring FASP."
}

proc RunIt { f command0 } {
    global status
    global build_dir
    if {![info exists build_dir]} {
	InstallHelper "Build Dir does not exist or is not a directory; Run Config First!"
	InitOpts
	update idletasks
    } elseif { ![regexp {^BUILD_FASP$}  $build_dir]} {
	InstallHelper "Build Dir has incorrect value=\"$build_dir\"; Run Config first!"
	InitOpts
	update idletasks
    }	else {	
	set flag 1
	if {![file exists $build_dir/Makefile]} {
	    InstallHelper "File $build_dir/Makefile does not exist. Run Config first!"
	    update idletasks
	} else {
	    InstallHelper "Installing ... please be patient..."
	    switch -exact -- [CheckInpData $command0] {
		docs { set xcommand(cmake) "|make -C $build_dir $command0 2>@ stdout"}		
		headers { set xcommand(cmake) "|make -C $build_dir $command0 2>@ stdout"}
		install { set xcommand(cmake) "|make -C $build_dir $command0 2>@ stdout"}
		default { set flag 0 } 
	    }
	    if { $flag == 1 } {
		catch "open [list $xcommand(cmake)] r" in
		$f configure -state normal
		while {[gets $in line] >=0}  {
		    $f insert end "$line\n"
		    update idletasks
		    $f see end
		}
		close $in
		$f configure -state disabled 
	    }
	}
	InstallHelper "DONE Installing FASP."
    }
}

proc GentleExit {} {
    exit
}

proc CheckInpData { x } {
    set y {}
    regsub "^\[ \t\]*" $x {} y
    regsub "\[ \t\]*\$" $y {} y
    regsub -all -- "\[ \t\]\[ \t\]*" $y { } y
    return $y
}

proc InitOpts { } {
#Initialize options
    global c_options
    global option_lbl
    global compiler_flags
    global status
    global maxlen
    global build_type
    global build_dir
    global extra_lib_path
    set build_dir {}
    set option_lbl(verbose) "Verbose output \?"
    set option_lbl(shared) "Build shared library \?"
    set option_lbl(doxywizard) "Use Doxygen GUI (if found)\?"
    set option_lbl(openmp) "Build with OpenMP support \?"
    set option_lbl(umfpack) "Build with UMFPACK support \?"
    set build_type "RELEASE"    
    set compiler_flags(c) {}
    set compiler_flags(cxx) {}
    set compiler_flags(fortran) {}
    set extra_lib_path(suitesparse_dir) {}
    foreach el { verbose shared doxywizard openmp umfpack} {
	set c_options($el) 0
	set maxlen -1
	set len [expr [string length $option_lbl($el)] + 2]
	if { $maxlen < $len } {
	    set maxlen $len
	}
    }
}

proc ShowOpts { f } {
    global c_options
    global maxlen
    global build_type
    global option_lbl
    foreach el { verbose shared doxywizard openmp umfpack} {
	frame $f.$el
	label $f.$el.l -relief flat -width $maxlen -anchor w \
	    -text $option_lbl($el) -border 0 \
	    -font  -*-Helvetica-bold-r-*--*-140-* 
	radiobutton $f.$el.y -variable c_options($el) \
	    -text "Yes" -value 1 -state normal
	radiobutton $f.$el.n -variable c_options($el) \
	    -text "No" -value 0 -state normal 
	pack $f.$el.l -side left 
	pack $f.$el.n -side right -fill both
	pack $f.$el.y -side right -fill both
	pack $f.$el -side top -fill both
    }
## buttons for the Build type
    frame $f.build_type
    label $f.build_type.l -relief flat   -width $maxlen -anchor w \
	-text "Build type" -border 0 \
	-font  -*-Helvetica-bold-r-*--*-140-*
    radiobutton $f.build_type.debug -variable build_type \
	-text "Debug" -value DEBUG -state normal
    radiobutton $f.build_type.release -variable build_type \
	-text "Release" -value RELEASE -state normal 
    pack $f.build_type.l -side left 
    pack $f.build_type.release -side right
    pack $f.build_type.debug -side right
    pack $f.build_type -side top
}

proc EntryList { parent } {
    global compiler_flags
    global extra_lib_path

    text $parent.text 
    $parent.text insert end "The Debug build type assumes -g as a compiler flag. Release build type assumes -O3 flag. Aditional compiler flags could be added below."
    $parent.text configure -width 45 -height 4  -wrap word \
	-font -*-Helvetica-normal-r-*--*-140-* \
	-background white -foreground black -relief flat -state disabled 	
	
	text $parent.text1 
    $parent.text1 insert end "Specify the path to SUITESPARSE (optional, use it when you want to use specific UMFPACK)."
    $parent.text1 configure -width 45 -height 2  -wrap word \
	-font -*-Helvetica-normal-r-*--*-140-* \
	-background white -foreground black -relief flat -state disabled 

    foreach n { c cxx fortran } {
	frame $parent.e_$n 
	label $parent.e_$n.l -width 20 -anchor w \
		-foreground black \
		-font -*-Helvetica-bold-r-*--*-140-*
	entry $parent.e_$n.e \
		-relief sunken \
		-foreground black -width 25 \
		-font -*-Helvetica-bold-r-*--*-140-* 
    
    }
    $parent.e_c.l configure -text "C flags"
    $parent.e_c.e configure -text [list $compiler_flags(c)] \
	    -textvariable compiler_flags(c)  -relief sunken 
    
    $parent.e_cxx.l configure -text "CXX flags"
    $parent.e_cxx.e configure -text $compiler_flags(cxx) \
	-textvariable compiler_flags(cxx)  -relief sunken 
    
    $parent.e_fortran.l configure -text "F flags"
    $parent.e_fortran.e configure -text  $compiler_flags(fortran) \
	    -textvariable compiler_flags(fortran) -relief sunken 
	    
	frame $parent.suitesparse
	label $parent.suitesparse.l -width 20 -anchor w \
		-foreground black \
		-font -*-Helvetica-bold-r-*--*-140-*
	entry $parent.suitesparse.e \
		-relief sunken \
		-foreground black -width 25 \
		-font -*-Helvetica-bold-r-*--*-140-* 
		
	$parent.suitesparse.l configure -text "SUITESPARSE"
    $parent.suitesparse.e configure -text [list $extra_lib_path(suitesparse_dir)] \
	    -textvariable extra_lib_path(suitesparse_dir) -relief sunken

	pack $parent.suitesparse.l -side left
	pack $parent.suitesparse.e -side left
	pack $parent.suitesparse -side bottom
	
	pack $parent.text1 -side bottom -fill both 
		$parent.text1 configure -state disabled

    foreach n {cxx fortran c } {
	pack $parent.e_$n.l -side left 
	pack $parent.e_$n.e -side left 
	pack $parent.e_$n -side bottom 
    }
	pack $parent.text -side bottom -fill x
	
} 

global build_dir
global c_options
global maxlen
global status
global helper

#set dname [join [list BUILD [eval [list exec uname -m]]-[eval  [list exec uname -s]]] "_"]
#set build_dir [join [list [pwd] {/} $dname] ""]  

frame .f
frame .f.buttons -borderwidth 10 

frame .f.images -background "#00008B"
pack .f.images -side top -fill x
#set im [ image create photo -file doc/fasp_small.gif ]
#label .f.images.logo -background "#FFFFFF" -image $im 
label .f.images.logo -background "#00008B" -text "FASP" \
    -background "#00008B" -foreground yellow -font { Helvetica -24 bold } 
pack .f.images.logo -side left
label .f.images.text -text  "FAST AUXILIARY SPACE PRECONDITIONING\nInstallation script" \
    -background "#00008B" -foreground "#FFFFFF" -font { Helvetica -18 bold } 
pack .f.images.text 

pack .f  -expand true -fill both
update idletasks

# GUI WINDOWS SETTINGS

pack .f.buttons -side top -fill x 

##########main 

## Fortran Style: 72 Chars per line :)
ScrolledText1 .f 72 24
pack .f.text  -side right -fill both
    .f.text configure -state disabled
    
#   grab release .f.text

button .f.buttons.help -text "Help" -command {Help}
button .f.buttons.quit -text "Quit" -command {GentleExit}
button .f.buttons.config -text "Config" -command {Config .f.text}
button .f.buttons.install -text "Install" -command {RunIt .f.text install}
button .f.buttons.uninstall -text "Uninstall" -command {Uninstall}
button .f.buttons.docs -text "HTML Docs" -command {RunIt .f.text docs}
button .f.buttons.headers -text "Headers" -command {RunIt .f.text headers}

pack .f.buttons.quit .f.buttons.config .f.buttons.install .f.buttons.docs \
     .f.buttons.headers .f.buttons.uninstall -side left

pack .f.buttons.help -side right

InitOpts 

frame .f.opts -borderwidth 10
label .f.opts.lhelper -textvar helper -relief flat -width 55 \
    -background "#00008B" -foreground "yellow" -font { Helvetica -14 bold } 

##label .f.opts.lstatus -textvar status -relief raised -width 35  

pack .f.opts.lhelper -side bottom -fill both
##pack .f.opts.lstatus -side bottom 

ShowOpts .f.opts

EntryList .f.opts

pack .f.opts -side left -fill y
