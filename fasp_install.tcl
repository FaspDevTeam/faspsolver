#!/usr/bin/wish
proc InstallStatus { s } {
global status
set status $s
update idletasks
}

proc Help {} {
    global prompthelp
    if {[info exists prompthelp(ok)]} {unset prompthelp(ok)}
    set f [toplevel .prompthelp -borderwidth 10 -width 30]
    set b [frame $f.buttons -bd 10]
    ScrolledText $f 80 20
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
    if {[catch {open "FASP_GUI.help.txt" r} in1]} {
	puts stderr "Can not read help file Install_FASP.help.txt" 
    } else {
	while { [gets $in1 line] >=0 } {
#	    puts $line
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

proc DistClean { } {
    global build_dir
    catch "exec /bin/rm -r $build_dir" c3
#    puts [list "C3 is =" $c3]
#    set in [list {open} {|find . -name '*~'} r]
#    puts $in
#    catch [eval $in] result
#    puts $result
#    while {[gets $in line] >=0}  {
#	catch "exec ls -l $line" result
#	puts $line 
#	puts $result
#    }
}

proc  Uninstall { } {
global build_dir
    catch "exec xargs /bin/rm < $build_dir/install_manifest.txt" c1
    catch "exec /bin/rm -rf $build_dir/install_manifest.txt doc/htdocs" c2
    DistClean
}


proc Config { f } {
    global build_dir
    global c_options
    global option_lbl
    global maxlen
## Build dir name should be set here. 
    DistClean

## Create Build_Dir now
    if {[catch {exec /bin/mkdir -p "$build_dir"}]} {
	InstallError "Cannot create $dname. Exiting "
    }

    foreach el { verbose shared doxywizard openmp } {
	puts [list $el $c_options($el)]
    }

    set config_flags [list "-DCMAKE_VERBOSE_MAKEFILE=OFF" "-DCMAKE_RULE_MESSAGES=ON"]
    if {$c_options(verbose)} {
	set  config_flags [list "-DCMAKE_VERBOSE_MAKEFILE=ON" "-DCMAKE_RULE_MESSAGES=ON"]
    }
    if {$c_options(shared)} {
	set config_flags [lappend config_flags "-DSHARED=$c_options(shared)"]
    }
    if {$c_options(openmp)} {
	set config_flags [lappend config_flags "-DUSE_OPENMP=$c_options(openmp)"]
    }
    if {$c_options(doxywizard)} {
	set config_flags [lappend config_flags "-DDOXYWIZARD=$c_options(doxywizard)"]
    }	

    #CONFIG_FLAGS+=-DCMAKE_C_FLAGS=$(cflags)
    #CONFIG_FLAGS+=-DCMAKE_CXX_FLAGS=$(cxxflags)
    #CONFIG_FLAGS+=-DCMAKE_Fortran_FLAGS=$(fflags)

#### SHELL command 

    set currdir [pwd]
    catch "cd $build_dir" c2
    puts $c2
    set xcommand(cmake)  "|cmake .. $config_flags"
    puts [list $xcommand(cmake) [pwd]]
    set in  [eval {open} [list $xcommand(cmake) r]]
    $f configure -state normal
    while {[gets $in line] >=0}  {
	$f insert end "$line\n"
	update idletasks
	$f see end
    }
    $f configure -state disabled 
    catch "exec cd $currdir" c2
}
proc RunIt { f command0 } {
    global build_dir
    set flag 1
    if {![file exists $build_dir/Makefile]} {
     puts "Run Config first"
    } else {
	switch -exact -- $command0 {
	    docs { set xcommand(cmake) "|make -C $build_dir $command0" }
	    headers { set xcommand(cmake) "|make -C $build_dir $command0" }
	    install { set xcommand(cmake) "|make -C $build_dir $command0" }
	    default { set flag 0 } 
	}
	puts [list $xcommand(cmake) $flag]
	if { $flag == 1 } {
	    set in  [eval {open} [list $xcommand(cmake) r]]
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
proc InstallError { x } {
	    puts stderr {*****************************************************************}
	    puts stderr "ERROR: $x"
	    puts stderr {*****************************************************************}
	    exit
}
proc GetValue { whattofocus prompt01 } {
    global prompt0
    set f [toplevel .prompt0 -borderwidth 10 -width 200]
    message $f.msg -aspect 1000 -text $prompt01
    set b [frame $f.buttons -bd 10]
    button $b.ok -text Yes -command {set prompt0(ok) 1}
    button $b.cancel -text No -command {set prompt0(ok) 0}
    pack $b.ok $b.cancel -side left
    bind $b.ok <Tab> [list focus $b.cancel]
    bind $b.cancel <Tab> [list focus $b.ok]
    focus $b.$whattofocus
    pack $f.msg -side top
    pack $f.buttons -side top -fill x
    grab $f
    tkwait variable prompt0(ok)
    grab release $f
    destroy $f
}

proc ShowOpts { f } {
    global c_options
    global maxlen
    set option_lbl(verbose) "Verbose output \?"
    set option_lbl(shared) "Build shared library \?"
    set option_lbl(doxywizard) "Use Doxygen GUI (if found)\?"
    set option_lbl(openmp) "Build with OpenMP support \?"
    
    foreach el { verbose shared doxywizard openmp } {
	set c_options($el) 0
	set maxlen -1
	set len [expr [string length $option_lbl($el)] + 2]
	if { $maxlen < $len } {
	    set maxlen $len
	}
    }
    foreach el { verbose shared doxywizard openmp } {
	frame $f.$el
	label $f.$el.l -relief flat -width $maxlen -anchor w \
	    -text $option_lbl($el) -border 0 \
	    -font  -*-Helvetica-bold-r-*--*-140-* 
	radiobutton $f.$el.y -variable c_options($el) \
	    -text "Yes" -value 1 -state normal
	radiobutton $f.$el.n -variable c_options($el) \
	    -text "No" -value 0 -state normal 
	pack $f.$el.l -side left
	pack $f.$el.n -side right
	pack $f.$el.y -side right
	pack $f.$el -side top
    }
}
    #proc SortItems {} {
    #   global items0 files0
    #  set c [lsort -command mycompare $items0]
    # set d {}
    # foreach el $c {
    #	set ix [lsearch -exact $items0 $el]
    #	if { $ix >= 0 } { 
#	    set items0 [eval [list lreplace $items0 $ix $ix "\*$el"]]
#	    lappend d [lindex $files0 $ix]
#	} else {
#	    puts stderror "Could not find item in {SortItems}"
#	    exit
#	} 
#  }
#    set items0 $c
#    set files0 $d
#}
# end data Chapter 
global aposition
global prompt0
global build_dir
global c_options
global maxlen


set prompt0(ok) 0
set aposition -1

set dname [join [list BUILD FASP [eval [list exec uname -m]]-[eval  [list exec uname -s]]] "_"]
set build_dir [join [list [pwd] {/} $dname] ""]  

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

#frame .f.images.text 
#    .f.images.text insert end "Whatever Whatever Whatever Whatever\n"
#pack .f.images.text -side right

pack .f  -expand true -fill both
update idletasks

# GUI WINDOWS SETTINGS

pack .f.buttons -side top -fill x 

##########main 


ScrolledText1 .f 80 20
pack .f.text  -side right -fill both
    .f.text configure -state disabled
#   grab release .f.text


button .f.buttons.help -text Help -command {Help}
button .f.buttons.quit -text Quit -command {GentleExit}
button .f.buttons.config -text "Config" -command {Config .f.text}
button .f.buttons.install -text "Install" -command {RunIt .f.text install}
button .f.buttons.uninstall -text "Uninstall" -command {Uninstall}
button .f.buttons.docs -text "HTML docs" -command {RunIt .f.text docs}
button .f.buttons.headers -text "Headers" -command {RunIt .f.text headers}

pack .f.buttons.quit .f.buttons.config .f.buttons.install .f.buttons.docs \
    .f.buttons.headers .f.buttons.uninstall  -side left

pack .f.buttons.help -side right

frame .f.opts -borderwidth 10
ShowOpts .f.opts
pack .f.opts -side left -fill y



#pack .f.f1 -expand true -side left -anchor w -fill both
#entry .f.e.estatus -textvar status -relief flat -width 56 -background "#00008B" -foreground "#FFFFFF" -font { Helvetica -12 bold } -state disabled
#frame .f.e.dirs
#pack .f.e.estatus -side bottom -fill x



