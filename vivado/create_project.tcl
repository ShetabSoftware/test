# =====================================================================
#  create_project.tcl  -  Vivado project for the ASP datapath (XC7Z020)
#
#  Usage (from the repo root, with Vivado on PATH):
#
#    vivado -mode batch -source vivado/create_project.tcl
#    vivado -mode gui   -source vivado/create_project.tcl
#
#  This script creates the project, adds every RTL source and the timing
#  XDC, and prints the block-design IP checklist.  It does NOT invent a
#  board-specific pinout: the AD9361 interface and Zynq PS wiring belong
#  in a block design that matches the carrier in use (ZC702, custom, …).
#
#  After the project exists, open it and run vivado/bd_checklist.tcl for
#  the exact IP instances and settings required by the RTL headers.
# =====================================================================

set repo_root [file normalize [file join [file dirname [info script]] ..]]
set proj_dir  [file join $repo_root build vivado_asp]
set proj_name asp_xc7z020

file mkdir $proj_dir

create_project $proj_name $proj_dir -part xc7z020clg484-1 -force
set_property target_language VHDL [current_project]
set_property simulator_language Mixed [current_project]

# ------------------------------------------------------------------ RTL
set pkg_files [list \
  [file join $repo_root rtl pkg asp_pkg.vhd] \
  [file join $repo_root rtl pkg asp_coef_pkg.vhd] \
]

set core_files [glob -nocomplain [file join $repo_root rtl core *.vhd]]
set top_files  [list \
  [file join $repo_root rtl top asp_datapath.vhd] \
  [file join $repo_root rtl top asp_axi_lite_regs.vhd] \
  [file join $repo_root rtl top asp_top.vhd] \
]

add_files -norecurse [concat $pkg_files $core_files $top_files]
set_property library work [get_files $pkg_files]
set_property file_type {VHDL} [get_files [concat $pkg_files $core_files $top_files]]
set_property top asp_top [current_fileset]
update_compile_order -fileset sources_1

# Packages first: Vivado's default order is not always topological for
# VHDL-93 packages that every core file uses.
reorder_files -fileset sources_1 \
  [file join $repo_root rtl pkg asp_pkg.vhd] \
  [file join $repo_root rtl pkg asp_coef_pkg.vhd]

# ------------------------------------------------------------------ XDC
add_files -fileset constrs_1 -norecurse \
  [file join $repo_root constraints asp_timing.xdc]

puts "=============================================================="
puts " ASP Vivado project created: $proj_dir/$proj_name.xpr"
puts " Part: xc7z020clg484-1 (edit -part for your package/speed)"
puts " Top:  asp_top"
puts ""
puts " Next:"
puts "   1. Open the project in the GUI"
puts "   2. Build the Zynq PS + AD9361 block design for your board"
puts "   3. Source vivado/bd_checklist.tcl and tick every REQUIRED IP"
puts "   4. Run synth_design / implementation; review WNS on the three"
puts "      formerly-long paths (now pipelined — see docs/09 §8)"
puts "=============================================================="
