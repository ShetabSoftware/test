# =====================================================================
#  bd_checklist.tcl  -  required Vivado IPs for asp_top
#
#  Source inside an open project after the block design exists:
#    source vivado/bd_checklist.tcl
#
#  Every entry matches a "-- REQUIRED VIVADO IP CORE" comment in the RTL.
#  This script only REPORTS; it does not create the BD, because board
#  pinout and AD9361 LVDS mapping are carrier-specific.
# =====================================================================

puts "===== ASP block-design IP checklist ====="
puts ""
puts "1. Clocking Wizard (clk_wiz)"
puts "   - Input:  AD9361 DATA_CLK (board-specific rate)"
puts "   - Output: clk_dsp = 130.944 MHz EXACT (4*FS_ADC = 8*FS_WORK)"
puts "   - Jitter optimized; buffered output"
puts "   - create_generated_clock comes from the IP XDC — do not duplicate"
puts ""
puts "2. AXI Clock Converter (axi_clock_converter)"
puts "   - PROTOCOL      = AXI4LITE"
puts "   - ADDR_WIDTH    = 12"
puts "   - DATA_WIDTH    = 32"
puts "   - ID_WIDTH      = 0"
puts "   - ASYNC_CLK     = 1"
puts "   - SI clock: FCLK_CLK0 (PS, typically 100 MHz)"
puts "   - MI clock: clk_dsp (130.944 MHz)"
puts "   - S_AXI <- PS GP0 via AXI Interconnect"
puts "   - M_AXI -> asp_top / asp_axi_lite_regs"
puts ""
puts "3. Zynq-7000 PS (processing_system7)"
puts "   - Enable M_AXI_GP0"
puts "   - Enable IRQ_F2P (level, active-high) for the dwell interrupt"
puts "   - FCLK_CLK0 for the AXI interconnect / clock converter SI side"
puts "   - Do NOT use FCLK as clk_dsp (rate plan requires AD9361 derivation)"
puts ""
puts "4. AXI Interconnect"
puts "   - One master (PS GP0) to the clock-converter SI"
puts ""
puts "5. Two axi_ad9361 instances (Analog Devices)"
puts "   - Shared LO, MCS-synchronised"
puts "   - RX LO 1577.466 MHz, TX LO 1573.374 MHz"
puts "   - RX bandwidth ~10 MHz; sample rate -> FS_ADC = 32.736 MHz"
puts "   - Freeze DC tracking after calibration"
puts "   - Four antenna IQ streams + common valid into asp_top"
puts ""
puts "6. Optional: AXI DMA"
puts "   - Capture path for bring-up / logging only"
puts "   - Not required for the null-steering datapath"
puts ""
puts "7. Do NOT instantiate FIR Compiler for the shaping filter"
puts "   - asp_fir_shape is bit-exact against the golden model; replacing"
puts "     it with FIR Compiler breaks co-simulation unless coefficients,"
puts "     rounding and saturation are proven identical"
puts ""
puts "Register map: see rtl/top/asp_axi_lite_regs.vhd header and sw/include/asp_regs.h"
puts "Software bring-up: sw/src/asp_main.c"
puts "=========================================="
