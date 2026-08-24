# ASP PS software skeleton

Bare-metal C for the Zynq-7000 PS side of the anti-spoofing datapath.

| File | Role |
|------|------|
| `include/asp_regs.h` | AXI4-Lite register map (mirrors `asp_axi_lite_regs.vhd`) |
| `include/asp_driver.h` | Driver API |
| `src/asp_driver.c` | MMIO helpers, telemetry, M-of-N rank-2 policy |
| `src/asp_main.c` | Bring-up order: AD9361 → MCS → SHIFT_INIT → IRQ → enable |

The AD9361 / GIC stubs in `asp_main.c` must be wired to the Analog Devices
no-OS (or your BSP) for the carrier in use.  Register offsets and the
1 kHz ISR contract are fixed by the PL and are documented in
`docs/09-vhdl-implementation.md` §9.
