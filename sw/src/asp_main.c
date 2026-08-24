/* =====================================================================
 *  asp_main.c  -  bring-up skeleton for the ASP on Zynq PS
 *
 *  Fill in the board-specific pieces (AD9361 SPI, MCS, GIC) from the
 *  Analog Devices HDL / no-OS tree for your carrier.  The sequence below
 *  matches docs/09-vhdl-implementation.md §9.
 *
 *  Build (example, bare-metal):
 *    arm-none-eabi-gcc -I sw/include -c sw/src/asp_driver.c
 *    arm-none-eabi-gcc -I sw/include -c sw/src/asp_main.c
 * ===================================================================== */
#include "asp_driver.h"

/* Stubs — replace with ADI no-OS / your BSP. */
extern int  ad9361_init_shared_lo(void);
extern int  ad9361_set_rx_lo_hz(unsigned long long hz);
extern int  ad9361_set_tx_lo_hz(unsigned long long hz);
extern int  ad9361_run_calibration_and_freeze_dc(void);
extern int  ad9361_mcs_sync(void);
extern void platform_enable_pl_irq(void (*isr)(void *), void *arg);
extern int8_t measure_rx_shift_init(void);

static asp_dev_t g_asp;
static volatile uint32_t g_dwell_seen;

static void asp_dwell_isr(void *arg)
{
    asp_dev_t *dev = (asp_dev_t *)arg;
    asp_telemetry_t t;

    asp_read_telemetry(dev, &t);
    asp_dwell_policy(dev, &t);
    asp_irq_clear(dev);
    g_dwell_seen++;
}

int main(void)
{
    int8_t shift0;

    /* 1. RF front end: shared LO, MCS, freeze DC tracking. */
    if (ad9361_init_shared_lo() != 0)
        return 1;
    if (ad9361_set_rx_lo_hz(1577466000ull) != 0)  /* 1577.466 MHz */
        return 2;
    if (ad9361_set_tx_lo_hz(1573374000ull) != 0)  /* 1573.374 MHz */
        return 3;
    if (ad9361_run_calibration_and_freeze_dc() != 0)
        return 4;
    if (ad9361_mcs_sync() != 0)
        return 5;

    /* 2. PL register file. */
    if (asp_init(&g_asp, ASP_BASE_ADDR) != 0)
        return 6;

    shift0 = measure_rx_shift_init();
    asp_set_shift_init(&g_asp, shift0);

    /* 3. Interrupt + enable datapath. */
    platform_enable_pl_irq(asp_dwell_isr, &g_asp);
    asp_irq_enable(&g_asp, true);
    asp_enable(&g_asp, true);

    for (;;) {
        /* Userspace logger must not starve the 1 kHz ISR. */
    }
}
