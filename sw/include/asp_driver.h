/* =====================================================================
 *  asp_driver.h  -  minimal PS driver for the ASP PL datapath
 * ===================================================================== */
#ifndef ASP_DRIVER_H
#define ASP_DRIVER_H

#include <stdint.h>
#include <stdbool.h>
#include "asp_regs.h"

typedef struct {
    uintptr_t base;
    /* M-of-N hysteresis for rank-2 enable (PS policy, not PL). */
    unsigned  m_of_n_m;
    unsigned  m_of_n_n;
    unsigned  rank2_votes;
    unsigned  vote_window;
} asp_dev_t;

typedef struct {
    uint32_t status;
    uint32_t dwell_cnt;
    int32_t  lambda[4];
    uint64_t det_lhs;
    uint64_t det_rhs;
    int32_t  w_re[4];
    int32_t  w_im[4];
    int8_t   tx_shift;
    uint32_t clip_cnt;
} asp_telemetry_t;

/* Returns 0 on success, -1 if ID register does not match. */
int  asp_init(asp_dev_t *dev, uintptr_t base);

void asp_enable(asp_dev_t *dev, bool on);
void asp_set_bypass(asp_dev_t *dev, bool on);
void asp_set_rank2(asp_dev_t *dev, bool on);
void asp_set_shift_init(asp_dev_t *dev, int8_t shift);
void asp_soft_reset(asp_dev_t *dev);

void asp_irq_enable(asp_dev_t *dev, bool on);
void asp_irq_clear(asp_dev_t *dev);

void asp_read_telemetry(const asp_dev_t *dev, asp_telemetry_t *out);

/* Call from the 1 kHz dwell ISR after reading telemetry.
 * Updates CONTROL.rank2_en from an M-of-N vote on STATUS.rank2 eligible.
 * Defaults: M=3, N=5 (set in asp_init). */
void asp_dwell_policy(asp_dev_t *dev, const asp_telemetry_t *t);

#endif /* ASP_DRIVER_H */
