/* =====================================================================
 *  asp_driver.c  -  PS-side control for the ASP PL datapath
 *
 *  Bare-metal / freestanding C.  No libc heap.  Replace the MMIO macros
 *  in asp_regs.h if your BSP uses a cache-bypass accessor.
 * ===================================================================== */
#include "asp_driver.h"

int asp_init(asp_dev_t *dev, uintptr_t base)
{
    uint32_t id;

    if (dev == 0)
        return -1;

    dev->base        = base;
    dev->m_of_n_m    = 3;
    dev->m_of_n_n    = 5;
    dev->rank2_votes = 0;
    dev->vote_window = 0;

    id = asp_read(base, ASP_REG_ID);
    if (id != ASP_ID_EXPECT)
        return -1;

    /* Safe defaults: enabled off, bypass off, rank-2 gated off. */
    asp_write(base, ASP_REG_CONTROL, 0);
    asp_write(base, ASP_REG_IRQ_ENABLE, 0);
    asp_write(base, ASP_REG_IRQ_STATUS, ASP_IRQ_DWELL); /* W1C */
    asp_write(base, ASP_REG_SHIFT_INIT, 0);
    return 0;
}

void asp_enable(asp_dev_t *dev, bool on)
{
    uint32_t c = asp_read(dev->base, ASP_REG_CONTROL);
    if (on) c |= ASP_CTRL_ENABLE; else c &= ~ASP_CTRL_ENABLE;
    asp_write(dev->base, ASP_REG_CONTROL, c);
}

void asp_set_bypass(asp_dev_t *dev, bool on)
{
    uint32_t c = asp_read(dev->base, ASP_REG_CONTROL);
    if (on) c |= ASP_CTRL_BYPASS; else c &= ~ASP_CTRL_BYPASS;
    asp_write(dev->base, ASP_REG_CONTROL, c);
}

void asp_set_rank2(asp_dev_t *dev, bool on)
{
    uint32_t c = asp_read(dev->base, ASP_REG_CONTROL);
    if (on) c |= ASP_CTRL_RANK2_EN; else c &= ~ASP_CTRL_RANK2_EN;
    asp_write(dev->base, ASP_REG_CONTROL, c);
}

void asp_set_shift_init(asp_dev_t *dev, int8_t shift)
{
    asp_write(dev->base, ASP_REG_SHIFT_INIT, (uint32_t)(uint8_t)shift);
}

void asp_soft_reset(asp_dev_t *dev)
{
    uint32_t c = asp_read(dev->base, ASP_REG_CONTROL);
    asp_write(dev->base, ASP_REG_CONTROL, c | ASP_CTRL_SOFT_RST);
    /* Bit is self-clearing in the PL. */
}

void asp_irq_enable(asp_dev_t *dev, bool on)
{
    asp_write(dev->base, ASP_REG_IRQ_ENABLE, on ? ASP_IRQ_DWELL : 0u);
}

void asp_irq_clear(asp_dev_t *dev)
{
    asp_write(dev->base, ASP_REG_IRQ_STATUS, ASP_IRQ_DWELL);
}

void asp_read_telemetry(const asp_dev_t *dev, asp_telemetry_t *out)
{
    uint32_t lo, hi;
    unsigned i;

    out->status    = asp_read(dev->base, ASP_REG_STATUS);
    out->dwell_cnt = asp_read(dev->base, ASP_REG_DWELL_CNT);
    out->lambda[0] = (int32_t)asp_read(dev->base, ASP_REG_LAMBDA0);
    out->lambda[1] = (int32_t)asp_read(dev->base, ASP_REG_LAMBDA1);
    out->lambda[2] = (int32_t)asp_read(dev->base, ASP_REG_LAMBDA2);
    out->lambda[3] = (int32_t)asp_read(dev->base, ASP_REG_LAMBDA3);

    lo = asp_read(dev->base, ASP_REG_DET_LHS_LO);
    hi = asp_read(dev->base, ASP_REG_DET_LHS_HI);
    out->det_lhs = ((uint64_t)(hi & 0xffffu) << 32) | lo;

    lo = asp_read(dev->base, ASP_REG_DET_RHS_LO);
    hi = asp_read(dev->base, ASP_REG_DET_RHS_HI);
    out->det_rhs = ((uint64_t)(hi & 0xffffu) << 32) | lo;

    for (i = 0; i < 4; i++) {
        out->w_re[i] = (int32_t)asp_read(dev->base, ASP_REG_W0_RE + 8u * i);
        out->w_im[i] = (int32_t)asp_read(dev->base, ASP_REG_W0_IM + 8u * i);
    }

    out->tx_shift = (int8_t)(asp_read(dev->base, ASP_REG_TX_SHIFT) & 0xffu);
    out->clip_cnt = asp_read(dev->base, ASP_REG_CLIP_CNT);
}

void asp_dwell_policy(asp_dev_t *dev, const asp_telemetry_t *t)
{
    bool eligible = (t->status & ASP_STAT_RANK2_EL) != 0;

    /* Sliding M-of-N: count eligible votes in the last N dwells. */
    if (eligible)
        dev->rank2_votes++;
    dev->vote_window++;

    if (dev->vote_window >= dev->m_of_n_n) {
        asp_set_rank2(dev, dev->rank2_votes >= dev->m_of_n_m);
        dev->rank2_votes = 0;
        dev->vote_window = 0;
    }
}
