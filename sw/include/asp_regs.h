/* =====================================================================
 *  asp_regs.h  -  ASP AXI4-Lite register map (byte offsets)
 *
 *  Mirrors rtl/top/asp_axi_lite_regs.vhd.  Keep the two in lockstep: the
 *  PL is the source of truth; this header is the PS view of it.
 * ===================================================================== */
#ifndef ASP_REGS_H
#define ASP_REGS_H

#include <stdint.h>

/* Base address: set to the address map entry for asp_axi_lite_regs in
 * the Vivado address editor.  Override at compile time if needed. */
#ifndef ASP_BASE_ADDR
#define ASP_BASE_ADDR  0x43C00000u
#endif

#define ASP_REG_ID          0x00u
#define ASP_REG_CONTROL     0x04u
#define ASP_REG_STATUS      0x08u
#define ASP_REG_DWELL_CNT   0x0Cu
#define ASP_REG_LAMBDA0     0x10u
#define ASP_REG_LAMBDA1     0x14u
#define ASP_REG_LAMBDA2     0x18u
#define ASP_REG_LAMBDA3     0x1Cu
#define ASP_REG_DET_LHS_LO  0x20u
#define ASP_REG_DET_LHS_HI  0x24u
#define ASP_REG_DET_RHS_LO  0x28u
#define ASP_REG_DET_RHS_HI  0x2Cu
#define ASP_REG_W0_RE       0x30u
#define ASP_REG_W0_IM       0x34u
#define ASP_REG_W1_RE       0x38u
#define ASP_REG_W1_IM       0x3Cu
#define ASP_REG_W2_RE       0x40u
#define ASP_REG_W2_IM       0x44u
#define ASP_REG_W3_RE       0x48u
#define ASP_REG_W3_IM       0x4Cu
#define ASP_REG_TX_SHIFT    0x50u
#define ASP_REG_CLIP_CNT    0x54u
#define ASP_REG_SHIFT_INIT  0x58u
#define ASP_REG_IRQ_STATUS  0x5Cu
#define ASP_REG_IRQ_ENABLE  0x60u

#define ASP_ID_EXPECT       0x41535002u  /* "AS" + revision 2 */

/* CONTROL */
#define ASP_CTRL_ENABLE     (1u << 0)
#define ASP_CTRL_BYPASS     (1u << 1)
#define ASP_CTRL_RANK2_EN   (1u << 2)
#define ASP_CTRL_SOFT_RST   (1u << 3)

/* STATUS */
#define ASP_STAT_DETECTED   (1u << 0)
#define ASP_STAT_RANK_SHIFT 1
#define ASP_STAT_RANK_MASK  (3u << 1)
#define ASP_STAT_RANK2_EL   (1u << 3)
#define ASP_STAT_FIR_OVF    (1u << 4)

/* IRQ */
#define ASP_IRQ_DWELL       (1u << 0)

static inline uint32_t asp_read(uintptr_t base, uint32_t off)
{
    return *(volatile uint32_t *)(base + off);
}

static inline void asp_write(uintptr_t base, uint32_t off, uint32_t val)
{
    *(volatile uint32_t *)(base + off) = val;
}

#endif /* ASP_REGS_H */
